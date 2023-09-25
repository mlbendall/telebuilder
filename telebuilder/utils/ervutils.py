import sys
import os
from tempfile import NamedTemporaryFile

from .igvutils import IGV
from .utils import raw_input_stderr, overlap_length

_QUIET = False

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2023 Matthew L. Bendall"

PROMPT_HELP = '''Options are:
    ignore                     - Ignore conflict. Annotation is unchanged.
    reject loc1[,loc2,...]     - Remove loci from annotation.
    merge loc1+loc2[+loc3+...] - Merge loci into single locus. The original loci are
                                 removed from annotation and replaced with merged locus.
                                 Attributes are inherited from the first locus.
    diff loc1%loc2             - Remove the part of loc1 that overlaps with loc2. Original
                                 loc1 is removed from annotation and loc2 remains
                                 unchanged. In case loc2 is contained by loc1, two loci
                                 are created.
'''

RESOLVE_ACTIONS = ['ignore','reject','diff','merge','y']

def prompt_cmd(can_accept = True):
    ''' Get command from user input '''
    operation = 'begin'
    while operation.lower() not in RESOLVE_ACTIONS:
        curprompt = ''
        if operation not in ['begin','?']:
            curprompt += 'INPUT IS INVALID'.ljust(20)
        else:
            curprompt += ' '.ljust(20)
        if can_accept:
            curprompt += '***Action to take (? for help, Y to accept): '
        else:
            curprompt += '***Action to take (REQUIRED, ? for help): '
        z = raw_input_stderr(curprompt).strip()
        if z == '?':
            print(PROMPT_HELP, file=sys.stderr)
        inputcmd = z.strip().split()
        if len(inputcmd) == 1 or len(inputcmd) == 2:
            operation = inputcmd[0]
        else:
            operation = ''
    
    assert operation.lower() in RESOLVE_ACTIONS
    return inputcmd


def resolve_merge(congtf, locstr, quiet=_QUIET):
    if locstr is None:
        locids = [g.attr['locus'] for g in congtf]
    else:
        locids = locstr.split('+')
    
    to_remove = [g for g in congtf if g.attr['locus'] in locids]

    merged = to_remove[0].copy()
    for h in to_remove[1:]:
        for _m in h.members:
            hm = _m.copy()
            if 'locus' in hm.attr:
                hm.attr['locus'] = merged.attr['locus']
            merged.add(hm)
    merged.cleanup()
    return to_remove, [merged,]


class LocusConflict(object):
    def __init__(self, gtf):
        self.gtf = gtf
        self.reason = 'unknown'
        self.action = None
        self.locstr = ''
        self.is_auto = False
        # self.locids = [g.attr['locus'] for g in self.gtf]
    
    def region_str(self, padding=1000):
        rstart = max(1, min([h.start for h in self.gtf]) - padding)
        rend = max([h.end for h in self.gtf]) + padding
        return '%s:%d-%d' % (self.gtf[0].chrom, rstart, rend)
    
    def igv_gtf(self):
        ret = ''
        for g in self.gtf:
            _g = g.copy()
            _g.set_attr('transcript_id', '%s' % _g.attr['locus'])
            ret += '\n'.join(str(m) for m in _g.members) + '\n'
        return ret.strip('\n')
    
    def display_conflict_text(self):
        ret = '#%s#\n' % ('-' * 80)
        ret += '### %s%s%s%s\n' % ('Locus'.ljust(20), 'Category'.ljust(10), 'Model'.ljust(15), 'Model Cov.')
        for g in self.gtf:
            ret += '### %s%s%s%s\n' % (g.attr['locus'].ljust(20), g.attr['category'].ljust(10),
                                     g.attr['intModel'].ljust(15), g.attr['model_pct'])
        ret += '#%s#\n' % ('-' * 80)
        for g in self.gtf:
            ret += '%s\n' % str(g)
        ret += '#%s#\n' % ('-' * 80)        
        return ret
    
    def display_solution_text(self, auto=False):
        if self.action is not None:
            ret = 'AutoInspect identified "%s". ' % self.reason
            ret += 'Solution: ' if auto else 'Suggested solution: '
            ret += '%s %s' % (self.action, self.locstr)
        else:
            ret = 'AutoInspect failed to find solution.'
        return ret

    def display_review_text(self):
        if self.action is not None:
            ret = 'Reason: %s' % self.reason
            if self.is_auto:
                ret += ' (auto). '
            else:
                ret += '. '
            ret += 'Solution: '
            ret += '%s %s' % (self.action, self.locstr)
        else:
            ret = 'AutoInspect failed to find solution.'
        return ret

    def inspect_auto(self):
        def _is_tandem(g1, g2):
            left, right = (g1, g2) if g1.start < g2.start else (g2, g1)
            left_last = left.members[-1]
            right_first = right.members[0]
            if left_last.attr['geneRegion'] == 'ltr' and right_first.attr['geneRegion'] == 'ltr':
                if left_last.start == right_first.start and left_last.end == right_first.end:
                    return ('merge', '%s+%s' % (left.attr['locus'], right.attr['locus']))
            return False
        
        def _is_insert(g1, g2):
            larger, smaller = (g1, g2) if g1.length() > g2.length() else (g2, g1)    
            if larger.contains(smaller, stranded=False):
                if smaller.contains(larger, stranded=False) is False:
                    return ('merge', '%s+%s' % (larger.attr['locus'], smaller.attr['locus']))
            return False

        def _is_mergeable(*gtfs):
            _gtfs =  sorted(gtfs, key=lambda x:x.length(), reverse=True)
            strandset = set(g.strand for g in _gtfs)
            if len(strandset) == 1:
                return ('merge', '+'.join(g.attr['locus'] for g in _gtfs))
            return False

        def _is_opposite(g1, g2):
            g1, g2 = (g1, g2) if g1.attr['model_cov'] > g2.attr['model_cov'] else (g2, g1)
            if g1.strand != g2.strand and g1.overlap(g2, False) < 100:
                return ('merge', '%s+%s' % (g2.attr['locus'], g1.attr['locus']))
                # return ('diff', '%s%%%s' % (g2.attr['locus'], g1.attr['locus']))
            return False

        def _is_largegap(g1, g2):
            gaps = [g1.members[i+1].start - m1.end for i,m1 in enumerate(g1.members[:-1])]
            if max(gaps) > 1000:
                return ('diff', '%s%%%s' % (g1.attr['locus'], g2.attr['locus']))
            gaps = [g2.members[i + 1].start - m1.end for i, m1 in enumerate(g2.members[:-1])]
            if max(gaps) > 1000:
                return ('diff', '%s%%%s' % (g2.attr['locus'], g1.attr['locus']))
            return False

        if len(self.gtf) == 2:
            _funs = [('tandem', _is_tandem),
                     ('insert', _is_insert),
                     ('mergeable', _is_mergeable),
                     ('opposite', _is_opposite),
                     ('largegap', _is_largegap),
                     ]
        elif len(self.gtf) > 2:
            _funs = [('mergeable', _is_mergeable), ]

        for reason, fun in _funs:
            if fun(*self.gtf):
                self.reason, (self.action, self.locstr)  = reason, fun(*self.gtf)
                self.is_auto = True
                break

    def update_action(self, input):
        if input[0] in 'Yy':
            if self.action is None: # There was no auto assignment
                return self.update_action(prompt_cmd(can_accept = False))
            return

        assert input[0] in RESOLVE_ACTIONS, 'ERROR: unknown input: %s' % str(input)
        self.reason = 'user_curated'
        self.is_auto = False
        self.action = input[0]
        if len(input) == 2:
            self.locstr = input[1]
        else:
            sep = {'reject': ',', 'merge': '+', 'diff': '%', }        
            if self.action == 'ignore':
                self.locstr = ''
            elif self.action in sep:
                self.locstr = sep[self.action].join(g.attr['locus'] for g in self.gtf)
            else:
                assert False, "unknown input: %s" % str(input)
    
    def _resolve_ignore(self):
        return ([], [])

    def _resolve_reject(self):
        locs = self.locstr.split(',')
        to_remove = []
        for loc in locs:
            f = [g for g in self.gtf if g.attr['locus'] == loc]
            assert len(f) == 1, 'ERROR: locus %s was found %d times.' % (loc, len(f))
            to_remove.append(f[0])
        
        return (to_remove, [])

    def _resolve_merge(self):
        locs = self.locstr.split('+')
        assert len(locs) >= 2, 'ERROR: found %d loci, must have 2 or more.' % len(locs)
        to_remove = []
        for loc in locs:
            f = [g for g in self.gtf if g.attr['locus'] == loc]
            assert len(f) == 1, 'ERROR: locus %s was found %d times.' % (loc, len(f))
            to_remove.append(f[0])

        to_remove.sort(key=lambda x:x.attr['model_pct'], reverse=True)
        merged = to_remove[0].copy()
        for other in to_remove[1:]:
            for m in other.members:
                m_copy = m.copy()
                if 'locus' in m_copy.attr:
                    m_copy.attr['locus'] = merged.attr['locus']
                if m_copy.strand != merged.strand:
                    m_copy.attr['prev_strand'] = m_copy.strand
                    m_copy.strand = merged.strand
                merged.add(m_copy)
        merged.cleanup()
        return (to_remove, [merged,])
    
    def _resolve_diff(self):
        ''' Remove the part of A that overlaps with B. B remains unchanged '''
        locs = self.locstr.split('%')
        assert len(locs) == 2, 'ERROR: found %d loci, must have exactly 2.' % len(locs)
        locA = [g for g in self.gtf if g.attr['locus'] == locs[0]]
        locB = [g for g in self.gtf if g.attr['locus'] == locs[1]]
        assert len(locA) == 1 and len(locB) == 1, 'ERROR: loci could not be identified:\n%s\n%s' % (str(locA), str(locB))
        locA,locB = (locA[0], locB[0])
        
        if not locA.overlap(locB, stranded=False):
            return ([], [])            # No overlap. Do not change locA or B
        if locB.contains(locA, stranded=False):
            print('locB contains locA', file=sys.stderr)
            return ([locA, ], [])      # locA is contained by locB. Remove locA.
        
        if locA.contains(locB, stranded=False):
            print('locA contains locB', file=sys.stderr)
            copyA_L = locA.copy()            
            tmpB = locB.span()
            tmpB.end = copyA_L.end + 1
            copyA_L.subtract(tmpB, stranded=False)
            copyA_L.set_attr('locus', '%sL' % copyA_L.attr['locus'])

            copyA_R = locA.copy()
            tmpB = locB.span()
            tmpB.start = copyA_R.start - 1
            copyA_R.subtract(tmpB, stranded=False)           
            copyA_R.set_attr('locus', '%sR' % copyA_R.attr['locus'])
            
            return ([locA, ], [copyA_L, copyA_R, ])
        else:
            print('another arrangement', file=sys.stderr)        
            copyA = locA.copy()
            copyA.subtract(locB, stranded=False)
            return ([locA, ], [copyA, ])
    
    def resolve(self):
        ''' Return a tuple with list of objects to remove and append '''
        assert self.action in RESOLVE_ACTIONS, 'ERROR: unknown action: %s' % str(self.action)
        if self.action == 'ignore':
            return self._resolve_ignore()
        elif self.action == 'reject':
            return self._resolve_reject()
        elif self.action == 'merge':
            return self._resolve_merge()
        elif self.action == 'diff':
            return self._resolve_diff()
            
    def __str__(self):
        return self.display_conflict_text()
    
def inspect_conflicts(conflicts, auto, igv, dest=''):
    lcons = [LocusConflict(cgroup) for cgroup in conflicts]
    
    if igv:
        tmp_gtf = os.path.join(dest, 'tmp.conflict.gtf')
        with open(tmp_gtf, 'w') as outh:
            for lcon in lcons:
                print(lcon.igv_gtf(), file=outh)
        
        igv.load(os.path.abspath(tmp_gtf))
        igv.collapse().expand('tmp.conflict.gtf')
    
    for i,lcon in enumerate(lcons):
        print(str(lcon), file=sys.stderr)
        if igv: igv.goto(lcon.region_str())
        lcon.inspect_auto()
        print(lcon.display_solution_text(auto), file=sys.stderr)

        if not auto:
            lcon.update_action(prompt_cmd())
        else:
            if lcon.action is None:
                lcon.update_action(prompt_cmd(can_accept=False))

    return lcons
