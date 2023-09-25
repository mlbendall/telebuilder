# -*- coding: utf-8 -*-

from .utils import numstr, attrstr, overlap_length

from . import gtfutils

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2023 Matthew L. Bendall"


class GTFLine(object):
    COLS = [
        ('chrom', str, '.'),
        ('source', str, '.'),
        ('feature', str, '.'),
        ('start', int, -1),
        ('end', int, -1),
        ('score', numstr, -1),
        ('strand', str, '.'),
        ('frame', str, '.'),
        ('attr', attrstr, {}),
    ]
    # Attributes that should be at the front of attribute string
    ATTRORDER = ['gene_id', 'transcript_id', 'locus', 'repName']

    def __init__(self, val=None):
        if val is None:
            for (n, t, d) in self.COLS:
                setattr(self, n, d)
            self.attr = {}
        elif type(val) is list:
            for (n, t, d), v in zip(self.COLS, val):
                setattr(self, n, t(v))
        elif type(val) is GTFLine:
            for (n, t, d) in self.COLS:
                setattr(self, n, getattr(val, n))

    def copy(self):
        return GTFLine(self.fmt())

    def length(self):
        return self.end - self.start

    def _fmt_attrs(self):
        if not self.attr:
            return '.'
        keyorder = [k for k in self.ATTRORDER if k in self.attr]
        keyorder += sorted(k for k in self.attr.keys() if k not in self.ATTRORDER)
        ret = ['%s "%s";' % (k, self.attr[k]) for k in keyorder]
        return ' '.join(ret)

    def fmt(self):
        ret = [str(getattr(self, n)) for (n, t, d) in self.COLS[:8]]
        return ret + [self._fmt_attrs()]

    def sequence(self, gdict):
        if self.strand == '+':
            return gdict[self.chrom][self.start:self.end]
        else:
            return gdict[self.chrom][self.start:self.end].reverse_complement()

    def asdict(self):
        ret = {k: v for k, v in self.attr.items()}
        ret.update({n: getattr(self, n) for (n, t, d) in self.COLS[:8]})
        return ret

    def overlap(self, other, stranded=True):
        if self.chrom != other.chrom:
            return 0
        if stranded and self.strand != other.strand:
            return 0
        return overlap_length((self.start, self.end), (other.start, other.end))

    def contains(self, other, stranded=True):
        if self.chrom != other.chrom:
            return False
        if stranded and self.strand != other.strand:
            return False
        return self.start <= other.start and self.end >= other.end

    def __str__(self):
        return '\t'.join(self.fmt())


class GTFCluster(object):
    def __init__(self, val=None):
        self.members = []
        self.chrom = '.'
        self.start, self.end = (-1, -1)
        self.strand = '.'
        self.attr = {}
        if val is not None:
            self.add(val)

    def copy(self):
        ret = GTFCluster([g.copy() for g in self.members])
        assert ret.chrom == self.chrom
        assert ret.start == self.start
        assert ret.end == self.end
        assert ret.strand == self.strand
        ret.attr = self.attr.copy()
        return ret

    def length(self):
        return self.end - self.start

    def add(self, other):
        """

        Args:
            other:

        Returns:

        """
        if type(other) is list:  # Provided list of GTFLine
            if type(other[0]) is not GTFLine:
                raise TypeError("List must contain GTFLine objects")
            self.members += other
        elif type(other) is GTFLine:
            self.members.append(other)
        elif type(other) is type(self):
            self.members += other.members
        else:
            raise TypeError("Argument type must be GTFLine or GTFCluster")

        # Sort members
        self.members = sorted(set(self.members), key=lambda x: x.start)

        # Update chrom, start, end
        chromset = set(h.chrom for h in self.members)
        if len(chromset) == 1:
            # All on same chrom
            self.chrom = chromset.pop()
            self.start = min(h.start for h in self.members)
            self.end = max(h.end for h in self.members)
        else:
            self.chrom = '.'
            self.start, self.end = (-1, -1)
        # Update strand
        strandset = set(h.strand for h in self.members)
        self.strand = strandset.pop() if len(strandset) == 1 else '.'

    def subtract(self, other, stranded=True):
        if self.chrom != other.chrom:
            return
        if stranded and self.strand != other.strand:
            return

        ospan = other if type(other) is GTFLine else other.span()
        new_members = []
        for m in self.members:
            l, r = gtfutils.subtract_gtflines(m, ospan)
            if l.length(): new_members.append(l)
            if r.length(): new_members.append(r)
        self.members = sorted(new_members, key=lambda x: x.start)
        self.start = min(m.start for m in self.members)
        self.end = max(m.end for m in self.members)

    def cleanup(self):
        new_members = [self.members[0].copy()]
        for m in self.members[1:]:
            l, r = gtfutils.subtract_gtflines(m, new_members[-1])
            if l.length(): new_members.append(l)
            if r.length(): new_members.append(r)
        self.members = sorted(new_members, key=lambda x: x.start)

    def overlap(self, other, stranded=True):
        if type(other) is GTFLine:
            return self.span().overlap(other, stranded)
        elif type(other) is type(self):
            return self.span().overlap(other.span(), stranded)

    def contains(self, other, stranded=True):
        if type(other) is GTFLine:
            return self.span().contains(other, stranded)
        elif type(other) is type(self):
            return self.span().contains(other.span(), stranded)
        assert False

    def span(self):
        return GTFLine([
            self.chrom,
            ','.join(sorted(set(h.source for h in self.members))),
            'gene',
            self.start,
            self.end,
            '.',
            self.strand,
            '.',
            self.attr,
        ])

    def set_attr(self, k, v, set_members=True):
        self.attr[k] = v
        if set_members:
            for h in self.members:
                h.attr[k] = v

    def del_attr(self, k):
        return self.attr.pop(k, None)

    def set_attr_from_members(self, k, newk=None, operation='unique'):
        newk = newk if newk is not None else k
        all_v = [h.attr[k] for h in self.members if k in h.attr]
        if operation == 'unique':
            newv = ','.join(sorted(set(all_v)))
        elif operation == 'first':
            newv = all_v[0]
        elif operation == 'most':
            newv = Counter(all_v).most_common()[0][0]
        elif operation == 'simplify':
            newv = ','.join(simplify_list(all_v))
        self.attr[newk] = newv

    def region_str(self, padding=0):
        ''' Return a region string, i.e. <chrom>:<start>-<end> '''
        return '%s:%d-%d' % (
        self.chrom, max(1, self.start - padding), self.end + padding)

    def display_str(self, comment=True, span=True):
        ret = []
        if comment:
            if 'locus' in self.attr:
                ret.append('### %s ###\t\t\t\t\t\t\t\t' % self.attr['locus'])
            else:
                ret.append('###---###\t\t\t\t\t\t\t\t')
        if span:
            ret.append(str(self.span()))
        ret += [str(m) for m in self.members]
        return '\n'.join(ret)

    def __str__(self):
        return self.display_str()

