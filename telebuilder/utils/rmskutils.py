# -*- coding: utf-8 -*-

import sys
import re
from collections import defaultdict, Counter
from subprocess import Popen, PIPE

from utils import numstr
from gtfutils import GTFLine

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


class RMSKLine(object):
    COLS = [
        ('bin', int),
        ('swScore', numstr),
        ('milliDiv', int),
        ('milliDel', int),
        ('milliIns', int),
        ('genoName', str),
        ('genoStart', int),
        ('genoEnd', int),
        ('genoLeft', int),
        ('strand', str),
        ('repName', str), 
        ('repClass', str),
        ('repFamily', str),
        ('repStart', int),
        ('repEnd', int),
        ('repLeft', int),
        ('id', int), 
    ]
    
    def __init__(self, row):
        for (n, t), v in zip(self.COLS, row):
            setattr(self, n, t(v))

    def fix_model_coords(self, modellen):
        """ Fix the alignment positions on the model
        
        This adjusts repStart, repEnd, and repLeft to be consistent with the
        model length.
        """
        if self.strand == '+':
            trueend = modellen + self.repLeft
            if trueend != self.repEnd:
                # print str(self)            
                replen = self.repEnd - self.repStart
                self.repEnd = trueend
                self.repStart = trueend - replen
        else:
            trueend = modellen + self.repStart
            if trueend != self.repEnd:
                # print str(self)
                replen = self.repEnd - self.repLeft
                self.repEnd = trueend
                self.repLeft = trueend - replen
        
    def to_gtf(self):
        _g = GTFLine()
        _g.chrom = self.genoName        # chrom
        _g.source = 'rmsk'              # source
        _g.feature = 'exon'             # feature
        _g.start = self.genoStart + 1   # start (add 1 since GTF is 1-based)
        _g.end = self.genoEnd           # end
        _g.score = self.swScore         # score
        _g.strand = self.strand         # strand
        _g.frame = '.'                  # frame
        _g.attr = {
            'repStart': self.repStart,
            'repEnd': self.repEnd,
            'repLeft': self.repLeft,
            'repName': self.repName,
            'repClass': self.repClass,
            'repFamily': self.repFamily,
        }
        return _g
    
    def fmt(self):
        return [str(getattr(self, n)) for n,t in self.COLS]
    
    def __str__(self):
        return '\t'.join(self.fmt())        


def guess_rmsk_model_lengths(rls):
    ''' Calculate the length of models from attributes '''
    ret = defaultdict(list)
    for rl in rls:
        sz = (rl.repEnd - rl.repLeft) if rl.strand == '+' else (rl.repEnd - rl.repStart)
        ret[rl.repName].append(sz)
    return {k:Counter(v).most_common()[0][0] for k,v in ret.iteritems()}


RMSK_MODEL_QUERY = '''SELECT * FROM rmsk WHERE repName = '%s';'''

def ucsc_download(build, query, outh=PIPE):
    cmd = 'mysql -h genome-mysql.cse.ucsc.edu -u genome -D %s -A -e "%s"' % (build, query)
    print >>sys.stderr, cmd
    p1 = Popen(cmd, shell=True, stdout=outh, stderr=PIPE)
    o, e = p1.communicate()
    if e: print >>sys.stderr, e
    return o
