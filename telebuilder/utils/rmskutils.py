# -*- coding: utf-8 -*-

import sys
import re
from collections import defaultdict, Counter
from subprocess import Popen, PIPE

from .utils import numstr
from .gtfutils import GTFLine

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
    return {k:Counter(v).most_common()[0][0] for k,v in ret.items()}


RMSK_MODEL_QUERY = '''SELECT * FROM rmsk WHERE repName = '%s';'''

def ucsc_download(build, query, outh=PIPE):
    cmd = 'mysql -h genome-mysql.cse.ucsc.edu -u genome -D %s -A -e "%s"' % (build, query)
    print(cmd, file=sys.stderr)
    p1 = Popen(cmd, shell=True, stdout=outh, stderr=PIPE)
    o, e = p1.communicate()
    if e: print(e, file=sys.stderr)
    return o

class RepeatMaskerOut(object):

    def __init__(self, l):
        self.line = l.strip('\n')
        fields = self.line.strip().split()
        self.bitscore = int(fields[0])
        self.pctdiv = float(fields[1])
        self.pctdel = float(fields[2])
        self.pctint = float(fields[3])
        self.chrom = fields[4]
        self.start = int(fields[5])
        self.end = int(fields[6])
        self.left = int(fields[7].strip('()'))
        self.strand = '-' if fields[8] == 'C' else fields[8]
        self.repName = fields[9]
        if '/' in fields[10]:
            self.repClass, self.repFamily = fields[10].split('/')
        else:
            self.repClass = self.repFamily = fields[10]
        if self.strand == '+':
            self.repStart = int(fields[11])
            self.repEnd = int(fields[12])
            self.repLeft = int(fields[13].strip('()'))
        else:
            self.repStart = int(fields[13])
            self.repEnd = int(fields[12])
            self.repLeft = int(fields[11].strip('()'))

    def to_gtf(self):
        _g = GTFLine()
        _g.chrom = self.chrom  # chrom
        _g.source = 'RepeatMasker'  # source
        _g.feature = 'similarity'  # feature
        _g.start = self.start
        _g.end = self.end  # end
        _g.score = self.bitscore  # score
        _g.strand = self.strand  # strand
        _g.frame = '.'  # frame
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
        return self.line

    def __str__(self):
        return self.fmt()






