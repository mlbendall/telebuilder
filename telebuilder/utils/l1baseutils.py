# -*- coding: utf-8 -*-


from .utils import numstr, tsv
from .gtfutils import GTFLine

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


class L1BaseLine(object):
    COLS = [
        ('chrom', str),
        ('chromStart', int),
        ('chromEnd', int),
        ('name', str),
        ('score', numstr),
        ('strand', str),
        ('thickStart', int),
        ('thickEnd', int),
        ('itemRgb', str),
    ]

    def __init__(self, row):
        for (n, t), v in zip(self.COLS, row):
            setattr(self, n, t(v))

    def to_gtf(self):
        _g = GTFLine()
        _g.chrom = self.chrom  # chrom
        _g.source = 'l1base'  # source
        _g.feature = 'exon'  # feature
        _g.start = self.chromStart + 1  # start (add 1 since GTF is 1-based)
        _g.end = self.chromEnd  # end
        _g.score = self.score  # score
        _g.strand = self.strand  # strand
        _g.frame = '.'  # frame
        _g.attr = {
            'exon_s': self.thickStart+1,
            'exon_e': self.thickEnd,
            'l1base_id': self.name,
        }
        return _g

    def fmt(self):
        return [str(getattr(self, n)) for n, t in self.COLS]

    def __str__(self):
        return '\t'.join(self.fmt())


def read_l1base_file(infile, comment='#'):
    for r in tsv(infile, comment):
        if r[0].startswith('track') or r[0].startswith('browser'):
            continue
        yield(L1BaseLine(r))

