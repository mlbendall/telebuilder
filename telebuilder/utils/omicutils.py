# -*- coding: utf-8 -*-

from collections import OrderedDict

from utils import tsv

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


class ChromosomeDict(object):
    def __init__(self, fh=None):
        if fh is None:
            self.reflen = OrderedDict()
            self.reforder = None
        else:
            self.reflen = OrderedDict()
            for row in tsv(fh):
                self.reflen[row[0]] = int(row[1])

            self.reforder = list(self.reflen.keys())

    def __str__(self):
        return '\n'.join('%s%d' % (k.ljust(30), v) for k,v in self.d.iteritems())
