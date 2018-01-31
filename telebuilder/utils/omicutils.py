# -*- coding: utf-8 -*-

import re
from collections import OrderedDict

import urllib2

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

def get_sequence_ucsc(build, region):
    d = {
        'baseurl': 'http://genome.ucsc.edu/cgi-bin/das',
        'build': build,
        'region': '%s:%d,%d' % region,
    }
    url = '%(baseurl)s/%(build)s/dna?segment=%(region)s' % d

    response = urllib2.urlopen(url)
    xml = response.read()
    m = re.search('<DNA length="(\d+)">([a-zA-Z\n]+)</DNA>', xml)
    assert m is not None
    reported_length = int(m.group(1))
    seq = str(m.group(2)).replace('\n', '')
    assert len(seq) == reported_length, 'Sequence length != reported length'
    assert reported_length == region[2]-region[1]+1, 'Calculated length != reported length'
    return seq
