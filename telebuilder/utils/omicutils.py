# -*- coding: utf-8 -*-

import sys
import re
from collections import OrderedDict

import urllib.request, urllib.error, urllib.parse

from .utils import tsv


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
        return '\n'.join('%s%d' % (k.ljust(30), v) for k,v in self.d.items())

def get_sequence_ucsc(build, region):
    d = {
        'baseurl': 'http://genome.ucsc.edu/cgi-bin/das',
        'build': build,
        'region': '%s:%d,%d' % region,
    }
    url = '%(baseurl)s/%(build)s/dna?segment=%(region)s' % d

    response = urllib.request.urlopen(url)
    xml = response.read()
    m = re.search('<DNA length="(\d+)">([a-zA-Z\n]+)</DNA>', xml)
    assert m is not None
    reported_length = int(m.group(1))
    seq = str(m.group(2)).replace('\n', '')
    assert len(seq) == reported_length, 'Sequence length != reported length'
    assert reported_length == region[2]-region[1]+1, 'Calculated length != reported length'
    return seq

def get_sequence_fasta(seqdict, region):
    py_start = int(region[1]) - 1
    py_end = int(region[2])
    if region[0] not in seqdict:
        print('Skipping, sequence %s not loaded' % region[0], file=sys.stderr)
        return None
    seq = str(seqdict[region[0]][py_start : py_end].seq)
    assert len(seq) == region[2]-region[1]+1, 'Calculated length (%d) != reported length (%d)' % (len(seq), region[2]-region[1]+1)
    return seq
