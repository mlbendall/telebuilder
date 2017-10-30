#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import string
from glob import glob
from collections import defaultdict

from utils import _get_build_from_file
from utils.gtfutils import sort_gtf, intersect_gtf
from utils.gtfutils import read_gtf_file, write_gtf_file
from utils.omicutils import ChromosomeDict
from utils import l1baseutils

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

# Alphabet suffixes
SUFFIXES = list(string.letters[:26])
# In case there are more than 26:
SUFFIXES += [a+b for b in SUFFIXES for a in SUFFIXES]

def namelocs(locs, cytogtf, cdict):
    locs = sort_gtf(locs, cdict.reforder)
    byband = defaultdict(list)
    for g, bands in intersect_gtf(locs, [cytogtf,], stranded=False):
        chrom = g.chrom[3:] if g.chrom[:3] == 'chr' else g.chrom
        if bands:
            band = bands[0][1].attr['gene_id']
        else:
            band = ''
        byband[(chrom, band)].append(g)

    for (chrom, band), gtfs in byband.iteritems():
        if len(gtfs) == 1:
            suffixes = ['']
        else:
            suffixes = SUFFIXES[:len(gtfs)]
            # Add parentheses if band is not present
            if band == '':
                suffixes = ['(%s)' % _ for _ in suffixes]

        # Set the locus name, transcript_id, and gene_id
        for g, suf in zip(gtfs, suffixes):
            name = '%s_%s%s%s' % (g.attr['category'], chrom, band, suf)
            g.attr['locus'] = name
            g.attr['transcript_id'] = name
            g.attr['gene_id'] = name
    return locs

def main(args):
    # Load chromosome sizes
    if os.path.exists(args.chrom_sizes):
        cdict = ChromosomeDict(args.chrom_sizes)
    else:
        print >>sys.stderr, '[WARNING] Chromosome sizes were not found'
        cdict = ChromosomeDict()

    # Load cytoband
    if os.path.exists(args.cytoband):
        cytoband = list(read_gtf_file(args.cytoband))
    else:
        cytoband = None
        print >> sys.stderr, '[WARNING] Cytoband was not found'

    ''' Stage 1: Convert L1Base tracks to GTF.'''
    print >>sys.stderr, '*** Stage 1: Converting L1Base tracks to GTF.'
    all_locs = defaultdict(list)
    groups = [('L1FLI', 'flil1'),
              ('L1ORF2', 'orf2l1'),
              ('L1FLnI', 'flnil1'),
    ]
    bedfiles = glob(os.path.join(args.track_dir, '*.bed'))
    for bf in bedfiles:
        group = [gn for gn, abbr in groups if re.search(abbr, bf)]
        assert len(group) == 1
        group = group[0]
        for l1line in l1baseutils.read_l1base_file(bf):
            g = l1line.to_gtf()
            g.attr['category'] = group
            all_locs[group].append(g)

    for group, locs in all_locs.iteritems():
        print '%s: %d' % (group, len(locs))

    ''' Stage 2: Remove Redundant annotations'''
    print >>sys.stderr, '*** Stage 2: Removing redundant annotations.'
    for idx_a in range(len(groups)-1, 0, -1):
        groupA = groups[idx_a][0]
        if groupA not in all_locs: continue
        to_remove = []
        gtfA = all_locs[groupA] #sort_gtf(all_locs[groupA], cdict.reforder)
        groupBs = [groups[b][0] for b in range(idx_a)]
        gtfBs = [all_locs[gB] for gB in groupBs if gB in all_locs]
        for g, overlaps in intersect_gtf(gtfA, gtfBs, stranded=False):
            if overlaps:
                to_remove.append(g)
        for g in to_remove:
            all_locs[groupA].remove(g)

    for group, locs in all_locs.iteritems():
        print '%s: %d' % (group, len(locs))

    print >>sys.stderr, '*** Stage 8: Naming loci'
    final_locs = {}
    for group, locs in all_locs.iteritems():
        final_locs[group] = namelocs(locs, cytoband, cdict)


    for group, locs in final_locs.iteritems():
        with open(os.path.join(args.outdir, '%s.gtf' % group), 'w') as outh:
            write_gtf_file(sort_gtf(locs, cdict.reforder), outh)


def console():
    import argparse
    parser = argparse.ArgumentParser(
        description='''Construct L1 annotations from L1 Base'''
    )
    group1 = parser.add_argument_group('Input')
    group1.add_argument('--genome_build', default=_get_build_from_file(),
                        help='''Genome build to use, i.e. hg19, hg38.
                                Required. Can be configured by creating a text
                                file called "build.txt" with the build ID.''')
    group1.add_argument('--chrom_sizes', default='chrom.sizes',
                        help='''File with chromosome sizes and order. Default
                                is "chrom.sizes" (in working directory).''')
    group1.add_argument('--cytoband', default='cytoband.gtf',
                        help='''File with cytoband information. Required for
                                naming loci according to chromosome band.
                                Default is "cytoband.gtf" (in working
                                directory).''')
    group2 = parser.add_argument_group('Output')
    group2.add_argument('--outdir', default='',
                        help='''Parent output directory. A new directory with
                                the family name will be created in this
                                directory. Default is working directory.''')
    group2.add_argument('--track_dir', default='l1base_tracks',
                        help='''Directory to store track files downloaded from
                                L1base.''')

    try:
        main(parser.parse_args())
    except KeyboardInterrupt:
        sys.exit('\nCancelled by user')

if __name__ == '__main__':
    console()