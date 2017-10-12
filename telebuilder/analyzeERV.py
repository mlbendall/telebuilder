#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip

from utils.utils import tsv
from utils.omicutils import ChromosomeDict
from utils import rmskutils
from utils import gtfutils
from collections import defaultdict, Counter
from utils.gtfclasses import GTFCluster

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def flanking_ltr_distribution(args, ltr_files):
    # Load chromosome sizes
    if os.path.exists(args.chrom_sizes):
        cdict = ChromosomeDict(args.chrom_sizes)
    else:
        print >>sys.stderr, '[WARNING] Chromosome sizes were not found'
        cdict = ChromosomeDict()

    # Read and filter by chromosome
    clusters = gtfutils.read_gtf_clusters(args.internal_gtf)
    if args.filter_chrom:
        chroms = set(args.filter_chrom.split(','))
        clusters = [c for c in clusters if c.chrom in chroms]

    iclusters = []
    for clust in clusters:
        newclust = GTFCluster([m for m in clust.members if m.attr['geneRegion'] == 'internal'])
        newclust.attr = clust.attr
        iclusters.append(newclust)

    locids = [c.attr['locus'] for c in iclusters]
    print >> sys.stderr, "[VERBOSE] Loaded %d internal clusters" % len(locids)

    # Create flanking annotations
    flanks = []
    slop = gtfutils.slop_gtf(iclusters, 500, cdict.reflen)
    for s, c in zip(slop, iclusters):
        left, right = gtfutils.subtract_gtflines(s, c.span())
        left.attr['side'] = "left"
        right.attr['side'] = "right"
        if left.length():
            flanks.append(left)
        else:
            print left
        if right.length():
            flanks.append(right)
        else:
            print right

    gtfBs = []
    for f in ltr_files:
        assert os.path.exists(f)
        gtfBs.append(list(gtfutils.read_gtf_file(gzip.open(f, 'rb'))))
    print >> sys.stderr, "[VERBOSE] Loaded %d LTR records from %d files." % (sum(len(_) for _ in gtfBs), len(gtfBs))

    # Intersect
    isect = gtfutils.intersect_gtf(flanks, gtfBs)
    # Summarize hits
    hits = defaultdict(lambda: {'left': set(), 'right': set()})
    for k, vlist in isect.iteritems():
        for lidx, h in vlist:
            hits[k.attr['locus']][k.attr['side']].add(h.attr['repName'])

    # Hits by LTR model
    totalhits = Counter()
    onesided = Counter()
    twosided = Counter()
    for locid in locids:
        if locid not in hits:
            continue  # nomatch += 1
        else:
            hasone = []
            hastwo = []
            hit = hits[locid]
            for h in list(hit['right']) + list(hit['left']):
                totalhits[h] += 1
            uhits = list(hit['right'] | hit['left'])
            for u in uhits:
                if u in hit['right'] and u in hit['left']:
                    twosided[u] += 1
                elif u in hit['right'] or u in hit['left']:
                    onesided[u] += 1

    MAXPRINT = 20

    print >> sys.stderr, 'Total number of flanking hits by model:'
    print >> sys.stderr, '\n'.join(
        ['\t%s%d' % (k.ljust(16), v) for k, v in totalhits.most_common() if
         v > 1][:MAXPRINT])

    print >> sys.stderr, 'Number of loci with hits on both sides:'
    print >> sys.stderr, '\n'.join(
        ['\t%s%d' % (k.ljust(16), v) for k, v in twosided.most_common() if
         v > 1][:MAXPRINT])

    print >> sys.stderr, 'Number of loci with hit on one side, but not both:'
    print >> sys.stderr, '\n'.join(
        ['\t%s%d' % (k.ljust(16), v) for k, v in onesided.most_common() if
         v > 1][:MAXPRINT])

    consider = [m for m, v in totalhits.most_common() if v > 1][:MAXPRINT]
    for i, m in enumerate(consider):
        print >> sys.stderr, 'Using LTRs:\n\t%s' % '+'.join(consider[:(i + 1)])
        chosen = set(consider[:(i + 1)])
        nlocboth = nlocone = nlocnone = 0
        for locid in locids:
            if locid not in hits:
                nlocnone += 1
            else:
                if len(hits[locid]['left'] & chosen) > 0 and len(
                                hits[locid]['right'] & chosen) > 0:
                    nlocboth += 1
                elif len(hits[locid]['left'] & chosen) > 0 or len(
                                hits[locid]['right'] & chosen) > 0:
                    nlocone += 1
                else:
                    nlocnone += 1

        print >> sys.stderr, '\t\t%d have both, %d have one, and %d have none' % (
        nlocboth, nlocone, nlocnone)
        pcts = [100 * float(_) / len(locids) for _ in
                (nlocboth, nlocone, nlocnone)]
        print >> sys.stderr, '\t\t(%.1f have both, %.1f have one, and %.1f have none)' % tuple(
            pcts)


def download_ltr(args):
    # Genome Build
    if args.genome_build is None:
        sys.exit("ERROR: --genome_build is required.")
    # print >>sys.stderr, "[VERBOSE] Genome build: %s" % args.genome_build

    # Load chromosome sizes
    if os.path.exists(args.chrom_sizes):
        cdict = ChromosomeDict(args.chrom_sizes)
    else:
        print >>sys.stderr, '[WARNING] Chromosome sizes were not found'
        cdict = ChromosomeDict()

    if args.filter_chrom:
        chroms = set(args.filter_chrom.split(','))
    else:
        chroms = cdict.reforder

    # Make output directory
    dest = os.path.abspath(args.ltr_dir)
    if not os.path.exists(dest):
        os.makedirs(dest)

    LTRQUERY = '''SELECT * FROM rmsk WHERE genoName = '%s' AND repClass = 'LTR' AND repName NOT REGEXP '-[Ii](nt)?$';'''
    rmskfiles = []
    for chrom in chroms:
        # print >> sys.stderr, "[VERBOSE] Chromosome: %s" % chrom
        rmskfile = os.path.join(dest, '%s.ltr.txt.gz' % chrom)
        if os.path.exists(rmskfile) and not args.overwrite:
            # print >> sys.stderr, "[VERBOSE] Found existing: %s" % rmskfile
            rmskfiles.append(rmskfile)
        else:
            response = rmskutils.ucsc_download(args.genome_build, LTRQUERY % chrom)
            if response:
                print >>sys.stderr, '[VERBOSE] Downloaded %d rows for %s' % (len(response.strip('\n').split('\n')), chrom)
                with gzip.open(rmskfile, 'wb') as outh:
                    print >>outh, response.strip('\n')
                rmskfiles.append(rmskfile)
            else:
                print >>sys.stderr, '[VERBOSE] No rows for %s' % chrom

    # Convert to GTF
    gtffiles = []
    for rf in rmskfiles:
        gf = '%s.gtf.gz' % '.'.join(rf.split('.')[:-2])
        if os.path.exists(gf) and not args.overwrite:
            # print >> sys.stderr, "[VERBOSE] Found existing: %s" % gf
            gtffiles.append(gf)
        else:
            rows = tsv(gzip.open(rf, 'rb'))
            header = rows.next()
            gtf = []
            for r in rows:
                g = rmskutils.RMSKLine(r).to_gtf()
                g.attr['geneRegion'] = 'ltr'
                gtf.append(g)
            gtf.sort(key=lambda x:x.start)
            with gzip.open(gf, 'wb') as outh:
                gtfutils.write_gtf_file(gtf, outh)
            print >> sys.stderr, "[VERBOSE] Wrote %d rows to %s" % (len(gtf), gf)
            gtffiles.append(gf)

    return gtffiles


def _get_build_from_file():
    if os.path.exists('build.txt'):
        return open('build.txt', 'rU').read().strip()



def console():
    import argparse
    parser = argparse.ArgumentParser(
        description='''Examine the distribution of various LTR models that flank internal
                       regions''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--genome_build', default=_get_build_from_file(),
                        help='''Genome build to use, i.e. hg19, hg38.
                                Required. Can be configured by creating a text
                                file called "build.txt" with the build ID.''')
    parser.add_argument('--chrom_sizes', default='chrom.sizes',
                        help='''File with chromosome sizes and order. Default
                                is "chrom.sizes" (in working directory).''')
    parser.add_argument('--ltr_dir', default='all_ltr',
                        help='Directory containing LTR annotations, one file per chromosome.')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite downloaded LTR annotations.')
    parser.add_argument('--filter_chrom',
                        help='''A comma-separated list of chromosomes to include. Default
                                is to examine all chromosomes.''')
    parser.add_argument('internal_gtf',
                        help="GTF file with internal annotations. Clustering is recommended.")

    args = parser.parse_args()
    ltrfiles = download_ltr(args)
    flanking_ltr_distribution(args, ltrfiles)


if __name__ == '__main__':
    console()
