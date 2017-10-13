#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

from utils.omicutils import ChromosomeDict
# from utils.gtfutils import cluster_gtf, slop_gtf, intersect_gtf, conflict_gtf
from utils.gtfutils import sort_gtf, intersect_gtf
from utils.gtfutils import read_gtf_file, write_gtf_file, read_gtf_clusters

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def gtftools_sort(args):
    # Load chromosome sizes
    if args.chrom_sizes and os.path.exists(args.chrom_sizes):
        cdict = ChromosomeDict(args.chrom_sizes)
    else:
        cdict = ChromosomeDict()

    giter = sort_gtf(read_gtf_file(args.infile), cdict.reforder)
    write_gtf_file(giter, args.outfile)

def gtftools_sortclust(args):
    # Load chromosome sizes
    if args.chrom_sizes and os.path.exists(args.chrom_sizes):
        cdict = ChromosomeDict(args.chrom_sizes)
    else:
        cdict = ChromosomeDict()

    giter = sort_gtf(read_gtf_clusters(args.infile), cdict.reforder)
    write_gtf_file(giter, args.outfile)

def gtftools_intersect(args):
    gtfA = read_gtf_file(args.a)
    gtfBs = [read_gtf_file(bname) for bname in args.b]
    for gA, gBs in intersect_gtf(gtfA, gtfBs):
        if not gBs:
            if args.v or args.all:
                l = "{}\t{}\t{}".format(gA, ".", ("\t"*9))
                print >>sys.stdout, l
        else:
            for i,gB in gBs:
                if not args.v:
                    l = "{}\t{}\t{}".format(gA, i, gB)
                    print >> sys.stdout, l

def console():
    import argparse
    parser = argparse.ArgumentParser(
        description='''GTF tools'''
    )
    subparsers = parser.add_subparsers()

    ''' sort '''
    parser_sort = subparsers.add_parser('sort',
                                        help='Sort GTF file.')
    parser_sort.add_argument('--chrom_sizes', default='chrom.sizes',
                        help='''File with chromosome sizes and order. Default
                                is "chrom.sizes" (in working directory).''')
    parser_sort.add_argument('infile',
                             nargs='?', type=argparse.FileType('rU'),
                             default=sys.stdin,
                             help="Input GTF file. Default: stdin.")
    parser_sort.add_argument('outfile',
                             nargs='?', type=argparse.FileType('w'),
                             default=sys.stdout,
                             help="Output GTF file. Default: stdout.")
    parser_sort.set_defaults(func=gtftools_sort)

    ''' sortclust '''
    parser_sortclust = subparsers.add_parser('sortclust',
                                             help='Sort clustered GTF file.')
    parser_sortclust.add_argument('--chrom_sizes', default='chrom.sizes',
                        help='''File with chromosome sizes and order. Default
                                is "chrom.sizes" (in working directory).''')
    parser_sortclust.add_argument('infile',
                             nargs='?', type=argparse.FileType('rU'),
                             default=sys.stdin,
                             help="Input GTF file. Default: stdin.")
    parser_sortclust.add_argument('outfile',
                             nargs='?', type=argparse.FileType('w'),
                             default=sys.stdout,
                             help="Output GTF file. Default: stdout.")
    parser_sortclust.set_defaults(func=gtftools_sortclust)

    ''' intersect '''
    parser_intersect = subparsers.add_parser('intersect',
                                             help='Intersect GTF files.')
    parser_intersect.add_argument('-v', action='store_true',
                                  help='''Only report entries in A that have no
                                          overlap with B. (Default is to report
                                          only entries in A that have
                                          overlap).''')
    parser_intersect.add_argument('--all', action='store_true',
                                  help='''Report all entries in A whether or
                                          not there is overlap with B. (Default
                                          is to report only entries in A that
                                          have overlap).''')
    parser_intersect.add_argument('-a',
                                  help='''GTF A.''')
    parser_intersect.add_argument('-b', action='append',
                                  help='''GTF B. Can be specified multiple
                                          times''')
    parser_intersect.set_defaults(func=gtftools_intersect)

    try:
        args = parser.parse_args()
        args.func(args)
    except KeyboardInterrupt:
        sys.exit('\nCancelled by user')


if __name__ == '__main__':
    console()

