#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

from utils.omicutils import ChromosomeDict
# from utils.gtfutils import cluster_gtf, slop_gtf, intersect_gtf, conflict_gtf
from utils.gtfutils import sort_gtf
from utils.gtfutils import read_gtf_file, write_gtf_file

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


def console():
    import argparse
    parser = argparse.ArgumentParser(
        description='''GTF tools'''
    )
    subparsers = parser.add_subparsers()
    parser_sort = subparsers.add_parser('sort', help='Sort GTF file.')

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

    try:
        args = parser.parse_args()
        args.func(args)
    except KeyboardInterrupt:
        sys.exit('\nCancelled by user')


if __name__ == '__main__':
    console()

