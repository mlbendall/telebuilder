#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re

from Bio import SeqIO
from Bio.Seq import Seq

from .utils.omicutils import ChromosomeDict
from .utils.omicutils import get_sequence_ucsc, get_sequence_fasta
# from utils.gtfutils import cluster_gtf, slop_gtf, intersect_gtf, conflict_gtf
from .utils import _get_build_from_file
from .utils.utils import wraplines
from .utils.gtfutils import sort_gtf, intersect_gtf
from .utils.gtfutils import read_gtf_file, write_gtf_file, read_gtf_clusters
from .utils.gtfclasses import GTFLine, GTFCluster

from . import __version__

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
                print(l, file=sys.stdout)
        else:
            for i,gB in gBs:
                if not args.v:
                    l = "{}\t{}\t{}".format(gA, i, gB)
                    print(l, file=sys.stdout)

def gtftools_tsv(args):
    gtfA = [g for g in read_gtf_file(args.infile) if g.feature == args.feat]
    all_attrs = set()
    for gA in gtfA:
        all_attrs |= set(gA.attr.keys())

    cols = [args.key, 'chrom', 'start', 'end', 'strand', 'score', ]
    all_attrs = [a for a in sorted(all_attrs) if a != args.key]
    print('\t'.join(cols + all_attrs), file=sys.stdout)
    for i, gA in enumerate(gtfA):
        r = [
            gA.attr[args.key] if args.key in gA.attr else 'LOC%06d' % (i+1),
            gA.chrom,
            str(gA.start),
            str(gA.end),
            gA.strand,
            gA.score if gA.score != -1 else '',
        ]
        for a in all_attrs:
            r.append(gA.attr[a] if a in gA.attr else '')
        print('\t'.join(map(str, r)), file=sys.stdout)

from Bio.Seq import Seq
import re

def gtftools_extract(args):
    if args.gtfout:
        outgtf = open(args.gtfout, 'w')

    if args.genome:
        genome_dict = {s.id:s for s in SeqIO.parse(args.genome, 'fasta')}
        print('Loading genome complete.', file=sys.stderr)

    # Create set of selected locus search strings
    selected_loci = set()
    if args.locus:
        selected_loci |= set([r'^%s$' % _ for _ in args.locus])
    if args.locus_regex:
        selected_loci |= set([r'%s' % _ for _ in args.locus_regex])
    if args.locus_file:
        selected_loci |= set([r'^%s$' % l.strip('\n') for l in args.locus_file])

    selected_regex = '|'.join('(?:{0})'.format(x) for x in selected_loci)
    #print selected_regex
    selected_regex = re.compile(selected_regex)

    clusters = read_gtf_clusters(args.infile)
    for gc in clusters:
        if not selected_regex.search(gc.attr['locus']):
            continue

        reg = (gc.chrom, gc.start, gc.end)

        # Adjust annotations to locus coordinates
        adjusted = GTFCluster()
        adjusted.attr = dict(gc.attr)
        adjusted.attr['reg'] = '%s:%d-%d(%s)' % (gc.chrom, gc.start, gc.end, gc.strand)
        for g in gc.members:
            h = g.copy()
            h.attr['reg'] = '%s:%d-%d(%s)' % (g.chrom, g.start, g.end, g.strand)
            h.chrom = gc.attr['locus']
            h.start = g.start - gc.start + 1
            h.end = g.end - gc.start + 1
            if h.attr['geneRegion'] == 'ltr':
                h.feature = 'LTR'
                # h.attr = {'name': h.attr['repName']}
                del h.attr['transcript_id']
                del h.attr['gene_id']
                h.attr['name'] = h.attr['repName']
            adjusted.add(h)

        if args.genome:
            rawseq = get_sequence_fasta(genome_dict, reg)
            if rawseq is None:
                continue
        else:
            rawseq = get_sequence_ucsc(args.genome_build, reg)

        if gc.strand == '-':
            # Reverse complement the sequence
            rawseq = str(Seq(rawseq).reverse_complement())
            # Reverse the annotations
            rev = GTFCluster()
            rev.attr = dict(adjusted.attr)
            for g in adjusted.members:
                h = g.copy()
                h.start = adjusted.end - g.end + 1
                h.end = adjusted.end - g.start + 1
                h.strand = '+'
                rev.add(h)
            adjusted = rev

        # Make regions
        allexons = rawseq.lower()
        # selexons = rawseq.lower()
        for n in adjusted.members:
            allexons = allexons[0:(n.start-1)] + allexons[(n.start-1):n.end].upper() + allexons[n.end:]
            # if 'all' in upregions or n.attr['geneRegion'] in upregions:
            #     selexons = selexons[0:(n.start-1)] + selexons[(n.start-1):n.end].upper() + selexons[n.end:]

        assert len(allexons) == len(rawseq)
        assert allexons.lower() == rawseq.lower(), '%s\n%s' % (allexons.lower(), rawseq)

        # Create annotations for introns
        for m in re.finditer('[a-z]+', allexons):
            ig = GTFLine([
                gc.attr['locus'], gc.members[0].source, 'intron',
                m.start() + 1, m.end(),
                '.', '+', '.', {'geneRegion':'nohit'}
            ])
            # ig.chrom = gc.attr['locus']
            # ig.start, ig.end = (m.start()+1, m.end())
            adjusted.add(ig)

        # Create internal annotation
        internals = [g for g in adjusted.members if g.attr['geneRegion'] == 'internal']
        misc_g = GTFLine([
            gc.attr['locus'], gc.members[0].source, 'misc_feature',
            min(g.start for g in internals), max(g.end for g in internals),
            '.', '+', '.', {'name':'internal',
                            'geneRegion': 'internal'},
        ])
        adjusted.add(misc_g)

        # Print
        print('>%s' % adjusted.attr['locus'], file=args.outfile)
        print(wraplines(allexons), file=args.outfile)
        if args.gtfout:
            print(adjusted, file=outgtf)
    #
    if args.gtfout:
        outgtf.close()
    return


def console():
    import argparse
    parser = argparse.ArgumentParser(
        description='''GTF tools'''
    )
    parser.add_argument(
        '-V', '--version',
        action='version',
        help='Show the version number and exit.',
        version=f"telebuilder {__version__}",
    )
    subparsers = parser.add_subparsers()

    ''' sort '''
    parser_sort = subparsers.add_parser('sort',
                                        help='Sort GTF file.')
    parser_sort.add_argument('--chrom_sizes', default='chrom.sizes',
                        help='''File with chromosome sizes and order. Default
                                is "chrom.sizes" (in working directory).''')
    parser_sort.add_argument('infile',
                             nargs='?', type=argparse.FileType('r'),
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
                             nargs='?', type=argparse.FileType('r'),
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

    ''' tsv '''
    parser_tsv = subparsers.add_parser('tsv',
                                        help='GTF to TSV.')
    parser_tsv.add_argument('--feat', default='gene',
                            help='''Only include records with feature. Default
                                    is "gene".''')
    parser_tsv.add_argument('--key', default='locus',
                            help='''Key used for record ID. Default is
                                    "locus".''')
    parser_tsv.add_argument('infile',
                             nargs='?', type=argparse.FileType('r'),
                             default=sys.stdin,
                             help="Input GTF file. Default: stdin.")
    parser_tsv.add_argument('outfile',
                             nargs='?', type=argparse.FileType('w'),
                             default=sys.stdout,
                             help="Output GTF file. Default: stdout.")
    parser_tsv.set_defaults(func=gtftools_tsv)

    ''' extract '''
    parser_extract = subparsers.add_parser('extract',
                                             help='Extract sequences.')
    parser_extract.add_argument('--genome_build',
                                default=_get_build_from_file(),
                                help='''Genome build to use, i.e. hg19, hg38.
                                        Can be configured by creating a text
                                        file called "build.txt" with the
                                        build ID.''')
    # parser_extract.add_argument('--nohit',
    #                             default='lowercase',
    #                             choices=['lowercase', 'gap', 'exclude', 'none'],
    #                             help='''How to report intervals that are not
    #                                     covered by the model.'''
    #                             )
    parser_extract.add_argument('--genome',
                                help='''FASTA file containing reference
                                        sequence.''')
    # parser_extract.add_argument('--regions',
    #                             default='all',
    #                             help='''Gene region(s) to include as uppercase,
    #                                     as a comma-separated list. Default is
    #                                     to have all exon annotations in
    #                                     uppercase ("all")''')
    parser_extract.add_argument('--locus',
                                action='append',
                                help='''Name of locus to be extracted. Multiple
                                        loci can be specified by providing this
                                        argument multiple times.''')
    parser_extract.add_argument('--locus_regex',
                                action='append',
                                help='''Regular expression matching locus to be
                                        extracted. Multiple loci can be
                                        specified by providing this argument
                                        multiple times.''')
    parser_extract.add_argument('--locus_file',
                                type=argparse.FileType('r'),
                                help='''File containing locus names to be
                                        extracted, one per line.''')
    parser_extract.add_argument('--gtfout',
                                help='Output GTF file (optional).')
    parser_extract.add_argument('infile',
                                nargs='?', type=argparse.FileType('r'),
                                default=sys.stdin,
                                help="Input GTF file. Default: stdin.")
    parser_extract.add_argument('outfile',
                                nargs='?', type=argparse.FileType('w'),
                                default=sys.stdout,
                                help="Output FASTA file. Default: stdout.")
    parser_extract.set_defaults(func=gtftools_extract)

    try:
        args = parser.parse_args()
        args.func(args)
    except KeyboardInterrupt:
        sys.exit('\nCancelled by user')


if __name__ == '__main__':
    console()

