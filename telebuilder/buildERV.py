#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import string
from glob import glob
from collections import defaultdict, Counter
from itertools import chain
import atexit

from .utils import rmskutils
from .utils import ervutils
from .utils import _get_build_from_file

from .utils.utils import tsv, collapse_list, raw_input_stderr
from .utils.gtfutils import cluster_gtf, sort_gtf, slop_gtf, intersect_gtf, conflict_gtf, region_gtf
from .utils.gtfutils import read_gtf_file, write_gtf_file
from .utils.omicutils import ChromosomeDict
from .utils.igvutils import igv_init

from . import __version__

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2023 Matthew L. Bendall"


def stage1(models, trackdir, gbuild, overwrite=False, logh=sys.stderr):
    ''' Stage 1: Download RMSK tracks from UCSC
    '''
    print('*** Stage 1: Downloading RepeatMasker tracks from UCSC', file=logh)    
    track_files = []
    for m in models:
        tf = os.path.join(trackdir, '%s.%s.txt' % (gbuild, m))
        if not os.path.exists(tf) or overwrite:
            print('\tDownloading %s RMSK track for %s' % (gbuild, m), file=sys.stderr)
            with open(tf, 'w') as outh:
                qry = rmskutils.RMSK_MODEL_QUERY % m
                response = rmskutils.ucsc_download(gbuild, qry, outh=outh)
        else:
            print('\tFound track for %s (%s)' % (m, tf), file=sys.stderr)        
        track_files.append((m, tf))
    
    return track_files

def stage2(track_files, int_model, logh=sys.stderr):
    ''' Stage 2: Convert RMSK tracks to GTF
    '''
    print('*** Stage 2: Converting RepeatMasker tracks to GTF', file=logh)    
    gtfs = defaultdict(list)
    mlens = {}
    for m,f in track_files:
        lines = tsv(f)
        try:
            header = next(lines)
            assert all(h==n for h,(n,t) in zip(header, rmskutils.RMSKLine.COLS))
            rmsk = [rmskutils.RMSKLine(l) for l in lines]
            mlens.update(rmskutils.guess_rmsk_model_lengths(rmsk))
            for rl in rmsk:
                rl.fix_model_coords(mlens[rl.repName])
                g = rl.to_gtf()
                g.attr['geneRegion'] = 'internal' if m in int_model else 'ltr'
                gtfs[m].append(g)
        except StopIteration:
            print('[WARNING] %s has no records (%s is empty).' % (m, f), file=logh)

    print("\tInternal model lengths:", file=logh)
    for m in int_model:
        print('\t\t%s%d' % (m.ljust(20), mlens[m]), file=logh)

    print("\tLTR model lengths:", file=logh)
    lmods = sorted(k for k in list(mlens.keys()) if k not in int_model)
    for m in lmods:
        print('\t\t%s%d' % (m.ljust(20), mlens[m]), file=logh)
    
    return gtfs, mlens

def stage3(igtfs, cdict, shortdist, longdist, logh=sys.stderr):
    ''' Stage 3: Merge internal annotations '''
    print('*** Stage 3: Merging internal annotations', file=logh)        
    print('\tFound %d internal annotations' % len(igtfs), file=logh)        
    # Merge the internal annotations that are very close (<10 bp)
    print('\tMerging annotations < %d bp apart' % shortdist, file=logh)
    iclusters = cluster_gtf(igtfs, dist=shortdist)
    print('\t%d merged annotations after first merge' % len(iclusters), file=logh)
    
    # Merge the internal annotations that are fairly close (<10kb) and have consecutive models
    def consecutive_rmsk_model(a, b):
        left, right = (a, b) if a.start < b.start else (b, a)
        if left.strand != right.strand: return False
        if left.strand == '+':
            # Determine if the right cluster is a continuation of the left cluster
            return left.members[-1].attr['repLeft'] < right.members[0].attr['repLeft']
        else:
            # Determine if the left cluster is a continuation of the right cluster
            return right.members[0].attr['repStart'] < left.members[-1].attr['repStart']

    print('\tMerging annotations < %d bp apart that have consecutive models' % (longdist), file=logh)    
    iclusters = cluster_gtf(iclusters, dist=longdist, criteria=consecutive_rmsk_model)
    print('\t%d merged annotations after second merge' % len(iclusters), file=logh)
    
    # Sort clusters and add attributes for locus and internal model
    iclusters = sort_gtf(iclusters, cdict.reforder)
    for i, g in enumerate(iclusters):
        g.set_attr_from_members('repName', 'intModel')
        g.set_attr('locus', '%s_%04d' % (g.attr['intModel'], i+1))
    
    return iclusters

def stage4(iclusters, ltr_gtfs, flanksize, cdict, logh=sys.stderr):
    ''' Stage 4: Find LTRs flanking internal clusters '''
    print('*** Stage 4: Finding LTRs that flank internal regions', file=logh)    
    print('\tUsing flanksize = %d' % flanksize, file=logh)
    print('\tFound %d LTR annotations' % sum(map(len, ltr_gtfs)), file=logh)
    
    # Create annotations with flanking regions (slop)
    slop = slop_gtf(iclusters, flanksize, cdict.reflen)
    # Find LTR annotations that overlap with the "slop" annotations
    isect = intersect_gtf(slop, ltr_gtfs)
    
    # Create annotations combining internal and LTR
    mclusters = []
    iclust_d = {g.attr['locus']:g for g in iclusters} # Quickly retrieve cluster by locus ID
    for islop, flanks in isect:
        # Retrieve icluster with locus matching islop
        g = iclust_d.pop(islop.attr['locus'])
        for lidx, h in flanks:
            g.add(h.copy())
        mclusters.append(g)
    
    print('\t%d merged clusters' % len(mclusters), file=logh)
    print('\t%d unmerged internal clusters' % len(list(iclust_d.values())), file=logh)
    
    # return sort_gtf(mclusters + iclust_d.values(), cdict.reforder)
    return sort_gtf(mclusters, cdict.reforder)


def stage5(mclusters, mlens, logh=sys.stderr):
    ''' Add cluster attributes '''
    print('*** Stage 5: Adding cluster attributes', file=logh)
    def calculate_internal_coverage(_clust):
        ''' Calculate number of internal model bases covered in locus
            This is the number of "query" bases represented, not reference bases.
        '''
        if _clust.strand == '+':
            return sum(a.attr['repEnd'] - a.attr['repStart'] for a in _clust.members if a.attr['geneRegion'] == 'internal')
        else:
            return sum(a.attr['repEnd'] - a.attr['repLeft'] for a in _clust.members if a.attr['geneRegion'] == 'internal')
    
    # Include category and model coverage attributes
    for clust in mclusters:
        # Determine category
        geneord = collapse_list([g.attr['geneRegion'] for g in clust.members])
        if geneord == ['ltr', 'internal', 'ltr', ]:
            clust.attr['category'] = "prototype"
        elif geneord == ['internal',]:
            clust.attr['category'] = "internal"
        elif geneord == ['ltr', 'internal',] or geneord == ['internal', 'ltr',]:
            clust.attr['category'] = "oneside"
        else:
            clust.attr['category'] = 'unknown'
        # Add model coverage and percent
        model_cov = calculate_internal_coverage(clust)
        model_pct = min(100, (float(model_cov) / mlens[clust.attr['intModel']]) * 100)
        clust.attr['model_cov'] =  model_cov
        clust.attr['model_pct'] =  '%.1f' % model_pct
        # Add "locus" attribute to all members
        assert 'locus' in clust.attr
        clust.set_attr('locus', clust.attr['locus'])
        
    return mclusters


def stage6(mclusters, min_model_pct, logh=sys.stderr):
    print('*** Stage 6: Filtering short clusters', file=logh) 
    rej = [c for c in mclusters if float(c.attr['model_pct']) < (min_model_pct * 100)]
    mc = [c for c in mclusters if float(c.attr['model_pct']) >= (min_model_pct * 100)]
    print('\t%d rejected clusters' % len(rej), file=logh)
    print('\t%d clusters' % len(mc), file=logh)
    return mc, rej

def stage7(mclusters, auto, igv, dest, logh=sys.stderr):
    print('*** Stage 7: Resolving conflicts', file=logh)
    conflicts = conflict_gtf(mclusters)
    if len(conflicts) == 0:
        print('\tNo conflicts found', file=sys.stderr)
        lcons = []
    else:
        print('\t%d conflict(s) to resolve' % len(conflicts), file=sys.stderr)
        lcons = ervutils.inspect_conflicts(conflicts, auto, igv, dest)
        
        for lcon in lcons:
            rem, app = lcon.resolve()
            for clust in rem:
                mclusters.remove(clust)

            for clust in app:
                clust.attr['category'] = '%s*' % lcon.reason
                mclusters.append(clust)
    
    return mclusters, lcons

# Alphabet suffixes
SUFFIXES = list(string.ascii_letters[:26])
# In case there are more than 26:
SUFFIXES += [a+b for b in SUFFIXES for a in SUFFIXES]

def stage8(mclusters, cytogtf, fam, logh=sys.stderr):
    print('*** Stage 8: Naming loci', file=logh)
    if cytogtf is None:
        for g in mclusters:
            g.set_attr('transcript_id', g.attr['locus'])
            g.set_attr('gene_id', g.attr['locus'])
    else:
        isect = intersect_gtf(mclusters, [cytogtf, ], stranded=False)
        byband = defaultdict(list)
        for g, bands in isect:
            chrom = g.chrom[3:] if g.chrom[:3] == 'chr' else g.chrom
            if bands:
                if len(bands) > 1:
                    _bands = sort_gtf([_[1] for _ in bands])
                    band = _bands[0].attr['gene_id']
                else:
                    band = bands[-1][1].attr['gene_id']
            else:
                band = ''
            byband[(chrom, band)].append(g)

        for (chrom, band), gtfs in byband.items():
            if len(gtfs) == 1:
                suffixes = ['']
            else:
                suffixes = SUFFIXES[:len(gtfs)]
                # Add parentheses if band is not present
                if band == '':
                    suffixes = ['(%s)' % _ for _ in suffixes]

            # Set the locus name, transcript_id, and gene_id
            for g, suf in zip(gtfs, suffixes):
                g.set_attr('locid', g.attr['locus'])
                name = '%s_%s%s%s' % (fam, chrom, band, suf)
                g.set_attr('locus', name)
                g.set_attr('transcript_id', name)
                g.set_attr('gene_id', name)
    return mclusters


def main(args):
    ''' Setup '''
    logh = sys.stderr if args.noisy else open(os.devnull, 'w')

    if args.genome_build is None:
        sys.exit("ERROR: --genome_build is required.")
    print("[VERBOSE] Genome build: %s" % args.genome_build, file=logh)

    # Create track directory, if necessary
    if not os.path.isdir(args.track_dir):
        print("[VERBOSE] Creating track directory: %s" % args.track_dir, file=logh)
        os.makedirs(args.track_dir)
    else:
        print("[VERBOSE] Using track directory: %s" % args.track_dir, file=logh)

    # Create output directory, if necessary
    dest = os.path.join(args.outdir, args.fam)
    if not os.path.isdir(dest):
        print("[VERBOSE] Creating output directory: %s" % dest, file=logh)
        os.makedirs(dest)
    else:
        print("[VERBOSE] Using output directory: %s" % dest, file=logh)

    
    # Setup IGV and snapshot directory
    igv = None if args.no_igv else igv_init(args.genome_build)
    if igv is None:
        if not args.no_igv:
            print("[WARNING] Could not connect to IGV.", file=sys.stderr)
        else:
            print("[VERBOSE] Not using IGV.", file=logh)
        snapshot_dir = None
        snapshot_final = False
    else:
        print("[VERBOSE] Using IGV.", file=logh)
        atexit.register(igv.close)
        if args.compare_gtfs and os.path.isdir(args.compare_gtfs):
            other_gtfs = glob(os.path.join(args.compare_gtfs, '*.gtf'))
            other_gtfs += glob(os.path.join(args.compare_gtfs, '*.gtf.gz'))
            for og in other_gtfs:
                igv.load(os.path.abspath(og))

        if args.no_snapshot is False:
            snapshot_dir = os.path.join(dest, 'snapshots')
            if not os.path.isdir(snapshot_dir):
                print("[VERBOSE] Creating snapshot directory: %s" % snapshot_dir, file=logh)
                os.makedirs(snapshot_dir)
            igv.snapshotDirectory(snapshot_dir)
            print("[VERBOSE] Snapshots will be saved to %s" % snapshot_dir, file=logh)
            snapshot_final = args.snapshot_final
        else:
            snapshot_dir = None
            snapshot_final = False

    if snapshot_final:
        print("[VERBOSE] Taking snapshots of final loci.", file=logh)

    # Keep intermediate files?
    save_intermediate = args.save_intermediate
    if save_intermediate:
        print("[VERBOSE] Saving intermediate GTFs", file=logh)

    # Load chromosome sizes
    if os.path.exists(args.chrom_sizes):
        cdict = ChromosomeDict(args.chrom_sizes)
    else:
        print('[WARNING] Chromosome sizes were not found', file=sys.stderr)
        cdict = ChromosomeDict()


    int_model = args.intmodel.split(',')
    ltr_model = args.ltrmodel.split(',')

    ''' Stage 1: Download RMSK tracks from UCSC '''
    track_files = stage1(int_model + ltr_model, args.track_dir, args.genome_build)

    ''' Stage 2: Convert RMSK tracks to GTF '''
    gtfs, mlens = stage2(track_files, int_model)
    ltr_model = [m for m in ltr_model if m in gtfs] # Revise LTR model

    if save_intermediate:
        for m, gtf in gtfs.items():
            with open(os.path.join(dest, '%s.gtf' % m), 'w') as outh:
                write_gtf_file(sort_gtf(gtf, cdict.reforder), outh)

    ''' Stage 3: Merge internal annotations '''
    igtfs = list(chain.from_iterable(gtfs[im] for im in int_model))
    iclusters = stage3(igtfs, cdict,
                       shortdist=args.short_dist,
                       longdist=args.long_dist)

    # Remove records that are not in the chromosome list
    if cdict.reforder is not None:
        iclusters = region_gtf(iclusters, cdict.reforder)

    if save_intermediate:
        with open(os.path.join(dest, 'internal.gtf'), 'w') as outh:
            write_gtf_file(sort_gtf(iclusters, cdict.reforder), outh)

    ''' Stage 4: Find flanking LTRs '''
    if args.flank_size:
        flanksize = args.flank_size
    else:
        flanksize = max([mlens[lm] for lm in ltr_model])
        flanksize = int(round(flanksize / 100.) * 100) # Make it a round 100

    mclusters = stage4(iclusters, [gtfs[lm] for lm in ltr_model], flanksize, cdict)

    ''' Stage 5: Add cluster attributes'''
    mclusters = stage5(mclusters, mlens)

    ''' Stage 6: Filter short '''
    mclusters, rejected = stage6(mclusters, args.min_model_pct)
    if save_intermediate:
        with open(os.path.join(dest, 'rejected.gtf'), 'w') as outh:
            write_gtf_file(sort_gtf(rejected, cdict.reforder), outh)
        with open(os.path.join(dest, 'merged.gtf'), 'w') as outh:
            write_gtf_file(sort_gtf(mclusters, cdict.reforder), outh)
    
    ''' Stage 7: Resolve conflicts'''
    mclusters, lcons = stage7(mclusters, args.auto, igv, dest)
    
    ''' Stage 8: Name loci '''
    # Load cytoband
    if os.path.exists(args.cytoband):
        cytoband = read_gtf_file(args.cytoband)
    else:
        cytoband = None
        print('[WARNING] Cytoband was not found', file=sys.stderr)

    mclusters = stage8(mclusters, cytoband, args.fam)
    
    ''' Stage 9: Output '''
    print('*** Stage 9. Final Output', file=sys.stderr)
    final_gtf_file = os.path.join(dest, '%s.gtf' % args.fam)
    final_gtf = sort_gtf(mclusters, cdict.reforder)
    with open(final_gtf_file, 'w') as outh:
        write_gtf_file(final_gtf, outh)

    if igv:
        igv.load(os.path.abspath(final_gtf_file))
        igv.collapse().expand(f'{args.fam}.gtf')

    # Review conflicts and display or snapshot
    if len(lcons) == 0:
        print("[VERBOSE] No conflicts.", file=logh)
    else:
        for i, lcon in enumerate(lcons):
            print('[REVIEW] Conflict %02d.' % (i+1), file=sys.stderr)
            print('\t%s' % lcon.display_review_text(), file=sys.stderr)
            if snapshot_dir is not None or args.review:
                igv.goto(lcon.region_str())
            if snapshot_dir is not None:
                igv.snapshot('conflict.%02d.%s.png' % ((i+1), lcon.action))
            if args.review:
                z = raw_input_stderr('\tPress [ENTER] to continue. ')

    # Create final snapshots
    if snapshot_final and snapshot_dir is not None:
        for g in final_gtf:
            igv.goto(g.attr['transcript_id'])
            igv.snapshot('%s.png' % g.attr['transcript_id'])

    ''' Summary '''
    print('*** Stage 10. Summary', file=sys.stderr)
    categories = Counter()
    ltr_usage = defaultdict(Counter)
    for g in final_gtf:
        categories[g.attr['category']] += 1
        if g.attr['category'] == 'prototype':
            rn = collapse_list([h.attr['repName'] for h in g.members])
            if rn[0] == rn[-1]:
                ltr_usage['prototype'][rn[0]] += 1
            else:
                lr = '%s/%s' % tuple(sorted([rn[0],rn[-1]]))
                ltr_usage['prototype'][lr] += 1
        if g.attr['category'] == 'oneside':
            rn = collapse_list([h.attr['repName'] for h in g.members])
            if rn[-1] in ltr_model:
                ltr_usage['oneside'][rn[-1]] += 1
            else:
                ltr_usage['oneside'][rn[0]] += 1

    print('\n\n', file=sys.stderr)
    print('%s %s summary %s' % ('*'*20, args.fam, '*'*20), file=sys.stderr)
    print('Locus types:', file=sys.stderr)
    for cat in ['prototype', 'oneside', 'internal',]:
        print('\t%s%d' % (cat.ljust(20), categories[cat]), file=sys.stderr)
    for cat,v in categories.most_common():
        if cat not in ['prototype', 'oneside', 'internal',]:
            print('\t%s%d' % (cat.ljust(20), categories[cat]), file=sys.stderr)
    print('LTR usage (prototype):', file=sys.stderr)
    for k,v in ltr_usage['prototype'].most_common():
        print('\t%s%d' % (k.ljust(20), v), file=sys.stderr)

    print('LTR usage (oneside):', file=sys.stderr)
    for k, v in ltr_usage['oneside'].most_common():
        print('\t%s%d' % (k.ljust(20), v), file=sys.stderr)



def console():
    import argparse
    parser = argparse.ArgumentParser(
        description='''Construct ERV annotations for an ERV family'''
    )
    parser.add_argument(
        '-V', '--version',
        action='version',
        help='Show the version number and exit.',
        version=f"telebuilder {__version__}",
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
    group1.add_argument('--compare_gtfs', default='compare',
                        help='''Directory containing GTFs to load for visual
                                comparison. Default is "compare" (in working
                                directory).''')

    group2 = parser.add_argument_group('Output')
    group2.add_argument('--outdir', default='',
                        help='''Parent output directory. A new directory with
                                the family name will be created in this
                                directory. Default is working directory.''')
    group2.add_argument('--track_dir', default='ucsc_tracks',
                        help='''Directory to store track files downloaded from
                                UCSC.''')
    group2.add_argument('--save_intermediate', action='store_true',
                        help='Write intermediate GTF files.')
    group2.add_argument('--snapshot_final', action='store_true',
                        help='''Take snapshots of all final loci. Requires IGV.
                                Default: Snapshots are only taken for final
                                conflicting loci (not all loci).''')
    group2.add_argument('--noisy', action='store_true',
                        help='Noisy output')

    group4 = parser.add_argument_group('Conflicts')
    group4.add_argument('--auto', action='store_true',
                        help='Attempt to automatically resolve conflicts.')
    group4.add_argument('--no_igv', action='store_true',
                        help='Do not use IGV. Also disables snapshots.')
    group4.add_argument('--no_snapshot', action='store_true',
                        help='''Do not take snapshots of conflicting loci.
                                Default: Snapshots are taken of conflicting
                                loci after conflicts are resolved, showing the
                                initial and final annotations.''')
    group4.add_argument('--review', action='store_true',
                        help='Review conflicts after resolving.')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--min_model_pct', type=float, default=0.1,
                        help='''Minimum percentage of internal model covered by
                                locus. Loci that do not meet this criteria are
                                rejected. Default: 0.1.''')
    group3.add_argument('--flank_size', type=int,
                        help='''Flanking distance to search for LTRs.
                                Default: length of longest LTR model, rounded
                                to nearest 100.''')
    group3.add_argument('--short_dist', type=int, default=10,
                        help='''Automatically merge internal annotations less
                                than short_dist apart. In this case, we assume
                                that the two annotations are part of the same
                                locus but were created as two annotations
                                because of failed extension. Default: 10.''')
    group3.add_argument('--long_dist', type=int, default=2500,
                        help='''Merge internal annotations less than long_dist
                                apart if the model alignment positions agree.
                                In this case, we assume that a single locus was
                                split due to an insertion, often by another
                                retroelement. Default: 2500.''')

    parser.add_argument('fam',
                        help="Name for family")
    parser.add_argument('intmodel',
                        help="Name(s) for internal models")
    parser.add_argument('ltrmodel',
                        help="Name(s) for LTR models, comma separated.")

    try:
        main(parser.parse_args())
    except KeyboardInterrupt:
        sys.exit('\nCancelled by user')


if __name__ == '__main__':
    console()

