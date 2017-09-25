def console():
    import argparse
    parser = argparse.ArgumentParser(
        description='''Construct ERV annotations for an ERV family'''
    )
    parser.add_argument('--genome_build', default=_get_build_from_file(),
                        help='''Genome build to use, i.e. hg19, hg38.
                                Required. Can be configured by creating a text
                                file called "build.txt" with the build ID.''')
    parser.add_argument('--chrom_sizes', default='chrom.sizes',
                        help='File with chromosome sizes and order.')
    parser.add_argument('--cytoband', default='cytoband.gtf',
                        help='''File with cytoband information. Required for
                                naming loci according to chromosome band.''')
    parser.add_argument('--min_model_pct', type=float, default=0.1,
                        help='''Minimum percentage of internal model covered by
                                locus. Loci that do not meet this criteria are
                                rejected.''')
    parser.add_argument('--flank_size', type=int,
                        help='''Flanking distance to search for LTRs.''')

    parser.add_argument('--outdir', default='',
                        help='Output directory for family')

    parser.add_argument('--track_dir', default='ucsc_tracks',
                        help='Directory to store UCSC track output.')

    parser.add_argument('--save_intermediate', action='store_true',
                        help='Write intermediate GTF files.')
    parser.add_argument('--noisy', action='store_true',
                        help='Noisy output')

    parser.add_argument('--auto', action='store_true',
                        help='Attempt to automatically resolve conflicts.')
    parser.add_argument('--review', action='store_true',
                        help='Review conflicts after resolving.')
    parser.add_argument('--no_igv', action='store_true',
                        help='Do not use IGV. Also disables snapshots.')
    parser.add_argument('--no_snapshot', action='store_true',
                        help='''Do not take snapshots of conflicting loci.
                                Default: Snapshots are taken of conflicting
                                loci after conflicts are resolved, showing the
                                initial and final annotations.''')
    parser.add_argument('--snapshot_final', action='store_true',
                        help='''Take snapshots of all final loci. Requires IGV.
                                Default: Snapshots are only taken for final
                                conflicting loci (not all loci).''')

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