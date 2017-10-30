# -*- coding: utf-8 -*-


from collections import defaultdict

from intervaltree import Interval, IntervalTree

from utils import tsv, overlap_length

from gtfclasses import GTFLine, GTFCluster

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


def subtract_gtflines(gA, gB):
    """ Subtract overlapping portions from GTF record

    Returns two copies of A with portion(s) overlapping B removed. Depending
    on the nature of the overlap, the operation can result in 0, 1, or 2 new
    records.

    Args:
        gA:
        gB:

    Returns:

    """
    r1 = gA.copy()
    r2 = gA.copy()
    r1.end = max(gA.start, min(gA.end, gB.start))
    r2.start = min(max(gA.start, gB.end), gA.end)
    return r1, r2


def sort_gtf(giter, chrom_order=None):
    """ Sort GTF

    Sort GTF records by chromosome and start position

    Args:
        giter (iterable): List of GTFLine or GTFCluster objects
        chrom_order (list): Chromosome sort order. Default is alphabetical.

    Returns:
        list: Sorted GTF records.

    """
    ret = sorted(giter, key=lambda x:x.start)
    if chrom_order is None:
        ret.sort(key=lambda x:x.chrom)
    else:
        c_d = {k:i for i,k in enumerate(chrom_order)}
        ret.sort(key=lambda x:c_d[x.chrom] if x.chrom in c_d else x.chrom)
    return ret


def region_gtf(giter, regions):
    """ Filter GTF records according to regions

    Args:
        giter:
        regions (list): List of regions to include. Regions can be specified
            as chrom[:start-end]. If start and end are not specified, the whole
            chromosome is used as the region.

    Yields:
        GTF records that a

    """
    _regions = defaultdict(list)
    for r in regions:
        if ':' in r:
            rchrom, t = r.split(':')
            rs, re = map(int, t.split('-'))
            assert rs <= re
            _regions[rchrom].append((rs,re))
        else: # Only chromosome was given
            _regions[r].append((float('-inf'), float('inf')))

    ret = []
    for g in giter:
        if g.chrom in _regions:
            for iv in _regions[g.chrom]:
                if overlap_length((g.start, g.end), iv) > 0:
                    ret.append(g)
                    break
    return ret


def slop_gtf(giter, b=0, chrom_lengths=None):
    """ Increase the size of each GTF record

    Args:
        giter (iterable): List of GTFLine or GTFCluster objects
        b (int): Number of base pairs to increase. Default is 0.
        chrom_lengths (dict): Chromosome name (str) -> chromosome length (int)

    Yields:
        GTFLine: New object with increased size.

    """
    chrom_lengths = {} if chrom_lengths is None else chrom_lengths
    for g in giter:
        yield GTFLine([
            g.chrom, '.', 'slop',
            max(1, g.start - b), # start
            min(chrom_lengths.setdefault(g.chrom, int(3e9)), g.end + b), # end
            '.', g.strand, '.', g.attr,
        ])


def cluster_gtf(gtf, dist=0, stranded=True, criteria=None):
    bychrom = defaultdict(list)
    for g in gtf:
        if stranded:
            bychrom['%s#%s' % (g.chrom, g.strand)].append(g)
        else:
            bychrom[g.chrom].append(g)
    
    ret = []
    for cchrom, glist in bychrom.iteritems():
        glist.sort(key=lambda x:x.start)
        cur = GTFCluster(glist[0]) if type(glist[0]) is GTFLine else glist[0]
        for g1 in glist[1:]:
            domerge = (g1.start - cur.end) <= dist
            if criteria is not None:
                domerge &= criteria(cur, g1)
            if domerge:
                cur.add(g1)
            else:
                ret.append(cur)
                cur = GTFCluster(g1) if type(g1) is GTFLine else g1
        ret.append(cur)
    return ret


def _chromstrand(g, stranded):
    return '{}#{}'.format(g.chrom, g.strand) if stranded else g.chrom

def intersect_gtf(gtfA, gtfBs, stranded=True):
    bychrom = defaultdict(IntervalTree)
    for i,gtfB in enumerate(gtfBs):
        for h in gtfB:
            bychrom[_chromstrand(h, stranded)].addi(h.start, h.end, (i, h))
    for g in gtfA:
        cchrom = _chromstrand(g, stranded)
        if cchrom in bychrom:
            m = [iv.data for iv in bychrom[cchrom].search(g.start, g.end)]
        else:
            m = []
        yield (g, m)


def conflict_gtf(gtf, dist=0, stranded=False):
    bychrom = defaultdict(list)
    for g in gtf:
        if stranded:
            bychrom['%s#%s' % (g.chrom, g.strand)].append(g)
        else:
            bychrom[g.chrom].append(g)
    
    ret = []
    for cchrom, glist in bychrom.iteritems():
        glist.sort(key=lambda x:x.start)
        tmp = [glist[0],]
        for g1 in glist[1:]:
            isconflict = (g1.start - tmp[-1].end) <= dist
            if isconflict:
                tmp.append(g1)
            else:
                if len(tmp) > 1:
                    ret.append(tmp)
                tmp = [g1, ]            
        if len(tmp) > 1:
            ret.append(tmp)
    
    return ret # return None if len(ret) == 0 else ret

"""
def subtract(gA, gB, both=False):
    r1 = gA.copy()
    r2 = gA.copy()
    r1.end = max(gA.start, min(gA.end, gB.start))
    r2.start = min(max(gA.start, gB.end), gA.end)
    return r1, r2
"""


def write_gtf_file(gtf, outfile, comment=True, span=True):
    outh = open(outfile, 'w') if type(outfile) is str else outfile
    for g in gtf:
        if type(g) is GTFCluster:
            print >>outh, g.display_str(comment, span)
        else:
            print >>outh, str(g)
    if type(outfile) is str: outh.close()

def read_gtf_file(infile, comment='#'):
    return (GTFLine(r) for r in tsv(infile, comment))

def read_gtf_clusters(infile, group_by=None, comment='#'):
    if group_by is None:
        groups = []
        for g in read_gtf_file(infile, comment):
            if g.feature == 'gene':
                groups.append([g])
            else:
                groups[-1].append(g)
    else:
        groups = defaultdict(list)
        for g in read_gtf_file(infile, comment):
            assert group_by in g.attr, 'ERROR: All row must contain group_by attribute "%s"' % group_by
            groups[g.attr[group_by]].append(g)
        groups = groups.values()
    
    clusters = []
    for grp in groups:
        spn = [g for g in grp if g.feature == 'gene']
        spn_attr = spn[0].attr
        clusters.append(GTFCluster([g for g in grp if g.feature != 'gene']))
        clusters[-1].attr = spn_attr
    return clusters
