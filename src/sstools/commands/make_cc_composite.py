#!/usr/bin/env python#!/usr/bin/env python
import sys
import os
import gzip
import optparse
import multiprocessing as mp
import pysam


def _load_bam_file_config(tsvfile):
    bamfiles = dict()
    for line in open(tsvfile):
        sample, tsvfile = line.strip("\n").split("\t")[:2]
        bamfiles[sample] = tsvfile
    return bamfiles


def _get_chromosomes(bedfile):
    chroms = []
    for line in gzip.open(bedfile, "rt"):
        chrom = line.strip("\n").split("\t")[0]
        chroms.append(chrom)
    chroms = list(sorted(set(chroms)))
    return chroms


def _get_strand(segment, region_strand):
    strand = "-" if segment.is_reverse else "+"
    if region_strand == "-":
        strand = "-" if strand == "+" else "+"
    return strand


def _get_color(strand):
    color = "0,0,0"
    if strand == "+":
        color = "107,137,138"
    elif strand == "-":
        color = "248,173,97"
    return color


def _worker(bedfile, chrom, bamfiles, outfile):
        with pysam.TabixFile(bedfile) as f, open(outfile, "w+") as fw:
            for line in f.fetch(chrom):
                chrom, start_r, end_r, cell, score_r, strand_r = line.strip("\n").split("\t")
                start_r = int(start_r)
                end_r = int(end_r)
                bamfile = bamfiles[cell]
                with pysam.AlignmentFile(bamfile) as bam:
                    for segment in bam.fetch(chrom, start_r, end_r):
                        if segment.is_duplicate:
                            continue
                        start = segment.reference_start
                        end = segment.reference_end
                        score = segment.mapping_quality
                        strand = _get_strand(segment, strand_r)
                        color = _get_color(strand)
                        row = [chrom, start, end, cell, score, strand, start, end, color]
                        fw.write("\t".join(map(str, row)) + "\n")
        return outfile


def make_cc_composite(args=None):
    parser = optparse.OptionParser(usage="%prog [options] <regions.bed.gz> <bamlist.tsv> <outdir>")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1,
                      help="Threads to run in parallel. [%default]")
    options, args = parser.parse_args(args)
    threads = options.threads
    bedfile, tsvfile, outdir = args
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    bamfiles = _load_bam_file_config(tsvfile)
    chroms = _get_chromosomes(bedfile)

    pool = mp.Pool(threads)
    for chrom in chroms:
        params = (bedfile, chrom, bamfiles, os.path.join(outdir, "%s.bed" % chrom))
        pool.apply_async(_worker, params)
    pool.close()
    pool.join()
            
