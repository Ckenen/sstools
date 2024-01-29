#!/usr/bin/env python
import sys
import os
import re
import multiprocessing
import optparse
from collections import Counter
import numpy as np
import pysam

REMOVE_DUPLICATES = True
CHECK_PATTERN = True
PATTERN = "^chr([0-9]+|[XY])$"


def _worker(bamfile, fafile, chrom):
    data = []
    with pysam.AlignmentFile(bamfile) as bam, pysam.FastaFile(fafile) as fa:
        for s in bam.fetch(chrom):
            if REMOVE_DUPLICATES and s.is_duplicate:
                continue
            start = s.reference_start
            end = s.reference_end
            counter = Counter(fa.fetch(chrom, start, end).upper())
            gc = counter["G"] + counter["C"]
            at = counter["A"] + counter["T"]
            ratio = gc / (gc + at)
            data.append([s.query_name, end - start, at, gc, ratio])
    return data


def stat_gc(args=None):
    parser = optparse.OptionParser(usage="%prog [options] input.bam genome.fasta output.tsv")
    parser.add_option("-t", "--threads", dest="threads", metavar="INT", default=1, type="int", 
                      help="Threads for running. [%default]")
    # parser.add_option("-l", "--layout", dest="layout", metavar="STR", default="SE", 
    #                   help="Layout of input file. [%default]")
    parser.add_option("-s", "--summary", dest="summary", metavar="PATH", default=None, 
                      help="Output summary to PATH. [%default]")
    options, args = parser.parse_args(args)
    bamfile, fafile, outfile = args
    threads = options.threads
    # layout = options.layout
    smrfile = options.summary
    name = os.path.splitext(os.path.basename(bamfile))[0]
    
    # running
    
    results = []
    pool = multiprocessing.Pool(int(threads))
    with pysam.AlignmentFile(bamfile) as bam:
        for chrom in bam.references:
            if CHECK_PATTERN and re.match(PATTERN, chrom) is None:
                continue
            r = pool.apply_async(_worker, (bamfile, fafile, chrom))
            results.append(r)
    pool.close()
    pool.join()
    results = [r.get() for r in results]
    
    # output
    
    gc_ratios = []
    fw = sys.stdout if outfile == "-" else open(outfile, "w+")
    fw.write("Read\tLength\tAT\tGC\tGC.Ratio\n")
    for data in results:
        for row in data:
            gc_ratios.append(row[4])
            fw.write("\t".join(map(str, row)) + "\n")
    fw.close()
    
    # summary
    
    mean, median, std = np.nan, np.nan, np.nan
    if len(gc_ratios) > 0:
        mean = np.mean(gc_ratios)
        median = np.median(gc_ratios)
        std = np.std(gc_ratios)
    fw = sys.stderr if smrfile is None else open(smrfile, "w+")
    fw.write("Name\tGC.Mean\tGC.Median\tGC.Std\n")
    fw.write("\t".join(map(str, [name, mean, median, std])) + "\n")
    fw.close()
        
        
if __name__ == "__main__":
    stat_gc()
    