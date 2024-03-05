#!/usr/bin/env python
import sys
import os
import re
from collections import defaultdict
from optparse import OptionParser
import multiprocessing
import numpy as np
import pysam


REMOVE_DUPLICATES = True
CHECK_PATTERN = True
PATTERN = "^chr([0-9]+|[XY])$"


def _worker(bamfile, chrom, layout):
    array = []
    with pysam.AlignmentFile(bamfile) as bam:
        if layout.upper() == "SE":
            for s in bam.fetch(chrom):
                if REMOVE_DUPLICATES and s.is_duplicate:
                    continue
                if s.is_supplementary or s.is_secondary:
                    continue
                start = s.reference_start
                end = s.reference_end
                array.append([s.query_name, end - start])
        elif layout == "PE":
            segments = defaultdict(list)
            for s in bam.fetch(chrom):
                if REMOVE_DUPLICATES and s.is_duplicate:
                    continue
                if s.is_supplementary or s.is_secondary:
                    continue
                if s.is_proper_pair:
                    segments[s.query_name].append(s)
            for query_name, items in segments.items():
                if len(items) != 2:
                    continue
                read1, read2 = items
                start = min(read1.reference_start, read2.reference_start)
                end = max(read1.reference_end, read2.reference_end)
                array.append([query_name, end - start])
        else:
            raise RuntimeError()
    return array


def stat_length(args=None):
    parser = OptionParser(usage="%prog [options] input.bam output.tsv")
    parser.add_option("-t", "--threads", dest="threads", metavar="INT", default=1, type="int", 
                      help="Threads for running. [%default]")
    parser.add_option("-l", "--layout", dest="layout", metavar="STR", default="SE", 
                      help="Layout of input file. [%default]")
    parser.add_option("-s", "--summary", dest="summary", metavar="PATH", default=None, 
                      help="Output summary to PATH. [%default]")
    options, args = parser.parse_args(args)
    bamfile, outfile = args
    threads = options.threads
    layout = options.layout
    smrfile = options.summary
    name = os.path.splitext(os.path.basename(bamfile))[0]
    
    results = []    
    pool = multiprocessing.Pool(threads)
    with pysam.AlignmentFile(bamfile) as f:
        for chrom in f.references:
            if CHECK_PATTERN and re.match(PATTERN, chrom) is None:
                continue
            r = pool.apply_async(_worker, (bamfile, chrom, layout))
            results.append(r)
    pool.close()
    pool.join()
    results = [r.get() for r in results]
    
    lengths = list()
    fw = sys.stdout if (outfile is None or outfile == "-") else open(outfile, "w+")
    fw.write("Read\tLength\n")
    for r in results:
        for query_name, length in r:
            lengths.append(length)
            fw.write("%s\t%d\n" % (query_name, length))
    fw.close() 
    
    fw = sys.stderr if smrfile is None else open(smrfile, "w+")
    mean, median, std, vmin, vmax = np.nan, np.nan, np.nan, np.nan, np.nan
    if len(lengths) > 0:
        mean = np.mean(lengths) 
        median = np.median(lengths)
        std = np.std(lengths)
        vmin = np.min(lengths)
        vmax = np.max(lengths)
    fw.write("Name\tLength.Mean\tLength.Median\tLength.Std\tLength.Min\tLength.Max\n")
    fw.write("\t".join(map(str, [name, mean, median, std, vmin, vmax])) + "\n")
    fw.close() 
    

if __name__ == "__main__":
    stat_length()
    
    
    