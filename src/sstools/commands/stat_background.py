#!/usr/bin/env python
import sys
import os
import re
import optparse
import multiprocessing as mp
import numpy as np
import pandas as pd
import pysam

LOWEST_N = 4
CHECK_PATTERN = True
PATTERN = "^chr[0-9]+$"
ANCHOR_START = 0
ANCHOR_CENTER = 1
ANCHOR_END = 2


def get_reads(bamfile, chrom, rm_dup, rm_low_conf):
    crick, watson = 0, 0
    with pysam.AlignmentFile(bamfile) as f:
        for s in f.fetch(chrom):
            if s.is_supplementary or s.is_secondary:
                continue
            if rm_dup and s.is_duplicate:
                continue
            if rm_low_conf:
                if s.has_tag("XH") and s.get_tag("XH") == "N":
                    continue   
            strand = "-" if s.is_reverse else "+"
            if s.is_paired and s.is_read2: # Paired-end
                strand = "-" if strand == "+" else "+"
            if strand == "+":
                crick += 1
            else:
                watson += 1
    return chrom, crick, watson


def stat_background(args=None):
    parser = optparse.OptionParser(usage="%prog [options] input.bam output.tsv")
    parser.add_option("-t", "--threads", dest="threads", metavar="INT", default=1, type="int", 
                      help="Threads for running. [%default]")
    parser.add_option("-s", "--summary", dest="summary", metavar="PATH", default=None, 
                      help="Output summary to PATH. [%default]")
    options, args = parser.parse_args(args)
    bamfile, outfile = args
    threads = options.threads
    smrfile = options.summary
    name = os.path.basename(os.path.splitext(bamfile)[0])

    results = []
    pool = mp.Pool(threads)
    with pysam.AlignmentFile(bamfile) as f:
        for chrom in f.references:
            if CHECK_PATTERN and re.match(PATTERN, chrom) is None:
                continue
            r = pool.apply_async(get_reads, (bamfile, chrom, True, False))
            results.append(r)
    pool.close()
    pool.join()
    
    rows = []
    for r in results:
        chrom, crick, watson = r.get()
        background = np.divide(min(crick, watson), (crick + watson))
        rows.append([chrom, crick, watson, background])
    df = pd.DataFrame(rows, columns=["Chrom", "Crick", "Watson", "Background"])
    assert LOWEST_N < len(df)
    df.to_csv(outfile, sep="\t", index=False)
    
    if smrfile:
        vs = list(sorted(df["Background"]))
        background = np.mean(vs[:LOWEST_N])
        with open(smrfile, "w+") as fw:
            fw.write("Name\tBackground\n")
            fw.write("%s\t%f\n" % (name, background))
    
    
if __name__ == "__main__":
    stat_background()