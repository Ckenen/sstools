#!/usr/bin/env python
import os
import re
import optparse
import multiprocessing as mp
import numpy as np
import pysam


CHECK_PATTERN = True
CHROM_PATTERN = "^chr[0-9]+$"
ANCHOR_START = 0
ANCHOR_CENTER = 1
ANCHOR_END = 2


def _worker(bamfile, chrom, width, rm_dup, rm_low_conf, anchor):
    with pysam.AlignmentFile(bamfile) as f:
        length = f.get_reference_length(chrom)
        nbin = int(length / width)
        if length % width > 0:
            nbin += 1
        counts = np.zeros(nbin)
        for s in f.fetch(chrom):
            if s.is_supplementary or s.is_secondary:
                continue
            if rm_dup and s.is_duplicate:
                continue
            if rm_low_conf and s.has_tag("XH") and s.get_tag("XH") == "N":
                continue
            start = s.reference_start
            end = s.reference_end
            if anchor == ANCHOR_START:
                a = start
            elif anchor == ANCHOR_CENTER:
                a = int((start + end) / 2)
            elif anchor == ANCHOR_END:
                a = end - 1
            else:
                assert False
            idx = int(a / width)
            counts[idx] += 1
        return counts


def calculate_spikiness(counts):
    v1 = 0
    for i in range(1, len(counts)):
        v1 += abs(counts[i] - counts[i - 1])
    v2 = sum(counts)
    spikiness = 0
    if v2 > 0:
        spikiness = v1 / v2
    return spikiness
    

def stat_spikiness(args=None):
    parser = optparse.OptionParser(usage="%prog [options] input.bam output.tsv")
    parser.add_option("-w", "--width", dest="width", type="int", default=1000000, metavar="INT", 
                      help="Width of bin. [%default]")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1, metavar="INT", 
                      help="Threads. [%default]")
    options, args = parser.parse_args(args)
    bamfile, outfile = args
    width = options.width
    threads = options.threads
    name = os.path.basename(os.path.splitext(bamfile)[0])

    pool = mp.Pool(threads)
    results = []
    with pysam.AlignmentFile(bamfile) as f:
        for chrom in f.references:
            if CHECK_PATTERN and re.match(CHROM_PATTERN, chrom) is None:
                continue
            params = (bamfile, chrom, width, True, False, ANCHOR_START)
            r = pool.apply_async(_worker, params)
            results.append(r)
    pool.close()
    pool.join()

    counts = []
    for r in results:
        counts.extend(r.get())
    spikiness = calculate_spikiness(counts)
    
    with open(outfile, "w+") as fw:
        fw.write("\t".join(["Name", "Spikiness"]) + "\n")
        fw.write("\t".join(map(str, [name, spikiness])) + "\n")
        
    
if __name__ == "__main__":
    stat_spikiness()