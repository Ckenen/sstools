#!/usr/bin/env python
import re
import optparse
from collections import Counter
import pickle
import numpy as np
import pandas as pd
import multiprocessing as mp
import pysam


WIDTH = 100
PATTERN = "^chr([0-9]+|[XY])$"
REMOVE_DUPLICATES = True


def load_window_gc_from_fasta(fafile, chrom, width):
    with pysam.FastaFile(fafile) as fa:
        length = fa.get_reference_length(chrom)
        seq = fa.fetch(chrom, 0, length)
        nbin = int(length / width) # Number of window
        if length % width > 0:
            nbin += 1
        vs = [None] * nbin # windows
        for i in range(nbin):
            start = i * width
            end = min(start + width, length)
            if end - start < width:
                continue
            counter = Counter(seq[start:end].upper())
            gc = counter["G"] + counter["C"]
            at = counter["A"] + counter["T"]
            if gc + at != width:
                continue
            vs[i] = gc
        return vs
    

def _worker(fafile, window_gc_list, bamfile, chrom, width):
    m = np.zeros((101, 2), dtype=np.int)
    
    if window_gc_list:
        vs1 = window_gc_list
    else:
        vs1 = load_window_gc_from_fasta(fafile, chrom, width)
    
    with pysam.AlignmentFile(bamfile) as bam:
        vs2 = [0] * len(vs1) # reads
        for s in bam.fetch(chrom):
            if REMOVE_DUPLICATES and s.is_duplicate:
                continue
            if s.is_supplementary or s.is_secondary:
                continue
            start = s.reference_start
            end = s.reference_end
            i1 = int(start / width)
            i2 = int((end - 1) / width)
            for i3 in range(i1, i2 + 1):
                vs2[i3] += 1
        for gc, reads in zip(vs1, vs2):
            if gc is None:
                continue
            m[gc][0] += 1
            m[gc][1] += reads
    return m


usage = """%prog [options] input.bam genome.fasta output.tsv

For a large number of samples, we should build a prepare file (genome.pkl) to save time.
Step 1: %prog --prepare -k genome.pkl genome.fasta
Step 2: %prog -t 4 -k genome.pkl input.bam genome.fasta output.tsv  
"""

def stat_gc_density(args=None):
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-p", "--prepare", dest="prepare", action="store_true", 
                      help="Generate genome.pkl file.")
    parser.add_option("-t", "--threads", dest="threads", metavar="INT", default=1, type="int", 
                      help="Threads for running. [%default]")
    parser.add_option("-k", "--pickle", dest="pklfile", metavar="STR",
                      help="Pre-calculated genome GC in pickle format. [%default]")
    options, args = parser.parse_args(args)
    pklfile = options.pklfile
    threads = options.threads
    
    if options.prepare:
        fafile = args[0]
        
        data = dict()
        pool = mp.Pool(threads)
        with pysam.FastaFile(fafile) as fa:
            for chrom in fa.references:
                params = (fafile, chrom, WIDTH)
                r = pool.apply_async(load_window_gc_from_fasta, params)
                d = dict()
                d["Chrom"] = chrom
                d["Length"] = fa.get_reference_length(chrom)
                d["Width"] = WIDTH
                d["GC"] = r
                data[chrom] = d
        pool.close()
        pool.join()
        
        for d in data.values():
            d["GC"] = d["GC"].get()
        
        with open(pklfile, "wb") as fw:
            pickle.dump(data, fw)
            
    else:
        bamfile, fafile, outfile = args

        data = None if pklfile is None else pickle.load(open(pklfile, "rb"))
        
        results = []
        pool = mp.Pool(threads)
        with pysam.FastaFile(fafile) as fa:
            for chrom in fa.references:
                if re.match(PATTERN, chrom) is None:
                    continue
                window_gc_list = None if data is None else data[chrom]["GC"]
                params = (fafile, window_gc_list, bamfile, chrom, WIDTH)
                r = pool.apply_async(_worker, params)
                results.append(r)
        pool.close()
        pool.join()
        
        m = np.zeros((101, 2), dtype=np.int)
        for r in results:
            m += r.get()
        d = pd.DataFrame(m)
        d.index.name = "GC"
        d.columns = ["Windows", "Reads"]
        d["Density"] = d["Reads"] / d["Windows"]
        d.to_csv(outfile, sep="\t")
            

if __name__ == "__main__":
    stat_gc_density()