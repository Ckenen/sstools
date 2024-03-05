#!/usr/bin/env python
import optparse
import pandas as pd
import multiprocessing as mp
import pysam

ANCHOR_START = 0
ANCHOR_CENTER = 1
ANCHOR_END = 2

def worker(bamfile, chrom, width, rmdup, rmlc, anchor):
    with pysam.AlignmentFile(bamfile) as f:
        length = f.get_reference_length(chrom)
        nbin = int(length / width)
        if length % width > 0:
            nbin += 1
        
        rows = []
        offset = 5
        for i in range(nbin):
            start = i * width
            end = min(start + width, length)
            row = [chrom, i, start, end, end - start, 0, 0, 0, 0, 0, 0]
            rows.append(row)
            
        for s in f.fetch(chrom):
            if s.is_supplementary:
                continue
            if s.is_secondary:
                continue
            if rmdup and s.is_duplicate:
                continue
            if rmlc and s.has_tag("XH") and s.get_tag("XH") == "N":
                continue
                
            start = s.reference_start
            end = s.reference_end
            strand = "-" if s.is_reverse else "+"
            if s.is_paired and s.is_read2: # Paired-end
                strand = "-" if strand == "+" else "+"
                
            hp = "U"
            if s.has_tag("XP"):
                hp = s.get_tag("XP")
                
            if anchor == ANCHOR_START:
                a = start
            elif anchor == ANCHOR_CENTER:
                a = int((start + end) / 2)
            elif anchor == ANCHOR_END:
                a = end - 1
            else:
                assert False
            idx = int(a / width)
            if strand == "+":
                rows[idx][offset] += 1
                if hp == "1":
                    rows[idx][offset + 1] += 1
                elif hp == "2":
                    rows[idx][offset + 2] += 1
            else:
                rows[idx][offset + 3] += 1
                if hp == "1":
                    rows[idx][offset + 4] += 1
                elif hp == "2":
                    rows[idx][offset + 5] += 1
                
        columns = ["Chrom", "Bin", "Start", "End", "Width", 
               "Crick", "Crick.HP1", "Crick.HP2", 
               "Watson", "Watson.HP1", "Watson.HP2"]
        df = pd.DataFrame(rows, columns=columns)
        
        return df


def run_pipeline(inbam, outfile, threads=1, width=1000000, rmdup=False, rmlc=False, anchor=ANCHOR_START):
    results = []
    pool = mp.Pool(threads)
    with pysam.AlignmentFile(inbam) as f:
        for chrom in f.references:
            args = (inbam, chrom, width, rmdup, rmlc, anchor)
            r = pool.apply_async(worker, args)
            results.append(r)
    pool.close()
    pool.join()
    
    dats = [r.get() for r in results]
    dat = pd.concat(dats, axis=0)
    dat.to_csv(outfile, sep="\t", index=False)
    
usage = """
    sstools StatBinRead [options] <input.bam> <outfile>
"""

def stat_bin_read(args):
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-d", "--remove-duplicates", dest="rmdup", action="store_true", default=False,
                      help="Remove duplicate reads.")
    parser.add_option("-c", "--remove-low-confidence", dest="rmlc", action="store_true", default=False, 
                      help="Remove low confidence reads.")
    parser.add_option("-w", "--width", dest="width", type="int", default=1000000, metavar="INT", 
                      help="")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1, metavar="INT", 
                      help="")
    parser.add_option("-a", "--anchor", dest="anchor", type="int", default=ANCHOR_START, metavar="INT", 
                      help="Anchor to assign bin. (0: start, 1: center, 2: end) [%default]")
    options, args = parser.parse_args(args)
    inbam, outfile = args
    
    run_pipeline(
        inbam=inbam, 
        outfile=outfile, 
        threads=options.threads, 
        width=options.width, 
        rmdup=options.rmdup, 
        rmlc=options.rmlc,
        anchor=options.anchor
    )


if __name__ == "__main__":
    stat_bin_read()