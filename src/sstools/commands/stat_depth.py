#!/usr/bin/env python
import optparse
import re
import multiprocessing as mp
import pandas as pd
import pysam


DEFAULT_CHROM_PATTERN = "^chr[0-9]+$"


def _worker(bamfile, chrom, rm_dup):
    data = dict()
    with pysam.AlignmentFile(bamfile) as bam:
        length = bam.get_reference_length(chrom)
        bases = 0
        for s in bam.fetch(chrom):
            if rm_dup and s.is_duplicate:
                continue
            start = s.reference_start
            end = s.reference_end
            bases += (end - start)
        data["Chrom"] = chrom
        data["Length"] = length
        data["Base"] = bases
        data["Depth"] = bases / length
    return data


def stat_depth(args=None):
    parser = optparse.OptionParser(usage="%prog [options] input.bam output.tsv")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1, metavar="INT", 
                      help="The number of parallel threads. [%default]")
    parser.add_option("-r", "--remove-duplicate", dest="remove_duplicate", default=False, action="store_true", 
                      help="Remove duplicate reads. [%default]")
    options, args = parser.parse_args(args)
    bamfile, outfile = args
    threads = options.threads
    rmdup = options.remove_duplicate
    
    results = []
    pool = mp.Pool(threads)
    with pysam.AlignmentFile(bamfile) as f:
        for chrom in f.references:
            if re.match(DEFAULT_CHROM_PATTERN, chrom) is None:
                continue
            r = pool.apply_async(_worker, (bamfile, chrom, rmdup))
            results.append(r)
    pool.close()
    pool.join()
    
    rows = []
    for r in results:
        data = r.get()
        rows.append([data["Chrom"], data["Length"], data["Base"], data["Depth"]])

    total_length = sum([row[1] for row in rows])
    total_base = sum([row[2] for row in rows])
    total_depth = total_base / total_length
    rows.append(["Overall", total_length, total_base, total_depth])

    d = pd.DataFrame(rows)
    d.columns = ["Name", "Length", "Base", "Depth"]
    d.to_csv(outfile, sep="\t", index=False)


if __name__ == "__main__":
    stat_depth()