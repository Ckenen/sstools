#!/usr/bin/env python
import sys
import optparse
import re
import json
import multiprocessing as mp
import pysam


DEFAULT_CHROM_PATTERN = "^chr[0-9]+$"


def worker(bamfile, chrom):
    data = dict()
    with pysam.AlignmentFile(bamfile) as bam:
        length = bam.get_reference_length(chrom)
        bases = 0
        for s in bam.fetch(chrom):
            start = s.reference_start
            end = s.reference_end
            bases += (end - start)
        data["Chrom"] = chrom
        data["Length"] = length
        data["Base"] = bases
        data["Depth"] = bases / length
    return data


def cal_genome_depth(args):
    parser = optparse.OptionParser(usage="%prog CalGenomeDepth [options] input.bam")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1, metavar="INT", 
                      help="The number of parallel threads. [%default]")
    parser.add_option("-o", "--output", dest="output", metavar="PATH", 
                      help="Output to PATH. [%default]")
    options, args = parser.parse_args(args)
    
    bamfile = args[0]
    
    results = []
    pool = mp.Pool(options.threads)
    with pysam.AlignmentFile(bamfile) as f:
        for chrom in f.references:
            r = pool.apply_async(worker, (bamfile, chrom))
            results.append(r)
    pool.close()
    pool.join()
    
    all_results = [r.get() for r in results]
    results = list(filter(lambda r: re.match(DEFAULT_CHROM_PATTERN, r["Chrom"]), all_results))
    
    data = dict()
    data["Length"] = sum([r["Length"] for r in results])
    data["Base"] = sum([r["Base"] for r in results])
    data["Depth"] = data["Base"] / data["Length"]
    data["Chromosomes"] = results
    
    if options.output:
        with open(options.output, "w+") as fw:
            json.dump(data, fw, indent=4)
    else:
        json.dump(data, sys.stdout, indent=4)
        