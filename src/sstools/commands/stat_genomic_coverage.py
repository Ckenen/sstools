#!/usr/bin/env python
import sys
import random
import optparse
import re
import multiprocessing
import gzip
from collections import defaultdict
import pandas as pd
import pysam


usage = """

    sstools StatGenomicCoverage [options] <input.bam> <output.tsv>
"""

class StatGenomicCoverage(object):
    def __init__(self):
        self.infile = None # input.bam
        self.outfile = None # output.tsv
        self.options = None
        self.init_parameters()
        self.percentages = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        self.execute()

    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-p", "--processors", dest="threads", type="int", default=1,
                          help="Processors to run in parallel. [%default]")
        parser.add_option("-b", "--blank-region", dest="blank", default=None,
                          help="BED format Ns regions in genome. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        self.options = options
        self.infile = args[0]
        self.outfile = args[1]

    @staticmethod
    def stat_by_chrom(path, chrom, percs, seed):
        regions = []
        with pysam.AlignmentFile(path) as f:
            for segment in f.fetch(chrom):
                if segment.is_unmapped:
                    continue
                if segment.is_secondary:
                    continue
                if segment.is_supplementary:
                    continue
                if segment.is_paired and (not segment.is_proper_pair):
                    continue
                start = segment.reference_start
                end = segment.reference_end
                regions.append([start, end])
        random.seed(seed)
        n1 = len(regions) # total reads
        result = []
        for p in percs:
            b = 0 # bases
            n2 = int(n1 * p) # sample reads
            if n2 > 0:
                tmp = [[x, y] for x, y in random.sample(regions, n2)]
                tmp.sort()
                i = 0
                while i < len(tmp) - 1:
                    r1 = tmp[i]
                    r2 = tmp[i + 1]
                    if r1[1] >= r2[0]:
                        r1[1] = max(r1[1], r2[1])
                        tmp.pop(i + 1)
                    else:
                        i += 1
                for x, y in tmp:
                    b += y - x
            result.append([chrom, p, n1, n2, b])
        return result

    def execute(self):
        pool = multiprocessing.Pool(self.options.threads)
        results = []
        chroms = []
        lengths = dict()
        with pysam.AlignmentFile(self.infile) as f:
            for i, chrom in enumerate(f.references):
                chroms.append(chrom)
                lengths[chrom] = f.get_reference_length(chrom)
                res = pool.apply_async(self.stat_by_chrom, (self.infile, chrom, self.percentages, i))
                results.append(res)
        pool.close()
        pool.join()

        rows = []
        for chrom, item in zip(chroms, results):
            assert item.successful()
            res = item.get()
            for row in res:
                rows.append(row)
                
        dat = pd.DataFrame(rows)
        dat.columns = ["Chrom", "Percentage", "Reads", "Sample", "Base"]
        
        # Ns
        blanks = defaultdict(int)
        if self.options.blank is not None:
            with gzip.open(self.options.blank, "rt") as f:
                for line in f:
                    chrom, start, end = line.strip("\n").split("\t")[:3]
                    blanks[chrom] += int(end) - int(start)
                
        dat["Length"] = [lengths[c] for c in dat["Chrom"]]
        dat["Blank"] = [blanks[c] for c in dat["Chrom"]]
        dat["Fill"] = dat["Length"] - dat["Blank"]
        dat["Coverage"] = dat["Base"] / dat["Fill"]
        
        chroms = list(filter(lambda c: re.match("^chr[0-9]+$", c) is not None, chroms))

        rows = []
        for p in self.percentages:
            d = dat[dat["Percentage"] == p]
            d = d[[x in chroms for x in d["Chrom"]]]
            reads = d["Reads"].sum()
            sample = d["Sample"].sum()
            base = d["Base"].sum()
            length = d["Length"].sum()
            blank = d["Blank"].sum()
            row = [p, reads, sample, base, length, blank]
            rows.append(row)
        dat = pd.DataFrame(rows)
        dat.columns = ["Percentage", "Reads", "Sample", "Base", "Length", "Blank"]
        dat["Fill"] = dat["Length"] - dat["Blank"]
        dat["Coverage"] = dat["Base"] / dat["Fill"]
        dat.to_csv(self.outfile, sep="\t", index=False)

