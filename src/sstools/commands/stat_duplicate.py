#!/usr/bin/env python
import sys
import optparse
import multiprocessing
import random
import pysam
import pandas as pd
import re

usage = """

    sstools StatDuplicate [options] <input.bam> <output.tsv>
"""

class StatDuplicate(object):
    def __init__(self):
        self.infile = None # input.bam
        self.outfile = None # output.tsv
        self.options = None
        self.init_parameters()
        self.percentages = [1/128, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1.0]
        self.execute()

    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-p", "--processors", dest="threads", type="int", default=1,
                          help="Processors to run in parallel. [%default]")
        parser.add_option("-d", "--max-diff", dest="max_diff", default=20, type="int",
                          help="Maximal difference. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        self.options = options
        self.infile = args[0]
        self.outfile = args[1]

    @staticmethod
    def cluster(regions, max_diff):
        array1 = [] # clustering by start
        tmp = None
        for item in sorted(regions, key=lambda item: item[0]):
            if tmp is None:
                tmp = [item]
            elif item[0] - tmp[-1][0] <= max_diff:
                tmp.append(item)
            else:
                array1.append(tmp)
                tmp = [item]
        if tmp is not None:
            array1.append(tmp)
                
        array2 = [] # clustering by end
        for items in array1:
            tmp = None
            for item in sorted(items, key=lambda item: item[1]):
                if tmp is None:
                    tmp = [item]
                elif item[1] - tmp[-1][1] <= max_diff:
                    tmp.append(item)
                else:
                    array2.append(tmp)
                    tmp = [item]
            if tmp is not None:
                array2.append(tmp)        
        return array2

    @staticmethod
    def stat_by_chrom(path, chrom, percs, max_diff, seed):
        random.seed(seed)
        se = None # single-end
        regions = []
        segments = []
        with pysam.AlignmentFile(path) as f:
            for segment in f.fetch(chrom):
                if segment.is_unmapped:
                    continue
                if se is None:
                    if segment.is_paired:
                        se = False
                    else:
                        se = True
                if se:
                    start = segment.reference_start
                    end = segment.reference_end
                    regions.append([start, end])
                else:
                    if not segment.is_proper_pair:
                        continue
                    segments.append(segment)
        if not se:
            segments = list(sorted(segments, key=lambda item: item.query_name))
            tmp1 = []
            tmp2 = None
            for segment in segments:
                if tmp2 is None:
                    tmp2 = [segment]
                elif segment.query_name == tmp2[0].query_name:
                    tmp2.append(segment)
                else:
                    tmp1.append(tmp2)
                    tmp2 = [segment]
            if tmp2 is not None:
                tmp1.append(tmp2)
            for items in tmp1:
                if len(items) != 2:
                    continue
                r1, r2 = items
                if r1.next_reference_start != r2.reference_start:
                    continue
                if r1.reference_start != r2.next_reference_start:
                    continue
                start = min(r1.reference_start, r2.reference_start)
                end = max(r1.reference_end, r2.reference_end)
                regions.append([start, end])
        
        result = []
        n1 = len(regions) # total reads
        for p in percs:
            n2 = int(n1 * p) # sample reads
            n3 = 0 # uniq reads
            tmp = None
            if n2 == 0:
                tmp = list()
            else:
                tmp = [[x, y] for x, y in random.sample(regions, n2)]
                tmp.sort()
            tmp1 = StatDuplicate.cluster(tmp, max_diff)
            n3 = len(tmp1)
            result.append([chrom, p, n1, n2, n3])
        return result

    def execute(self):
        pool = multiprocessing.Pool(self.options.threads)
        results = []
        chroms = []
        with pysam.AlignmentFile(self.infile) as f:
            for i, chrom in enumerate(f.references):
                chroms.append(chrom)
                results.append(pool.apply_async(self.stat_by_chrom, \
                    (self.infile, chrom, self.percentages, self.options.max_diff, i)))
        pool.close()
        pool.join()
        
        rows = []
        for chrom, item in zip(chroms, results):
            assert item.successful()
            res = item.get()
            for row in res:
                rows.append(row)
                
        dat = pd.DataFrame(rows)
        dat.columns = ["Chrom", "Percentage", "Reads", "Sample", "UniqReads"]
        dat["UniqRatio"] = dat["UniqReads"] / dat["Sample"]
        dat = dat.fillna(0)
        
        chroms = list(filter(lambda c: re.match("^chr[0-9]+$", c) is not None, chroms))
                
        rows = []
        for p in self.percentages:
            d = dat[dat["Percentage"] == p]
            d = d[[x in chroms for x in d["Chrom"]]]
            reads = d["Reads"].sum()
            sample = d["Sample"].sum()
            uniq = d["UniqReads"].sum()
            row = [p, reads, sample, uniq,]
            rows.append(row)
        dat = pd.DataFrame(rows)
        dat.columns = ["Percentage", "Reads", "Sample", "UniqReads"]
        dat["UniqRatio"] = dat["UniqReads"] / dat["Sample"]
        dat = dat.fillna(0)
        dat.to_csv(self.outfile, sep="\t", index=False)      
