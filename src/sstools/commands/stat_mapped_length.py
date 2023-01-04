#!/usr/bin/env python
import sys
import os
import optparse
from collections import defaultdict
import pysam
import multiprocessing
import numpy as np
from sstools.utils import load_fragments

usage = """

    sstools StatMappedLen [options] <input.bam> <out.tsv>
"""

class StatMappedLen(object):
    def __init__(self):
        self.infile = None # input.bam
        self.outfile = None # outfile.tsv
        self.options = None
        self.init_parameters()
        self.execute()

    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-P", "--paired-end", dest="is_paired", action="store_true", default=False,
                          help="Input is paired-end reads. [%default]")
        parser.add_option("-p", "--processors", dest="threads", type="int", default=1,
                          help="Processors to run in parallel. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        self.options = options
        self.infile = args[0]
        self.outfile = args[1]

    def stat_chrom_mapped_length(cls, bamfile, chrom, is_paired=False):
        data = dict()
        detail = []
        counter = defaultdict(int)
        for fragment in load_fragments(bamfile, chrom, is_paired):
            chrom, start, end = fragment[:3]
            name = fragment[3].query_name
            length = end - start
            detail.append([name, length])
            counter[length] += 1
        data["Detail"] = detail
        data["Summary"] = counter
        return data
    
    def execute(self):
        pool = None
        if self.options.threads > 1:
            pool = multiprocessing.Pool(self.options.threads)
        results = []
        with pysam.AlignmentFile(self.infile) as f:
            for chrom in f.references:
                args = (self.infile, chrom, self.options.is_paired)
                if pool:
                    r = pool.apply_async(self.stat_chrom_mapped_length, args)
                else:
                    r = self.stat_chrom_mapped_length(*args)
                results.append(r)
        if pool:
            pool.close()
            pool.join()
            results = [r.get() for r in results]
            
        lengths = []
        with open(self.outfile, "w+") as fw:
            fw.write("Read\tLength\n")
            for data in results:
                for name, length in data["Detail"]:
                    lengths.append(length)
                    fw.write("%s\t%d\n" % (name, length))
        with open(self.outfile + ".counter", "w+") as fw:
            fw.write("Length\tCount\n")
            counter = defaultdict(int)
            for data in results:
                for k, v in data["Summary"].items():
                    counter[k] += v
            for k, v in sorted(counter.items()):
                fw.write("%s\t%s\n" % (k, v))
        with open(self.outfile + ".summary", "w+") as fw:
            fw.write("Name\tCount\tMean\tMedian\tStd\n")
            fw.write("\t".join(map(str, [
                os.path.basename(os.path.splitext(self.infile)[0]), 
                len(lengths), 
                np.mean(lengths), 
                np.median(lengths), 
                np.std(lengths)])) + "\n")