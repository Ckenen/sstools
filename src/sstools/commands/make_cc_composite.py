#!/usr/bin/env python#!/usr/bin/env python
import sys
import os
import re
import json
import gzip
from collections import Counter
import optparse
import multiprocessing
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "arial"
import matplotlib.pyplot as plt
import pysam
from sstools import utils


usage = """

    sstools MakeCCComposite [options] <regions.bed.gz> <outdir>
"""


class MakeCCComposite(object):
    def __init__(self):
        self.infile = None
        self.outdir = None
        self.init_parameters()
        self.execute()

    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-p", "--processors", dest="threads", type="int", default=1,
                          help="Processors to run in parallel. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        self.options = options
        self.infile = args[0]
        self.outdir = args[1]

    @classmethod
    def process_chromosome(cls, infile, chrom, outfile):
        with pysam.TabixFile(infile) as f, open(outfile, "w+") as fw:
            for line in f.fetch(chrom):
                chrom, start0, end0, cell, score0, strand0 = line.strip("\n").split("\t")
                start0 = int(start0)
                end0 = int(end0)
                run = cell.split(".")[0]
                bamfile = "results/mapping/mark_parental/%s/%s.bam" % (run, cell)
                with pysam.AlignmentFile(bamfile) as bam:
                    for segment in bam.fetch(chrom, start0, end0):
                        if segment.is_duplicate:
                            continue
                        start = segment.reference_start
                        end = segment.reference_end
                        score = segment.mapping_quality
                        strand = "-" if segment.is_reverse else "+"
                        if strand0 == "-":
                            strand = "-" if strand == "+" else "+"
                        color = "107,137,138"
                        if strand == "-":
                            color = "248,173,97"
                        line = "\t".join(map(str, [chrom, start, end, cell, \
                            score, strand, start, end, color]))
                        fw.write(line + "\n")
        return outfile

    def execute(self):
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        with gzip.open(self.infile, "rt") as f:
            chroms = [line.strip("\n").split("\t")[0] for line in f]
            chroms = list(sorted(set(chroms)))

        pool = None
        if self.options.threads > 1:
            pool = multiprocessing.Pool(self.options.threads)
        results = []
        for chrom in chroms:
            args = (self.infile, chrom, "%s/%s.bed" % (self.outdir, chrom))
            if pool is None:
                r = self.process_chromosome(*args)
            else:
                r = pool.apply_async(self.process_chromosome, args)
            results.append(r)
            # break
        if pool is not None:
            pool.close()
            pool.join()
            results = [r.get() for r in results]
        