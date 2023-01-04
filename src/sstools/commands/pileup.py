#!/usr/bin/env python
import sys
import os
import re
import json
import optparse
import multiprocessing
import gzip
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
import pysam
from pyBioInfo.IO.File import Alignment
from pyBioInfo.Utils import SegmentTools, BundleBuilder
from sstools import utils

usage = """

    sstools Pileup [options] <input.bam> <outfile|outdir>
"""



class Pileup(object):
    def __init__(self):
        self.bamfile = None
        self.outfile = None
        self.init_parameters()

        self.execute()

    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-f", "--fasta", dest="fasta", help="[%default]")
        # parser.add_option("-c", "--chrom", dest="chrom", help="[%default]")
        # parser.add_option("-s", "--start", dest="start", type="int", help="[%default]")
        # parser.add_option("-e", "--end", dest="end", type="int", help="[%default]")
        # parser.add_option("-s", "--strand", dest="strand", default=None, help="[%default]")
        #parser.add_option("-m", "--multi-file", dest="multi_file", action="store_true", default=False, help="[%default]")
        #parser.add_option("-g", "--split-RG", dest="split_rg", action="store_true", default=False, help="[%default]")
        parser.add_option("-p", "--processors", dest="threads", type="int", default=1, \
            help="Processors to run in parallel. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        self.options = options
        self.bamfile = args[0]
        self.outfile = args[1]

    @classmethod
    def load_alignments(cls, bamfile, chrom, start, end, strand=None):
        with pysam.AlignmentFile(bamfile) as f:
            for segment in f.fetch(chrom, start, end):
                obj = Alignment(segment)
                if strand is not None and obj.strand != strand:
                    continue
                obj.events = SegmentTools.get_events(segment)
                yield obj

    @classmethod
    def counter_list_to_str(cls, counter_list):
        items1 = []
        for counter in counter_list:
            items2 = []
            for k, v in sorted(counter.items()):
                if v == 0:
                    continue
                items2.append("%s:%d" % (k, v))
            if len(items2) > 0:
                items1.append(",".join(items2))
        if len(items1) > 0:
            return ";".join(items1)
        else:
            return ""

    @staticmethod
    def execute_bin(bamfile, fasta, chrom, bin_start, bin_end, strand=None, outfile_or_fw=None):
        fh_fasta = pysam.FastaFile(fasta)

        fw = None
        need_close = False
        if isinstance(outfile_or_fw, str):
            if outfile_or_fw.endswith(".gz"):
                fw = gzip.open(outfile_or_fw, "wt")
            else:
                fw = open(outfile_or_fw, "w+")
            need_close = True
        else:
            fw = outfile_or_fw
            need_close = False
        
        alignments = Pileup.load_alignments(bamfile, chrom, bin_start, bin_end, strand)
        bundles = BundleBuilder(alignments, keep=True)
        for n, bundle in enumerate(bundles):
            # if n % 100 == 0:
            #     pass
                # print(n)
            bundle_start = bundle.start_min
            bundle_end = bundle.end_max
            width = bundle_end - bundle_start
            sequence = fh_fasta.fetch(chrom, bundle_start, bundle_end)

            data = defaultdict(list)
            for obj in bundle.data:
                data[obj.segment.get_tag("RG")].append(obj)

            array1 = [] # read groups
            array2 = [] # alignments
            array3 = [] # coverages
            array4 = [] # events
            for read_group, reads in data.items():
                array1.append(read_group)
                array2.append(reads)
                covs = np.zeros(width, dtype=np.int)
                for read in reads:
                    for idx in np.arange(read.start - bundle_start, read.end - bundle_start):
                        covs[idx] += 1
                array3.append(covs)

                events = defaultdict(list)
                for read in reads:
                    for e in read.events:
                        if e[1] == "-":
                            continue
                        for bi in range(len(e[1])):
                            events[e[0] + bi].append(e[2])
                array4.append(events)

            for pos in range(bundle_start, bundle_end):
                if pos < bin_start:
                    continue
                elif pos >= bin_end:
                    break
                idx = pos - bundle_start # offset to bundle start
                covs_all = [covs[idx] for covs in array3]
                cell_count = 0
                read_count = 0
                ref = sequence[idx]
                counter_list = [] # base counter
                for i, cov in enumerate(covs_all):
                    if cov == 0:
                        continue
                    cell_count += 1
                    read_count += cov
                    base_counter = Counter(array4[i][pos])
                    tot = sum(base_counter.values())
                    ref_count = cov - tot
                    if ref_count > 0:
                        base_counter[ref] = int(ref_count)
                    counter_list.append(base_counter)
                line = "\t".join(map(str, [pos, ref, read_count, cell_count, Pileup.counter_list_to_str(counter_list)]))
                fw.write(line + "\n")
        if need_close:
            fw.close()
        fh_fasta.close()
        
            
    def execute(self):
        bamfile = self.bamfile
        fasta = self.options.fasta
        outdir = self.outfile
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        pool = None
        results = []
        if self.options.threads > 1:
            pool = multiprocessing.Pool(self.options.threads)
        with pysam.AlignmentFile(self.bamfile) as f:
            for chrom in f.references:
                length = f.get_reference_length(chrom)
                for strand in ["+", "-"]:
                    outfile = os.path.join(outdir, "%s.%s.tsv.gz" % (chrom, strand))
                    if pool is None:
                        Pileup.execute_bin(bamfile, fasta, chrom, 0, length, strand, outfile)
                    else:
                        r = pool.apply_async(Pileup.execute_bin, (bamfile, fasta, chrom, 0, length, strand, outfile))
                        results.append(r)
        if pool is not None:
            pool.close()
            pool.join()
            for r in results:
                assert r.successful()