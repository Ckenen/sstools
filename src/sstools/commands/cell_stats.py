#!/usr/bin/env python#!/usr/bin/env python
import sys
import os
import re
import json
from collections import Counter
import optparse
import multiprocessing
import numpy as np
import pysam
from sstools import utils


usage = """

    sstools CellStats [options] <input.bam> <genome.fasta> <outfile>
"""


class CellStats(object):
    def __init__(self):
        self.bamfile = None
        self.fastafile = None
        self.outfile = None
        self.init_parameters()
        self.execute()

    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-b", "--bin-width", dest="bin_width", type="int", default=1000000,
                          help="Bin width. [%default]")
        parser.add_option("-p", "--processors", dest="threads", type="int", default=1,
                          help="Processors to run in parallel. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 3:
            parser.print_help()
            exit(1)
        self.options = options
        self.bamfile = args[0]
        self.fastafile = args[1]
        self.outfile = args[2]

    @classmethod
    def process_chromosome(cls, bamfile, fastafile, chrom, is_paired, bin_width):
        with pysam.AlignmentFile(bamfile) as bam, pysam.FastaFile(fastafile) as fasta:
            length = bam.get_reference_length(chrom)
            bin_count = int(length / bin_width)
            if length % bin_width > 0:
                bin_count += 1
            counts1 = [0] * bin_count
            counts2 = [0] * bin_count
            gc_values = []
            total_bases = 0
            total_reads = 0
            for fragment in utils.load_fragments(bamfile, chrom, is_paired):
                chrom, start, end, segment = fragment[:4]
                idx = int(start / bin_width)
                if segment.is_reverse:
                    counts2[idx] += 1
                else:
                    counts1[idx] += 1
                total_reads += 1
                seq = fasta.fetch(chrom, start, end).upper()
                counter = Counter(seq)
                gc = (counter["G"] + counter["C"]) / len(seq)
                gc_values.append(gc)
                total_bases += len(seq)
        c1, c2 = sum(counts1), sum(counts2)
        bg = 0
        if c1 + c2 > 0:
            bg = min(c1, c2) / (c1 + c2)

        data = dict()
        data["Chrom"] = chrom
        data["Is_Paired"] = is_paired
        data["Bin_Width"] = bin_width
        data["Length"] = length
        data["Bin_Count"] = bin_count
        data["Bin_Values"] = [counts1, counts2]
        data["GC_Values"] = gc_values
        data["Total_Bases"] = total_bases
        data["Total_Reads"] = total_reads
        data["Background"] = bg
        return data

    def execute(self):
        pool = None
        if self.options.threads > 1:
            pool = multiprocessing.Pool(self.options.threads)
        results = []

        is_paired = utils.infer_library_layout(self.bamfile)
        with pysam.AlignmentFile(self.bamfile) as f:
            for chrom in f.references:
                if re.match("^chr([0-9]+|[XY])$", chrom) is None:
                    continue
                args = (self.bamfile, self.fastafile, chrom,
                        is_paired, self.options.bin_width)
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

        # Total reads
        total_reads = sum([r["Total_Reads"] for r in results])

        # Total base
        total_bases = sum([r["Total_Bases"] for r in results])

        # Total size
        total_size = sum([r["Length"] for r in results])
        depth = 0
        if total_size > 0:
            depth = total_bases / total_size

        # Background
        bgs = list(sorted([r["Background"] for r in results]))
        background = np.mean(bgs[:4])

        # GC content
        gc_values = []
        for r in results:
            gc_values.extend(r["GC_Values"])
            r.pop("GC_Values")
        gc_mean = np.mean(gc_values)
        gc_median = np.median(gc_values)
        gc_std = np.std(gc_values)

        # Spikiness
        counts = []
        for r in results:
            for c1, c2 in zip(r["Bin_Values"][0], r["Bin_Values"][1]):
                counts.append(c1 + c2)
        v1 = 0
        for i in range(1, len(counts)):
            v1 += abs(counts[i] - counts[i - 1])
        v2 = sum(counts)
        spikiness = 0
        if v2 > 0:
            spikiness = v1 / v2

        # report
        data = dict()
        data["Total_Reads"] = total_reads
        data["Total_Bases"] = total_bases
        data["Total_Size"] = total_size
        data["Depth"] = depth
        data["GC_Mean"] = gc_mean
        data["GC_Median"] = gc_median
        data["GC_Std"] = gc_std
        data["Background"] = background
        data["Spikiness"] = spikiness
        data["Detail"] = results
        with open(self.outfile, "w+") as fw:
            json.dump(data, fw, indent="\t")

        print("File\tReads\tDepth\tGC_Mean\tGC_Median\tGC_Std\tBackground\tSpikiness")
        print(os.path.basename(self.bamfile), total_reads, depth,
              gc_mean, gc_median, gc_std, background, spikiness, sep="\t")
