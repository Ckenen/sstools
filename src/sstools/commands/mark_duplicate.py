#!/usr/bin/env python
import sys
import optparse
import pysam
from sstools import utils

usage = """

    sstools MarkDuplicate [options] <input.bam> <output.bam>

This command will added three TAGs to the BAM file:
    DN: duplicate set name
    DS: duplicate set size
    DI: duplicate set index
"""


class MarkDuplicate(object):
    def __init__(self):
        self.options = None
        self.infile = None
        self.outfile = None
        self.init_parameters()

        self.total = 0
        self.uniq = 0
        self.duplicate_set_name = 0

        self.execute()

        print("Input: %s" % self.infile)
        print("Output: %s" % self.outfile)
        print("Total: %d" % self.total)
        print("Uniq: %d (%.2f%%)" % (self.uniq, utils.divide_zero(self.uniq * 100, self.total)))

    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-d", "--max-diff", dest="max_diff", default=20, type="int",
                          help="Maximal difference. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        self.options = options
        self.infile = args[0]
        self.outfile = args[1]

    def clustering_by_start(self, segments):
        array = None
        for segment in sorted(segments, key=lambda item: item.reference_start):
            if array is None:
                array = [segment]
            elif segment.reference_start - array[-1].reference_start <= self.options.max_diff:
                array.append(segment)
            else:
                yield array
                array = [segment]
        if array is not None:
            yield array

    def clustering_by_end(self, segments):
        array = None
        for segment in sorted(segments, key=lambda item: item.reference_end):
            if array is None:
                array = [segment]
            elif segment.reference_end - array[-1].reference_end <= self.options.max_diff:
                array.append(segment)
            else:
                yield array
                array = [segment]
        if array is not None:
            yield array

    def mark_duplicates(self, segments):
        segments = list(sorted(
            segments, key=lambda item: item.reference_end - item.reference_start, reverse=True))
        duplicate_set_size = len(segments)
        for i, segment in enumerate(segments):
            if i == 0:
                segment.flag = segment.flag & 0b101111111111
            else:
                segment.flag = segment.flag | 0b010000000000
            segment.set_tag("DS", duplicate_set_size)
            segment.set_tag("DI", i)  # duplicate index

    def execute(self):
        f = pysam.AlignmentFile(self.infile)
        header = utils.add_pg_to_header(f.header.as_dict(), cl=" ".join(sys.argv))
        fw = pysam.AlignmentFile(self.outfile, "wb", header=header)
        for chrom in f.references:
            segments = []
            for segment in f.fetch(chrom):
                if segment.is_secondary:
                    continue
                if segment.is_supplementary:
                    continue
                segments.append(segment)
            self.total += len(segments)
            for items1 in self.clustering_by_start(segments):
                for items2 in self.clustering_by_end(items1):
                    self.uniq += 1
                    for segment in items2:
                        segment.set_tag("DN", self.duplicate_set_name)
                    self.mark_duplicates(items2)
                    self.duplicate_set_name += 1
            for segment in segments:
                fw.write(segment)
        f.close()
        fw.close()
        
