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
        
        self.fh_ds = None
        self.fh_read = None
        if self.options.report:
            self.fh_ds = open(self.outfile + ".ds.txt", "w+")
            self.fh_ds.write("DuplicateSetName\tDuplicatSetSize\tReadNames\n")
            self.fh_read = open(self.outfile + ".read.txt", "w+")
            self.fh_read.write("ReadName\tDuplicateSetName\tDuplicateSetSize\tDuplicateSetIndex\n")
            
        self.execute()
        
        if self.fh_ds is not None:
            self.fh_ds.close()
        if self.fh_read is not None:
            self.fh_read.close()
        
        print("Input: %s" % self.infile)
        print("Output: %s" % self.outfile)
        print("Total: %d" % self.total)
        print("Uniq: %d (%.2f%%)" % (self.uniq, utils.divide_zero(self.uniq * 100, self.total)))

    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-d", "--max-diff", dest="max_diff", default=20, type="int",
                          help="Maximal difference. [%default]")
        parser.add_option("-r", "--report", dest="report", action="store_true", default=False, 
                          help="Report DuplicateSet details. (output.bam.ds.txt and output.bam.read.txt) [%default]")
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
            for segments_clustered_by_start in self.clustering_by_start(segments):
                for segments_duplicate_set in self.clustering_by_end(segments_clustered_by_start):
                    self.uniq += 1
                    for segment in segments_duplicate_set:
                        segment.set_tag("DN", self.duplicate_set_name)
                    self.mark_duplicates(segments_duplicate_set)
                    if self.fh_ds is not None:
                        read_names = ",".join([item.query_name for item in segments_duplicate_set])
                        line = "\t".join(map(str, [self.duplicate_set_name, len(segments_duplicate_set), read_names]))
                        self.fh_ds.write(line + "\n")
                    self.duplicate_set_name += 1
            for segment in segments:
                fw.write(segment)
                if self.fh_read is not None:
                    line = "\t".join(map(str, [segment.query_name, 
                                            segment.get_tag("DN"), 
                                            segment.get_tag("DS"), 
                                            segment.get_tag("DI")]))
                    self.fh_read.write(line + "\n")
        f.close()
        fw.close()
        
