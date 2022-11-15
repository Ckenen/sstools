#!/usr/bin/env python
import sys
import optparse
import re
from collections import defaultdict
import pysam
from sstools import utils

usage = """

    sstools FilterBam [options] <input.bam> <output.bam>
"""

class FilterBam(object):
    def __init__(self):
        self.options = None
        self.infile = None
        self.outfile = None
        self.init_parameters()
        
        # summary
        self.total = 0
        self.unmapped = 0
        self.secondary = 0
        self.supplementary = 0
        self.primary = 0
        self.improper_paired = 0
        self.outer_length = 0
        self.outer_mapping_quality = 0
        self.too_many_clipping = 0
        self.exclude_seqname = 0
        self.final = 0
        
        # counter
        self.mapq_counter = defaultdict(int)
        self.seqname_counter = defaultdict(int)
        self.length_counter = defaultdict(int)
        self.clip_counter = defaultdict(int)
        
        # run
        if self.options.require_proper_pair:
            raise NotImplementedError()
        else:
            self.filter_single_end_bam()
            
        self.report_summary()
        
    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        group = optparse.OptionGroup(parser, title="Mapping status")
        group.add_option("-u", "--remove-unmapped", dest="remove_unmapped", action="store_true", default=False, 
                          help="Remove the unmapped alignments. [%default]")
        group.add_option("-s", "--remove-secondary", dest="remove_secondary", action="store_true", default=False, 
                          help="Remove secondary alignments. [%default]")
        group.add_option("-S", "--remove-supplementary", dest="remove_supplementary", action="store_true", default=False, 
                          help="Remove supplementary alignments. [%default]")
        group.add_option("-n", "--seqname-pattern", dest="seqname_pattern", default=None, 
                         help="Only output the chromosomes that match the regular expression. [%default]")
        parser.add_option_group(group)
        group = optparse.OptionGroup(parser, title="Mapping quality")
        group.add_option("-q", "--min-mapq", dest="min_mapq", default=0, type="int", 
                         help="Minimal mapping quality. Only available for primary alignments. [%default]")
        group.add_option("-Q", "--max-mapq", dest="max_mapq", default=int(1e9), type="int", 
                         help="Maximal mapping quality. Only available for primary alignments. [%default]")
        parser.add_option_group(group)
        group = optparse.OptionGroup(parser, title="Mapped region length")
        group.add_option("-m", "--min-length", dest="min_length", default=0, type="int", 
                         help="Minimal mapped region length. Only available for primary alignments. Including deletion and excluding insertion and clipping. [%default]")
        group.add_option("-M", "--max-length", dest="max_length", default=int(1e9), type="int", 
                         help="Maximal mapped region length. Only available for primary alignments. Including deletion and excluding insertion and clipping. [%default]")
        parser.add_option_group(group)
        group = optparse.OptionGroup(parser, title="Soft/hard clipping")
        group.add_option("-c", "--max-clipping", dest="max_clipping", default=int(1e9), type="int", 
                         help="Maximal soft-clipping and hard-clipping. Only available for primary alignments. [%default]")
        parser.add_option_group(group)
        group = optparse.OptionGroup(parser, title="Paired-end")
        group.add_option("-r", "--proper-pair", dest="require_proper_pair", action="store_true", default=False, 
                         help="Only output the proper-paired alignments. [%default]")
        group.add_option("-b", "--both-pass", dest="both_pass", action="store_true", default=False, 
                         help="Only output the proper-paired alignments that both mates satisfy all threshold simultaneously. [%default]")
        parser.add_option_group(group)
                  
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        self.options = options
        self.infile = args[0]
        self.outfile = args[1]
        
    def report_summary(self):
        print("Input: %s" % self.infile)
        print("Output: %s" % self.outfile)
        print("Total: %d" % self.total)
        print("Removing details:")
        print(" - Unmapped: %d (%.2f%%)" % (self.unmapped, utils.divide_zero(self.unmapped * 100, self.total)))
        print(" - Secondary: %d (%.2f%%)" % (self.secondary, utils.divide_zero(self.secondary * 100, self.total)))
        print(" - Supplementary: %d (%.2f%%)" % (self.supplementary, utils.divide_zero(self.supplementary * 100, self.total)))
        print(" - Seqname: %d (%.2f%%)" % (self.exclude_seqname, utils.divide_zero(self.exclude_seqname * 100, self.total)))
        print(" - MapQ: %d (%.2f%%)" % (self.outer_mapping_quality, utils.divide_zero(self.outer_mapping_quality * 100, self.total)))
        print(" - Clipping: %d (%.2f%%)" % (self.too_many_clipping, utils.divide_zero(self.too_many_clipping * 100, self.total)))
        print(" - Length: %d (%.2f%%)" % (self.outer_length, utils.divide_zero(self.outer_length * 100, self.total)))
        print("Pass: %d (%.2f%%)" % (self.final, utils.divide_zero(self.final * 100, self.total)))
        
    def check_mapping_quality(self, segment):
        mapq = segment.mapping_quality
        self.mapq_counter[mapq] += 1
        return mapq >= self.options.min_mapq and mapq <= self.options.max_mapq
    
    def check_clipping(self, segment):
        five_clipping = 0
        three_clipping = 0
        cigars = segment.cigartuples
        if cigars[0][0] == pysam.CHARD_CLIP or cigars[0][0] == pysam.CSOFT_CLIP:
            five_clipping = cigars[0][1]
        if cigars[-1][0] == pysam.CHARD_CLIP or cigars[-1][0] == pysam.CSOFT_CLIP:
            three_clipping = cigars[1][1]
        clipping = max(five_clipping, three_clipping)
        self.clip_counter[clipping] += 1
        return clipping <= self.options.max_clipping
      
    def check_mapped_length(self, segment):
        length = 0
        for cigar in segment.cigartuples:
            if cigar[0] == pysam.CMATCH or cigar[0] == pysam.CDEL:
                length += cigar[1]
        self.length_counter[length] += 1
        return length >= self.options.min_length and length <= self.options.max_length
                
    def check_seqname_pattern(self, segment):
        pattern = self.options.seqname_pattern
        seqname = segment.reference_name
        self.seqname_counter[seqname]
        return pattern is None or re.match(pattern, seqname) is not None
    
    def filter_single_end_bam(self, require_proper_paired=False):
        with pysam.AlignmentFile(self.infile) as f:
            header = utils.add_pg_to_header(f.header.as_dict(), cl=" ".join(sys.argv))
            with pysam.AlignmentFile(self.outfile, "wb", header=header) as fw:
                for segment in f.fetch(until_eof=True):
                    self.total += 1
                    if segment.is_unmapped:
                        if self.options.remove_unmapped:
                            self.unmapped += 1
                            continue
                    else:
                        if not self.check_seqname_pattern(segment):
                            self.exclude_seqname += 1
                            continue
                        if segment.is_secondary:
                            if self.options.remove_secondary:
                                self.secondary += 1
                                continue
                        elif segment.is_supplementary:
                            if self.options.remove_supplementary:
                                self.supplementary += 1
                                continue
                        else:
                            self.primary += 1
                            if require_proper_paired and not segment.is_proper_paired:
                                self.improper_paired += 1
                                continue
                            if not self.check_mapping_quality(segment):
                                self.outer_mapping_quality += 1
                                continue
                            if not self.check_clipping(segment):
                                self.too_many_clipping += 1
                                continue
                            if not self.check_mapped_length(segment):
                                self.outer_length += 1
                                continue               
                    self.final += 1
                    fw.write(segment)
                    
    def filter_paired_end_bam(self):
        if self.options.both_pass:
            raise NotImplementedError()
        self.filter_single_end_bam(require_proper_paired=True)
        
        
if __name__ == "__main__":
    FilterBam()