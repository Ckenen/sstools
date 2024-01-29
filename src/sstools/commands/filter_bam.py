#!/usr/bin/env python
import sys
import optparse
import re
from collections import defaultdict
import pysam
from pyBioInfo.Utils import SegmentTools

usage = """

    sstools FilterBam [options] <input.bam> <output.bam>
"""

    
def run_pipeline_for_paired_end(inbam, outbam, 
                                pattern=None, 
                                min_mapq=0):
    
    n_total = 0
    n_filtered = 0
    with pysam.AlignmentFile(inbam) as f, pysam.AlignmentFile(outbam, "wb", f) as fw:
        
        for chrom in f.references:
            
            if pattern and re.match("^chr([0-9]+|[XY])$", chrom) is None:
                continue
            
            ss = defaultdict(list)
            for s in f.fetch(chrom):
                n_total += 1
                if s.is_secondary:
                    n_secondary += 1
                    continue
                if s.is_supplementary:
                    n_supplementary += 1
                    continue
                if not s.is_proper_pair:
                    continue
                if s.mapping_quality < min_mapq:
                    n_invalid_mapping_quality += 1
                    continue
                ss[s.query_name].append(s)
    
            paired_ss = []
            for k, v in ss.items():
                if len(v) != 2:
                    continue
                r1, r2 = v
                if r1.is_read2:
                    r1, r2 = r2, r1
                assert r1.is_read1 and r2.is_read2
                paired_ss.append(r1)
                paired_ss.append(r2)
                n_filtered += 2
                
            for s in sorted(paired_ss, key=lambda s: s.reference_start):
                fw.write(s)
        

def run_pipeline_for_single_end(inbam, outbam, 
                            pattern=None,
                            min_mapq=0, 
                            max_mapq=1e8,
                            min_len=0,
                            max_len=1e8,
                            max_clip=1e8):
    
    n_total = 0
    n_unmapped = 0
    n_seqname_pattern = 0
    n_secondary = 0
    n_supplementary = 0
    n_invalid_mapq = 0
    n_invalid_len = 0
    n_invalid_clip = 0
    n_pass = 0
    
    with pysam.AlignmentFile(inbam) as f, pysam.AlignmentFile(outbam, "wb", f) as fw:
        for s in f.fetch(until_eof=True):
            n_total += 1
            if s.is_unmapped:
                n_unmapped += 1
                continue
            else:
                if pattern and re.match(pattern, s.reference_name) is None:
                    n_seqname_pattern += 1
                    continue
                if s.is_secondary:
                    n_secondary += 1
                    continue
                elif s.is_supplementary:
                    n_supplementary += 1
                    continue
                else:
                    q = s.mapping_quality
                    if q < min_mapq or q > max_mapq:
                        n_invalid_mapq += 1
                        continue
                    if max(SegmentTools.get_clipped(s)) > max_clip:
                        n_invalid_clip += 1
                        continue
                    mapped_len = SegmentTools.get_mapped_length(s)
                    if mapped_len < min_len or mapped_len > max_len:
                        n_invalid_len += 1
                        continue
            n_pass += 1    
            fw.write(s)


def filter_bam(args=None):
    usage = "sstools FilterBam [options] <input.bam> <output.bam>"
    parser = optparse.OptionParser(usage=usage)
        
    group = optparse.OptionGroup(parser, title="Mapping status")
    group.add_option("-n", "--seqname-pattern", dest="seqname_pattern", default=None, metavar="STR",
                        help="Only output the chromosomes that match the regular expression. [%default]")
    parser.add_option_group(group)
    
    group = optparse.OptionGroup(parser, title="Mapping quality")
    group.add_option("-q", "--min-mapq", dest="min_mapq", default=0, type="int", metavar="INT",
                        help="Minimal mapping quality. Only available for primary alignments. [%default]")
    group.add_option("-Q", "--max-mapq", dest="max_mapq", default=int(1e8), type="int",  metavar="INT",
                        help="Maximal mapping quality. Only available for primary alignments. [%default]")
    parser.add_option_group(group)
    
    group = optparse.OptionGroup(parser, title="Mapped region length")
    group.add_option("-m", "--min-length", dest="min_length", default=0, type="int",  metavar="INT",
                        help="Minimal mapped region length. Only available for primary alignments. Including deletion and excluding insertion and clipping. [%default]")
    group.add_option("-M", "--max-length", dest="max_length", default=int(1e8), type="int",  metavar="INT",
                        help="Maximal mapped region length. Only available for primary alignments. Including deletion and excluding insertion and clipping. [%default]")
    parser.add_option_group(group)
    
    group = optparse.OptionGroup(parser, title="Soft/hard clipping")
    group.add_option("-c", "--max-clip", dest="max_clip", default=int(1e8), type="int",  metavar="INT",
                        help="Maximal soft-clipping and hard-clipping. Only available for primary alignments. [%default]")
    parser.add_option_group(group)
    
    group = optparse.OptionGroup(parser, title="Paired-end")
    group.add_option("-p", "--paired", dest="is_paired", action="store_true", default=False, 
                        help="Only output the proper-paired alignments. [%default]")
    parser.add_option_group(group)
                
    options, args = parser.parse_args(args)
    print(options, args)
    exit(0)
    
    if len(args) != 2:
        parser.print_help()
        exit(1)

    if options.is_paired:
        run_pipeline_for_paired_end(
            inbam=args[0], 
            outbam=args[1], 
            pattern=options.seqname_pattern,
            min_mapq=options.min_mapq)
    else:
        run_pipeline_for_single_end(
            inbam=args[0], 
            outbam=args[1], 
            pattern=options.seqname_pattern, 
            min_mapq=options.min_mapq, 
            max_mapq=options.max_mapq, 
            min_len=options.min_length, 
            max_len=options.max_length, 
            max_clip=options.max_clip)
        
        
if __name__ == "__main__":
    filter_bam()