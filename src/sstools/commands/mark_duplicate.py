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


def cluster_by_start(segments, diff=0, strand=False):
    
    if strand:
        
        segments_forward = list()
        segments_reverse = list()
        for s in segments:
            if s.is_reverse:
                segments_reverse.append(s)
            else:
                segments_forward.append(s)
                
        for segments in [segments_forward, segments_reverse]:
            array = None
            for s in sorted(segments, key=lambda item: item.reference_start):
                if array is None:
                    array = [s]
                elif s.reference_start - array[-1].reference_start <= diff:
                    array.append(s)
                else:
                    yield array
                    array = [s]
            if array is not None:
                yield array
            
    else:
        
        array = None
        for s in sorted(segments, key=lambda item: item.reference_start):
            if array is None:
                array = [s]
            elif s.reference_start - array[-1].reference_start <= diff:
                array.append(s)
            else:
                yield array
                array = [s]
        if array is not None:
            yield array
                

def cluster_by_end(segments, diff=0):
    array = None
    for segment in sorted(segments, key=lambda item: item.reference_end):
        if array is None:
            array = [segment]
        elif segment.reference_end - array[-1].reference_end <= diff:
            array.append(segment)
        else:
            yield array
            array = [segment]
    if array is not None:
        yield array
            
            
def set_duplicate_flags(segments):
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
            
              
def run_pipeline(inbam, outbam,                
                 diff=20,
                 strand=False,
                 f_read=None, 
                 f_group=None):
    
    n_total = 0
    n_uniq = 0
    index = 0 # group index
    
    h_read = None
    if f_read:
        h_read = open(f_read, "w+")
        h_read.write("Read\tGroup\tSize\tIndex\n")
        
    h_group = None
    if f_group:
        h_group = open(f_group, "w+")
        h_group.write("Group\tSize\tReads\n")
        
    with pysam.AlignmentFile(inbam) as f, pysam.AlignmentFile(outbam, "wb", f) as fw:
        for chrom in f.references:
            
            segments = []
            for segment in f.fetch(chrom):
                if segment.is_secondary:
                    continue
                if segment.is_supplementary:
                    continue
                segments.append(segment)
            n_total += len(segments)
            
            for segments1 in cluster_by_start(segments, diff, strand):
                for segments2 in cluster_by_end(segments1, diff):
                    n_uniq += 1
                    for segment in segments2:
                        segment.set_tag("DN", index)
                    set_duplicate_flags(segments2)
                    if h_group is not None:
                        read_names = ",".join([item.query_name for item in segments2])
                        line = "\t".join(map(str, [index, len(segments2), read_names]))
                        h_group.write(line + "\n")
                    index += 1
                    
            for segment in segments:
                fw.write(segment)
                if h_read is not None:
                    line = "\t".join(map(str, [
                        segment.query_name, 
                        segment.get_tag("DN"), 
                        segment.get_tag("DS"), 
                        segment.get_tag("DI")]))
                    h_read.write(line + "\n")
                    
    if h_group:
        h_group.close()
        
    print("Input: %s" % inbam)
    print("Output: %s" % outbam)
    print("Total: %d" % n_total)
    print("Uniq: %d (%.2f%%)" % (n_uniq, utils.divide_zero(n_uniq * 100, n_total)))
        
    
def mark_duplicate(args):
    
    parser = optparse.OptionParser(usage=usage)
    
    parser.add_option("-d", "--max-distance", dest="max_distance", type="int", default=20, metavar="INT",
                        help="Maximal distance between duplicates edge. [%default]")
    parser.add_option("-s", "--strand-sense", dest="strand_sense", action="store_true", default=False, 
                      help="")
    parser.add_option("-r", "--read-matrix", dest="read_matrix", metavar="PATH", 
                      help="")
    parser.add_option("-g", "--group-matrix", dest="group_matrix", metavar="PATH", 
                      help="")
    
    options, args = parser.parse_args(args)
    if len(args) != 2:
        parser.print_help()
        exit(1)
    inbam, outbam = args
    
    run_pipeline(inbam=inbam, 
                 outbam=outbam, 
                 diff=options.max_distance, 
                 strand=options.strand_sense, 
                 f_read=options.read_matrix, 
                 f_group=options.group_matrix)
