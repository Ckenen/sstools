#!/usr/bin/env python
import optparse
import re
from collections import defaultdict
import pysam
from pyBioInfo.Utils import SegmentTools

   
def pipeline_pe(inbam, outbam, pattern, min_mapq):
    n_total = 0
    n_unmapped = 0
    n_seqname = 0
    n_secondary = 0
    n_proper_pair = 0
    n_mapq = 0
    n_other = 0
    n_pass = 0
    with pysam.AlignmentFile(inbam) as f, pysam.AlignmentFile(outbam, "wb", f) as fw:
        data = defaultdict(list)
        for s in f.fetch(until_eof=True):
            n_total += 1
            if s.is_unmapped:
                n_unmapped += 1
                continue
            if pattern and re.match(pattern, s.referene_name) is None:
                n_seqname += 1
                continue
            if s.is_secondary:
                n_secondary += 1
                continue
            if not s.is_proper_pair:
                n_proper_pair += 1
                continue
            data[s.query_name].append(s)
            
        segments = []
        for name, ss in data.items():
            ss1 = [] # primary
            ss2 = [] # supplementary
            for s in ss:
                if s.is_supplementary:
                    ss2.append(s)
                else:
                    ss1.append(s)
            if len(ss1) == 2:
                r1, r2 = ss1
                if r1.is_read2:
                    r1, r2 = r2, r1
                assert r1.is_read1 and r2.is_read2
                q = min(r1.mapping_quality, r2.mapping_quality)
                if q < min_mapq:
                    n_mapq += len(ss)
                    continue
                n_pass += len(ss)
                segments.extend(ss)
            else:
                n_other += len(ss)
                continue
        for s in SegmentTools.sort_segments(segments):
            fw.write(s)
    
    print("Total:", n_total)
    print("Unmapped:", n_unmapped)
    print("InvalidSeqname:", n_seqname)
    print("SecondaryMapped:", n_secondary)
    print("UnProperPair:", n_proper_pair)
    print("LowMapQuality:", n_mapq)
    print("Other:", n_other)
    print("Pass:", n_pass)

def pipeline_se(inbam, outbam, pattern, mapq):
    n_total = 0
    n_unmapped = 0
    n_seqname = 0
    n_secondary = 0
    n_mapq = 0
    n_other = 0
    n_pass = 0
    
    with pysam.AlignmentFile(inbam) as f, pysam.AlignmentFile(outbam, "wb", f) as fw:
        data = defaultdict(list)
        for s in f.fetch(until_eof=True):
            n_total += 1
            if s.is_unmapped:
                n_unmapped += 1
                continue
            if pattern and re.match(pattern, s.reference_name) is None:
                n_seqname += 1
                continue
            if s.is_secondary:
                n_secondary += 1
                continue
            data[s.query_name].append(s)
            
        segments = []
        for name, ss in data.items():
            ss1 = [] # primary
            ss2 = [] # supplementary
            for s in ss:
                if s.is_supplementary:
                    ss2.append(s)
                else:
                    ss1.append(s)
            if len(ss1) == 1:
                s = ss1[0]
                q = s.mapping_quality
                if q < mapq:
                    n_mapq += len(ss)
                    continue
                n_pass += len(ss)
                segments.extend(ss)
            else:
                n_other += len(ss)
                continue
        for s in SegmentTools.sort_segments(segments):
            fw.write(s)

    print("Total:", n_total)
    print("Unmapped:", n_unmapped)
    print("InvalidSeqname:", n_seqname)
    print("SecondaryMapped:", n_secondary)
    print("LowMapQuality:", n_mapq)
    print("Other:", n_other)
    print("Pass:", n_pass)
    

def filter_bam(args=None):
    usage = "sstools FilterBam [options] <input.bam> <output.bam>"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-n", "--seqname-pattern", dest="pattern", default=None, metavar="STR",
                        help="Only output the chromosomes that match the regular expression. [%default]")
    parser.add_option("-q", "--min-mapq", dest="min_mapq", default=0, type="int", metavar="INT",
                        help="Minimal mapping quality. Only available for primary alignments. [%default]")
    parser.add_option("-p", "--paired", dest="is_paired", action="store_true", default=False, 
                        help="Only output the proper-paired alignments. [%default]")
    options, args = parser.parse_args(args)
    inbam, outbam = args
    is_paired = options.is_paired
    pattern = options.pattern
    mapq = options.min_mapq
    if is_paired:
        print("Running at paired-end mode.")
        pipeline_pe(inbam, outbam, pattern, mapq)
    else:
        print("Running at single-end mode.")
        pipeline_se(inbam, outbam, pattern, mapq)
        
        
if __name__ == "__main__":
    filter_bam()