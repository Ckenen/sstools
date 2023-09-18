#!/usr/bin/env python
import optparse
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import ShiftLoader
from sstools import utils

usage = """

    sstools MarkRegion [options] <input.bam> <output.bam>
    
This command will add a TAG to the BAM file.
"""

def load_regions(f, chrom):
    try:
        for line in f.fetch(chrom):
            chrom, start, end = line.split("\t")[:3]
            start, end = int(start), int(end)
            yield GRange(chrom=chrom, start=start, end=end)
    except ValueError:
        pass


def run_pipeline(inbam, outbam, inbed, f_matrix=None, 
                 tag_name="XR", rm_hit=False, rm_not_hit=False):
    
    n_total = 0
    n_hit = 0
    n_not_hit = 0
    
    h_matrix = None
    if f_matrix:
        h_matrix = open(f_matrix, "w+")
        h_matrix.write("Read\tHit\n")
    
    with pysam.AlignmentFile(inbam) as f, \
        pysam.TabixFile(inbed) as f2, \
        pysam.AlignmentFile(outbam, "wb", f) as fw:
            
        for chrom in f.references:
            regions = load_regions(f2, chrom)
            loader = ShiftLoader(regions)
            for s in f.fetch(chrom):
                n_total += 1
                start, end = s.reference_start, s.reference_end
                is_hit = len(list(loader.fetch(chrom=chrom, start=start, end=end))) > 0
                if is_hit:
                    n_hit += 1
                else:
                    n_not_hit += 1
                v = "Y" if is_hit else "N"
                if h_matrix:
                    h_matrix.write("%s\t%s\n" % (s.query_name, v))
                s.set_tag(tag_name, v)
                if rm_hit and is_hit:
                    continue
                if rm_not_hit and not is_hit:
                    continue
                fw.write(s)
                
    if h_matrix:
        h_matrix.close()
        
    print("Infile: %s" % inbam)
    print("Region: %s" % inbed)
    print("Outfile: %s" % outbam)
    print("Total: %d" % n_total)
    print("Hit: %d (%.2f%%)" % (n_hit, utils.divide_zero(n_hit * 100, n_total)))
    print("NotHit: %d (%.2f%%)" % (n_not_hit, utils.divide_zero(n_not_hit * 100, n_total)))


def mark_region(args):
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-f", "--region-file", dest="region", metavar="PATH",
                        help="PATH of regions in BED format. [%default]")
    parser.add_option("-m", "--matrix-file", dest="matrix", metavar="PATH",
                        help="PATH of output matrix file. [%default]")
    parser.add_option("-n", "--tag-name", dest="tag_name", default="XR", metavar="STR",
                        help="Tag name record hit region or not. [%default]")
    parser.add_option("-r", "--remove-hit", dest="remove_hit", action="store_true", default=False, 
                        help="Remove alignments that hit regions. [%default]")
    parser.add_option("-R", "--remove-not-hit", dest="remove_not_hit", action="store_true", default=False, 
                        help="Remove alignments that not hit any region. [%default]")
    options, args = parser.parse_args(args)
    if len(args) != 2:
        parser.print_help()
        exit(1)
    
    run_pipeline(
        inbam=args[0], 
        outbam=args[1], 
        inbed=options.region, 
        f_matrix=options.matrix, 
        tag_name=options.tag_name,
        rm_hit=options.remove_hit, 
        rm_not_hit=options.remove_not_hit)
    