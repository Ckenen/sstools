#!/usr/bin/env python
import sys
import optparse
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import ShiftLoader
from sstools import utils

usage = """

    sstools MarkRegion [options] <input.bam> <output.bam>
    
This command will add a TAG to the BAM file.
"""


class MarkRegion(object):
    def __init__(self):
        self.options = None
        self.infile = None
        self.outfile = None
        self.init_parameters()
        
        self.total = 0
        self.hit = 0
        self.not_hit = 0
        
        self.fh = None
        self.fh = open(self.outfile + ".hits.txt", "w+")
        self.fh.write("ReadName\tHit\n")
        
        self.execute()
        
        if self.fh is not None:
            self.fh.close()
                    
        print("Infile: %s" % self.infile)
        print("Region: %s" % self.options.region)
        print("Outfile: %s" % self.outfile)
        print("Total: %d" % self.total)
        print("Hit: %d (%.2f%%)" % (self.hit, utils.divide_zero(self.hit * 100, self.total)))
        print("NotHit: %d (%.2f%%)" % (self.not_hit, utils.divide_zero(self.not_hit * 100, self.total)))
        
        
    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-f", "--region-file", dest="region", 
                          help="Region file in BED format. [%default]")
        parser.add_option("-n", "--tag-name", dest="tag_name", default="XR", 
                          help="Tag name record hit region or not. [%default]")
        parser.add_option("-r", "--remove-hit", dest="remove_hit", action="store_true", default=False, 
                          help="Remove alignments that hit regions. [%default]")
        parser.add_option("-R", "--remove-not-hit", dest="remove_not_hit", action="store_true", default=False, 
                          help="Remove alignments that not hit any region. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        self.options = options
        self.infile = args[0]
        self.outfile = args[1]
    
    def load_regions(self, f, chrom):
        try:
            for line in f.fetch(chrom):
                chrom, start, end = line.split("\t")[:3]
                start, end = int(start), int(end)
                yield GRange(chrom=chrom, start=start, end=end)
        except ValueError:
            pass
    
    def execute(self):
        f1 = pysam.AlignmentFile(self.infile)
        f2 = pysam.TabixFile(self.options.region)
        header = utils.add_pg_to_header(f1.header.as_dict(), cl=" ".join(sys.argv))
        fw = pysam.AlignmentFile(self.outfile, "wb", header=header)
        for chrom in f1.references:
            regions = self.load_regions(f2, chrom)
            loader = ShiftLoader(regions)
            for segment in f1.fetch(chrom):
                self.total += 1
                start, end = segment.reference_start, segment.reference_end
                is_hit = len(list(loader.fetch(chrom=chrom, start=start, end=end))) > 0
                if is_hit:
                    self.hit += 1
                else:
                    self.not_hit += 1
                s = "Y" if is_hit else "N"
                if self.fh is not None:
                    self.fh.write("%s\t%s\n" % (segment.query_name, s))
                segment.set_tag(self.options.tag_name, s)
                if self.options.remove_hit and is_hit:
                    continue
                if self.options.remove_not_hit and not is_hit:
                    continue
                fw.write(segment)
        f1.close()
        f2.close()
        fw.close()