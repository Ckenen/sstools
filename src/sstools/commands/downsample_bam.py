#!/usr/bin/env python
import sys
import optparse
import random
import pysam


def downsample_bam(args):
    parser = optparse.OptionParser(usage="%prog [options] input.bam downsample.bam")
    parser.add_option("-p", "--proportion", dest="proportion", type="float", metavar="FLOAT", 
                      help="Proportion of downsample. [%default]")
    parser.add_option("-s", "--seed", dest="seed", type="int", default=0, metavar="INT", 
                      help="Random seed. [%default]")
    
    options, args = parser.parse_args(args)
    bamfile, outfile = args
    seed = options.seed
    proportion = options.proportion

    random.seed(seed)

    with pysam.AlignmentFile(bamfile) as f, pysam.AlignmentFile(outfile, "wb", f) as fw:
        for segment in f:
            v = random.random()
            if v < proportion:
                fw.write(segment)
