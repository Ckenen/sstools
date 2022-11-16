#!/usr/bin/env python
import sys
import optparse
import multiprocessing
import numpy as np
import pandas as pd
import pysam

usage = """

    sstools CalculateCWR [options] <input.bam> <output.tsv>
"""

class CalculateCWR(object):
    def __init__(self):
        self.infile = None # input.bam
        self.outfile = None # output.tsv
        self.options = None # 
        self.init_parameters()
        self.execute()

    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-p", "--processors", dest="threads", type="int", default=1,
                          help="Processors to run in parallel. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        self.options = options
        self.infile = args[0]
        self.outfile = args[1]

    @staticmethod
    def cal_cwr_by_chrom(infile, chrom):
        # 0: all
        # 1: high confidence
        # 2: not duplicates
        # 3: high confidence and not duplicates
        m = np.zeros((4, 2), dtype=np.int)
        with pysam.AlignmentFile(infile) as f:
            for segment in f.fetch(chrom):
                if segment.mapping_quality < 0:
                    continue
                strand = "-" if segment.is_reverse else "+"
                if segment.is_paired and segment.is_read2:
                    strand = "+" if strand == "-" else "-"
                    
                if strand == "+": # 0: all
                    m[0][0] += 1
                else:
                    m[0][1] += 1
                xh = segment.get_tag("XH") == "Y"
                dup = segment.is_duplicate
                if xh:
                    if strand == "+": # 1: high confidence
                        m[1][0] += 1
                    else:
                        m[1][1] += 1
                if not dup:
                    if strand == "+": # 2: not duplicate
                        m[2][0] += 1
                    else:
                        m[2][1] += 1
                    if xh:
                        if strand == "+": # 3: high confidence and not duplicate
                            m[3][0] += 1
                        else:
                            m[3][1] += 1
        return m

    def execute(self):
        array = []
        pool = multiprocessing.Pool(self.options.threads)
        with pysam.AlignmentFile(self.infile) as f:
            for chrom in f.references:
                res = pool.apply_async(self.cal_cwr_by_chrom, (self.infile, chrom))
                array.append([chrom, res])
        pool.close()
        pool.join()
        
        rows = []
        for chrom, res in array:
            assert res.successful()
            m = res.get()
            c1, w1 = m[0]
            c2, w2 = m[1]
            c3, w3 = m[2]
            c4, w4 = m[3]
            row = [chrom, c1, w1, c2, w2, c3, w3, c4, w4]
            rows.append(row)
            
        dat = pd.DataFrame(rows)
        dat.columns = ["Chrom", 
                    "Crick","Watson", 
                    "Crick.HC", "Watson.HC", 
                    "Crick.RmDup", "Watson.RmDup", 
                    "Crick.HC.RmDup", "Watson.HC.RmDup"]
        dat["Log2CWR"] = np.log2(dat["Crick"] / dat["Watson"])
        dat["Log2CWR.HC"] = np.log2(dat["Crick.HC"] / dat["Watson.HC"])
        dat["Log2CWR.RmDup"] = np.log2(dat["Crick.RmDup"] / dat["Watson.RmDup"])
        dat["Log2CWR.HC.RmDup"] = np.log2(dat["Crick.HC.RmDup"] / dat["Watson.HC.RmDup"])
        dat.to_csv(self.outfile, sep="\t", index=False)

