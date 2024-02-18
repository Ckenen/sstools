#!/usr/bin/env python#!/usr/bin/env python
import sys
import optparse
import multiprocessing as mp
import numpy as np
import pysam

usage = """

    sstools CallInversion [options] <composite.bed.gz> <outfile.bed>
"""

def _get_chromosomes(bedfile):
    with pysam.TabixFile(bedfile) as f:
        return f.contigs
    

def _get_counts(regions, peaks):
    counts = []
    idx = 0
    for i1, i2 in peaks:
        count = 0
        j = 0
        while True:
            k = idx + j
            if k >= len(regions):
                break
            item = regions[k]
            if item[1] <= i1:
                if j == 0:
                    idx += 1
                    continue
                else:
                    j += 1
                    continue
            if item[0] >= i2:
                break
            count += 1
            j += 1
        counts.append(count)
    return counts


def _worker(bedfile, chrom):
    regions = []
    max_end = 0
    with pysam.TabixFile(bedfile) as f:
        for line in f.fetch(chrom):
            row = line.strip("\n").split("\t")
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            max_end = max(max_end, end)
            strand = row[5]
            regions.append([start, end, strand])
    # max_end = max([x[1] for x in regions])
            
    regions1 = list(filter(lambda item: item[2] == "+", regions))
    regions2 = list(filter(lambda item: item[2] == "-", regions))

    covs1 = np.zeros(max_end + 1000, dtype=np.int)
    covs2 = np.zeros(max_end + 1000, dtype=np.int)
    for region in regions:
        if region[2] == "+":
            covs = covs1
        else:
            covs = covs2
        for i in range(region[0], region[1]):
            covs[i] += 1
            
    i1 = None
    tmp1 = []
    for i2, c in enumerate(covs2):
        if c > 0:
            if i1 is None:
                i1 = i2
        else:
            if i1 is not None:
                tmp1.append([i1, i2])
                i1 = None
    if i1 is not None:
        tmp1.append([i1, len(covs2)])

    crick = _get_counts(regions1, tmp1)
    watson = _get_counts(regions2, tmp1)
    
    peaks = []
    for peak, c, w in zip(tmp1, crick, watson):
        peaks.append([chrom, peak[0], peak[1], c, w])
    return peaks
    
    
def call_inversion(args=None):
    parser = optparse.OptionParser(usage="%prog [options] <composite.bed.gz> <outfile.bed>")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1,
                        help="Threads to run in parallel. [%default]")
    options, args = parser.parse_args(args)
    threads = options.threads
    compfile, outfile = args
            
    pool = mp.Pool(threads)
    results = []
    for chrom in _get_chromosomes(compfile):
        params = (compfile, chrom)
        results.append(pool.apply_async(_worker, params))
    pool.close()
    pool.join()
    results = [r.get() for r in results]

    with open(outfile, "w+") as fw:
        for peaks in results:
            for peak in peaks:
                chrom, start, end, crick, watson = peak
                if watson < 3:
                    continue
                ratio = watson / (crick + watson)
                name = "%d;%d;%.2f" % (crick, watson, ratio)
                if ratio < 0.1:
                    continue
                color = "0,0,255"
                if ratio >= 0.9:
                    color = "0,245,255"
                line = "\t".join(map(str, [chrom, start, end, name, \
                    ".", "+", start, end, color]))
                fw.write(line + "\n")
    

# class CallInversion(object):
#     def __init__(self):
#         self.infile = None
#         self.outfile = None
#         self.init_parameters()
#         self.execute()

#     def init_parameters(self):
#         parser = optparse.OptionParser(usage=usage)
#         parser.add_option("-p", "--processors", dest="threads", type="int", default=1,
#                           help="Processors to run in parallel. [%default]")
#         options, args = parser.parse_args(sys.argv[2:])
#         if len(args) != 2:
#             parser.print_help()
#             exit(1)
#         self.options = options
#         self.infile = args[0]
#         self.outfile = args[1]

#     @classmethod
#     def get_counts(cls, regions, peaks):
#         counts = []
#         idx = 0
#         for i1, i2 in peaks:
#             count = 0
#             j = 0
#             while True:
#                 k = idx + j
#                 if k >= len(regions):
#                     break
#                 item = regions[k]
#                 if item[1] <= i1:
#                     if j == 0:
#                         idx += 1
#                         continue
#                     else:
#                         j += 1
#                         continue
#                 if item[0] >= i2:
#                     break
#                 count += 1
#                 j += 1
#             counts.append(count)
#         return counts

#     @classmethod
#     def process_chromosome(cls, bedfile, chrom):
#         regions = []
#         max_end = 0
#         with pysam.TabixFile(bedfile) as f:
#             for line in f.fetch(chrom):
#                 row = line.strip("\n").split("\t")
#                 chrom = row[0]
#                 start = int(row[1])
#                 end = int(row[2])
#                 max_end = max(max_end, end)
#                 strand = row[5]
#                 regions.append([start, end, strand])
#         # max_end = max([x[1] for x in regions])
                
#         regions1 = list(filter(lambda item: item[2] == "+", regions))
#         regions2 = list(filter(lambda item: item[2] == "-", regions))

#         covs1 = np.zeros(max_end + 1000, dtype=np.int)
#         covs2 = np.zeros(max_end + 1000, dtype=np.int)
#         for region in regions:
#             if region[2] == "+":
#                 covs = covs1
#             else:
#                 covs = covs2
#             for i in range(region[0], region[1]):
#                 covs[i] += 1
                
#         i1 = None
#         tmp1 = []
#         for i2, c in enumerate(covs2):
#             if c > 0:
#                 if i1 is None:
#                     i1 = i2
#             else:
#                 if i1 is not None:
#                     tmp1.append([i1, i2])
#                     i1 = None
#         if i1 is not None:
#             tmp1.append([i1, len(covs2)])

#         crick = cls.get_counts(regions1, tmp1)
#         watson = cls.get_counts(regions2, tmp1)
        
#         peaks = []
#         for peak, c, w in zip(tmp1, crick, watson):
#             peaks.append([chrom, peak[0], peak[1], c, w])
#         return peaks
         

#     def execute(self):
#         with pysam.TabixFile(self.infile) as f:
#             chroms = f.contigs

#         pool = None
#         if self.options.threads > 1:
#             pool = multiprocessing.Pool(self.options.threads)
#         results = []
#         for chrom in chroms:
#             args = (self.infile, chrom)
#             if pool is None:
#                 r = self.process_chromosome(*args)
#             else:
#                 r = pool.apply_async(self.process_chromosome, args)
#             results.append(r)
#         if pool:
#             pool.close()
#             pool.join()
#             results = [r.get() for r in results]

#         with open(self.outfile, "w+") as fw:
#             for peaks in results:
#                 for peak in peaks:
#                     chrom, start, end, crick, watson = peak
#                     if watson < 3:
#                         continue
#                     r = watson / (crick + watson)
#                     name = "%d;%d;%.2f" % (crick, watson, r)
#                     if r < 0.1:
#                         continue
#                     color = "0,0,255"
#                     if r >= 0.9:
#                         color = "0,245,255"
#                     line = "\t".join(map(str, [chrom, start, end, name, \
#                         ".", "+", start, end, color]))
#                     fw.write(line + "\n")

            

