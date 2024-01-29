#!/usr/bin/env python
import re
import optparse
import multiprocessing
import pysam


PATTERN = "^chr[0-9]+$"


def _worker(bamfile, chrom):
    with pysam.AlignmentFile(bamfile) as bam:
        length = bam.get_reference_length(chrom)
        regions = []
        for s in bam.fetch(chrom):
            start = s.reference_start
            end = s.reference_end
            regions.append([start, end])
        regions.sort()
        i = 0
        while i < len(regions) - 1:
            r1 = regions[i]
            r2 = regions[i + 1]
            if r1[1] >= r2[0]:
                r1[1] = max(r1[1], r2[1])
                regions.pop(i + 1)
            else:
                i += 1
        cov = 0
        for start, end in regions:
            cov += (end - start)
        return [chrom, length, cov]


def stat_coverage(args=None):
    parser = optparse.OptionParser(usage="%prog [options] input.bam output.tsv")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1, metavar="INT", 
                      help="Running threads. [%default]")
    options, args = parser.parse_args(args)
    bamfile, outfile = args
    threads = options.threads

    pool = multiprocessing.Pool(threads)
    results = []
    with pysam.AlignmentFile(bamfile) as bam:
        for chrom in bam.references:
            if re.match(PATTERN, chrom) is None:
                continue
            args = (bamfile, chrom)
            r = pool.apply_async(_worker, (bamfile, chrom))
            results.append(r)
    pool.close()
    pool.join()
    
    rows = []
    for r in results:
        chrom, length, coverage = r.get()
        ratio = coverage / length
        rows.append([chrom, length, coverage, ratio])
    total_length = sum([row[1] for row in rows])
    total_coverage = sum([row[2] for row in rows])
    total_ratio = total_coverage / total_length
    rows.append(["Overall", total_length, total_coverage, total_ratio])
    
    with open(outfile, "w+") as fw:
        fw.write("Name\tLength\tCoverage\tRatio\n")
        for row in rows:
            fw.write("\t".join(map(str, row)) + "\n")
        

if __name__ == "__main__":
    stat_coverage()