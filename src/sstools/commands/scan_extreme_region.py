#!/usr/bin/env python
import optparse
import re
import multiprocessing as mp
import numpy as np
import pyBigWig


def _worker(bwfile, chrom, width):
    width2 = width * 10000
    with pyBigWig.open(bwfile) as f:
        length = f.chroms()[chrom]
        nbin = int(length / width) + (1 if length % width > 0 else 0)
        vs = np.zeros(nbin)
        i = 0
        for i1 in range(0, length, width2):
            i2 = min(i1 + width2, length)
            covs = f.values(chrom, i1, i2)
            covs = np.nan_to_num(covs, 0)
            for j1 in range(i1, i2, width):
                j2 = min(j1 + width, length)
                vs[i] = np.mean(covs[j1-i1:j2-i1])
                i += 1
    return vs


def stat_coverage(bwfile, width, threads):
    pool = mp.Pool(threads)
    results = []
    with pyBigWig.open(bwfile) as f:
        chroms = f.chroms()
        for chrom in f.chroms():
            results.append(pool.apply_async(_worker, (bwfile, chrom, width)))
    pool.close()
    pool.join()
    data = dict()
    for chrom, r in zip(chroms, results):
        data[chrom] = r.get()
    return chroms, data


def get_regions(chroms, data, width, p1, p2):
    pattern1 = "^chr[0-9]+$"
    pattern2 = "^chrX$"
    pattern3 = "^chrY$"
    
    regions = []
    for pattern in [pattern1, pattern2, pattern3]:
        # cutoff
        vs = []
        for chrom, covs in data.items():
            if re.match(pattern, chrom):
                vs.extend(data[chrom])
        vs = np.array(vs)
        vs1 = vs[vs>0]
        vs2 = np.sort(vs1)
        if p1 == 0:
            vmin = 0
        else:
            vs2[int(len(vs2) * p1)]
        vmax = vs2[min(int(len(vs2) * p2), len(vs2) - 1)]
        
        # extreme region
        for chrom, covs in data.items():
            tmp = []
            if re.match(pattern, chrom):
                length = chroms[chrom]
                for i, v in enumerate(covs):
                    if v < vmin or v > vmax:
                        start = i * width
                        end = min(start + width, length)
                        r = [chrom, start, end]
                        if len(tmp) == 0:
                            tmp.append(r)
                        elif start <= tmp[-1][2]:
                            tmp[-1][2] = end
                        else:
                            tmp.append(r)
            regions.extend(tmp)
    return regions


def scan_extreme_region(args=None):
    parser = optparse.OptionParser(usage="%prog [options] <input.bw> <output.bed>")
    parser.add_option("-p", "--proportion", dest="proportion", default="0.01,0.99",help="[%default]")
    parser.add_option("-w", "--width", dest="width", type="int", default=500, help="[%default]")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1, help="[%default]")
    options, args = parser.parse_args(args)
    bwfile, outfile = args
    width = options.width
    threads = options.threads
    p1, p2 = options.proportion.split(",")
    p1, p2 = float(p1), float(p2)
    
    chroms, data = stat_coverage(bwfile, width, threads)
    regions = get_regions(chroms, data, width, p1, p2)

    with open(outfile, "w+") as fw:
        for r in sorted(regions):
            fw.write("\t".join(map(str, r)) + "\n")
            

if __name__ == "__main__":
    scan_extreme_region()