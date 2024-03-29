#!/usr/bin/env python
import os
import shutil
import optparse
import subprocess
import multiprocessing as mp
import pysam


REUSE = False

def _worker(bamfile, sizefile, chrom, outdir, 
            min_mapq=0, 
            keep_secondary=False, 
            keep_supplementary=False):
    f_bed = os.path.join(outdir, "%s.bed" % chrom)
    f_bg1 = os.path.join(outdir, "%s.bedgraph" % chrom)
    f_bg2 = os.path.join(outdir, "%s.sorted.bedgraph" % chrom)

    if not REUSE or not os.path.exists(f_bed):
        with pysam.AlignmentFile(bamfile) as f, open(f_bed, "w+") as fw:
            for s in f.fetch(chrom):
                if s.is_secondary:
                    if not keep_secondary:
                        continue
                elif s.is_supplementary:
                    if not keep_supplementary:
                        continue
                    if s.mapping_quality < min_mapq:
                        continue
                else:
                    if s.mapping_quality < min_mapq:
                        continue
                start = s.reference_start
                end = s.reference_end
                
                # name = s.query_name
                # score = s.mapping_quality
                # strand = "-" if s.is_reverse else "+"
                # items = [chrom, start, end, name, score, strand]
                
                items = [chrom, start, end]
                
                line = "\t".join(map(str, items))
                fw.write(line + "\n")

    cmd1 = "bedtools genomecov -bg -i %s -g %s > %s" % (f_bed, sizefile, f_bg1)
    if not REUSE or not os.path.exists(f_bg1):
        subprocess.check_call(cmd1, shell=True)

    cmd2 = "sort -k1,1 -k2,2n -k3,3n %s > %s" % (f_bg1, f_bg2)
    if not REUSE or not os.path.exists(f_bg2):
        subprocess.check_call(cmd2, shell=True)

    return f_bg2
    

def bam_to_bigwig(args=None):
    parser = optparse.OptionParser(usage="%prog [options] <input.bam> <output.bw>")
    parser.add_option("-q", "--mapq", dest="mapq", type="int", metavar="INT", default=20, 
                      help="[%default]")
    parser.add_option("--secondary", dest="keep_secondary", action="store_true")
    parser.add_option("--supplementary", dest="keep_supplementary", action="store_true")
    parser.add_option("-t", "--threads", dest="threads", type="int", metavar="INT",  default=1,
                      help="Threads. [%default]")
    options, args = parser.parse_args(args)
    bamfile, bwfile = args
    threads = options.threads
    mapq = options.mapq
    keep_sec = options.keep_secondary
    keep_sup = options.keep_supplementary

    tmpdir = bwfile + ".TMP"
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    sizefile = os.path.join(tmpdir, "genome.sizes")
    bgfile = os.path.join(tmpdir, "merged.bg")

    chroms = None
    with pysam.AlignmentFile(bamfile) as f, open(sizefile, "w+") as fw:
        chroms = list(sorted(f.references))
        for chrom in chroms:
            length = f.get_reference_length(chrom)
            fw.write("%s\t%s\n" % (chrom, length))
        
    pool = mp.Pool(threads)
    results = []
    for chrom in chroms:
        params = (bamfile, sizefile, chrom, tmpdir, mapq, keep_sec, keep_sup)
        results.append(pool.apply_async(_worker, params))
    pool.close()
    pool.join()

    with open(bgfile, "w+") as fw:
        for r in results:
            path = r.get()
            with open(path) as f:
                for line in f:
                    fw.write(line)


    cmd = "bedGraphToBigWig %s %s %s" % (bgfile, sizefile, bwfile)
    subprocess.check_call(cmd, shell=True)

    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    bam_to_bigwig()