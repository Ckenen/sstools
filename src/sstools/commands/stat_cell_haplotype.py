#!/usr/bin/env python
import sys
import os
import optparse
from collections import defaultdict
import multiprocessing
import numpy as np
import pandas as pd
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import ShiftLoader

usage = """

    sstools StatCellHap [options] <input.bam> <outdir>
"""

class StatCellHaplotype(object):
    def __init__(self):
        self.infile = None
        self.outdir = None
        self.options = None
        self.init_patameters()
        self.execute()

    def init_patameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-p", "--processors", dest="threads", type="int", default=1,
                          help="Processors to run in parallel. [%default]")
        parser.add_option("-g", "--genome", dest="genome", default=None,
                          help="Genome sequences in FASTA format. [%default]")
        parser.add_option("-a", "--phased", dest="phased", default=None,
                          help="Phased alleles in BED/VCF format. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        self.options = options
        self.infile = args[0]
        self.outdir = args[1]
    
    @staticmethod
    def load_snvs(path, chrom):
        assert path.endswith(".bed.gz")
        with pysam.TabixFile(path) as f:
            try:
                for line in f.fetch(chrom):
                    line = line.strip("\n")
                    if line.startswith("#"):
                        continue
                    row = line.split("\t")
                    chrom, start, end, name = row
                    if "*" in name:
                        continue
                    start, end = int(start), int(end)
                    hp1_base, hp2_base = name.split("|")
                    obj = GRange(chrom=chrom, start=start, end=end, name=name)
                    obj.hp1_base = hp1_base
                    obj.hp2_base = hp2_base
                    yield obj
            except ValueError:
                pass

    @staticmethod
    def assign_ref_base(snvs, fasta):
        for snv in snvs:
            ref = fasta.fetch(snv.chrom, snv.start, snv.end)
            snv.ref = ref

    @staticmethod
    def load_alignments(path, chrom):
        with pysam.AlignmentFile(path) as f:
            for segment in f.fetch(chrom):
                if segment.is_unmapped:
                    continue
                start = segment.reference_start
                end = segment.reference_end
                name = segment.query_name
                strand = "-" if segment.is_reverse else "+"
                me = segment.get_tag("ME")
                events = dict()
                if len(me) > 0:
                    for s in me.split(";"):
                        items = s.split(",")
                        pos = int(items[0])
                        ref = items[1]
                        alt = items[2]
                        if ref == "-": # insertion
                            continue
                        elif alt == "-": # deletion
                            qua = int(items[3])
                            for i, base in enumerate(ref):
                                events[pos + i] = [base, alt, qua]
                        else: # mismatch
                            qua = int(items[3])
                            events[pos] = [ref, alt, qua]
                obj = GRange(chrom=chrom, start=start, end=end, 
                            name=name, strand=strand)
                obj.events = events
                yield obj

    @staticmethod
    def get_base_counter(snv, alignments):
        counter = defaultdict(int)
        for obj in alignments:
            event = obj.events.get(snv.start)
            if event is None:
                continue
            counter[event[1]] += 1
        ref_count = len(alignments) - sum(counter.values())
        assert ref_count >= 0
        if ref_count > 0:
            counter[snv.ref] = ref_count
        return counter

    @staticmethod
    def counter_to_str(counter):
        s = []
        for k, v in sorted(counter.items()):
            if v > 0:
                s.append("%s:%d" % (k, v))
        s = ",".join(s)
        return s

    @staticmethod
    def determine_base(counter):
        total = sum(counter.values())
        if total > 0:
            items = list(sorted(counter.items(), key=lambda item: item[1]))
            k, v = items[-1]
            if v >= total * 0.75:
                return k
        return ""

    @staticmethod
    def process_one_chrom(kwargs):
        bam = kwargs["bam"]
        phased = kwargs["phased"]
        fasta = kwargs["fasta"]
        outdir = kwargs["outdir"]
        chrom = kwargs["chrom"]
        
        if not os.path.exists(outdir + "/%s" % chrom):
            os.mkdir(outdir + "/%s" % chrom)
            
        fw = open(outdir + "/%s/details.tsv" % chrom, "w+")
        fw.write("Chrom\tPosition\tRefBase\tPatBase\tMatBase\tCrickNum\tWatsonNum\tCrickBase\tWatsonBase\tCrickDetail\tWatsonDetail\n")
        
        snvs = list(StatCellHaplotype.load_snvs(phased, chrom))
        with pysam.FastaFile(fasta) as f:
            StatCellHaplotype.assign_ref_base(snvs, f)
        
        summary = dict()
        for read_number in range(1, 11):
            summary[read_number] = np.zeros((2, 5), dtype=np.int)
                    
        alignments = StatCellHaplotype.load_alignments(bam, chrom)
        loader = ShiftLoader(alignments)
        for snv in sorted(snvs):
            reads1 = []
            reads2 = []
            for read in loader.fetch(obj=snv):
                if read.strand == "+":
                    reads1.append(read)
                else:
                    reads2.append(read)
            count1 = len(reads1)
            count2 = len(reads2)
            if count1 == 0 and count2 == 0:
                continue
            counter1 = StatCellHaplotype.get_base_counter(snv, reads1)
            counter2 = StatCellHaplotype.get_base_counter(snv, reads2)
            s1 = StatCellHaplotype.counter_to_str(counter1)
            s2 = StatCellHaplotype.counter_to_str(counter2)
            base1 = StatCellHaplotype.determine_base(counter1)
            base2 = StatCellHaplotype.determine_base(counter2)
            row = [snv.chrom, snv.start, snv.ref, snv.hp1_base, snv.hp2_base, 
                    count1, count2, base1, base2, s1, s2]
            line = "\t".join(map(str, row))
            fw.write(line + "\n")
            
            if count1 in summary.keys():
                matrix = summary[count1]
                if base1 == snv.hp1_base:
                    matrix[0][0] += 1
                elif base1 == snv.hp2_base:
                    matrix[0][1] += 1
                elif base1 == "-":
                    matrix[0][2] += 1
                elif base1 == "": # Ambiguity
                    matrix[0][4] += 1
                else:
                    matrix[0][3] += 1
                    
            if count2 in summary.keys():
                matrix = summary[count2]
                if base2 == snv.hp1_base:
                    matrix[1][0] += 1
                elif base2 == snv.hp2_base:
                    matrix[1][1] += 1
                elif base2 == "-":
                    matrix[1][2] += 1
                elif base2 == "": # Ambiguity
                    matrix[1][4] += 1
                else:
                    matrix[1][3] += 1
                    
        fw.close()
                    
        fw = open(outdir + "/%s/summary.log" % chrom, "w+")
        for k in summary.keys():
            matrix = summary[k]
            v1 = matrix[0][0] + matrix[1][1]
            v2 = matrix[0][1] + matrix[1][0]
            if v1 + v2 > 0:
                precision = max(v1, v2) / (v1 + v2)
            else:
                precision = 0
            fw.write("=" * 40 + "\n")
            fw.write("Read number: %d\n" % k)
            fw.write("%s\n" % str(matrix))
            fw.write("Precision: %f\n" % precision)
            dat = pd.DataFrame(matrix)
            dat.columns = ["Paternal", "Maternal", "Deletion", "OtherBase", "Ambiguity"]
            dat.index = ["Crick", "Watson"]
            dat.to_csv(outdir + "/%s/summary.%d.tsv" % (chrom, k), sep="\t")
        fw.close()

    def execute(self):       
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        
        pool = None
        if self.options.threads > 1:
            pool = multiprocessing.Pool(self.options.threads)
        results = []
        with pysam.AlignmentFile(self.infile) as f:
            for chrom in f.references:
                kwargs = {"bam": self.infile, 
                        "phased": self.options.phased,
                        "fasta": self.options.genome,
                        "chrom": chrom,
                        "outdir": self.outdir}
                if pool is None:
                    self.process_one_chrom(kwargs)
                else:
                    res = pool.apply_async(self.process_one_chrom, (kwargs,))
                    results.append(res)
        if pool is not None:
            pool.close()
            pool.join()
            for res in results:
                assert res.successful()

