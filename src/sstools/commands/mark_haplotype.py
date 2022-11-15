#!/usr/bin/env python
import sys
import optparse
from collections import Counter, defaultdict
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import SegmentTools
from pyBioInfo.Utils import ShiftLoader
from sstools import utils

usage = """

    sstools MarkHaplotype [options] <input.bam> <output.bam>
"""


class MarkHaplotype(object):
    def __init__(self):
        self.options = None
        self.infile = None
        self.outfile = None
        self.init_parameters()
        self.counter = defaultdict(int)
        self.execute()
        total = sum(self.counter.values())
        print("Total: %d" % total)
        print("Paternal: %d (%.2f%%)" % (self.counter["P"], utils.divide_zero(self.counter["P"] * 100, total)))
        print("Maternal: %d (%.2f%%)" % (self.counter["M"], utils.divide_zero(self.counter["M"] * 100, total)))
        print("Other: %d (%.2f%%)" % (self.counter["O"], utils.divide_zero(self.counter["O"] * 100, total)))
        print("Ambigous: %d (%.2f%%)" % (self.counter["A"], utils.divide_zero(self.counter["A"] * 100, total)))
        print("Missing: %d (%.2f%%)" % (self.counter["-"], utils.divide_zero(self.counter["-"] * 100, total)))
        print("Unknown: %d (%.2f%%)" % (self.counter["U"], utils.divide_zero(self.counter["U"] * 100, total)))

    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-p", "--phased-file", dest="phased",
                          help="Region file in indexed BED/VCF format. [%default]")
        parser.add_option("-n", "--tag-name", dest="tag_name", default="HP",
                          help="Tag name of haplotype. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        assert options.phased.endswith(
            ".bed.gz") or options.phased.endswith(".vcf.gz")
        self.options = options
        self.infile = args[0]
        self.outfile = args[1]

    def load_phased_snps(self, f, chrom):
        if isinstance(f, pysam.VariantFile):
            name = list(f.header.samples)[0]
            try:
                for record in f.fetch(chrom):
                    start = record.pos - 1
                    ps = record.samples[name]["PS"]
                    if ps == ".":
                        continue
                    if ps != "PATMAT" or ps != "HOMVAR": # Only available for GIAB VCF
                        continue
                    gt = record.samples[name]["GT"]
                    allele1 = record.alleles[gt[0]]
                    allele2 = record.alleles[gt[1]]
                    if len(allele1) > 1 or len(allele2) > 1:
                        continue
                    if allele1 == allele2:
                        continue
                    snp = GRange(chrom=chrom, start=start, end=start + 1)
                    snp.base_p = allele1
                    snp.base_m = allele2
                    yield snp
            except ValueError:
                pass
        else:
            try:
                for line in f.fetch(chrom):
                    chrom, start, end, name = line.split("\t")[:4]
                    start, end = int(start), int(end)
                    base_p, base_m = name.split("|")
                    if base_p == base_m or len(base_p) > 1 or len(base_m) > 1:
                        continue
                    snp = GRange(chrom=chrom, start=start, end=end)
                    snp.base_p = base_p
                    snp.base_m = base_m
                    yield snp
            except ValueError:
                pass

    def get_parentals(self, segment, snps):
        parentals = []
        if len(snps) > 0:
            parsed_cigar = SegmentTools.parse_cigar(segment)
            sequence = segment.query_sequence
            for item in snps:
                start = item.start
                base = None
                for block in parsed_cigar:
                    if block[3][0] <= start < block[3][1]:
                        if block[0] == "M":
                            offset = start - block[3][0]
                            read_idx = block[2][0] + offset
                            base = sequence[read_idx]
                        elif block[0] == "D":
                            base = "-"
                        else:
                            assert False
                        break
                if base == "-":
                    parentals.append("-")
                elif base == item.base_p:
                    parentals.append("P")  # Paternal
                elif base == item.base_m:
                    parentals.append("M")  # Maternal
                else:
                    parentals.append("O")  # Other
        return parentals

    def determine_parental(self, parentals):
        parental = "U"  # Unknown
        if len(parentals) > 0:
            items = list(sorted(Counter(parentals).items(),
                         key=lambda item: item[1]))
            if len(items) > 0:
                parental = items[-1][0]
                if (len(parental) > 1) and (items[-2][1] == items[-1][1]):
                    parental = "A"  # Ambigous
            else:
                assert False
        return parental

    def execute(self):
        f1 = pysam.AlignmentFile(self.infile)
        header = utils.add_pg_to_header(
            f1.header.as_dict(), cl=" ".join(sys.argv))
        if self.options.phased.endswith(".vcf.gz"):
            f2 = pysam.VariantFile(self.options.phased)
        else:
            f2 = pysam.TabixFile(self.options.phased)
        fw = pysam.AlignmentFile(self.outfile, "wb", header=header)
        for chrom in f1.references:
            loader = ShiftLoader(self.load_phased_snps(f2, chrom))
            for segment in f1.fetch(chrom):
                start = segment.reference_start
                end = segment.reference_end
                snps = list(loader.fetch(chrom=chrom, start=start, end=end))
                parentals = self.get_parentals(segment, snps)
                parental = self.determine_parental(parentals)
                self.counter[parental] += 1
                segment.set_tag(self.options.tag_name, parental)
                fw.write(segment)
        f1.close()
        f2.close()
        fw.close()
