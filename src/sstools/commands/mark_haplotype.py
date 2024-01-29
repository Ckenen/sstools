#!/usr/bin/env python
import sys
import optparse
import json
from collections import Counter, defaultdict
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import SegmentTools
from pyBioInfo.Utils import ShiftLoader
from sstools import utils

usage = """

    sstools MarkHaplotype [options] <input.bam> <output.bam>
"""


def load_phased_snps(f, chrom):
    if isinstance(f, pysam.VariantFile):
        sample = list(f.header.samples)[0]
        try:
            for record in f.fetch(chrom):
                start = record.pos - 1
                ps = record.samples[sample]["PS"]
                if ps == ".":
                    continue
                if ps != "PATMAT" and ps != "0": # GIAB VCF
                    continue
                gt = record.samples[sample]["GT"]
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
    elif isinstance(f, pysam.TabixFile):
        try:
            for line in f.fetch(chrom):
                chrom, start, end, sample = line.strip("\n").split("\t")[:4]
                start, end = int(start), int(end)
                base_p, base_m = sample.split("|")
                if base_p == base_m or len(base_p) > 1 or len(base_m) > 1:
                    continue
                snp = GRange(chrom=chrom, start=start, end=end)
                snp.base_p = base_p
                snp.base_m = base_m
                yield snp
        except ValueError:
            pass
    else:
        raise RuntimeError()
    

def get_parentals(segment, snps):
    parentals = []
    parsed_cigar = SegmentTools.parse_cigar(segment)
    if len(snps) > 0:
        for snp in snps:
            base = SegmentTools.get_query_base(segment, snp.start, parsed_cigar)
            if base is None:
                continue
            if base == "-":
                parentals.append("-")
            elif base == snp.base_p:
                parentals.append("P")  # Paternal
            elif base == snp.base_m:
                parentals.append("M")  # Maternal
            else:
                parentals.append("O")  # Other
    return parentals


def determine_parental(parentals):
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


def run_pipeline(inbam, outbam, f_phased, f_matrix=None, tag_name="HP"):
    
    counter = defaultdict(int)
    
    h_matrix = None
    if f_matrix:
        h_matrix = open(f_matrix, "w+")
            
    with pysam.AlignmentFile(inbam) as f, pysam.AlignmentFile(outbam, "wb", f) as fw:
        if f_phased.endswith(".vcf.gz"):
            f2 = pysam.VariantFile(f_phased)
        elif f_phased.endswith(".bed.gz"):
            f2 = pysam.TabixFile(f_phased)
        else:
            raise RuntimeError()
        
        for chrom in f.references:
            loader = ShiftLoader(load_phased_snps(f2, chrom))
            for segment in f.fetch(chrom):
                start = segment.reference_start
                end = segment.reference_end
                snps = list(loader.fetch(chrom=chrom, start=start, end=end))
                parentals = get_parentals(segment, snps)
                parental = determine_parental(parentals)
                if h_matrix is not None:
                    detail = json.dumps(Counter(parentals))
                    h_matrix.write("%s\t%s\t%s\n" % (segment.query_name, parental, detail))
                counter[parental] += 1
                segment.set_tag(tag_name, parental)
                fw.write(segment)
                
        f2.close()
        
    if h_matrix:
        h_matrix.close()
        
    total = sum(counter.values())
    print("Total: %d" % total)
    print("Paternal: %d (%.2f%%)" % (counter["P"], utils.divide_zero(counter["P"] * 100, total)))
    print("Maternal: %d (%.2f%%)" % (counter["M"], utils.divide_zero(counter["M"] * 100, total)))
    print("Other: %d (%.2f%%)" % (counter["O"], utils.divide_zero(counter["O"] * 100, total)))
    print("Ambiguous: %d (%.2f%%)" % (counter["A"], utils.divide_zero(counter["A"] * 100, total)))
    print("Missing: %d (%.2f%%)" % (counter["-"], utils.divide_zero(counter["-"] * 100, total)))
    print("Unknown: %d (%.2f%%)" % (counter["U"], utils.divide_zero(counter["U"] * 100, total)))


def mark_haplotype(args):
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-p", "--phased-file", dest="phased", metavar="PATH",
                        help="PATH of phased SNPs file in BED/VCF format. [%default]")
    parser.add_option("-m", "--matrix-file", dest="matrix", metavar="PATH",
                        help="PATH of output matrix file. [%default]")
    parser.add_option("-n", "--tag-name", dest="tag_name", default="HP",
                        help="Tag name of haplotype. [%default]")
    options, args = parser.parse_args(args)
    if len(args) != 2:
        parser.print_help()
        exit(1)
        
        
    assert options.phased.endswith(".bed.gz") or options.phased.endswith(".vcf.gz")
    run_pipeline(
        inbam=args[0],
        outbam=args[1],
        f_phased=options.phased,
        f_matrix=options.matrix,
        tag_name=options.tag_name
    )
    