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
                a1 = record.alleles[gt[0]]
                a2 = record.alleles[gt[1]]
                if len(a1) > 1 or len(a2) > 1:
                    continue
                if a1 == a2:
                    continue
                snp = GRange(chrom=chrom, start=start, end=start + 1)
                snp.allele1 = a1
                snp.allele2 = a2
                yield snp
        except ValueError:
            pass
    elif isinstance(f, pysam.TabixFile):
        try:
            for line in f.fetch(chrom):
                chrom, start, end, sample = line.strip("\n").split("\t")[:4]
                start, end = int(start), int(end)
                a1, a2 = sample.split("|")
                if a1 == a2 or len(a1) > 1 or len(a2) > 1:
                    continue
                snp = GRange(chrom=chrom, start=start, end=end)
                snp.allele1 = a1
                snp.allele2 = a2
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
            elif base == snp.allele1:
                parentals.append("P")  # Paternal
            elif base == snp.allele2:
                parentals.append("M")  # Maternal
            else:
                parentals.append("O")  # Other
    return parentals


def determine_parental(parentals):
    parental = "U"  # Unknown
    if len(parentals) > 0:
        items = list(sorted(Counter(parentals).items(), key=lambda item: item[1]))
        if len(items) > 0:
            parental = items[-1][0]
            if (len(parental) > 1) and (items[-2][1] == items[-1][1]):
                parental = "A"  # Ambigous
        else:
            assert False
    return parental


def run_pipeline(inbam, outbam, phased, detail, tag_name):
    counter = defaultdict(int)
    
    h_detail = None
    if detail:
        h_detail = open(detail, "w+")
            
    with pysam.AlignmentFile(inbam) as f, pysam.AlignmentFile(outbam, "wb", f) as fw:
        if phased.endswith(".vcf.gz"):
            f2 = pysam.VariantFile(phased)
        elif phased.endswith(".bed.gz"):
            f2 = pysam.TabixFile(phased)
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
                if h_detail is not None:
                    detail = json.dumps(Counter(parentals))
                    h_detail.write("%s\t%s\t%s\n" % (segment.query_name, parental, detail))
                counter[parental] += 1
                segment.set_tag(tag_name, parental)
                fw.write(segment)   
        f2.close()
    if h_detail:
        h_detail.close()
        
    total = sum(counter.values())
    print("Total: %d" % total)
    print("Paternal: %d (%.2f%%)" % (counter["P"], utils.divide_zero(counter["P"] * 100, total)))
    print("Maternal: %d (%.2f%%)" % (counter["M"], utils.divide_zero(counter["M"] * 100, total)))
    print("Other: %d (%.2f%%)" % (counter["O"], utils.divide_zero(counter["O"] * 100, total)))
    print("Ambiguous: %d (%.2f%%)" % (counter["A"], utils.divide_zero(counter["A"] * 100, total)))
    print("Missing: %d (%.2f%%)" % (counter["-"], utils.divide_zero(counter["-"] * 100, total)))
    print("Unknown: %d (%.2f%%)" % (counter["U"], utils.divide_zero(counter["U"] * 100, total)))


usage = """sstools MarkHaplotype [options] <input.bam> <phased.bed.gz|phased.vcf.gz> <output.bam>

<phased.bed.gz|phased.vcf.gz>: PATH of phased SNPs file in BED/VCF format.
For BED format, the name of record is expected to 'P|M' (eg. 'A|G').
For VCF format, the phase set (PS) is expected to 'PATMAT' or '0'.
"""

def mark_haplotype(args):
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-d", "--detail-file", dest="detail", metavar="PATH",
                        help="PATH of output detail file. [%default]")
    parser.add_option("-n", "--tag-name", dest="tag_name", default="HP",
                        help="Tag name of haplotype. [%default]")
    options, args = parser.parse_args(args)
    inbam, phased, outbam = args
    detail = options.detail
    tag_name = options.tag_name
    assert phased.endswith(".bed.gz") or phased.endswith(".vcf.gz")
    run_pipeline(inbam, outbam, phased, detail, tag_name)
    
    
if __name__ == "__main__":
    mark_haplotype()
    