#!/usr/bin/env python
# import sys
import re
import optparse
import multiprocessing as mp
import json
from collections import defaultdict
import numpy as np
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import ShiftLoader


def load_snps(path, chrom=None, het_only=False, hom_only=None):
    
    snps = []
    
    with pysam.VariantFile(path) as f:
        
        sample = f.header.samples[0]
        
        try:
            
            for record in f.fetch(chrom):
                
                gt = record.samples[sample]["GT"]
                try:
                    ps = record.samples[sample]["PS"]
                except Exception:
                    ps = None
                a1 = record.alleles[gt[0]]
                a2 = record.alleles[gt[1]]
                
                if len(a1) > 1 or len(a2) > 1:
                    continue
                
                if het_only and a1 == a1:
                    continue
                
                if hom_only and a1 != a2:
                    continue
                
                snp = GRange(chrom=record.chrom, start=record.start, end=record.stop, name="SNP")
                snp.allele1 = a1
                snp.allele2 = a2
                snp.ps = ps
                
                if len(snps) == 0:
                    snps.append(snp)
                elif snp.start > snps[-1].start:
                    snps.append(snp)
                else:
                    pass
                    # print(chrom, snps[-1].start, snp.start)
                
        except Exception:
            pass
        
    snps.sort()
    
    return snps


def load_regions(path, chrom):
    regions = []
    with pysam.TabixFile(path) as f:
        try:
            for line in f.fetch(chrom):
                row = line.strip("\t").split("\t")
                chrom, start, end = row[:3]
                obj = GRange(chrom=chrom, start=int(start), end=int(end))
                regions.append(obj)
        except Exception:
            pass
    return regions


def filter_snps(snps, regions):
    loader = ShiftLoader(regions)
    snps1 = []
    for het in snps:
        if len(list(loader.fetch(obj=het))) > 0:
            snps1.append(het)
    return snps1


def snp_to_set(snps, het_only=False, hom_only=False):
    array = []
    for snp in snps:
        a1, a2 = snp.allele1, snp.allele2
        if a1 > a2:
            a1, a2 = a2, a1
        if het_only and a1 == a2:
            continue
        if hom_only and a1 != a2:
            continue
        array.append((snp.chrom, snp.start, a1, a2))
    return set(array)


def worker(chrom, f_reference, f_query, f_bed=None, 
           het_only=False, hom_only=False, 
           reference_ps="PATMAT", query_ps="PATMAT"):
    data = dict()
    data["Chrom"] = chrom
    data["Reference_VCF"] = f_reference
    data["Query_VCF"] = f_query
    data["Region_BED"] = f_bed
    
    snps1 = load_snps(f_reference, chrom, het_only, hom_only)
    snps2 = load_snps(f_query, chrom, het_only, hom_only)
    data["Reference_SNPs"] = len(snps1)
    data["Query_SNPs"] = len(snps2)
    
    if f_bed:
        regions =load_regions(f_bed, chrom)
        snps1 = filter_snps(snps1, regions)
        snps2 = filter_snps(snps2, regions)
    data["Filtered_Reference_SNPs"] = len(snps1)
    data["Filtered_Query_SNPs"] = len(snps2)
        
    # All SNPs
    set1 = snp_to_set(snps1)
    set2 = snp_to_set(snps2)
    set3 = set1 & set2
    n1, n2, n3 = len(set1), len(set2), len(set3)
    recall = np.divide(n3, n1)
    precision = np.divide(n3, n2)
    f1 = np.divide(2 * recall * precision, (recall + precision))
    data["All_SNP_Reference"] = n1
    data["All_SNP_Query"] = n2
    data["All_SNP_Overlap"] = n3
    data["All_SNP_Recall"] = recall
    data["All_SNP_Precision"] = precision
    data["All_SNP_F1"] = f1
    
    # Het SNPs
    set1 = snp_to_set(snps1, het_only=True)
    set2 = snp_to_set(snps2, het_only=True)
    set3 = set1 & set2
    n1, n2, n3 = len(set1), len(set2), len(set3)
    recall = np.divide(n3, n1)
    precision = np.divide(n3, n2)
    f1 = np.divide(2 * recall * precision, (recall + precision))
    data["Het_SNP_Reference"] = n1
    data["Het_SNP_Query"] = n2
    data["Het_SNP_Overlap"] = n3
    data["Het_SNP_Recall"] = recall
    data["Het_SNP_Precision"] = precision
    data["Het_SNP_F1"] = f1
    
    # Hom SNPs
    set1 = snp_to_set(snps1, hom_only=True)
    set2 = snp_to_set(snps2, hom_only=True)
    set3 = set1 & set2
    n1, n2, n3 = len(set1), len(set2), len(set3)
    recall = np.divide(n3, n1)
    precision = np.divide(n3, n2)
    f1 = np.divide(2 * recall * precision, (recall + precision))
    data["Hom_SNP_Reference"] = n1
    data["Hom_SNP_Query"] = n2
    data["Hom_SNP_Overlap"] = n3
    data["Hom_SNP_Recall"] = recall
    data["Hom_SNP_Precision"] = precision
    data["Hom_SNP_F1"] = f1
        
    # Calling, genotypeing and phasing
    
    set1 = set([(snp.chrom, snp.start) for snp in snps1])
    set2 = set([(snp.chrom, snp.start) for snp in snps2])
    set3 = set1 & set2
    n1, n2, n3 = len(set1), len(set2), len(set3)
    recall = np.divide(n3, n1)
    precision = np.divide(n3, n2)
    f1 = np.divide(2 * recall * precision, (recall + precision))
    data["Calling_Reference"] = n1
    data["Calling_Query"] = n2
    data["Calling_Overlap"] = n3
    data["Calling_Recall"] = recall
    data["Calling_Precision"] = precision
    data["Calling_F1"] = f1
    
    snps1 = list(filter(lambda snp: (snp.chrom, snp.start) in set3, snps1))
    snps2 = list(filter(lambda snp: (snp.chrom, snp.start) in set3, snps2))
    counter = defaultdict(int)
    n1, n2 = 0, 0
    for snp1, snp2 in zip(snps1, snps2):
        assert snp1.chrom == snp2.chrom and snp1.start == snp2.start
        if snp1.allele1 == snp1.allele2:
            if snp2.allele1 == snp2.allele2:
                if snp1.allele1 == snp2.allele1:
                    counter["HOM-HOM"] += 1
                else:
                    counter["HOM-HOM.2"] += 1
            else:
                counter["HOM-HET"] += 1
        else:
            if snp2.allele1 == snp2.allele2:
                counter["HET-HOM"] += 1
            else:
                if snp1.allele1 == snp2.allele1 and snp1.allele2 == snp2.allele2:
                    counter["HET-HET"] += 1
                    if snp1.ps == reference_ps and snp2.ps == query_ps:
                        n1 += 1
                elif snp1.allele1 == snp2.allele2 and snp1.allele2 == snp2.allele1:
                    counter["HET-HET"] += 1
                    if snp1.ps == reference_ps and snp2.ps == query_ps:
                        n2 += 1
                else:
                    counter["HET-HET.2"] += 1
                    
    data["Genotyping_Detail"] = counter
    data["Genotyping_Total"] = sum(counter.values())
    data["Genotyping_Identical"] = counter["HOM-HOM"] + counter["HET-HET"]
    data["Genotyping_Precision"] = np.divide(data["Genotyping_Identical"], data["Genotyping_Total"])
                    
    data["Phasing_Total"] = n1 + n2
    data["Phasing_Identical"] = max(n1, n2)
    data["Phasing_Precision"] = np.divide(data["Phasing_Identical"], data["Phasing_Total"])
    
    # print(data)
    
    return data
    

def benchmark_snp(args):
    parser = optparse.OptionParser(usage="%prog BenchmarkSNP [options] reference.vcf.gz query.vcf.gz")
    parser.add_option("-b", "--bed", dest="bed", metavar="PATH", 
                      help="Ignore SNPs outside of regions. [%default]")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1, metavar="INT", 
                      help="The number of parallel threads. [%default]")
    parser.add_option("-o", "--output", dest="output", metavar="PATH", default="benchmark_overall.json",
                      help="Save the benchmark results (JSON) to PATH. [%default]")
    options, args = parser.parse_args(args)
    f_ref, f_que = args
    
    chroms = []
    with pysam.VariantFile(f_ref) as f:
        chroms.extend(f.header.contigs)
    with pysam.VariantFile(f_que) as f:
        chroms.extend(f.header.contigs)
    chroms = list(sorted(set(chroms)))
    
    pool = mp.Pool(options.threads)
    results = []
    for chrom in chroms:
        r = pool.apply_async(worker, (chrom, f_ref, f_que, options.bed))
        results.append(r)
    pool.close()
    pool.join()
    
    all_results = [r.get() for r in results]
    
    pattern = "^chr[0-9]+$"
    results = list(filter(lambda r: re.match(pattern, r["Chrom"]), all_results))
        
    data = dict()
    
    data["Chrom_Pattern"] = pattern
    data["Reference_SNPs"] = sum([r["Reference_SNPs"] for r in results])
    data["Query_SNPs"] = sum([r["Query_SNPs"] for r in results])
    data["Filtered_Reference_SNPs"] = sum([r["Filtered_Reference_SNPs"] for r in results])
    data["Filtered_Query_SNPs"] = sum([r["Filtered_Query_SNPs"] for r in results])
    
    data["All_SNP_Reference"] = sum([r["All_SNP_Reference"] for r in results])
    data["All_SNP_Query"] = sum([r["All_SNP_Query"] for r in results])
    data["All_SNP_Overlap"] = sum([r["All_SNP_Overlap"] for r in results])
    data["All_SNP_Recall"] = np.divide(data["All_SNP_Overlap"], data["All_SNP_Reference"])
    data["All_SNP_Precision"] = np.divide(data["All_SNP_Overlap"], data["All_SNP_Query"])
    data["All_SNP_F1"] = np.divide(2 * data["All_SNP_Recall"] * data["All_SNP_Precision"], data["All_SNP_Recall"] + data["All_SNP_Precision"])
    
    data["Het_SNP_Reference"] = sum([r["Het_SNP_Reference"] for r in results])
    data["Het_SNP_Query"] = sum([r["Het_SNP_Query"] for r in results])
    data["Het_SNP_Overlap"] = sum([r["Het_SNP_Overlap"] for r in results])
    data["Het_SNP_Recall"] = np.divide(data["Het_SNP_Overlap"], data["Het_SNP_Reference"])
    data["Het_SNP_Precision"] = np.divide(data["Het_SNP_Overlap"], data["Het_SNP_Query"])
    data["Het_SNP_F1"] = np.divide(2 * data["Het_SNP_Recall"] * data["Het_SNP_Precision"], data["Het_SNP_Recall"] + data["Het_SNP_Precision"])
    
    data["Hom_SNP_Reference"] = sum([r["Hom_SNP_Reference"] for r in results])
    data["Hom_SNP_Query"] = sum([r["Hom_SNP_Query"] for r in results])
    data["Hom_SNP_Overlap"] = sum([r["Hom_SNP_Overlap"] for r in results])
    data["Hom_SNP_Recall"] = np.divide(data["Hom_SNP_Overlap"], data["Hom_SNP_Reference"])
    data["Hom_SNP_Precision"] = np.divide(data["Hom_SNP_Overlap"], data["Hom_SNP_Query"])
    data["Hom_SNP_F1"] = np.divide(2 * data["Hom_SNP_Recall"] * data["Hom_SNP_Precision"], data["Hom_SNP_Recall"] + data["Hom_SNP_Precision"])
    
    data["Calling_Reference"] = sum([r["Calling_Reference"] for r in results])
    data["Calling_Query"] = sum([r["Calling_Query"] for r in results])
    data["Calling_Overlap"] = sum([r["Calling_Overlap"] for r in results])
    data["Calling_Recall"] = np.divide(data["Calling_Overlap"], data["Calling_Reference"])
    data["Calling_Precision"] = np.divide(data["Calling_Overlap"], data["Calling_Query"])
    data["Calling_F1"] = np.divide(2 * data["Calling_Recall"] * data["Calling_Precision"], data["Calling_Recall"] + data["Calling_Precision"])
    
    counter = defaultdict(int)
    for r in results:
        print(r["Chrom"])
        for k, v in r["Genotyping_Detail"].items():
            counter[k] += v
    data["Genotyping_Detail"] = counter
    data["Genotyping_Total"] = sum([r["Genotyping_Total"] for r in results])
    data["Genotyping_Identical"] = sum([r["Genotyping_Identical"] for r in results])
    data["Genotyping_Precision"] = np.divide(data["Genotyping_Identical"], data["Genotyping_Total"])
                    
    data["Phasing_Total"] = sum([r["Phasing_Total"] for r in results])
    data["Phasing_Identical"] = sum([r["Phasing_Identical"] for r in results])
    data["Phasing_Precision"] = np.divide(data["Phasing_Identical"], data["Phasing_Total"])   
    
    data["Chromosomes"] = all_results
            
    if options.output:  
        with open(options.output, "w+") as fw:
            json.dump(data, fw, indent=4)
            