#!/usr/bin/env python
import sys
import os
import json
import re
import shutil
import optparse
from collections import OrderedDict
import multiprocessing as mp
import pysam


def get_deletion_blocks(segment):
    start = segment.reference_start
    blocks = [] # [start, end]
    for flag, count in segment.cigartuples:
        if flag == pysam.CMATCH:
            start += count
        elif flag == pysam.CDEL:
            end = start + count
            blocks.append([start, end])
            start = end
        elif flag == pysam.CINS:
            continue
        elif flag == pysam.CSOFT_CLIP or flag == pysam.CHARD_CLIP:
            continue
        elif flag == pysam.CEQUAL:
            start += count
        elif flag == pysam.CDIFF:
            start += count
        else:
            assert False
    return blocks


def get_insertion_blocks(segment):
    start = segment.reference_start
    offset = 0
    ins = [] # [start, length]
    for flag, count in segment.cigartuples:
        if flag == pysam.CMATCH:
            offset += count
        elif flag == pysam.CDEL:
            offset += count
        elif flag == pysam.CINS:
            ins.append([offset + start, count])
        elif flag == pysam.CSOFT_CLIP or flag == pysam.CHARD_CLIP:
            continue
        elif flag == pysam.CEQUAL:
            offset += count
        elif flag == pysam.CDIFF:
            offset += count
        else:
            print(flag, count)
            assert False
    return ins


def check_agree_and_disagree(record, segments):
    
    cells_agree = []
    cells_disagree = []
    
    if record.info["SVTYPE"] == "DEL":
        
        del_start, del_end = record.start, record.stop
        svlen = abs(record.info["SVLEN"])
        
        for segment in segments:
            
            segment_start, segment_end = segment.reference_start, segment.reference_end
            overlap_start, overlap_end = max(del_start, segment_start), min(del_end, segment_end)
            overlap_len = overlap_end - overlap_start
            
            assert overlap_len > 0
            # if overlap_len <= 0:
            #     continue
            
            if segment_start < del_start:
                if segment_end > del_end:
                    n = 0
                    for x, y in get_deletion_blocks(segment):
                        if max(x, del_start) < min(y, del_end):
                            n += (y - x)
                    if min(svlen, n) > max(svlen, n) * 0.7:
                    # if svlen * 0.5 < n < svlen * 1.5:
                        cells_agree.append(segment)
                    else:
                        cells_disagree.append(segment)
                else:
                    if overlap_len >= svlen * 0.1 or overlap_len >= 50:
                        cells_disagree.append(segment)
                    else:
                        continue
            else:
                if segment_end > del_end:
                    if overlap_len >= svlen * 0.1 or overlap_len >= 50:
                        cells_disagree.append(segment)
                    else:
                        continue
                else:
                    cells_disagree.append(segment)
    elif record.info["SVTYPE"] == "INS":
        pos = record.pos
        svlen = record.info["SVLEN"]
        pmin, pmax = pos - svlen * 0.5, pos + svlen * 0.5
        for segment in segments:
            distance_to_edge = min(pos - segment.reference_start, segment.reference_end - pos)
            if distance_to_edge < svlen or distance_to_edge < 200:
                continue            
            n = 0
            for v1, v2 in get_insertion_blocks(segment):
                if pmin < v1 < pmax:
                    n += v2
            if svlen * 0.5 < n < svlen * 1.5:
                cells_agree.append(segment)
            else:
                cells_disagree.append(segment)
    
    return cells_agree, cells_disagree


def get_cell_number(items):
    cells = [item[1] for item in items]
    cells = list(filter(lambda c: c != "Unknwon", cells))
    return len(set(cells))


def get_read_group(segment):
    if segment.has_tag("RG"):
        return segment.get_tag("RG")
    else:
        return "Bulk"


def _worker(f_vcf, f_bam, chrom, f_out):
    with pysam.VariantFile(f_vcf) as h_vcf, \
        pysam.AlignmentFile(f_bam) as h_bam, \
        open(f_out, "w+") as fw:

        # header
        line = "\t".join(["Chrom", "Start", "End", "Name", "Length",
                        "SupportRead", "SupportCell", 
                        "OverlapRead", "OverlapCell", 
                        "AgreeRead", "AgreeCell", 
                        "DisagreeRead", "DisagreeCell", 
                        "Detail"])
        fw.write(line + "\n")

        for record in h_vcf.fetch(chrom):
            
            svtype = record.info["SVTYPE"]
            if svtype != "DEL" and svtype != "INS":
                continue
            
            chrom, start, end = record.chrom, record.start, record.stop
            name = record.id
            svlen = abs(record.info["SVLEN"])
            
            segments = [s for s in h_bam.fetch(chrom, start, end)]
            segments = list(filter(lambda s: s.mapping_quality >= 30, segments))
            rname_rgs = {s.query_name: get_read_group(s) for s in segments}
            agree, disagree = check_agree_and_disagree(record, segments)

            data = OrderedDict()
            data["Support"] = [[rname, rname_rgs.get(rname, "Unknown")] for rname in record.info["RNAMES"]]
            data["Overlap"] = [[s.query_name, get_read_group(s)] for s in segments]
            data["Agree"] = [[s.query_name, get_read_group(s)] for s in agree]
            data["Disagree"] = [[s.query_name, get_read_group(s)] for s in disagree]
            s = json.dumps(data)
            
            n_support_cell = get_cell_number(data["Support"])
            n_overlap_cell = get_cell_number(data["Overlap"])
            n_agree_cell = get_cell_number(data["Agree"])
            n_disagree_cell = get_cell_number(data["Disagree"])
            
            line = "\t".join(map(str, [
                chrom, start, end, name, svlen,
                record.info["RE"], n_support_cell, 
                len(segments), n_overlap_cell, 
                len(agree), n_agree_cell, 
                len(disagree), n_disagree_cell, 
                s]))
            fw.write(line + "\n")
    return f_out


def quantify_sv(args=None):
    parser = optparse.OptionParser(usage="%prog [options] <sv.vcf.gz> <aligned.bam> <output.tsv>")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1, metavar="INT")
    parser.add_option("-d", "--temp-directory", dest="tmpdir", metavar="DIR")
    options, args = parser.parse_args(args)
    threads = options.threads
    tmpdir = options.tmpdir
    f_vcf, f_bam, f_out = args
    
    if tmpdir is None:
        tmpdir = f_out + ".TMP"
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
        
    f_tsv_list = []
    pool = mp.Pool(threads)
    with pysam.VariantFile(f_vcf) as f:
        for chrom in f.header.contigs:
            if re.match("^chr([0-9]+|[XYM])$", chrom) is None:
                continue
            f_tsv = os.path.join(tmpdir, "%s.tsv" % chrom)
            pool.apply_async(_worker, (f_vcf, f_bam, chrom, f_tsv))
            f_tsv_list.append(f_tsv)
    pool.close()
    pool.join()

    with open(f_out, "w+") as fw:
        for i, path in enumerate(f_tsv_list):
            with open(path) as f:
                for j, line in enumerate(f):
                    if j == 0 and i > 0:
                        continue
                    fw.write(line)            
    
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    
    
if __name__ == "__main__":
    quantify_sv()