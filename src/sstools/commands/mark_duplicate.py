#!/usr/bin/env python
import optparse
import numpy as np
from collections import defaultdict
import pysam
from pyBioInfo.Range import GRange
# from pyBioInfo.IO.File import Alignment
from pyBioInfo.Utils import BundleBuilder


class Read(GRange):
    def __init__(self, primary, supplementaries, secondaries):
        chrom = primary.reference_name
        start = primary.reference_start
        end = primary.reference_end
        name = primary.query_name
        strand = "-" if primary.is_reverse else "+"
        super(Read, self).__init__(chrom=chrom, start=start, end=end, name=name, strand=strand)
        self.primary = primary
        self.supplementaries = supplementaries
        self.secondaries = secondaries
        self.query_length = primary.query_length
        

# def get_query_length(s):
#     return s.query_length


# def load_alignments(segments):
#     for s in segments:
#         a = Alignment(s)
#         a.query_length = get_query_length(s)
#         yield a


def cluster_by_strand(clusters):
    array = []
    for cluster in clusters:
        tmp1 = []
        tmp2 = []
        for a in cluster:
            if a.strand == "+":
                tmp1.append(a)
            elif a.strand == "-":
                tmp2.append(a)
            else:
                raise RuntimeError()
        if len(tmp1) > 0:
            array.append(tmp1)
        if len(tmp2) > 0:
            array.append(tmp2)
    return array


def cluster_by_length(clusters, max_foldchange):
    array = []
    for cluster in clusters:
        tmp = []
        for a in sorted(cluster, key=lambda x: x.query_length):
            if len(tmp) == 0:
                tmp.append([a])
            else:
                v1 = a.query_length
                v2 = tmp[-1][-1].query_length
                assert v1 >= v2
                if v1 <= v2 * max_foldchange:
                    tmp[-1].append(a)
                else:
                    tmp.append([a])
        array.extend(tmp)
    return array


def cluster_by_start(clusters, max_distance):
    array = []
    for cluster in clusters:
        tmp = []
        for a in sorted(cluster, key=lambda x: x.start):
            if len(tmp) == 0:
                tmp.append([a])
            else:
                if a.start - tmp[-1][-1].start <= max_distance:
                    tmp[-1].append(a)
                else:
                    tmp.append([a])
        array.extend(tmp)
    return array
                    

def cluster_by_end(clusters, max_distance):
    array = []
    for cluster in clusters:
        tmp = []
        for a in sorted(cluster, key=lambda x: x.end):
            if len(tmp) == 0:
                tmp.append([a])
            else:
                if a.end - tmp[-1][-1].end <= max_distance:
                    tmp[-1].append(a)
                else:
                    tmp.append([a])
        array.extend(tmp)
    return array
            
            
# def set_duplicate_flags(alignments, index):
#     alignments = list(sorted(alignments, key=lambda x: len(x), reverse=True))
#     for i, a in enumerate(alignments):
#         if i == 0:
#             a.segment.flag = a.segment.flag & 0b101111111111 # unique
#         else:
#             a.segment.flag = a.segment.flag | 0b010000000000 # dupplicate
#         a.segment.set_tag("DN", index)
#         a.segment.set_tag("DS", len(alignments))
#         a.segment.set_tag("DI", i)  # duplicate index


def set_duplicate_flags(reads, index):
    reads = list(sorted(reads, key=lambda x: len(x), reverse=True))
    for i, r in enumerate(reads):
        ss = [r.primary]
        if r.supplementaries:
            ss.extend(r.supplementaries)
        if r.secondaries:
            ss.extend(r.secondaries)
        for s in ss:
            if i == 0:
                s.flag = s.flag & 0b101111111111 # unique
            else:
                s.flag = s.flag | 0b010000000000 # dupplicate
            s.set_tag("DN", index)
            s.set_tag("DS", len(reads))
            s.set_tag("DI", i)  # duplicate index
            
      
# def run_pipeline(inbam, outbam, max_distance, strand_sense, read_matrix, group_matrix):
#     n_total = 0
#     n_uniq = 0
#     index = 0 # group index
    
#     h_read = None
#     if read_matrix:
#         h_read = open(read_matrix, "w+")
#         h_read.write("Read\tGroup\tSize\tIndex\n")
        
#     h_group = None
#     if group_matrix:
#         h_group = open(group_matrix, "w+")
#         h_group.write("Group\tSize\tReads\n")
        
#     with pysam.AlignmentFile(inbam) as f, pysam.AlignmentFile(outbam, "wb", f) as fw:
#         for chrom in f.references:
#             for bundle in BundleBuilder(load_alignments(f.fetch(chrom)), keep=True):
#                 alignments = bundle.data
#                 alignments_primary = []
#                 for a in alignments:
#                     if not a.segment.is_supplementary:
#                         alignments_primary.append(a)
#                 n_total += len(alignments_primary)
#                 clusters = [alignments_primary]
#                 if strand_sense:
#                     clusters = cluster_by_strand(clusters)
#                 while True:
#                     n = len(clusters)
#                     clusters = cluster_by_length(clusters, 1.2)
#                     clusters = cluster_by_start(clusters, max_distance)
#                     clusters = cluster_by_end(clusters, max_distance)
#                     if len(clusters) == n:
#                         break
#                 for cluster in clusters:
#                     set_duplicate_flags(cluster, index)                    
#                     if h_group:
#                         read_names = ",".join([a.name for a in cluster])
#                         line = "\t".join(map(str, [index, len(clusters), read_names]))
#                         h_group.write(line + "\n")
#                     index += 1
#                     n_uniq += 1
                
#                 for a in alignments:
#                     fw.write(a.segment)
#                     if h_read and not a.segment.is_supplementary:
#                         line = "\t".join(map(str, [
#                             a.segment.query_name, 
#                             a.segment.get_tag("DN"), 
#                             a.segment.get_tag("DS"), 
#                             a.segment.get_tag("DI")]))
#                         h_read.write(line + "\n")
#     if h_read:
#         h_read.close()        
#     if h_group:
#         h_group.close()
        
#     print("Total: %d" % n_total)
#     print("Uniq: %d (%.2f%%)" % (n_uniq, np.divide(n_uniq * 100, n_total)))
    

def run_pipeline(inbam, outbam, max_distance, strand_sense, read_matrix, group_matrix, mark_supplementary=False):
    n_total = 0
    n_uniq = 0
    index = 0 # group index
    
    h_read = None
    if read_matrix:
        h_read = open(read_matrix, "w+")
        h_read.write("Read\tGroup\tSize\tIndex\n")
        
    h_group = None
    if group_matrix:
        h_group = open(group_matrix, "w+")
        h_group.write("Group\tSize\tReads\n")
        
    with pysam.AlignmentFile(inbam) as f, \
        pysam.AlignmentFile(outbam, "wb", f) as fw:
            
        chroms = list(f.references)
        
        if mark_supplementary:
            all_reads = []
            all_segments = []
            data = defaultdict(list)
            for s in f.fetch(until_eof=True):
                all_segments.append(s)
                if s.is_unmapped:
                    continue
                data[s.query_name].append(s)
            for k, v in data.items():
                primaries = [] # primary
                supplementaries = [] # supplementaries
                secondaries = [] # secondaries
                for s in v:
                    if s.is_supplementary:
                        supplementaries.append(s)
                    elif s.is_secondary:
                        secondaries.append(s)
                    else:
                        primaries.append(s)
                if len(primaries) == 0:
                    raise RuntimeError()
                elif len(primaries) == 1:
                    r = Read(primaries[0], supplementaries, secondaries)
                    all_reads.append(r)
                else:
                    raise RuntimeError()
            chrom_reads = defaultdict(list)
            for r in all_reads:
                chrom_reads[r.chrom].append(r)
                
            for chrom in chroms:
                for bundle in BundleBuilder(sorted(chrom_reads[chrom]), keep=True):
                    n_total += len(bundle.data)
                    clusters = [bundle.data]
                    if strand_sense:
                        clusters = cluster_by_strand(clusters)
                    while True:
                        n = len(clusters)
                        clusters = cluster_by_length(clusters, 1.2)
                        clusters = cluster_by_start(clusters, max_distance)
                        clusters = cluster_by_end(clusters, max_distance)
                        if len(clusters) == n:
                            break
                    for cluster in clusters:
                        set_duplicate_flags(cluster, index)
                        
                        if h_group:
                            read_names = ",".join([a.name for a in cluster])
                            line = "\t".join(map(str, [index, len(clusters), read_names]))
                            h_group.write(line + "\n")
                        index += 1
                        n_uniq += 1
                        
            for s in all_segments:
                fw.write(s)
            
        else:
            raise NotImplementedError()
            
    if h_read:
        h_read.close()        
    if h_group:
        h_group.close()
        
    print("Total: %d" % n_total)
    print("Uniq: %d (%.2f%%)" % (n_uniq, np.divide(n_uniq * 100, n_total)))
        
    
usage = """sstools MarkDuplicate [options] <input.bam> <output.bam>

This command will added three TAGs to the BAM file:
    DN: name of duplicate set
    DS: size of duplicate set
    DI: index of read at duplicate set
"""


def mark_duplicate(args):
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-d", "--max-distance", dest="max_distance", type="int", default=20, metavar="INT",
                        help="Maximal distance between duplicates edge. [%default]")
    parser.add_option("-s", "--strand-sense", dest="strand_sense", action="store_true", default=False, 
                      help="")
    parser.add_option("-r", "--read-matrix", dest="read_matrix", metavar="PATH", 
                      help="")
    parser.add_option("-g", "--group-matrix", dest="group_matrix", metavar="PATH", 
                      help="")
    options, args = parser.parse_args(args)
    inbam, outbam = args
    read_matrix = options.read_matrix
    group_matrix = options.group_matrix
    strand_sense = options.strand_sense
    max_distance = options.max_distance
    
    run_pipeline(inbam, outbam, max_distance, strand_sense, read_matrix, group_matrix, mark_supplementary=True)


if __name__ == "__main__":
    mark_duplicate()