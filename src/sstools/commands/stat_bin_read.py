#!/usr/bin/env python
import optparse
import pandas as pd
import multiprocessing as mp
import pysam
import matplotlib
matplotlib.use("agg")

usage = """

    sstools StatBinRead [options] <input.bam> <outfile>
"""


def worker(bamfile, chrom, width, rmdup, rmlc):
    with pysam.AlignmentFile(bamfile) as f:
        length = f.get_reference_length(chrom)
        nbin = int(length / width)
        if length % width > 0:
            nbin += 1
        
        rows = []
        offset = 5
        for i in range(nbin):
            start = i * width
            end = min(start + width, length)
            row = [chrom, i, start, end, end - start, 0, 0, 0, 0, 0, 0]
            rows.append(row)
            
        for s in f.fetch(chrom):
            if rmdup and s.is_duplicate:
                continue
            
            if rmlc and s.has_tag("XH") and s.get_tag("XH") == "N":
                continue
                
            start = s.reference_start
            end = s.reference_end
            strand = "-" if s.is_reverse else "+"
            if s.is_paired and s.is_read2: # Paired-end
                strand = "-" if strand == "+" else "+"
                
            parental = "U"
            if s.has_tag("XP"):
                parental = s.get_tag("XP")
                
            idx = int(start / width)
            if strand == "+":
                rows[idx][offset] += 1
                if parental == "P":
                    rows[idx][offset + 1] += 1
                elif parental == "M":
                    rows[idx][offset + 2] += 1
            else:
                rows[idx][offset + 3] += 1
                if parental == "P":
                    rows[idx][offset + 4] += 1
                elif parental == "M":
                    rows[idx][offset + 5] += 1
                
        columns = ["Chrom", "Bin", "Start", "End", "Width", 
               "Crick", "Crick.P", "Crick.M", 
               "Watson", "Watson.P", "Watson.M"]
        df = pd.DataFrame(rows, columns=columns)
        
        return df


def run_pipeline(inbam, outfile, threads=1, width=1000000, rmdup=False, rmlc=False):
    results = []
    pool = mp.Pool(threads)
    with pysam.AlignmentFile(inbam) as f:
        for chrom in f.references:
            args = (inbam, chrom, width, rmdup, rmlc)
            r = pool.apply_async(worker, args)
            results.append(r)
    pool.close()
    pool.join()
    
    dats = [r.get() for r in results]
    dat = pd.concat(dats, axis=0)
    dat.to_csv(outfile, sep="\t", index=False)
    

def stat_bin_read(args):
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-d", "--remove-duplicates", dest="rmdup", action="store_true", default=False,
                      help="Remove duplicate reads.")
    parser.add_option("-c", "--remove-low-confidence", dest="rmlc", action="store_true", default=False, 
                      help="Remove low confidence reads.")
    parser.add_option("-w", "--width", dest="width", type="int", default=1000000, metavar="INT", 
                      help="")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1, metavar="INT", 
                      help="")
    options, args = parser.parse_args(args)
    if len(args) != 2:
        parser.print_help()
        exit(1)
    
    run_pipeline(
        inbam=args[0], 
        outfile=args[1], 
        threads=options.threads, 
        width=options.width, 
        rmdup=options.rmdup, 
        rmlc=options.rmlc
    )

    

# class StatBinRead(object):
#     def __init__(self):
#         self.infile = None # input.bam
#         self.outdir = None # output directory
#         self.options = None
#         self.init_parameters()
#         if not os.path.exists(self.outdir):
#             os.mkdir(self.outdir)
#         self.execute()

#     def init_parameters(self):
#         parser = optparse.OptionParser(usage=usage)
#         parser.add_option("-w", "--bin-width", dest="bin_width", default=int(1e6), type="int",
#                           help="Bin width. [%default]")
#         parser.add_option("-n", "--seqname", dest="seqname", default="^chr([0-9]+|[XY])$",
#                           help="Seqname pattern. [%default]")
#         parser.add_option("-q", "--min-mapq", dest="min_mapq", default=30, type="int",
#                           help="Minimal mapping quality. [%default]")
#         parser.add_option("-s", "--remove-secondary", dest="remove_secondary", action="store_true", default=False, 
#                           help="Remove secondary alignments. [%default]")
#         parser.add_option("-S", "--remove-supplementary", dest="remove_supplementary", action="store_true", default=False, 
#                           help="Remove supplementary alignments. [%default]")
#         options, args = parser.parse_args(sys.argv[2:])
#         if len(args) != 2:
#             parser.print_help()
#             exit(1)
#         self.options = options
#         self.infile = args[0]
#         self.outdir = args[1]

#     def load_bam(self, path):
#         chroms = []
#         lengths = dict()
#         reads = defaultdict(list)
#         with pysam.AlignmentFile(path) as f:
#             for chrom in f.references:
#                 length = f.get_reference_length(chrom)
#                 if self.options.seqname != "" and re.match(self.options.seqname, chrom) is None:
#                     continue
#                 chroms.append(chrom)
#                 lengths[chrom] = length
#                 for segment in f.fetch(chrom):
#                     if segment.is_unmapped:
#                         continue
#                     if self.options.remove_secondary and segment.is_secondary:
#                         continue
#                     if self.options.remove_supplementary and segment.is_supplementary:
#                         continue
#                     if segment.mapping_quality < self.options.min_mapq:
#                         continue
#                     start = segment.reference_start
#                     end = segment.reference_end
#                     strand = "-" if segment.is_reverse else "+"
#                     if segment.is_paired and segment.is_read2: # Paired-end
#                         strand = "-" if strand == "+" else "+"
#                     is_dup = segment.is_duplicate
#                     is_hc = True # in high confidence region
#                     if segment.has_tag("XH"):
#                         is_hc = segment.get_tag("XH") == "Y"
#                     parental = "U"
#                     if segment.has_tag("XP"):
#                         parental = segment.get_tag("XP")
#                     reads[chrom].append([start, end, strand, is_dup, is_hc, parental])
#         return chroms, lengths, reads

#     def make_bin_read_count_table(self, chroms, lengths, reads, outdir):
#         data1 = dict() # key = [chrom, require_uniq, require_hc]
#         data2 = dict() # key = [require_uniq, reguire_hc]

#         bin_width = self.options.bin_width
        
#         chrom_bin_dir = "%s/chrom_bins" % outdir
#         if not os.path.exists(chrom_bin_dir):
#             os.mkdir(chrom_bin_dir)
            
#         for require_uniq in [True, False]:
#             for require_hc in [True, False]:
#                 suffix = ""
#                 if require_uniq:
#                     suffix += ".RmDup"
#                 if require_hc:
#                     suffix += ".IsHC"
#                 suffix += ".tsv"
                    
#                 array = []
#                 for chrom in chroms:
#                     length = lengths[chrom]
                    
#                     nbin = int(length / bin_width)
#                     if length % bin_width > 0:
#                         nbin += 1
                        
#                     rows = []
#                     for bi in range(nbin):
#                         start = bi * bin_width
#                         end = min((bi + 1) * bin_width, length)
#                         rows.append([chrom, bi, start, end, end - start])
#                     d0 = pd.DataFrame(rows)
#                     d0.columns = ["Chrom", "Bin", "Start", "End", "Width"]

#                     matrix = np.zeros((nbin, 6), dtype=np.int)
#                     for read in reads[chrom]:
#                         start, end, strand, is_dup, is_hc, parental = read
#                         if require_uniq and is_dup:
#                             continue
#                         if require_hc and not is_hc:
#                             continue
#                         bi = int(start / bin_width)
#                         if strand == "+":
#                             matrix[bi][0] += 1
#                             if parental == "P":
#                                 matrix[bi][1] += 1
#                             elif parental == "M":
#                                 matrix[bi][2] += 1
#                         else:
#                             matrix[bi][3] += 1
#                             if parental == "P":
#                                 matrix[bi][4] += 1
#                             elif parental == "M":
#                                 matrix[bi][5] += 1
#                     d1 = pd.DataFrame(matrix)
#                     d1.columns = ["Crick", "Crick.Paternal", "Crick.Maternal", "Watson", "Watson.Paternal", "Watson.Maternal"]
#                     d2 = pd.concat([d0, d1], axis=1)
#                     d2.to_csv("%s/%s%s" % (chrom_bin_dir, chrom, suffix), sep="\t", index=False)
#                     array.append(d2)
#                     data1[(chrom, require_uniq, require_hc)] = d2
#                 d3 = pd.concat(array, axis=0, ignore_index=True)
#                 d3.to_csv("%s/bins%s" % (outdir, suffix), sep="\t", index=False)
#                 data2[(require_uniq, require_hc)] = d3
#         return data1, data2

#     def plot_chrom_barplot_axes(self, ax, d, xlim=None, ylim=None, title=None, xlabel=None, ylabel=None, ytick=True):
#         # data
#         ys = np.arange(len(d)) + 0.5
#         xs1 = d["Crick"]
#         xs2 = d["Crick.Paternal"]
#         xs3 = d["Crick.Maternal"]
#         xs4 = d["Watson"]
#         xs5 = d["Watson.Paternal"]
#         xs6 = d["Watson.Maternal"]
#         mean1 = np.mean(xs1)
#         mean2 = np.mean(xs4)
#         if xlim is None:
#             xlim = max(max(max(xs1), max(xs4)) * 1.2, 10)

#         ax0 = plt.gca()
#         plt.sca(ax)
#         if title is None:
#             v1 = d["Watson"].mean()
#             v2 = d["Crick"].mean()
#             v3 = np.divide(min(v1, v2) * 100, v1 + v2)
#             s = "%.2f, %.2f (%.4f%%)" % (v1, v2, v3)
#             plt.title(s)
#         else:
#             plt.title(title)
#         # barplot
#         plt.barh(ys, xs1, height=1, color="C0")
#         plt.barh(ys, xs2, height=1, color="blue")
#         plt.barh(ys, xs3, left=xs2, height=1, color="red")
#         plt.barh(ys, -xs4, height=1, color="C1")
#         plt.barh(ys, -xs5, left=-xs6, height=1, color="blue")
#         plt.barh(ys, -xs6, height=1, color="red")
#         plt.plot([0, 0], [0, len(xs1)], lw=1, color="black")
#         plt.plot([mean1, mean1], [0, len(xs1)], ls="--", lw=1, color="grey")
#         plt.plot([-mean2, -mean2], [0, len(xs1)], ls="--", lw=1, color="grey")
        
#         for x, y in zip(xs1, ys):
#             if x >= xlim:
#                 plt.scatter(xlim, y + 0.5, s=10, marker="o", color="black", clip_on=False)
#         for x, y in zip(xs4, ys):
#             if x >= xlim:
#                 plt.scatter(-xlim, y + 0.5, s=10, marker="o", color="black", clip_on=False)
#         if xlabel:
#             plt.xlabel(xlabel)
#         if ylabel:
#             plt.ylabel(ylabel)
#         if not ytick:
#             plt.yticks([])
#         plt.xlim(-xlim, xlim)
#         plt.ylim(0, ylim)
#         if True:
#             ax.spines["top"].set_visible(False)
#             ax.spines["left"].set_visible(False)
#             ax.spines["right"].set_visible(False)
#         plt.sca(ax0)
          
#     def plot_chrom_bin_read_count(self, chroms, data, max_bin_count, outdir, cell=None):
#         barplot_dir = "%s/chrom_barplot" % outdir
#         if not os.path.exists(barplot_dir):
#             os.mkdir(barplot_dir)
            
#         for require_uniq in [True, False]:
#             for require_hc in [True, False]:
#                 suffix = ""
#                 title = None
#                 if require_uniq:
#                     suffix += ".RmDup"
#                 if require_hc:
#                     suffix += ".IsHC"
#                 suffix += ".pdf"
#                 for chrom in chroms:
#                     dat = data[(chrom, require_uniq, require_hc)]
#                     plt.figure(figsize=(3, 6))
#                     self.plot_chrom_barplot_axes(plt.gca(), dat, 
#                                             xlim=None, 
#                                             ylim=max_bin_count, 
#                                             title=title, 
#                                             xlabel=chrom, 
#                                             ylabel=cell, 
#                                             ytick=False)
#                     plt.tight_layout()
#                     plt.savefig("%s/%s%s" % (barplot_dir, chrom, suffix), dpi=300)
#                     plt.close()

#     def plot_genome_bin_read_count(self, chroms, data, max_bin_count, outdir, cell=None):
#         for require_uniq in [True, False]:
#             for require_hc in [True, False]:
#                 for trim in [True, False]:
#                     suffix = ""
#                     if require_uniq:
#                         suffix += ".RmDup"
#                     if require_hc:
#                         suffix += ".IsHC"
#                     if trim:
#                         suffix += ".Trim"
#                     suffix += ".pdf"
                
#                     d = data[(require_uniq, require_hc)]
#                     vs = list(d["Crick"]) + list(d["Watson"])
#                     if trim:
#                         vs = list(filter(lambda item: item > 0, vs))
#                         if len(vs) > 10:
#                             vs = vs[int(len(vs) * 0.1):int(len(vs) * 0.9)]
#                         if len(vs) > 0:
#                             mean = np.mean(vs)
#                             std = np.std(vs)
#                         else:
#                             mean = 0
#                             std = 0
#                         xlim = max(mean + 2 * std, 10)
#                         xlim = mean * 3
#                     else:
#                         xlim = max(vs) * 1.2
#                     nrow = 2
#                     ncol = 12
#                     fig, axs = plt.subplots(nrow, ncol, figsize=(ncol * 3, nrow * 6))

#                     if cell is None:
#                         cell = "Unknown"
                        
#                     title = "Cell: %s, Suffix: %s" % (cell, suffix)
#                     plt.suptitle(title)
                        
#                     for ci, chrom in enumerate(chroms):
#                         ax = axs[int(ci / ncol)][ci % ncol]
#                         d2 = d[d["Chrom"] == chrom]
#                         self.plot_chrom_barplot_axes(ax, d2, xlim=xlim, ylim=max_bin_count, xlabel=chrom, ytick=False)

#                     for ai in range(nrow * ncol):
#                         if ai >= len(chroms):
#                             axs[int(ai/ncol)][ai%ncol].set_visible(False)
#                     # axs[0][5].text(0, max_bin_count * 0.9, title, fontsize=24, va="center", ha="center")
#                     plt.tight_layout(rect=(0, 0, 1, 0.95))
#                     plt.savefig("%s/barplot%s" % (outdir, suffix), dpi=300)
#                     plt.close()
                    
#     def execute(self):
#         cell = self.infile.split("/")[-1]
#         if cell.endswith(".bam"):
#             cell = cell[:-4]
        
#         print("Loading bam file.")
#         chroms, lengths, reads = self.load_bam(self.infile)
#         # for chrom in chroms:
#         #     length = lengths[chrom]
#         #     print(chrom, length, sep="\t")
        
#         print("Making bin read count table.")
#         data1, data2 = self.make_bin_read_count_table(chroms=chroms, 
#                                                 lengths=lengths,
#                                                 reads=reads, 
#                                                 outdir=self.outdir)

#         max_chrom_length = max(lengths.values())
#         max_bin_count = int(max_chrom_length / self.options.bin_width)
#         if max_chrom_length % self.options.bin_width > 0:
#             max_bin_count += 1
#         max_bin_count
        
#         print("Plotting chromosome barplot.")
#         self.plot_chrom_bin_read_count(chroms=chroms, 
#                                 data=data1, 
#                                 max_bin_count=max_bin_count, 
#                                 outdir=self.outdir, 
#                                 cell=cell)
        
#         print("Plotting genome barplot.")
#         self.plot_genome_bin_read_count(chroms=chroms, 
#                                 data=data2, 
#                                 max_bin_count=max_bin_count, 
#                                 outdir=self.outdir, cell=cell)
