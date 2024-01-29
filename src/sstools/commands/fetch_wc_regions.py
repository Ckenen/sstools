#!/usr/bin/env python#!/usr/bin/env python
import sys
import os
import re
import json
import optparse
from collections import defaultdict
from PyPDF2 import PdfFileMerger
import numpy as np
from scipy.stats import ttest_ind_from_stats
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "arial"
import matplotlib.pyplot as plt
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import BedFile
from pyBioInfo.Utils import ShiftLoader


usage = """

    sstools FetchWCRegion [options] <input.bam> <outdir>
"""


class FetchWCRegion(object):
    def __init__(self):
        self.infile = None
        self.outdir = None
        self.init_parameters()
        self.execute()

    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-b", "--blackholes", dest="blackhole",
                          help="Blackhole regions in BED format. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        self.options = options
        self.infile = args[0]
        self.outdir = args[1]

    def load_reads(self, bamfile):
        chroms = []
        lengths = dict()
        reads = dict()
        with pysam.AlignmentFile(bamfile) as f:
            for chrom in f.references:
                # if re.match("^chr([0-9]+|[XY])$", chrom) is None:
                #     continue
                length = f.get_reference_length(chrom)
                chroms.append(chrom)
                lengths[chrom] = length
            for chrom in chroms:
                reads1 = []
                for segment in f.fetch(chrom):
                    if segment.is_duplicate:
                        continue
                    chrom = segment.reference_name
                    start = segment.reference_start
                    end = segment.reference_end
                    strand = "-" if segment.is_reverse else "+"
                    read = GRange(chrom=chrom, start=start, end=end, strand=strand)
                    read.parental = "U"
                    if segment.has_tag("XP"):
                        read.parental = segment.get_tag("XP")
                    read.duplicate_set_name = segment.get_tag("DN")
                    reads1.append(read)
                reads[chrom] = reads1
        return chroms, lengths, reads

    def load_blackholes(self, bedfile):
        blackholes = defaultdict(list)
        if bedfile:
            with BedFile(bedfile) as f:
                for region in f:
                    blackholes[region.chrom].append(region)
        return blackholes

    def remove_blackhole_reads(self, chrom_reads, chrom_blackholes):
        new_chrom_reads = dict()
        for chrom, reads in chrom_reads.items():
            new_reads = [] # fitlered
            regions = chrom_blackholes[chrom]
            regions.sort()
            loader = ShiftLoader(regions)
            for read in reads:
                hits = list(loader.fetch(obj=read))
                if len(hits) == 0:
                    new_reads.append(read)
            new_chrom_reads[chrom] = new_reads
        return new_chrom_reads

    def make_matrix(self, reads, chrom_length, bin_width):
        # Matrix: [[crick, crick.p, crick.m, watson, watson.p, watson.m]]
        bin_count = int(chrom_length / bin_width)
        if chrom_length % bin_width > 0:
            bin_count += 1
        matrix = np.zeros((bin_count, 6), dtype=np.int)
        for read in reads:
            idx = int(read.start / bin_width)
            parental = read.parental
            if read.strand == "+":
                matrix[idx][0] += 1
                if parental == "P":
                    matrix[idx][1] += 1
                elif parental == "M":
                    matrix[idx][2] += 1
            else:
                matrix[idx][3] += 1
                if parental == "P":
                    matrix[idx][4] += 1
                elif parental == "M":
                    matrix[idx][5] += 1
        return matrix

    def plot_genome_bin_barplot(self, chroms, chrom_matrixes, outdir):
        ys = []
        colors = []
        seperate_lines = []
        offset = 0
        xticks = []
        ps = [] # xs of xticks
        for i, chrom in enumerate(chroms):
            if i != 0:
                seperate_lines.append(offset)
            matrix = chrom_matrixes[chrom]
            color = "C%d" % (i % 10)
            for crick, watson in zip(matrix[:, 0], matrix[:, 3]):
                ys.append(crick + watson)
                colors.append(color)
            s = chrom
            if s.startswith("chr"):
                s = s[3:]
            xticks.append(s)
            ps.append(offset + len(matrix) / 2)
            offset += len(matrix)
        ys = np.array(ys)
        
        ys1 = ys[ys>0]
        mean = np.mean(ys1)
        median = np.median(ys1)
        std = np.std(ys1, ddof=1)
            
        plt.figure(figsize=((len(ys) + 150) * 0.005, 3))
        plt.title("mean: %.2f, median: %d, std: %.2f" % (mean, median, std))
        plt.bar(np.arange(len(ys)) + 0.5, ys, width=1, color=colors)
        for x in seperate_lines:
            plt.axvline(x, lw=1, ls="--", color="grey")
        plt.axhline(median, ls="--", lw=1, color="red")
        plt.xticks(ps, xticks, rotation=0)
        plt.xlim(0, len(ys))
        plt.ylabel("Read count / bin")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "genome_bin_barplot.pdf"), dpi=300)

    def make_xs_ys(self, reads):
        x, y = 0, 0
        xs, ys = [x], [y]
        for read in reads:
            x += 1
            if read.strand == "+":
                y += 1
            else:
                y -= 1
            xs.append(x)
            ys.append(y)
        xs, ys = np.array(xs), np.array(ys)
        return xs, ys

    def make_ks(self, xs, ys, window_size):
        ks = [] 
        for i1 in range(0, len(xs)):
            i2 = i1 + window_size
            if i2 >= len(xs):
                break
            v1, v2 = ys[i1], ys[i2]
            k = abs((v2 - v1) / window_size)
            ks.append(k)
        ks = np.array(ks)
        return ks

    def filter_reads_by_k(self, chrom, reads, window_size, cutoff, outdir):
        reads = reads.copy()
        round_num = 1 # round

        draw_data = []
        
        while len(reads) >= window_size:
            xs, ys = self.make_xs_ys(reads)
            ks = self.make_ks(xs, ys, window_size)

            if False:
                d = dict()
                d["round"] = round_num
                d["xs"] = xs
                d["ys"] = ys
                d["ks"] = ks
                draw_data.append(d)
            
            i1 = None
            tmp1 = []
            for i2, k in enumerate(ks):
                if k >= cutoff:
                    if i1 is None:
                        i1 = i2
                else:
                    if i1 is not None:
                        tmp1.append([i1, i2])
                        i1 = None
            if i1 is not None:
                tmp1.append([i1, len(ks)])
            if len(tmp1) == 0:
                break

            tmp2 = []
            for i1, i2 in tmp1:
                i2 += window_size
                if len(tmp2) == 0:
                    tmp2.append([i1, i2])
                elif i1 <= tmp2[-1][1]:
                    tmp2[-1][1] = max(tmp2[-1][1], i2)
                else:
                    tmp2.append([i1, i2])

            e = 0
            tmp3 = []
            for i1, i2 in tmp2:
                if i1 > e:
                    tmp3.append([e, i1])
                e = i2
            if len(reads) > e:
                tmp3.append([e, len(reads)])

            tmp4 = []
            for i1, i2 in tmp3:
                if i2 - i1 < window_size:
                    continue
                tmp4.extend(reads[i1:i2])
            reads = tmp4

            round_num += 1

        if len(draw_data) > 0:
            fig, axs = plt.subplots(len(draw_data), 2, figsize=(8, 4 * len(draw_data)))
            for i, d in enumerate(draw_data):
                round_num = d["round"]
                xs, ys, ks = d["xs"], d["ys"], d["ks"]
                title = "%s, size: %d, k: %.2f, round: %d" % (chrom, window_size, cutoff, round_num)
                xlim1, xlim2 = 0, len(xs) - 1
                ylim1 = (max(ys) + min(ys)) / 2 - (xlim2 - xlim1) / 2
                ylim2 = (max(ys) + min(ys)) / 2 + (xlim2 - xlim1) / 2
                # xs, ys
                if len(draw_data) == 1:
                    plt.sca(axs[0])
                else:
                    plt.sca(axs[i][0])
                plt.title(title)
                plt.plot(xs, ys)
                plt.xlim(xlim1, xlim2)
                plt.ylim(ylim1, ylim2)
                plt.xlabel("Read index")
                # ks
                if len(draw_data) == 1:
                    plt.sca(axs[1])
                else:
                    plt.sca(axs[i][1])
                plt.title(title)
                plt.plot(np.arange(len(ks)), ks)
                plt.ylim(-0.1, 1.1)
                plt.axhline(cutoff, ls="--", lw=1, color="grey")
                plt.xlabel("Window index")
                plt.ylabel("|K|")
            plt.tight_layout()
            plt.savefig(outdir + "/%s.pdf" % chrom, dpi=300)
            plt.close()

        return reads

    def filter_ex(self, chroms, chrom_matrixes, chrom_reads_no_ccw, bin_width):
        read_counts = []
        for chrom in chroms:
            matrix = chrom_matrixes[chrom]
            read_counts.extend(matrix[:,0] + matrix[:,3])
        read_counts = np.array(read_counts)
        read_counts_no_zero = read_counts[read_counts > 0]
        bin_read_mean = np.mean(read_counts_no_zero)
        bin_read_std = np.std(read_counts_no_zero, ddof=1)
        chrom_reads_final = dict()
        for chrom in chroms:
            matrix = chrom_matrixes[chrom]
            vs = matrix[:, 0] + matrix[:, 3]
            flags = [1] * len(matrix)
            for i1 in range(0, len(vs)):
                i2 = i1 + 10
                if i2 > len(vs):
                    break
                mean = np.mean(vs[i1:i2])
                std = np.std(vs[i1:i2], ddof=1)
                fc =  mean / bin_read_mean
                p = ttest_ind_from_stats(bin_read_mean, bin_read_std, 10, mean, std, 10)[1]
                if fc > 2 and p < 0.01:
                    for i3 in range(i1, i2):
                        flags[i3] = 0
            new_reads = []
            for read in chrom_reads_no_ccw[chrom]:
                idx = int(read.start / bin_width)
                if flags[idx] == 0:
                    continue
                new_reads.append(read)
            chrom_reads_final[chrom] = new_reads
        return chrom_reads_final

    def plot_barplot(self, chroms, chrom_matrixes, chrom_matrixes_no_ccw, chrom_matrixes_final, bin_width, outdir):
        max_bin_read_count = 0
        max_bin_count = 0
        for chrom in chroms:
            matrix = chrom_matrixes[chrom]
            max_bin_read_count = max(max_bin_read_count, max(matrix[:, 0]))
            max_bin_read_count = max(max_bin_read_count, max(matrix[:, 3]))
            max_bin_count = max(max_bin_count, len(matrix))
        
        for chrom in chroms:
            matrix1 = chrom_matrixes[chrom]
            matrix2 = chrom_matrixes_no_ccw[chrom]
            matrix3 = chrom_matrixes_final[chrom]
            bin_count = len(matrix1)
            xs = np.arange(bin_count)
            ylim = max(max(matrix1[:, 0]), max(matrix1[:, 3])) * 1.2
            
            fig, axs = plt.subplots(1, 3, figsize=(18, 2), sharex=True, sharey=True)
            
            plt.sca(axs[0])
            plt.title("All reads")
            plt.bar(xs, matrix1[:, 0], width=1, color="C0")
            plt.bar(xs, matrix1[:, 1], width=1, color="blue")
            plt.bar(xs, matrix1[:, 2], bottom=matrix1[:, 1], width=1, color="red")
            plt.bar(xs, -matrix1[:, 3], width=1, color="C1")
            plt.bar(xs, -matrix1[:, 4], bottom=-matrix1[:, 5], width=1, color="blue")
            plt.bar(xs, -matrix1[:, 5], width=1, color="red")
            plt.plot([0, bin_count], [0, 0], lw=0.5, color="black")
            plt.xlim(0, max_bin_count)
            plt.ylim(-ylim, ylim)
            plt.xlabel("Bins (%d)" % bin_width)
            plt.ylabel(chrom)
            plt.tight_layout()

            plt.sca(axs[1])
            plt.title("WC reads")
            plt.bar(xs, matrix1[:, 0], width=1, color="grey")
            plt.bar(xs, -matrix1[:, 3], width=1, color="grey")
            plt.bar(xs, matrix2[:, 0], width=1, color="C0")
            plt.bar(xs, matrix2[:, 1], width=1, color="blue")
            plt.bar(xs, matrix2[:, 2], bottom=matrix2[:, 1], width=1, color="red")
            plt.bar(xs, -matrix2[:, 3], width=1, color="C1")
            plt.bar(xs, -matrix2[:, 4], bottom=-matrix2[:, 5], width=1, color="blue")
            plt.bar(xs, -matrix2[:, 5], width=1, color="red")
            plt.plot([0, bin_count], [0, 0], lw=0.5, color="black")
            plt.xlim(0, max_bin_count)
            plt.ylim(-ylim, ylim)
            plt.xlabel("Bins (%d)" % bin_width)
            plt.tight_layout()
            
            plt.sca(axs[2])
            plt.title("Final reads")
            plt.bar(xs, matrix1[:, 0], width=1, color="grey")
            plt.bar(xs, -matrix1[:, 3], width=1, color="grey")
            plt.bar(xs, matrix3[:, 0], width=1, color="C0")
            plt.bar(xs, matrix3[:, 1], width=1, color="blue")
            plt.bar(xs, matrix3[:, 2], bottom=matrix3[:, 1], width=1, color="red")
            plt.bar(xs, -matrix3[:, 3], width=1, color="C1")
            plt.bar(xs, -matrix3[:, 4], bottom=-matrix3[:, 5], width=1, color="blue")
            plt.bar(xs, -matrix3[:, 5], width=1, color="red")
            plt.plot([0, bin_count], [0, 0], lw=0.5, color="black")
            plt.xlim(0, max_bin_count)
            plt.ylim(-ylim, ylim)
            plt.xlabel("Bins (%d)" % bin_width)
            plt.tight_layout()
            
            plt.savefig(outdir + "/%s.pdf" % chrom, dpi=300)
            plt.close()

    def execute(self):
        matrix_dir = self.outdir + "/matrix"
        filter_cc_dir = self.outdir + "/filtered.cc"
        filter_ccw_dir = self.outdir + "/filtered.ccw"
        filter_final_dir = self.outdir + "/filtered.final"
        for d in [self.outdir, matrix_dir, filter_cc_dir, filter_ccw_dir, filter_final_dir]:
            if not os.path.exists(d):
                os.mkdir(d)

        bin_width = 1000000

        chroms, chrom_lengths, chrom_reads = self.load_reads(self.infile)
        chrom_blackholes = self.load_blackholes(self.options.blackhole)
        chrom_reads = self.remove_blackhole_reads(chrom_reads, chrom_blackholes)
        chrom_matrixes = dict()
        for chrom in chroms:
            reads = chrom_reads[chrom]
            length = chrom_lengths[chrom]
            matrix = self.make_matrix(reads, length, bin_width)
            chrom_matrixes[chrom] = matrix
        self.plot_genome_bin_barplot(chroms, chrom_matrixes, self.outdir)

        # threshold
        read_counts = []
        for chrom in chroms:
            matrix = chrom_matrixes[chrom]
            read_counts.extend(matrix[:,0] + matrix[:,3])
        read_counts = np.array(read_counts)
        read_counts_no_zero = read_counts[read_counts > 0]
        bin_read_median = np.median(read_counts_no_zero)

        # filter CC or WW reads
        window_size = max(int(bin_read_median * 2), 50)
        chrom_reads_no_cc = dict() # not CC or WW
        chrom_matrixes_no_cc = dict()
        for chrom in chroms:
            length = chrom_lengths[chrom]
            reads = chrom_reads[chrom]
            new_reads = self.filter_reads_by_k(chrom, reads, window_size, 0.75, filter_cc_dir)
            chrom_reads_no_cc[chrom] = new_reads
            chrom_matrixes_no_cc[chrom] = self.make_matrix(new_reads, length, bin_width)

        # fitler CCW or CWW reads
        window_size = max(int(bin_read_median * 30), 750)
        chrom_reads_no_ccw = dict() # not CC or WW
        chrom_matrixes_no_ccw = dict()
        for chrom in chroms:
            reads = chrom_reads_no_cc[chrom]
            length = chrom_lengths[chrom]
            new_reads = self.filter_reads_by_k(chrom, reads, window_size, 0.25, filter_ccw_dir)
            chrom_reads_no_ccw[chrom] = new_reads
            chrom_matrixes_no_ccw[chrom] = self.make_matrix(new_reads, length, bin_width)

        # filter extremely high read count bin
        chrom_reads_final = self.filter_ex(chroms, chrom_matrixes, chrom_reads_no_ccw, bin_width)
        chrom_matrixes_final = dict()
        for chrom in chroms:
            length = chrom_lengths[chrom]
            reads = chrom_reads_final[chrom]
            chrom_matrixes_final[chrom] = self.make_matrix(reads, length, bin_width)

        # chromosome barplot
        self.plot_barplot(chroms, chrom_matrixes, chrom_matrixes_no_ccw, chrom_matrixes_final, bin_width, filter_final_dir)
        pdf = PdfFileMerger()
        for chrom in chroms:
            path = filter_final_dir + "/%s.pdf" % chrom
            pdf.append(path)
        pdf.write(self.outdir + "/chromosome_bin_barplot.pdf")
            
        # WC reads duplicate set name
        data = dict()
        for chrom in chroms:
            # WC bin number more than half of chromosome
            matrix = chrom_matrixes_final[chrom]
            vs = matrix[:, 0] + matrix[:, 3]
            if sum(vs > 0) / len(vs) < 0.5:
                data[chrom] = []
                continue
            names = []
            reads = chrom_reads_final[chrom]
            for read in reads:
                names.append(read.duplicate_set_name)
            data[chrom] = names
        with open(self.outdir + "/duplicate_set_names.json", "w+") as fw:
            json.dump(data, fw)

