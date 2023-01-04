#!/usr/bin/env python#!/usr/bin/env python
import sys
import os
import json
import optparse
import multiprocessing
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "arial"
import matplotlib.pyplot as plt


usage = """

    sstools FetchCCRegion [options] <config.txt> <outdir>
"""


class FetchCCRegion(object):
    def __init__(self):
        self.txtfile = None
        self.outdir = None
        self.init_parameters()
        self.execute()

    def init_parameters(self):
        parser = optparse.OptionParser(usage=usage)
        parser.add_option("-p", "--processors", dest="threads", type="int", default=1,
                          help="Processors to run in parallel. [%default]")
        options, args = parser.parse_args(sys.argv[2:])
        if len(args) != 2:
            parser.print_help()
            exit(1)
        self.options = options
        self.txtfile = args[0]
        self.outdir = args[1]

    @classmethod
    def get_status(cls, v):
        if v > 2:
            return 1
        elif v < -2:
            return -1
        else:
            return 0

    @classmethod
    def cal_bin_log2_cwr(cls, v1, v2, lim=4):
        if v1 == 0:
            if v2 == 0:
                return 0
            else:
                return -lim
        else:
            if v2 == 0:
                return lim
            else:
                v = np.log2(np.divide(v1, v2))
                if v > 4:
                    v = 4
                elif v < -4:
                    v = -4
                return v

    @classmethod
    def make_data_frame(cls, data, chrom, outdir):
        cells = list(sorted(data.keys()))
        rows_crick = []
        rows_watson = []
        rows_cwr = []
        for cell in cells:
            cwr_list = []
            for item in data[cell]["Detail"]:
                if item["Chrom"] == chrom:
                    counts1, counts2 = item["Bin_Values"]
                    rows_crick.append(counts1)
                    rows_watson.append(counts2)
                    log2cwr = np.log2(np.divide(sum(counts1), sum(counts2)))
                    for c1, c2 in zip(counts1, counts2):
                        v = cls.cal_bin_log2_cwr(c1, c2, lim=4)
                        if log2cwr < 0:
                            v = v * -1
                        cwr_list.append(v)
                    break

            if len(cwr_list) == 0:
                # print(chrom)
                assert False
            rows_cwr.append(cwr_list)
        crick = pd.DataFrame(rows_crick)
        crick.index = cells
        crick.index.name = "Cell"
        crick.to_csv(outdir + "/crick.raw.tsv", sep="\t")
        watson = pd.DataFrame(rows_watson)
        watson.index = cells
        watson.index.name = "Cell"
        watson.to_csv(outdir + "/watson.raw.tsv", sep="\t")
        dat = pd.DataFrame(rows_cwr)
        dat.index = cells
        dat.index.name = "Cell"
        dat.to_csv(outdir + "/cwr.raw.tsv", sep="\t")
        return crick, watson, dat

    @classmethod
    def filter_and_cluster_cells(cls, crick, watson, dat, outdir):
        flags = np.log2(crick.sum(axis=1) / watson.sum(axis=1)).abs() > 2
        crick, watson, dat = crick[flags], watson[flags], dat[flags]

        ret = sns.clustermap(dat, cmap="bwr", figsize=(8, 10), col_cluster=False)
        plt.savefig(outdir + "/clustermap.pdf", dpi=300)
        plt.close()
        dat = ret.data2d
        crick = crick.loc[dat.index]
        watson = watson.loc[dat.index]
        dat.to_csv(outdir + "/cwr.filtered.tsv", sep="\t")
        crick.to_csv(outdir + "/crick.filtered.tsv", sep="\t")
        watson.to_csv(outdir + "/watson.filtered.tsv", sep="\t")
        # Heatmap
        plt.figure(figsize=(10, 10))
        sns.heatmap(dat, cmap="bwr", vmin=-4, vmax=4, cbar_kws={"label": "Log2CWR"})
        plt.xlabel("Bins")
        plt.ylabel("Cells (%d)" % len(dat))
        plt.tight_layout()
        plt.savefig(outdir + "/heatmap.pdf", dpi=300)
        plt.close()

        return crick, watson, dat

    @classmethod
    def get_bin_width(cls, data):
        for d in data.values():
            return d["Detail"][0]["Bin_Width"]

    @classmethod
    def get_cc_regions(cls, crick, watson, dat, chrom, bin_width, outdir):
        means = dat.mean(axis=0)
        regions = []
        masks = []
        for row_idx in range(len(dat)):
            c = crick.iloc[row_idx].sum()
            w = watson.iloc[row_idx].sum()
            chrom_log2cwr = np.log2(np.divide(c, w))
            # chrom_status = "CC" if chrom_log2cwr > 0 else "WW"

            values = dat.iloc[row_idx]
            # status = list(map(cls.get_status, dat.iloc[row_idx]))
            status = [] # 1: same, 0: diff
            for bin_idx, (m, v) in enumerate(zip(means, values)):
                if abs(v - m) < 2:
                    status.append(1)
                else:
                    c = crick.values[row_idx][bin_idx]
                    w = watson.values[row_idx][bin_idx]
                    if c + w < 10:
                        status.append(1)
                    else:
                        status.append(0)
            array = [] # [start, end]
            start = None
            for j, v in enumerate(status):
                if v == 0:
                    if start is None:
                        continue
                    else:
                        array.append([start, j])
                        start = None
                else:
                    if start is None:
                        start = j
                    else:
                        continue
            if start is not None:
                array.append([start, len(status)])
                
            # Filter too short
            array = list(filter(lambda r: r[1] - r[0] >= 2, array))
            # Fill small gap
            j = 0
            while j < len(array) - 1:
                if array[j + 1][0] - array[j][1] < max(int(len(dat.columns) * 0.03), 2):
                    array[j][1] = array[j + 1][1]
                    array.pop(j + 1)
                else:
                    j += 1
            # Fill edge
            if len(array) > 0:
                if array[0][0] <= 2:
                    array[0][0] = 0
                if array[-1][1] + 2 >= len(status):
                    array[-1][1] = len(status)
            # Keep longest region
            if len(array) > 1:
                array = list(sorted(array, key=lambda r: r[1] - r[0]))
                array = array[-1:]
            # Trim bound
            new_regions = []
            for start, end in array:
                if start != 0:
                    start += 2
                if end != len(status):
                    end -= 2
                if end - start > len(status) * 0.3:
                    if sum(status[start:end]) / (end - start) >= 0.9:
                        new_regions.append([start, end])
            array = new_regions
            # Report
            for r in array:
                regions.append([chrom, r[0] * bin_width, r[1] * bin_width, \
                    dat.index.values[row_idx], ".", "+" if chrom_log2cwr > 0 else "-"])
                
            # New status
            new_status = np.zeros(len(status))
            for r in array:
                for pos in np.arange(r[0], r[1]):
                    new_status[pos] = 1
            masks.append(new_status)
        masks = pd.DataFrame(masks)
        masks.index = dat.index
        plt.figure(figsize=(10, 10))
        sns.heatmap(masks, cmap="bwr", vmin=0, vmax=1)
        plt.tight_layout()
        plt.savefig(outdir + "/heatmap.cc.pdf", dpi=300)

        with open(outdir + "/cc_regions.bed", "w+") as fw:
            for row in regions:
                line = "\t".join(map(str, row))
                fw.write(line + "\n")

        return masks, regions

    @classmethod
    def barplot(cls, crick, watson, masks, outdir):
        ncol = 4
        nrow = int(len(masks) / ncol)
        if len(masks) % ncol > 0:
            nrow += 1
        fig, axs = plt.subplots(nrow, ncol, figsize=(ncol * 3, nrow * 1), sharex=True)
        for i in range(len(masks)):
            ax = axs[int(i / ncol)][i % ncol]
            plt.sca(ax)
            plt.title(masks.index.values[i])
            counts1 = crick.iloc[i]
            counts2 = watson.iloc[i]
            status = masks.iloc[i]
            if sum(counts1) < sum(counts2):
                counts1, counts2 = counts2, counts1
            colors1 = []
            colors2 = []
            for s in status:
                if s == 0:
                    colors1.append("red")
                    colors2.append("red")
                else:
                    colors1.append("C0")
                    colors2.append("C1")
            xs = np.arange(len(counts1))
            plt.bar(xs, counts1, width=1, color=colors1)
            plt.bar(xs, -counts2, width=1, color=colors2)
        plt.tight_layout()
        plt.savefig(outdir + "/cell_cc_regions.pdf", dpi=300)
        plt.close()

    @classmethod
    def process_chromosome(cls, data, chrom, outdir):
        outdir = os.path.join(outdir, chrom)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        bin_width = cls.get_bin_width(data)
        crick, watson, dat = cls.make_data_frame(data, chrom, outdir)
        crick, watson, dat = cls.filter_and_cluster_cells(crick, watson, dat, outdir)
        masks, regions = cls.get_cc_regions(crick, watson, dat, chrom, bin_width, outdir)
        cls.barplot(crick, watson, masks, outdir)
        return regions

    def execute(self):
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        config = json.load(open(self.txtfile))
        data = dict()
        for cell in config["Cells"]:
            run = cell.split(".")[0]
            path = "results/mapping/mark_parental/%s/%s.stats.json" % (run, cell)
            data[cell] = json.load(open(path))

        pool = None
        if self.options.threads > 1:
            pool = multiprocessing.Pool(self.options.threads)
        results = []
        for chrom in config["Chroms"]:
            # if chrom != "chrX":
            #     continue
            args = (data, chrom, self.outdir)
            if pool is None:
                r = self.process_chromosome(*args)
            else:
                r = pool.apply_async(self.process_chromosome, args)
            results.append(r)
            # break
        if pool is not None:
            pool.close()
            pool.join()
            results = [r.get() for r in results]
        
        with open(self.outdir + "/cc_regions.bed", "w+") as fw:
            for regions in results:
                for row in regions:
                    line = "\t".join(map(str, row))
                    fw.write(line + "\n")