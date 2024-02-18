#!/usr/bin/env python
import sys
import os
import optparse
import glob
import multiprocessing as mp
from collections import OrderedDict
import subprocess
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "arial"
import matplotlib.pyplot as plt


def load_config(configfile):
    config = OrderedDict()
    with open(configfile) as f:
        for line in f:
            sample, path = line.strip("\n").split("\t")[:2]
            assert os.path.exists(path)
            config[sample] = path
    return config


def load_read_counts(config):
    sample_counts = dict()
    for sample, path in config.items():
        sample_counts[sample] = pd.read_csv(path, sep="\t")
    return sample_counts


def get_bin_width(data):
    width = None
    for cell, d in data.items():
        for chrom, d1 in d.groupby(by="Chrom"):
            if len(d1) > 1:
                width = d["Width"].values[0]
                break
        if width is not None:
            break
    assert width is not None
    return width


def get_chrom_lengths(data):
    info = None
    for sample, d in data.items():
        tmp = dict()
        for chrom, d1 in d.groupby(by="Chrom"):
            length = max(d1["End"])
            tmp[chrom] = length
        if info is None:
            info = tmp
        else:
            assert len(tmp) == len(info)
            for k in info.keys():
                assert info[k] == tmp[k]
    return info


def make_proportion_matrix(Mc, Mw, min_reads=2):
    rows = []
    for i in range(len(Mc)):
        crick = Mc.iloc[i].values
        watson = Mw.iloc[i].values
        ps = []
        for c, w in zip(crick, watson):
            if c <= min_reads and w <= min_reads:
                p = 0.5
            else:
                p = c / (c + w)
            ps.append(p)
        rows.append(ps)
    Mp = pd.DataFrame(rows, index=Mc.index, columns=Mc.columns)
    return Mp


def make_matrix(sample_counts, chrom):
    cells = list(sorted(sample_counts.keys()))
    rows_c, rows_w = [], []
    for cell in cells:
        d = sample_counts[cell]
        d = d[d["Chrom"] == chrom]
        rows_c.append(d["Crick"].values)
        rows_w.append(d["Watson"].values)
    Mc = pd.DataFrame(rows_c, index=cells)
    Mw = pd.DataFrame(rows_w, index=cells)
    Mp = make_proportion_matrix(Mc, Mw)
    return Mc, Mw, Mp


def predict_chrom_pattern(Mc, Mw):
    data = dict()
    for i, cell in enumerate(Mc.index):
        vs1 = Mc.iloc[i].values
        vs2 = Mw.iloc[i].values
        s1, s2 = sum(vs1), sum(vs2)
        if s1 > s2:
            pattern = "CC"
        elif s1 < s2:
            pattern = "WW"
        else:
            pattern = "WC"
        data[cell] = pattern
    return data


def cluster_matrix(Mc, Mw, Mp, outfile=None):
    ret = sns.clustermap(Mp, cmap="bwr", col_cluster=False, figsize=(10, 10))
    ret.fig.tight_layout()
    if outfile:
        ret.fig.savefig(outfile, dpi=300)
    samples = ret.data2d.index
    return Mc.loc[samples], Mw.loc[samples], Mp.loc[samples]


def filter_matrix(Mc, Mw, Mp, min_reads=100, min_fold_change=4):
    samples = []
    for i, sample in enumerate(Mc.index):
        vs1 = Mc.iloc[i].values
        vs2 = Mw.iloc[i].values
        s1, s2 = sum(vs1), sum(vs2)
        if s1 + s2 >= min_reads:
            if s1 >= s2 * min_fold_change or s2 >= s1 * min_fold_change:
                samples.append(sample)
    return Mc.loc[samples], Mw.loc[samples], Mp.loc[samples]
                   

def reverse_matrix(Mc, Mw, Mp):
    rows_c, rows_w = [], []
    for i in range(len(Mc)):
        vs1 = Mc.iloc[i].values
        vs2 = Mw.iloc[i].values
        if sum(vs1) < sum(vs2):
            vs1, vs2 = vs2, vs1
        rows_c.append(vs1)
        rows_w.append(vs2)
    m1 = pd.DataFrame(rows_c, index=Mc.index, columns=Mc.columns)
    m2 = pd.DataFrame(rows_w, index=Mc.index, columns=Mc.columns)
    m3 = make_proportion_matrix(m1, m2)
    return m1, m2, m3


def plot_heatmap(M, outfile):
    plt.figure(figsize=(6, 10))
    sns.heatmap(M, cmap="bwr")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()


def get_expected_proportions(Mp):
    return Mp.mean(axis=0)

 
def mask_matrix(Mc, Mw, Mp, expected_proportions):
    Mm = np.zeros(Mc.shape, dtype=np.int)
    for i, sample in enumerate(Mc.index):
        cricks = Mc.iloc[i]
        watsons = Mw.iloc[i]
        proportions = Mp.iloc[i]
        assert len(proportions) == len(expected_proportions)
        nbin = len(cricks)

        status = np.zeros(nbin)
        for j in np.arange(nbin):
            c, w, p, ep = cricks[j], watsons[j], proportions[j], expected_proportions[j]
            if abs(p - ep) < 0.2:
                s = 1
            elif c + w < 10:
                s = 1
            else:
                s = 0
            status[j] = s
        
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
                
        
        if True:
            # Filter too short
            array = list(filter(lambda r: r[1] - r[0] >= 2, array))
            
        if True:
            # Fill small gap
            j = 0
            while j < len(array) - 1:
                if array[j + 1][0] - array[j][1] < max(int(nbin * 0.03), 2):
                    array[j][1] = array[j + 1][1]
                    array.pop(j + 1)
                else:
                    j += 1
                    
        if True:
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
            
        if True:
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
            
        for r in array:
            for j in range(r[0], r[1]):
                Mm[i][j] = 1
    
    Mm = pd.DataFrame(Mm, index=Mc.index, columns=Mc.columns)
    return Mm


def plot_bin_reads(Mc, Mw, Mm, outfile):
    ncol = 4
    nrow = int(len(Mm) / ncol)
    if len(Mm) % ncol > 0:
        nrow += 1
    fig, axs = plt.subplots(nrow, ncol, figsize=(ncol * 3, nrow * 1), sharex=True)
    for i in range(len(Mm)):
        ax = axs[int(i / ncol)][i % ncol]
        plt.sca(ax)
        plt.title(Mm.index.values[i])
        counts1 = Mc.iloc[i]
        counts2 = Mw.iloc[i]
        status = Mm.iloc[i]
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
    plt.savefig(outfile, dpi=300)
    plt.close()


def report_cc_regions(Mm, chrom, length, bin_width, patterns):
    regions = []
    for i, sample in enumerate(Mm.index):
        status = Mm.values[i,:]
        j1 = None
        array = []
        for j, s in enumerate(status):
            if s == 0:
                if j1 is None:
                    continue
                else:
                    array.append([j1, j])
                    j1 = None
            else:
                if j1 is None:
                    j1 = j
                else:
                    continue
        if j1 is not None:
            array.append([j1, len(status)])
        for j1, j2 in array:
            start = j1 * bin_width
            end = min(j2 * bin_width, length)
            r = [chrom, start, end, sample, ".", "+" if patterns[sample] == "CC" else "-"]
            regions.append(r)
    return regions


def process_chromosome(sample_counts, chrom, outdir):
    bin_width = get_bin_width(sample_counts)
    chrom_lengths = get_chrom_lengths(sample_counts)
    length = chrom_lengths[chrom]

    if length < 10000000:
        sys.stderr.write("%s shorter than 10,000,000.\n" % chrom)
        return None

    chrom_outdir = os.path.join(outdir, chrom)
    if not os.path.exists(chrom_outdir):
        os.mkdir(chrom_outdir)

    Mc, Mw, Mp = make_matrix(sample_counts, chrom)
    Mc.to_csv(os.path.join(chrom_outdir, "matrix_of_crick.raw.tsv"), sep="\t")
    Mw.to_csv(os.path.join(chrom_outdir, "matrix_of_watson.raw.tsv"), sep="\t")
    Mp.to_csv(os.path.join(chrom_outdir, "matrix_of_crick_proportion.raw.tsv"), sep="\t")
    plot_heatmap(Mp, os.path.join(chrom_outdir, "heatmap_of_crick_proportion.raw.pdf"))
    
    patterns = predict_chrom_pattern(Mc, Mw)

    Mc, Mw, Mp = cluster_matrix(Mc, Mw, Mp, os.path.join(chrom_outdir, "clustermap.pdf"))
    Mc.to_csv(os.path.join(chrom_outdir, "matrix_of_crick.clustered.tsv"), sep="\t")
    Mw.to_csv(os.path.join(chrom_outdir, "matrix_of_watson.clustered.tsv"), sep="\t")
    Mp.to_csv(os.path.join(chrom_outdir, "matrix_of_crick_proportion.clustered.tsv"), sep="\t")
    plot_heatmap(Mp, os.path.join(chrom_outdir, "heatmap_of_crick_proportion.clustered.pdf"))

    Mc, Mw, Mp = filter_matrix(Mc, Mw, Mp)
    if len(Mc) < 10:
        sys.stderr.write("%s less than 10 samples after filtering.\n" % chrom)
    Mc.to_csv(os.path.join(chrom_outdir, "matrix_of_crick.filtered.tsv"), sep="\t")
    Mw.to_csv(os.path.join(chrom_outdir, "matrix_of_watson.filtered.tsv"), sep="\t")
    Mp.to_csv(os.path.join(chrom_outdir, "matrix_of_crick_proportion.filtered.tsv"), sep="\t")
    plot_heatmap(Mp, os.path.join(chrom_outdir, "heatmap_of_crick_proportion.filtered.pdf"))

    Mc, Mw, Mp = reverse_matrix(Mc, Mw, Mp)
    Mc.to_csv(os.path.join(chrom_outdir, "matrix_of_crick.reversed.tsv"), sep="\t")
    Mw.to_csv(os.path.join(chrom_outdir, "matrix_of_watson.reversed.tsv"), sep="\t")
    Mp.to_csv(os.path.join(chrom_outdir, "matrix_of_crick_proportion.reversed.tsv"), sep="\t")
    plot_heatmap(Mp, os.path.join(chrom_outdir, "heatmap_of_crick_proportion.reversed.pdf"))

    expected_proportions = get_expected_proportions(Mp)

    Mm = mask_matrix(Mc, Mw, Mp, expected_proportions)
    plot_heatmap(Mm, os.path.join(chrom_outdir, "heatmap_of_masked.pdf"))
    Mm.to_csv(os.path.join(chrom_outdir, "matrix_of_masked.tsv"), sep="\t")

    plot_bin_reads(Mc, Mw, Mm, os.path.join(chrom_outdir, "cc_regions.pdf"))

    regions = report_cc_regions(Mm, chrom, length, bin_width, patterns)

    with open(os.path.join(chrom_outdir, "cc_regions.bed"), "w+") as fw:
        for r in sorted(regions):
            fw.write("\t".join(map(str, r)) + "\n")


def fetch_cc_regions(args=None):
    parser = optparse.OptionParser(usage="%prog [options] config.tsv outdir")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1, metavar="INT", 
                      help="")
    options, args = parser.parse_args(args)
    threads = options.threads
    configfile, outdir = args

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    config = load_config(configfile)
    sample_counts = load_read_counts(config)
    chrom_lengths = get_chrom_lengths(sample_counts)
    chroms = list(sorted(chrom_lengths.keys()))

    pool = mp.Pool(threads)
    for chrom in chroms:
        params = (sample_counts, chrom, outdir)
        pool.apply_async(process_chromosome, params)
    pool.close()
    pool.join()

    # merge all CC regions
    paths =  list(glob.glob("%s/*/cc_regions.bed" % outdir))
    if len(paths) > 0:
        cmd = "cat %s | sort -k1,1 -k2,2n > %s" % (" ".join(paths), os.path.join(outdir, "all_cc_regions.bed"))
        subprocess.check_call(cmd, shell=True)
    else:
        with open(os.path.join(outdir, "all_cc_regions.bed"), "w+") as fw:
            pass

    print("All completed!")
        

if __name__ == "__main__":
    fetch_cc_regions()