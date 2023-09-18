#!/usr/bin/env python
import os
import re
import optparse
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "arial"
import matplotlib.pyplot as plt

usage = """

    sstools PlotBinRead [options] <input.tsv> <outfile.pdf>
"""

def plot_bin_read(args):
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-t", "--trim", dest="trim", action="store_true", default=False)
    options, args = parser.parse_args(args)
    
    trim = options.trim
    tsvfile, pdffile = args
 
    dat = pd.read_csv(tsvfile, sep="\t")
    dat = dat[[re.match("^chr([0-9]+|[XY])$", c) is not None for c in dat["Chrom"]]]
    
    chroms = []
    for c in dat["Chrom"]:
        if c not in chroms:
            chroms.append(c)
            
    ymax = max(Counter(dat["Chrom"]).values())

    ncol = 12
    nrow = int(len(chroms) / ncol)
    if len(chroms) % ncol > 0:
        nrow += 1

    xmax = max(dat["Crick"])
    xmax = max(xmax, max(dat["Watson"]))
    xmax = max(xmax * 1.2, 5)
    if trim:
        vs = []
        vs.extend(dat["Crick"])
        vs.extend(dat["Watson"])
        xmax = np.mean(vs) + np.std(vs) * 3
        
    fig, axs = plt.subplots(nrow, ncol, figsize=(ncol * 2, nrow * 4), sharex=True, sharey=True)
    plt.suptitle(os.path.basename(tsvfile))
    for i, chrom in enumerate(chroms):
        if nrow == 1:
            ax = axs[i]
        else:
            ax = axs[int(i / ncol)][i % ncol]
        plt.sca(ax)
        d = dat[dat["Chrom"] == chrom]
        c = np.mean(d["Crick"])
        w = np.mean(d["Watson"])
        b = min(c, w) / (c + w) * 100
        plt.title("%.2f,%.2f(%.2f%%)" % (w, c, b))
        ys = np.arange(len(d)) + 0.5
        plt.barh(ys, d["Crick"], height=1)
        plt.barh(ys, d["Crick.P"], height=1, color="blue")
        plt.barh(ys, d["Crick.M"], left=d["Crick.P"], height=1, color="red")
        plt.barh(ys, -d["Watson"], height=1)
        plt.barh(ys, -d["Watson.M"], height=1, color="red")
        plt.barh(ys, -d["Watson.P"], left=-d["Watson.M"], height=1, color="blue")
        plt.xlim(-xmax, xmax)
        plt.ylim(0, ymax)
        plt.xlabel(chrom)
        plt.yticks([])
        for loc in ["top", "left", "right"]:
            ax.spines[loc].set_visible(False)
    plt.tight_layout()
    plt.savefig(pdffile, dpi=300)
    plt.close()
