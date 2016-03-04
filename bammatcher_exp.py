'''
Created on 4/3/2016

@author: Paul Wang (ppswang@gmail.com)

Experimental methods for bam-matcher

'''

import math
from matplotlib import pyplot as plt





def count_AF_bins(gt_data, n_bins):
    freq_list = []
    fin = open(gt_data, "r")
    fin.readline() # ignore first line
    for line in fin:
        bits = line.strip("\n").split("\t")
        dp_ = int(bits[5])
        vdp_ = bits[6]

        if "," in vdp_: # then this is probably GATK output
            if vdp_.startswith("[") == False:
                vdp_ = int(vdp_.split(",")[2])
            else:
                vdp_ = vdp_.replace("[", "").replace("]", "").split(",")[0]
        elif vdp_ == "NA":
            vdp_ = 0
        else:
            vdp_ = int(vdp_)
        if dp_ >0:
            vaf_ = float(vdp_)/dp_
            freq_list.append(vaf_)

    freq_bins = [0]*n_bins
    window = 1.0/n_bins
    for fq_ in freq_list:
        bin_ = int(math.floor(fq_/window))
        if bin_ == n_bins:
            bin_ = n_bins-1
        freq_bins[bin_] += 1
    return freq_list, freq_bins

def plot_VAF(freqlist, bin_count, fig_out, bam_path):
    plt.clf()
    plt.hist(freqlist, bin_count)
    plt.xlabel("Variant allele frequency", fontsize=8)
    plt.ylabel("Locus count", fontsize=8)
    plt.xlim(0, 1)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.title(bam_path, fontsize=8)
    fig = plt.gcf()
    fig.set_size_inches(4,4)
    fig.savefig(fig_out, dpi=200)

# def plot_ascii_VAF(freqlist)
