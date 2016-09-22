'''
Created on 4/3/2016

@author: Paul Wang (ppswang@gmail.com)

Experimental methods for bam-matcher

'''

import math
import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot as plt
import gzip

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






# ===================================
# Using SNP array data

# functions

# check file type

def check_file_type(filepath):
    # if ends in .bam, assume is BAM file
    if filepath.endswith(".bam"):
        return "bam"

    # otherwise, read a file lines and check
    if filepath.endswith(".gz"):
        fin = gzip.open(filepath, "r")
    else:
        fin = open(filepath, "r")

    line1 = fin.readline()
    line2 = fin.readline()
    # SNP array data
    if line1.strip() == "[Header]" and line2.startswith("GSGT Version"):
        return "snp"
    else:
        return "unknown"


# convert SNP array input data to bam-matcher genotype data format

# INPUT

# [Header]
# GSGT Version	1.9.4
# Processing Date	9/22/2016 10:40 AM
# Content		HumanCytoSNP-12v2-1_L.bpm
# Num SNPs	294602
# Total SNPs	294602
# Num Samples	60
# Total Samples	60
# File 	1 of 60
# [Data]
# Sample ID	SNP Name	SNP Index	Chr	Position	SNP	Allele1 - Top	Allele2 - Top	Allele1 - Plus	Allele2 - Plus	Allele1 - Forward	Allele2 - Forward	Allele1 - Design	Allele2 - Design	Allele1 - AB	Allele2 - AB	GC Score	ILMN Strand	Customer Strand
# 3999565005_R01C01	cnvi0111185	1	13	100631779	[N/A]	-	-	-	-	-	-	-	-	-	-	0.0000	PLUS	PLUS
# 3999565005_R01C01	cnvi0111186	2	2	177058958	[N/A]	-	-	-	-	-	-	-	-	-	-	0.0000	PLUS	PLUS
# 3999565005_R01C01	cnvi0111187	3	17	35295593	[N/A]	-	-	-	-	-	-	-	-	-	-	0.0000	PLUS	PLUS
# 3999565005_R01C01       rs1000016       740     2       235690982       [T/C]   A       A       -       -       A       A       T       T       A       A       0.9483  BOT     TOP
# 3999565005_R01C01       rs10000180      741     4       83899764        [A/G]   G       G       -       -       G       G       G       G       B       B       0.9423  TOP     TOP
# 3999565005_R01C01       rs10000272      742     4       189690383       [A/G]   A       A       -       -       T       T       A       A       A       A       0.9556  TOP     BOT
# 3999565005_R01C01       rs1000031       743     18      46361441        [T/C]   A       G       -       -       T       C       T       C       A       B       0.9623  BOT     BOT
# 3999565005_R01C01       rs10000432      744     4       47511781        [A/G]   A       G       -       -       T       C       A       G       A       B       0.8771  TOP     BOT


# OUTPUT
# CHROM	POS	REF	ALT	QUAL	DP	AO	GT
# 1	1650807	T	C	1756.18	146	72	T/C
# 1	2488153	A	G	5984.89	486	250	A/G
# 1	7887248	G	.	0	7	0	G/G
# 1	9097241	C	.	0	446	0	C/C


def convert_snp_to_genotype(snp_data_file, gt_data_file):
    if snp_data_file.endswith(".gz"):
        fin = gzip.open(snp_data_file, "r")
    else:
        fin = open(snp_data_file, "r")

    fout = open(gt_data_file+".unsorted", "w")
    fields_to_extract = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AO", "GT"]
    fout.write("%s\n" % "\t".join(fields_to_extract))

    # columns to use:
    # SNP:    this is REF/ALT, e.g. [T/C]
    # Allele1 - AB, Allele2 - AB    - these are genotype calls in A/B form
    #    where A = REF, B = ALT

    snp_data_cols = ["Chr", "Position", "SNP", "Allele1 - AB", "Allele2 - AB"]
    snp_data_col_indices = {}

    # read header and get column positions
    for line in fin:
        # read the line after [Data]
        if line.startswith("[Data]"):
            snp_columns_line = fin.readline()
            break
    snp_columns = snp_columns_line.strip().split("\t")
    for col_idx, col in enumerate(snp_columns):
        if col in snp_data_cols:
            snp_data_col_indices[col] = col_idx

    # verify that the SNP data file contains all required columns
    for col in snp_data_cols:
        if col not in snp_data_col_indices:
            print "SNP DATA ERROR:\nInput SNP array data file (%s) is missing expected column (%s)" % (snp_data_file, col)
            exit(1)

    # extract required information
    for line in fin:
        snp_bits = line.strip().split("\t")
        CHROM_    = snp_bits[  snp_data_col_indices["Chr"]  ]
        POSITION_ = snp_bits[  snp_data_col_indices["Position"]  ]
        REF_ALT   = snp_bits[  snp_data_col_indices["SNP"]  ][1:-1]
        ALLELE1   = snp_bits[  snp_data_col_indices["Allele1 - AB"]  ]
        ALLELE2   = snp_bits[  snp_data_col_indices["Allele2 - AB"]  ]

        if REF_ALT == "N/A":
            continue
        if ALLELE1 == "-" or ALLELE2 == "-":
            continue
        REF_, ALT_ = REF_ALT.split("/")
        GT_ = "%s/%s" % (ALLELE1.replace("A", REF_).replace("B", ALT_), ALLELE2.replace("A", REF_).replace("B", ALT_)    )
        fout.write("\t".join([ CHROM_, POSITION_, REF_, ALT_,  "1000", "1000", "500", GT_ ])  + "\n")
    fout.close()

    # but this file is unsorted
























#
