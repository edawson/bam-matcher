#!/usr/bin/env python
'''
Created on 01/12/2014
@author: paul

Extract SNPs from Thousand Genome VCF file
- exonic variants (using RefSeq or Ensembl annotations)
- randomly select variants across the whole genome (excluding Y and MT)
- has global allele frequency of between 0.55 and 0.65

Thousand Genome 
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/
Download vcf.gz files for all chromosomes and combine into a single file.
These VCF files don't use 'chr' for chromosome names

RefSeq GTF file:
Generate using http://genome.ucsc.edu/cgi-bin/hgTables
This file uses 'chr' for chromosome names

'''
import vcf
import HTSeq
import os
import random
import sys
import subprocess
from argparse import ArgumentParser

# ==============================================================================
def handle_args():
    parser = ArgumentParser(description="Extract a smaller VCF from 1KG \
database for BAM file comparison")
    parser.add_argument("--invcf", "-i", required=True,  
                        help="Input VCF file (should be 1KG vcf or vcf.gz)")
    parser.add_argument("--outvcf", "-o", required=True, help="Output VCF file")
    parser.add_argument("--number_of_snps", "-n", required=False, default=1500, 
                        help="Total number of SNPs to include. Default = 1500")
    parser.add_argument("--AF-range", "-a", required=False, default="0.45-0.55", 
                        help="Range of global allele frequencies (1KG_AF) to \
extract. Default = 0.45-0.55" )
    parser.add_argument("--exon-only", "-x", required=False, action="store_true",
                         help="Select only exon variants")
    parser.add_argument("--refseq", "-R", required=False, 
                        help="RefSeq GTF file, required if you want to \
extract exon-only variants")
    parser.add_argument("--include-Y","-Y", required=False, action="store_true",
                         help="Include SNPs from Y chromosome. Default = False")
    parser.add_argument("--verbose", "-v", required=False, default=False,
                        action="store_true", 
                        help="Turn on verbose status report. Default = False")
    parser.add_argument("--p-threshold-multiplier", "-N", required=False, 
                        default=1.5, type=float, help="P-threshold multiplier. \
Default = 1.5. Increase this value if not enough SNPs are selected.")
    return parser.parse_args()

# 1KG_AA=G;1KG_AC=2;1KG_AF=0.0009;AFR_AF=0.0041;1KG_AN=2184;AVGPOST=0.9999;ERATE=0.0003;LDAF=0.001;RSQ=0.9547;SNPSOURCE=EXOME;THETA=0.0005;VT=SNP
# => 0.0009
def extract_1KG_AF(info_str):
    info_bits = info_str.split(";")
    for bit in info_bits:
        if bit.startswith("1KG_AF="):
            return float(bit.split("=")[1])
    return None
# ==============================================================================

args = handle_args()


# ------------------------------------------------------------------------------
# SET BINARIES AND PATHS
BEDTOOLS = "bedtools"

# ------------------------------------------------------------------------------
# Some variables
invcf_has_chr = False
gtf_has_chr = False
invcf_path = os.path.abspath(args.invcf)
outvcf_path = os.path.abspath(args.outvcf)

# ------------------------------------------------------------------------------
# check input VCF and determine whether it uses chr for chromosome names
if os.path.isfile(args.invcf) == False:
    print "Cannot find input VCF file (%s)" % args.invcf
    sys.exit(1)
vcf_in = vcf.Reader(open(args.invcf, "r"))
variant = vcf_in.next()
invcf_has_chr = variant.CHROM.startswith("chr")

# ------------------------------------------------------------------------------
# parse AF-range
AF_low, AF_high = args.AF_range.split("-")
AF_low = float(AF_low)
AF_high = float(AF_high)

# ------------------------------------------------------------------------------
# exon-only and refseq
if args.exon_only:
    if args.refseq == None:
        print "Extracting exonic variants requires a RefSeq GTF (--refseq/-R)"
        sys.exit(1)
    else:
        # check if RefSeq GTF file exists
        if os.path.isfile(args.refseq) == False:
            print "Cannot find input RefSeq GTF file (%s)" % args.refseq
            sys.exit(1)
        gtf_path = os.path.abspath(args.refseq)
        # determine whether the GTF file has chr or not
        fin = open(args.refseq, "r")
        line = fin.readline().strip("\n")
        gtf_has_chr = line.startswith("chr")

# -----------------------------------------------------------------------------
if args.verbose:
    run_summary = """
RUN SUMMARY
---------------------------------------
Input 1000G VCF file:            %s
 - has 'chr' in chromosome name: %r
Output VCF file:                 %s
Number of SNPs to write:         %d
Global allele frequency range:   %f - %f
""" % (invcf_path, invcf_has_chr, outvcf_path, args.number_of_snps, AF_low, AF_high) 
    print run_summary

    if args.exon_only:
        run_sum2 = """Extracting exonic variants only
RefSeq GTF file:                 %s
 - has 'chr' in chromosome name: %r
----------------------------------------
""" % (gtf_path, gtf_has_chr) 
        print run_sum2
    else:
        print "----------------------------------------"

# ------------------------------------------------------------------------------
# Random selection parameters 

# These numbers were calculated using Nimblegen exome capture bed file, 
# so likely to be a little bit off from refseq CDS regions
#   Total SNPs (non-indel) ~ 38,000,000
#   Total exonic SNPs ~ 900,000
#   Number of exonic SNPs with global AF 0.45-0.55 ~ 7550
#   Number of all SNPs with global AF 0.45-0.55 ~ 500,000

random.seed()
p_thresh = args.number_of_snps/730000.0    # heuristically adjusted
if args.exon_only:
    p_thresh = args.number_of_snps/7550.0

# this is a normalisation step, 
# otherwise, the number of variants selected will be too low
p_thresh = p_thresh * args.p_threshold_multiplier   
if args.verbose:
    print "Random selection probability threshold: ", p_thresh

# this is for stdout printing
backspace_length = len(str(args.number_of_snps)) * 2 + 20
backspace_str = "\b"*backspace_length

# -----------------------------------------------------------------------------
# Extracting variants from VCF file
# write header lines
fin = open(invcf_path, "r")
fout = open(outvcf_path, "w")
for line in fin:
    if line.startswith("#"):
        fout.write(line)
    else:
        break

# if not extracting exon-only variants, we can just load the 1KG vcf file directly
# but if we want to use exon only, need to run bed-intersect first and write to a temporary file
# because pyvcf cannot generated a vcf.Record from a single VCF line alone, it needs to parse the header

if args.exon_only:
    # STILL NEED TO SORT OUT THE STUFF WITH 'chr'
    if args.verbose:
        print "Getting exonic variants, this could take a few minutes..."
    get_vcf_cmd = "awk -F \"\\t\" '$3==\"CDS\" {print $1\"\\t\"($4-1)\"\\t\"$5}' %s | sed 's|^chr||' | %s intersect  -a %s -b -" % (gtf_path, BEDTOOLS, invcf_path)
    get_vcf_proc = subprocess.Popen([get_vcf_cmd], shell=True, stdout=subprocess.PIPE)
    input_stream = get_vcf_proc.stdout
else:
    input_stream = open(invcf_path, "r")

NUCLEOTIDES = { "A":1, "C":1, "G":1, "T":1}
var_written = 0
for line in input_stream:
    if line.startswith("#"):
        continue
    bits = line.strip("\n").split("\t")
    chr_ = bits[0]
    ref_ = bits[3]
    alt_ = bits[4]
    info_str = bits[7]

    # accept only SNPs, no indels
    if ref_ not in NUCLEOTIDES or alt_ not in NUCLEOTIDES:
        continue
    # accept only within specified range of AF values
    global_AF = extract_1KG_AF(info_str)
    if global_AF < AF_low or global_AF > AF_high:
        continue
    # throw the die
    pval = random.random()
    if pval > p_thresh:
        continue
    if args.include_Y == False:
        if chr_ == "Y" or chr_ == "chrY":
            break
    if args.verbose:
        print """die throw successful! (%f < %f)
1KG_AF = %f
%s
""" % (pval, p_thresh, global_AF, line.strip("\n")) 
    else:
        sys.stdout.write("%s%d/%d (%s)" % (backspace_str, var_written, args.number_of_snps, chr_))
        sys.stdout.flush()
    fout.write(line)
    var_written += 1
    if var_written > args.number_of_snps:
        break
fout.close()    

if args.verbose:
    print "\n\nFINISHED!\n"
    






































