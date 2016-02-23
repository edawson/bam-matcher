'''
Created on 23/22/2016

@author: Paul Wang (ppswang@gmail.com)

Moved methods from main script to here


'''

import gzip
import vcf

# Sort VCF entries by chromosome order ordered by the reference index
# This is somewhat fugly... requires some unix commands
def sort_vcf_by_chrom_order(invcf, outvcf, ref_index):
    chrom_order = []
    fin = open(ref_index, "r")
    for line in fin:
        chrom_order.append(line.strip().split()[0])

    fin = open(invcf, "r")
    fout = open(outvcf, "w")
    # write header
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
        else:
            break
    # write out entries in for each chromosome in order, and sorted by position
    for chrom in chrom_order:
        grep_str = "grep -v ^# '%s' | egrep ^%s[[:space:]] | sort -k2n" % \
                    (invcf, chrom)
        grep_cmd = subprocess.Popen([grep_str], stdout=subprocess.PIPE,
                                    shell=True)
        for line in grep_cmd.stdout:
            fout.write(line)
    # DONE
    fout.close()
    return


# Convert variants VCF file to intervals for variant callers
def convert_vcf_to_intervals(invcf, output, window, ntries, format="gatk",
                             filtering=True):
    vcf_read = vcf.Reader(open(invcf, "r"))
    fout = open(output, "w")
    n_written = 0

    for var in vcf_read:
        # Filtering variants
        # 1) 1KG_AF between 0.45 - 0.55
        # 2) SNPs only, no indels
        if filtering:
            if "1KG_AF" in var.INFO and (var.INFO["1KG_AF"] < 0.45 or
                                         var.INFO["1KG_AF"] > 0.55):
                continue
            if var.var_type != "snp":
                continue
        # intervals format
        if format == "gatk" or format == "varscan":
            start_pos = var.POS - window
            end_pos   = start_pos + window
            fout.write("%s:%d-%d\n" % (var.CHROM, start_pos, end_pos))
        # BED format
        elif format == "bed" or format == "freebayes":
            start_pos = var.POS - window - 1
            end_pos   = start_pos + window + 1
            fout.write("%s\t%d\t%d\n" % (var.CHROM, start_pos, end_pos))
        n_written += 1
        if ntries != None:
            if n_written >= ntries:
                break
    fout.close()
    return

# Convert VCF to TSV file
# This is to replace GATK -T VariantsToTable, which is too slow...
def VCFtoTSV(invcf, outtsv, caller):
    fout = open(outtsv, "w")
    vcf_in = vcf.Reader(open(invcf, "r"))

    if caller == "gatk" or caller == "varscan":
        fields_to_extract = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AD", "GT"]
    elif caller == "freebayes":
        fields_to_extract = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AO", "GT"]
    fout.write("%s\n" % "\t".join(fields_to_extract))
    for var in vcf_in:
        if var.var_type != "snp":
            continue
        chrom_ = var.CHROM
        pos_   = str(var.POS)
        ref_   = var.REF
        alt_   = var.ALT[0]
        alt_str = "."
        if alt_ != None:
            alt_str = alt_.sequence
        qual_ = str(var.QUAL)
        dp_   = str(var.INFO["DP"])

        ad_or_ao = "NA"
        ad_str = "NA"
        gt_ = "NA"
        if var.samples[0].called:
            if caller == "gatk" or caller == "varscan":
                ad_ = var.samples[0]["AD"]
                for a_ in ad_:
                    ad_str += ",%d" % a_
                ad_str = ad_str[1:]
            else:
                ad_str = str(var.samples[0]["AO"])
            gt_ = var.samples[0].gt_bases

        data_bits = [chrom_, pos_, ref_, alt_str, qual_, dp_, ad_str, gt_]
        fout.write("%s\n" % "\t".join(data_bits))
    fout.close()


# GENOTYPE COMPARISON FUNCTIONS
def is_hom(gt):
    gt_ = gt.split("/")
    if gt_[0] == gt_[1]:
        return True
    else:
        return False

def same_gt(gt1, gt2):
    if gt1 == gt2:
        return True
    else:
        gt1_ = sorted(gt1.split("/"))
        gt2_ = sorted(gt2.split("/"))
        if gt1_ == gt2_:
            return True
        else:
            return False

def is_subset(hom_gt, het_gt):
    gt_hom = hom_gt.split("/")[0]
    if gt_hom in het_gt:
        return True
    else:
        return False

# Get list of chromosome names from the BAM file

def get_chrom_names(bam_file):
    chrom_list = []
    inbam = gzip.open(bam_file, "r")
    for line in inbam:
        if line.startswith("@") == False and line.startswith("BAM") == False:
            break
        elif line.startswith("@SQ"):
            bits = line.strip().split("\t")
            for chunk_ in bits:
                if chunk_.startswith("SN:"):
                    chrom_list.append(chunk_[3:])
    return chrom_list
