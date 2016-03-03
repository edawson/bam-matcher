'''
Created on 23/22/2016

@author: Paul Wang (ppswang@gmail.com)

Moved methods from main script to here


'''

import gzip
import vcf
import os
import ConfigParser
import subprocess

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
def convert_vcf_to_intervals(invcf, output, window, ntries, format="gatk", cmap=None):
    vcf_read = vcf.Reader(open(invcf, "r"))
    fout = open(output, "w")
    n_written = 0

    for var in vcf_read:
        # intervals format
        chrom_ = var.CHROM
        if cmap != None:
            if var.CHROM not in cmap:
                continue
            else:
                chrom_ = cmap[var.CHROM]

        if format == "gatk" or format == "varscan":
            start_pos = var.POS - window
            end_pos   = start_pos + window
            fout.write("%s:%d-%d\n" % (chrom_, start_pos, end_pos))
            n_written += 1

        # BED format
        elif format == "bed" or format == "freebayes":
            start_pos = var.POS - window - 1
            end_pos   = start_pos + window + 1
            fout.write("%s\t%d\t%d\n" % (chrom_, start_pos, end_pos))
            n_written += 1

        if ntries >0:
            if n_written >= ntries:
                break

    fout.close()
    return

# Convert VCF to TSV file
# This is to replace GATK -T VariantsToTable, which is too slow...
def VCFtoTSV(invcf, outtsv, caller):
    fout = open(outtsv, "w")
    vcf_in = vcf.Reader(open(invcf, "r"))
    var_ct = 0
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
        if caller == "varscan":
            qual_ = str(var.samples[0]["GQ"])
        if caller == "freebayes" or caller == "gatk":
            dp_ = str(var.INFO["DP"])
        else:
            dp_ = str(var.INFO["ADP"])

        ad_or_ao = "NA"
        ad_str = "NA"
        gt_ = "NA"
        if var.samples[0].called:
            if caller == "gatk":
                ad_ = var.samples[0]["AD"]
                for a_ in ad_:
                    ad_str += ",%d" % a_
                ad_str = ad_str[1:]
            elif caller == "varscan":
                ad_str = str(var.samples[0]["AD"])
            else:
                ad_str = str(var.samples[0]["AO"])
            gt_ = var.samples[0].gt_bases

        data_bits = [chrom_, pos_, ref_, alt_str, qual_, dp_, ad_str, gt_]
        fout.write("%s\n" % "\t".join(data_bits))
        var_ct += 1
    fout.close()
    return var_ct

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
def get_chrom_names_from_BAM(bam_file):
    chrom_list = []
    inbam = gzip.open(bam_file, "r")
    for line in get_bam_header(bam_file):
        if line.startswith("@SQ"):
            bits = line.strip().split("\t")
            for chunk_ in bits:
                if chunk_.startswith("SN:"):
                    chrom_list.append(chunk_[3:])
    return chrom_list



def check_fasta_index(fastafile):
    ref_idx = fastafile + ".fai"
    if not check_file_read(ref_idx, "", None, silent=True):
        print """%s
Specified reference FASTA file (%s) needs to be indexed by samtools:

samtools faidx %s
""" % (FILE_ERROR, fastafile, fastafile)
        return False
    else:
        return True


def get_chrom_names_from_REF(ref_fasta):
    # first check that the FASTA file is indexed
    if not check_fasta_index(ref_fasta):
        exit(1)
    ref_idx = ref_fasta + ".fai"

    chrom_list = []
    fin = open(ref_idx, "r")
    for line in fin:
        if line.strip() == "":
            continue
        chrom_list.append(line.strip().split()[0])
    return chrom_list

def get_chrom_data_from_map(chrom_map_file):
    chrom_ct = 0
    default_chroms = []
    alternate_chroms = []
    def_to_alt = {}
    alt_to_def = {}
    fin = open(chrom_map_file, 'r')
    header = fin.readline()
    for line in fin:
        if line.strip() == "":
            continue
        chrom_ct += 1
        bits = line.strip().split()
        default_chroms.append(bits[0])
        alternate_chroms.append(bits[1])
        def_to_alt[bits[0]] = bits[1]
        alt_to_def[bits[1]] = bits[0]
    return default_chroms, alternate_chroms, def_to_alt, alt_to_def








def fetch_config_value(config_obj, config_section, config_keyword):
    try:
        value_ = config_obj.get(config_section, config_keyword)
    except ConfigParser.NoOptionError as e:
        print """%s

Missing keyword '%s' in configuration file.
Do not remove or comment out keywords in the config file.
Just leave blank value if not using the parameter.

Python error message:
%s
""" % (CONFIG_ERROR, config_keyword, e)
        exit(1)
    except Exception as e:
        print """%s

Unknown config error.

Python error message:
%s
""" (CONFIG_ERROR, e)
        exit(1)
    # if all good, return value_
    return value_







def check_caller(caller, caller_binary, JAVA="java", verbose=False, SAMTL="samtools", logfile=None):
    logwrite = open(logfile, "w")
    check_exitcode = 0
    caller_cmd = []

    # ------------------------------------------
    # checking GATK
    if caller == "gatk":
        if caller_binary == "":
            print """%s
GATK path was not specified.
Do this in the configuration file""" % CONFIG_ERROR
            exit(1)
        if os.access(caller_binary, os.R_OK) == False:
            print """%s
Cannot access GATK jar file (%s).
It is either missing or not readable.""" % (CONFIG_ERROR, caller_binary)
            exit(1)
        caller_cmd = [JAVA, "-jar", caller_binary, "-version"]

    # -------------------------------------------
    # Checking Freebayes
    elif caller == "freebayes":
        if caller_binary == "":
            print """%s
Freebayes path was not specified.
Do this in the configuration file""" % CONFIG_ERROR
            exit(1)
        caller_cmd = [caller_binary, "--version"]
    # ------------------------------------------
    # checking VARSCAN
    elif caller == "varscan":
        # ===================================
        # Testing samtools
        if SAMTL == "":
            print """%s
SAMtools path was not specified.
Do this in the configuration file""" % CONFIG_ERROR
            exit(1)

        sam_cmd = [SAMTL, "--version"]
        try:
            sam_proc = subprocess.Popen(sam_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        except Exception as e:
            # for line in sam_proc.stdout:
            #     print line.strip()
            print """%s
Samtools error. Please check samtools path and installation.

Command tested:
%s --version

Python error msg:
%s

""" % (CALLER_ERROR, SAMTL, e)
            exit(1)

        # ===================================
        # Testing VarScan itself
        if caller_binary == "":
            print """%s
VarScan2 path was not specified.
Do this in the configuration file""" % CONFIG_ERROR
            exit(1)
        if os.access(caller_binary, os.R_OK) == False:
            print """%s
Cannot access VarScan2 jar file (%s).
It is either missing or not readable.""" % (CONFIG_ERROR, caller_binary)
            exit(1)
        caller_cmd = [JAVA, "-jar", caller_binary]

    # --------------------------------------------

    try:
        caller_proc = subprocess.check_call(caller_cmd, stderr=logwrite, stdout=logwrite)
    except Exception as e:
        logwrite.close()
        print """%s
Caller check failed.

Command tested:
%s

Python error msg:
%s

See log file for system error msg:
%s

Please check the caller command or path to the binary.

""" % (CALLER_ERROR, " ".join(caller_cmd), e, logfile)
        exit(1)

    # --------------------------------------------
    if verbose:
        fin = open(logfile, "r")
        version_line = fin.readline()
        fin.close()

        if caller == "gatk":
            print "GATK version", version_line
        elif caller == "freebayes":
            print "Freebayes", version_line
        elif caller == "varscan":
            print "\n"+version_line

    # --------------------------------------------
    return True









def get_bam_header(bam_file):
    header_lines = []
    fin = gzip.open(bam_file, "r")
    firstline = fin.readline().strip()
    header_lines.append("@HD"+firstline.split("@HD")[1])
    for line in fin:
        if line.startswith("@") == False:
            break
        else:
            header_lines.append(line.strip())
    return header_lines


def check_file_read(file_path, file_object_name, error_type, silent=False):
    if os.access(file_path, os.R_OK) == False:
        if not silent:
            print """%s
Cannot access %s file (%s).
It either does not exist or is not readable.
""" % (error_type, file_object_name, file_path)
        return False
    else:
        return True


def check_file_write(file_path, file_object_name, error_type, silent=False):
    if os.access(file_path, os.W_OK) == False:
        if not silent:
            print """%s
Specified path for %s (%s) is not writable.
""" % (error_type, file_object_name, file_path)
        return False
    else:
        return True



def print_caller_failure_message(call_cmd, logfile_):
    print """%s
Variant calling failed.

Caller command:
%s

Check caller log: %s

""" % (CALLER_ERROR, call_cmd, logfile_)

















CONFIG_TEMPLATE_STR = """# BAM-matcher configuration file
# If not setting a specific parameter, just leave it blank, rather than deleting or commenting out the line
# Missing parameter keywords will generate errors

[VariantCallers]
# file paths to variant callers and other binaries

# This is the default caller to use
caller:    gatk

# These are paths (or commands) to the caller executables
# full paths is always required for *.jar files (GATK and VarScan2)
# sometime you may need to specify full path to the binary (for freebayes, samtools and java)
GATK:      GenomeAnalysisTK.jar
freebayes: freebayes
samtools:  samtools
varscan:   VarScan.jar
java:      java

[ScriptOptions]
DP_threshold:   15
number_of_SNPs:

# fast_freebayes enables --targets option for Freebayes, faster but more prone to Freebayes errors
# set to False will use --region, each variant is called separately
fast_freebayes: True

VCF_file: variants.vcf

[VariantCallerParameters]
# GATK memory usage in GB
GATK_MEM: 4

# GATK threads (-nt)
GATK_nt:  1

# VarScan memory usage in GB
VARSCAN_MEM: 4

[GenomeReference]
# default reference fasta file
REFERENCE: hg19.fasta
REF_ALTERNATE:
# CHROM_MAP is required if using two different genome references that have different (but compatible) chromosome names
# this is mainly to deal with the hg19 "chr" issue
CHROM_MAP:

[BatchOperations]
# you MUST specify a cache directory, this directory should be read/write-able by all BAM-matcher users
CACHE_DIR:

[Miscellaneous]
"""


















#-------------------------------------------------------------------------------
# Error messages
CONFIG_ERROR = """+--------------+
| CONFIG ERROR |
+--------------+"""

FILE_ERROR = """+------------+
| FILE ERROR |
+------------+"""

CALLER_ERROR = """+--------------+
| CALLER ERROR |
+--------------+"""

WARNING_MSG = """+---------+
| WARNING |
+---------+"""

ARGUMENT_ERROR = """+----------------+
| ARGUMENT ERROR |
+----------------+"""

CALLER_ERROR = """+-----------------------+
| VARIANT CALLING ERROR |
+-----------------------+"""



# --about-alternate-ref message

ABOUT_ALTERNATE_REF_MSG = """
If input BAM files are mapped with different but compatible genome references
(which may also have different chromosome names), user can also supply an
alternate genome reference file (--alternate-ref/-A) either at run time or
in the configuration file. However, a chromosome map file (--chromsome-map/-M)
is also required when doing so. This file is required to tell BAM-matcher what
are the matching chromosome names in the two different genome references.
The format is two columns of chromosome names, one column for each genome
reference, with each line representing an unique chromosome entry identifiable
in each genome reference.

This function is mainly provided to address the problem with the two major
versions of hg19 being used, where one has "chr" and the other without.

Format for chromosome map:

DEFAULT   ALTERNATE
1           chr1
11          chr11
2           chr2
3           chr3
...etc

- fields can be separated by tab or spaces.
- header is necessary,
- chromosome names in 'DEFAULT' column  must match the chromosome names in the
  genome reference file, and 'ALTERNATE' match the alternate reference

Notes:
1. The input variants VCF file should be referencing the default genome reference
(--reference/-R). Efforts have been made to resolve the hg19 "chr"
issue, however, it is unclear whether it is sufficient for other types of
compatible-but-not-quite-the-same genome references.

2. When using alternate reference, only variants which are found in chromosomes
listed in the chromosome map will be used for comparison. For example, when
working with hg19, minor contigs and the mitochondrial chromosome (chrM/MT)
should not be included in the chromosome map, as they are different between the
two versions.

3. Cached genotype data will store genomic positions using the format in the
default reference genome (--reference/-R).

"""
