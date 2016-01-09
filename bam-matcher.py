#!/usr/bin/env python
'''
Created on 09/12/2014

@author: Paul Wang (ppswang@gmail.com)

This is the pretty-formatted version for general release.

Still need to:
- improve confliction-resolution, error-checking, etc
- have better output reporting and more report modes (console-only, and HTML output)
- need to deal with space in file paths ??? (not sure if this is resolved now)
- add batch mode-friendly stuff, e.g. simpler TSV output, cached genotype data

'''
from argparse import ArgumentParser
# from uuid import uuid4 as random_file_name
import os
import subprocess
import sys
import vcf
import HTSeq
import random
import string
import ConfigParser
import shutil
from Cheetah.Template import Template
from hashlib import md5
# from fisher import pvalue

#============================================================================
# Methods
#============================================================================
def handle_args():
    parser = ArgumentParser(description="Compare two BAM files to see if \
    they are from the same samples, using frequently occuring SNPs \
    reported in the 1000Genome database")

    # these are always required
    parser_grp1 = parser.add_argument_group("REQUIRED")
    parser_grp1.add_argument("--bam1",           "-B1", required=True,
                             help="First BAM file")
    parser_grp1.add_argument("--bam2",           "-B2", required=True,
                             help="Second BAM file")

    # if not specified, will always look into the same location
    parser_grp2 = parser.add_argument_group("CONFIGURATION")
    parser_grp2.add_argument("--config",          "-c",  required=False,
                             help="Specify configuration file (default = \
                             /dir/where/script/is/located/bam-matcher.conf)")
    parser_grp2.add_argument("--generate-config", "-G",  required=False,
                             help="Specify where to generate configuration \
                             file template")

    # Output needs to be specified at command line as it's usually per-run
    parser_grp3 = parser.add_argument_group("OUTPUT REPORT")
    parser_grp3.add_argument("--output",         "-o",  required=False,
                             help="Specify output report path (default = \
                             /current/dir/bam_matcher.SUBFIX)")
    parser_grp3.add_argument("--short-output",   "-so",  required=False,
                             default = False, action="store_true",
                             help="Short output mode (tab-separated).")
    parser_grp3.add_argument("--html",           "-H",  required=False,
                             action="store_true", help="Enable HTML output. HTML file name = report + '.html'")
    parser_grp3.add_argument("--no-report",      "-n",  required=False,
                             action="store_true", help="Don't write output to file. Results output to command line only.")
    parser_grp3.add_argument("--scratch-dir",    "-s",  required=False,
                             help="Scratch directory for temporary files. If not specified, the report output directory will be used (default = /tmp/[random_string])")

    # Usually defined in
    parser_grp4 = parser.add_argument_group("VARIANTS")
    parser_grp4.add_argument("--vcf",            "-V",  required=False,
                             help="VCF file containing SNPs to check (default can be specified in config file instead)")
    parser_grp4.add_argument("--filter-vcf",      "-FT", required=False,
                             action="store_true", help="Enable filtering of the input VCF file")

    # Specifying these values here will override config values
    parser_grp5 = parser.add_argument_group("CALLERS AND SETTINGS (will override config values)")
    parser_grp5.add_argument("--caller",         "-CL", required=False,
                             default="gatk",
                             choices=('gatk', 'freebayes', 'varscan'),
                             help="Specify which caller to use (default = 'gatk')")
    parser_grp5.add_argument("--dp-threshold",   "-DP", required=False,
                             type=int, help="Minimum required depth for comparing variants")
    parser_grp5.add_argument("--number_of_snps", "-N", required=False,
                             type=int, help="Number of SNPs to compare.")
    parser_grp5.add_argument("--fastfreebayes",  "-FF", required=False,
                             action="store_true", help="Use --targets option for Freebayes.")
    parser_grp5.add_argument("--gatk-mem-gb" ,   "-GM", required=False,
                             type=int, help="Specify Java heap size for GATK (GB, int)")
    parser_grp5.add_argument("--gatk-nt" ,       "-GT", required=False,
                             type=int, help="Specify number of threads for GATK UnifiedGenotyper (-nt option)")
    parser_grp5.add_argument("--varscan-mem-gb", "-VM", required=False,
                             type=int, help="Specify Java heap size for VarScan2 (GB, int)")


    # overriding reference matching
    parser_grp6 = parser.add_argument_group("REFERENCES")
    parser_grp6.add_argument("--reference",      "-R",  required=False,
                             help="Default reference fasta file. Needs to be indexed with samtools faidx")
    parser_grp6.add_argument("--ref_noChr",      "-Rn", required=False,
                             help="Reference fasta file, no 'chr' in chromosome names. Needs to be indexed with samtools faidx")
    parser_grp6.add_argument("--ref_wChr",       "-Rw",  required=False,
                             help="Reference fasta file, has 'chr' in chromosome names. Needs to be indexed with samtools faidx")
    parser_grp6.add_argument("--bam1-reference", "-B1R", required=False,
                             help="Reference fasta file for BAM1. Requires --bam2-reference/-B2R, overrides other settings")
    parser_grp6.add_argument("--bam2-reference", "-B2R", required=False,
                             help="Reference fasta file for BAM2. Requires --bam1-reference/-B1R, overrides other settings")

    # for batch operations
    parser_grp7 = parser.add_argument_group("BATCH OPERATIONS")
    parser_grp7.add_argument("--do-not-cache",    "-NC", required=False,
                             default=False, action="store_true",
                             help="Do not keep variant-calling output for future comparison. By default (False) data is written to /bam/filepath/without/dotbam.GT_compare_data")
    parser_grp7.add_argument("--recalculate",    "-RC", required=False,
                             default=False, action="store_true",
                             help="Don't use cached variant calling data, redo variant-calling. Will overwrite cached data unless told not to (-NC)")
    parser_grp7.add_argument("--cache-dir",      "-CD", required=False,
                             help="Specify directory for cached data. Overrides configuration")

    # optional, not in config
    parser.add_argument("--debug", "-d", required=False, action="store_true",
                        help="Debug mode. Temporary files are not removed")
    parser.add_argument("--verbose", "-v", required=False, action="store_true",
                        help="Verbose reporting. Default = False")
    return parser.parse_args()

#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
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

#===============================================================================
# End of methods
#===============================================================================











#===============================================================================
# Parsing configuration file and command line options
#===============================================================================



#-------------------------------------------------------------------------------
# Before calling handle_args(), need to check if --generate-config was called
GENERATE_CONFIG_SYMS = ["--generate-config", "--generate_config", "-G"]
for idx, arg in enumerate(sys.argv):
    if arg in GENERATE_CONFIG_SYMS:
        try:
            config_template_output = os.path.abspath(sys.argv[idx+1])
        except:
            config_template_output = os.path.abspath("bam-matcher.conf.template")
        print """
====================================
Generating configuration file template
====================================

Config template will be written to %s

""" % config_template_output

        # make sure that it doesn't overwrite anything!
        if os.path.isfile(config_template_output):
            print "The specified path ('%s') for config template exists already." % config_template_output
            print "Write to another file."
            sys.exit(1)

        fout = open(config_template_output, "w")
        fout.write("""[VariantCallers]
# file paths to variant callers and other binaries
GATK:      GenomeAnalysisTK.jar
freebayes: freebayes
samtools:  samtools
varscan:   VarScan.jar
java:      java

[ScriptOptions]
DP_threshold:   15
filter_VCF:      False
number_of_SNPs: 1500
# enable --targets option for Freebayes, faster but more prone to Freebayes errors
# set to False will use --region, each variant is called separately
fast_freebayes: True
VCF_file: variants_noX.vcf

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

# Reference fasta file, with no chr in chromosome name (e.g. Broad19.fasta)
REF_noChr: Broad19.fasta

# Reference fasta file with 'chr' in chromosome names
REF_wChr:  genome.fa

[BatchOperations]
CACHE_DIR:  cache_dir

[Miscellaneous]
""")
        fout.close()
        exit()



for arg in enumerate(sys.argv):
    if arg == "-H" or arg == "-html":
        print """

===================================
HTML output is not implemented yet.
===================================

"""


# okay now to parse the arguments
args = handle_args()

#-------------------------------------------------------------------------------
# some random stuff
random.seed()
# vcf_random_str = str(random_file_name())
random_str = ""
for _ in range(12):
    random_str += random.choice(string.ascii_uppercase + string.digits)
temp_files = []

#-------------------------------------------------------------------------------
# Error messages
CONFIG_ERROR = """
+--------------+
| CONFIG ERROR |
+--------------+"""

FILE_ERROR = """
+------------+
| FILE ERROR |
+------------+"""

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

#-------------------------------------------------------------------------------
# Read in configuration file

config_file = ""
if args.config == None:
    config_file = os.path.join(SCRIPT_DIR, "bam-matcher.conf")
else:
    config_file = args.config

# Test if the config file exists
if os.access(config_file, os.R_OK) == False:
    print """%s
The config file '%s' either does not exist or is not readable
"""% (CONFIG_ERROR, os.path.abspath(config_file))
    sys.exit(1)

# Load config and check if all necessary bits are available
config = ConfigParser.ConfigParser()
config.read(config_file)

# are all the sections there?
config_sections = config.sections()
REQUIRED_CONFIG_SECTIONS = ["GenomeReference", "VariantCallerParameters",
                            "ScriptOptions", "VariantCallers", "Miscellaneous",
                            "BatchOperations"]
for sect in REQUIRED_CONFIG_SECTIONS:
    if sect not in config_sections:
        print """%s
Missing required section in config file: %s
""" % (CONFIG_ERROR, sect)
        sys.exit(1)

#-------------------------------------------------------------------------------
# setting variables using the config file
GATK      = config.get("VariantCallers", "GATK")
FREEBAYES = config.get("VariantCallers", "freebayes")
SAMTOOLS  = config.get("VariantCallers", "samtools")
VARSCAN   = config.get("VariantCallers", "varscan")
JAVA      = config.get("VariantCallers", "java")

DP_THRESH      = config.get("ScriptOptions", "DP_threshold")
FILTER_VCF     = config.get("ScriptOptions", "filter_VCF")
NUMBER_OF_SNPS = config.get("ScriptOptions", "number_of_SNPs")
FAST_FREEBAYES = config.get("ScriptOptions", "fast_freebayes")
VCF_FILE       = config.get("ScriptOptions", "VCF_file")

REFERENCE = config.get("GenomeReference", "REFERENCE")
REF_noChr = config.get("GenomeReference", "REF_noCHR")
REF_wChr  = config.get("GenomeReference", "REF_wCHR")

GATK_MEM    = config.get("VariantCallerParameters", "GATK_MEM")
GATK_NT     = config.get("VariantCallerParameters", "GATK_nt")
VARSCAN_MEM = config.get("VariantCallerParameters", "VARSCAN_MEM")
CACHE_DIR   = config.get("BatchOperations", "CACHE_DIR")

BATCH_RECALCULATE = False
BATCH_USE_CACHED  = True
BATCH_WRITE_CACHE = True

JUDGE_THRESHOLD = 0.95
RNA_THRESHOLD = 0.9

STDERR_ = open("/dev/null", "w")
if args.verbose:
    STDERR_ = None

#-------------------------------------------------------------------------------
# Input and output files

if args.verbose:
    print """
=========================
Checking input and output
"""

# check input bams
for bam_file in [args.bam1, args.bam2]:
    if os.access(bam_file, os.R_OK) == False:
        print "Cannot access bam file '%s'. Either it doesn't exist or it's not \
readable." % bam_file
        sys.exit(1)
    # check bam files are index:
    bam_index1 = bam_file.rstrip(".bam") + ".bai"
    bam_index2 = bam_file + ".bai"
    if (os.access(bam_index1, os.R_OK) == False and
        os.access(bam_index2, os.R_OK) == False):
        print "The input bam files need to be indexed."
        sys.exit(1)

# resolving output report path
REPORT_PATH = ""
current_dir = os.path.abspath("./")
if args.output == None:
    REPORT_PATH = os.path.join(current_dir, "bam_matcher_report")
    bam1_name = os.path.basename(args.bam1).rstrip(".bam").replace(".", "_")
    bam2_name = os.path.basename(args.bam2).rstrip(".bam").replace(".", "_")
    REPORT_PATH = REPORT_PATH + ".%s_%s_%s" % (bam1_name, bam2_name, random_str)
else:
    REPORT_PATH = os.path.abspath(args.output)
if args.no_report:
    REPORT_PATH = "/dev/null"
    args.html = False

# check output directory is writable
REPORT_DIR = os.path.dirname(REPORT_PATH)
if REPORT_PATH != "/dev/null" and os.access(  REPORT_DIR, os.W_OK   ) == False:
    print """%s
Specified output directory is not writable: %s
""" % (FILE_ERROR, REPORT_DIR)
    sys.exit(1)

# HTML path
html_path = ""
if args.html:
    html_path = REPORT_PATH + ".html"

# scratch dir
SCRATCH_DIR = "/tmp/%s" % random_str
if args.scratch_dir == None:
    if args.verbose:
        print "Scratch directory: %s" % SCRATCH_DIR
    os.mkdir(SCRATCH_DIR)
else:
    SCRATCH_DIR = args.scratch_dir
    # if specified scratch exists
    if os.path.isdir(SCRATCH_DIR) and os.access(SCRATCH_DIR, os.W_OK) == False:
        print """%s
Specified scratchd directory is not writable: %s
""" % (FILE_ERROR, SCRATCH_DIR)
        sys.exit(1)
    # if not, make it
    else:
        if args.verbose:
            print "Creating scratch directory: %s" % SCRATCH_DIR
        try:
            os.mkdir(SCRATCH_DIR)
        except OSError, e:
            print """%s
Unable to create specified scratch directory: %s
Python error: %s
""" % (FILE_ERROR, SCRATCH_DIR, e)

if args.verbose:
    print """
Input and output seem okay
--------------------------
"""

#-------------------------------------------------------------------------------
# Test caller binaries

if args.verbose:
    print """
=========================
Checking caller

Caller to use: %s
---------------------------------------------
This is not really catching errors correctly.
Needs to catch errors and exit if fails
---------------------------------------------
""" % args.caller


if JAVA == "":
    print "Java command was not specified. \
Do this in the configuration file"
    sys.exit(1)

#-------------------------------------------
# Testing GATK
if args.caller == "gatk":
    if GATK == "":
        print "GATK path was not specified. \
    Do this in the configuration file"
        sys.exit(1)

    if os.access(GATK, os.R_OK) == False:
        print "Cannot access GATK jar file (%s)" % GATK
        sys.exit(1)
    gatk_cmd = [JAVA, "-jar", GATK, "-version"]
    try:
        gatk_proc = subprocess.Popen(gatk_cmd, stdout=subprocess.PIPE,
                                     stderr=STDERR_)
        if args.verbose:
            print "Testing GATK, running:\n   '%s'" % " ".join(gatk_cmd)
            print "This should get the GATK version number:"
            for line in gatk_proc.stdout:
                print line
        else:
            gatk_proc.communicate()
    except subprocess.CalledProcessError, e:
        print "%s\nSomething wrong with GATK settings" % CONFIG_ERROR
        print "Python error msg: ", e
        sys.exit(1)

#-------------------------------------------
# Testing Freebayes
elif args.caller == "freebayes":
    if FREEBAYES == "":
        print "Path to Freebayes binary was not specified. \
Do this in the configuration file"
        sys.exit(1)

    # test freebayes
    if os.path.isfile(FREEBAYES) == False:
        print "Cannot find Freebayes (%s)" % FREEBAYES
        sys.exit(1)
    if os.access(FREEBAYES, os.R_OK) == False:
        print "Freebayes binary file (%s) is not readable" % FREEBAYES
        sys.exit(1)

    # free_cmd = "%s --version" % FREEBAYES
    free_cmd = [FREEBAYES, "--version"]
    try:
        free_proc = subprocess.Popen(free_cmd, stdout=subprocess.PIPE,
                                     stderr=STDERR_)
        if args.verbose:
            print "Testing Freebayes, running:\n   '%s'" % " ".join(free_cmd)
            print "This should get the Freebayes version number:"
            for line in free_proc.stdout:
                print line
        else:
            free_proc.communicate()
    except subprocess.CalledProcessError, e:
        print "%s\nSomething wrong with Freebayes" % CONFIG_ERROR
        print "Python error msg: ", e
        sys.exit(1)

#-------------------------------------------
# Testing Varscan
elif args.caller == "varscan":
    if VARSCAN == "":
        print "Path to VarScan2 jar file was not specified. \
Do this in the configuration file."
        sys.exit(1)

    if os.access(VARSCAN, os.R_OK) == False:
        print "Cannot access VarScan2 jar file (%s)" % VARSCAN
        sys.exit(1)

    varscan_cmd = [JAVA, "-jar", VARSCAN]
    try:
        varscan_proc = subprocess.Popen(varscan_cmd, stdout=subprocess.PIPE,
                                        stderr=STDERR_)
        if args.verbose:
            print "Testing Varscan, running:\n   '%s'" % " ".join(varscan_cmd)
            print "This should generate the version number and command menu:\n"
            for line in varscan_proc.stdout:
                print line
        else:
            varscan_proc.communicate()
    except subprocess.CalledProcessError, e:
        print "%s\nSomething went wrong with Varscan\n%s" % (CONFIG_ERROR, e)
        sys.exit(1)

if args.verbose:
    print """
Caller settings seem okay
--------------------------
"""
#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# Checking and validating settings and parameters
if args.verbose:
    print """
=========================
Checking settings and parameters"""

#-------------------------------------------
# Variants file

# is it overridden by args?
if args.vcf != None:
    VCF_FILE = args.vcf

# is it specified?
if VCF_FILE == "":
    print """%s
No variants file (VCF) has been specified.
Use --vcf/-V at command line or VCF_FILE in the configuration file.
""" % CONFIG_ERROR
    sys.exit(1)

# is it readable?
VCF_FILE = os.path.abspath(VCF_FILE)
if os.access(VCF_FILE, os.R_OK) == False:
    print """%s
Cannot find or read the variants VCF file: %s
""" % (CONFIG_ERROR, VCF_FILE)
    sys.exit(1)
if args.verbose:
    print "VCF file:      ", VCF_FILE

#-------------------------------------------
# DP threshold

# get from command line
if args.dp_threshold != None:
    DP_THRESH = args.dp_threshold
# get from config file
else:
    if DP_THRESH == "":
        print """
+---------+
| WARNING |
+---------+
DP_threshold value was not specified in the config file or arguments.
Default value (15) will be used instead.
Setting DP_threshold = 15
"""
    else:
        try:
            DP_THRESH = int(DP_THRESH)
        except ValueError, e:
            print """%s
DP_threshold value ('%s') in config file is not a valid integer.
""" % (CONFIG_ERROR, DP_THRESH)
            sys.exit(1)

if args.verbose:
    print "DP_threshold:  ", DP_THRESH

#-------------------------------------------
# filter_VCF
# only override if the flag is enabled at commandline
if args.filter_vcf:
    FILTER_VCF = args.filter_vcf
# if not specified at commandline, use config values
else:
    if FILTER_VCF == "":
        print """
+---------+
| WARNING |
+---------+
filter_VCF value was not set in the configuration file.
Default value (False) will be used instead.
Setting filter_VCF = False
"""
        FILTER_VCF = False
    else:
        if FILTER_VCF in ["False", "false", "F", "FALSE"]:
            FILTER_VCF = False
        elif FILTER_VCF in ["True", "true", "T", "TRUE"]:
            FILTER_VCF = True
        else:
            print """%s
Invalid value ('%s') was specified for filter_VCF in the configuration file.
Use 'False' or 'True'""" % (CONFIG_ERROR, FILTER_VCF)
            sys.exit(1)

if args.verbose:
    print "filter_VCF:    ", FILTER_VCF

#-------------------------------------------
# Number of SNPs
# get from command line
if args.number_of_snps != None:
    NUMBER_OF_SNPS = args.number_of_snps
# get from config file
else:
    if NUMBER_OF_SNPS == "":
        print """
+---------+
| WARNING |
+---------+
number_of_SNPs was not specified in the configuration file.
Default value (1500) will be used instead.
Setting number_of_SNPs = 1500
"""
        NUMBER_OF_SNPS = 1500
    else:
        try:
            NUMBER_OF_SNPS = int(NUMBER_OF_SNPS)
        except ValueError, e:
            print """%s
number_of_SNPs value ('%s') in config file is not a valid integer.
""" % (CONFIG_ERROR, NUMBER_OF_SNPS)
            sys.exit(1)

if NUMBER_OF_SNPS <= 200:
    print """
+---------+
| WARNING |
+---------+
Using fewer than 200 SNPs is not recommended, may not be sufficient to
correctly discriminate between samples.
"""

if args.verbose:
    print "number_of_SNPs:", NUMBER_OF_SNPS

#-------------------------------------------
# Fast Freebayes

# only check if actually using Freebayes for variant calling
if args.caller == "freebayes":
    # get from command line
    if args.fastfreebayes:
        FAST_FREEBAYES = args.fastfreebayes

    # get from config file
    else:
        if FAST_FREEBAYES == "":
            print """
+---------+
| WARNING |
+---------+
fast_freebayes was not set in the configuration file.
Default value (False) will be used instead.

Setting fast_freebayes = False
"""
            FAST_FREEBAYES = False
        else:
            if FAST_FREEBAYES in ["False", "false", "F", "FALSE"]:
                FAST_FREEBAYES = False
            elif FAST_FREEBAYES in ["True", "true", "T", "TRUE"]:
                FAST_FREEBAYES = True
            else:
                print """%s
Invalid value ('%s') was specified for fast_freebayes in the configuration file.
Use 'False' or 'True'""" % (CONFIG_ERROR, FAST_FREEBAYES)
                sys.exit(1)

    if args.verbose:
        print "fast_freebayes:", FAST_FREEBAYES

#-------------------------------------------
# GATK parameters
if args.caller == "gatk":

    # GATK_MEM
    # get from command line?
    # get from config file
    if GATK_MEM == "":
        print """
+---------+
| WARNING |
+---------+
GATK_MEM was not specified in the configuration file.
Default value (4) will be used instead.

Setting GATK_MEM = 4
"""
        GATK_MEM = 4
    else:
        try:
            GATK_MEM = int(GATK_MEM)
        except ValueError, e:
            print """%s
GATK_MEM value ('%s') in the config file is not a valid integer.
""" % (CONFIG_ERROR, GATK_MEM)
            sys.exit(1)

    # GATK_nt
    if GATK_NT == "":
        print """
+---------+
| WARNING |
+---------+
GATK_nt was not specified in the configuration file.
Default value (1) will be used instead.

Setting GATK_nt = 1
"""
        GATK_NT = 1
    else:
        try:
            GATK_NT = int(GATK_NT)
        except ValueError, e:
            print """%s
GATK_nt value ('%s') in the config file is not a valid integer.
""" % (CONFIG_ERROR, GATK_NT)
            sys.exit(1)

#-------------------------------------------
# VarScan parameters
if args.caller == "varscan":
    # get from command line?

    # get from config file
    if VARSCAN_MEM == "":
        print """
+---------+
| WARNING |
+---------+
VARSCAN_MEM was not specified in the configuration file.
Default value (4) will be used instead.

Setting VARSCAN_MEM = 4
"""
        VARSCAN_MEM = 4
    else:
        try:
            VARSCAN_MEM = int(VARSCAN_MEM)
        except ValueError, e:
            print """%s
VARSCAN_MEM value ('%s') in the config file is not a valid integer.
""" % (CONFIG_ERROR, VARSCAN_MEM)
            sys.exit(1)

#-------------------------------------------
# References
bam1_ref = ""
bam2_ref = ""
bam1_haschr = False
bam2_haschr = False
AVAILABLE_REFERENCES = []


# if any of the reference options are used, disable all config REFERENCE settings
if args.reference != None or args.ref_noChr != None or args.ref_wChr != None or args.bam1_reference != None or args.bam2_reference != None:

    REFERENCE = ""
    REF_noChr = ""
    REF_wChr  = ""
    if args.reference != None:
        REFERENCE = args.reference
    if args.ref_noChr != None:
        REF_noChr = args.ref_noChr
    if args.ref_wChr != None:
        REF_wChr = args.ref_wChr


if args.bam1_reference == None and args.bam2_reference == None:
    if REFERENCE == "" and REF_noChr == "" and REF_wChr == "": # and bam1_ref == "" and bam2_ref == "":
        print "No genome reference files were specified!"
        sys.exit(1)
    for ref in [REFERENCE, REF_noChr, REF_wChr]:
        if ref == "":
            continue
        else:
            if os.access(ref, os.R_OK) == False:
                print "Specified reference fasta file ('%s') is either not present or readable" % ref
                sys.exit(1)
        # check that the references are indexed
        ref_idx = ref + ".fai"
        if os.access(ref_idx, os.R_OK) == False:
            print "Make sure that the reference file ('%s') has been indexed by samtools." % ref
            sys.exit(1)
        # if all checks pass, add to list
        AVAILABLE_REFERENCES.append(ref)
    if args.verbose:
        print "Available references:\n%s" % "\n  ".join(AVAILABLE_REFERENCES)

    # ----------
    # compare chromosomes between reference and bam files
    bam1_chrlist = []
    bam2_chrlist = []

    # get bam1 chromosomes
    bam_in = HTSeq.BAM_Reader(args.bam1)
    bam_header = bam_in.get_header_dict()["SQ"]
    for headline in bam_header:
        chr_name = headline["SN"]
        bam1_chrlist.append(chr_name)
        if "chr" in chr_name:
            bam1_haschr = True

    # get bam2 chromosomes
    bam_in = HTSeq.BAM_Reader(args.bam2)
    bam_header = bam_in.get_header_dict()["SQ"]
    for headline in bam_header:
        chr_name = headline["SN"]
        bam2_chrlist.append(headline["SN"])
        if "chr" in chr_name:
            bam2_haschr = True

    # ---------
    # generate chr list for references
    ref_chrlist = {}
    bam1_ref_chr_ct_max = -1
    bam2_ref_chr_ct_max = -1

    for ref in AVAILABLE_REFERENCES:
        ref_chrlist[ref] = []
        ref_idx = ref + ".fai"
        for line in open(ref_idx, "r"):
            ref_chrlist[ref].append(line.strip("\n").split("\t")[0])

        bam1_ref_chr_ct = len(set(bam1_chrlist).intersection(set(ref_chrlist[ref])))
        if bam1_ref_chr_ct > bam1_ref_chr_ct_max:
            bam1_ref = ref
            bam1_ref_chr_ct_max = bam1_ref_chr_ct

        bam2_ref_chr_ct = len(set(bam2_chrlist).intersection(set(ref_chrlist[ref])))
        if bam2_ref_chr_ct > bam2_ref_chr_ct_max:
            bam2_ref = ref
            bam2_ref_chr_ct_max = bam2_ref_chr_ct

    if args.verbose:
        print "BAM1 ('%s') is matched to reference ('%s')" % (args.bam1, bam1_ref)
        print "BAM2 ('%s') is matched to reference ('%s')" % (args.bam2, bam2_ref)
    # ------------
    # check whether the matching was appropriate
    if bam1_ref_chr_ct_max < 5:
        print "Fewer than 5 chromosome names in BAM1 ('%s') are found in the \
matched reference ('%s')" % (args.bam1, bam1_ref)
        print "Please check that correct reference files are supplied"
        sys.exit(1)
    if bam2_ref_chr_ct_max < 5:
        print "Fewer than 5 chromosome names in BAM2 ('%s') are found in the \
matched reference ('%s')" % (args.bam2, bam2_ref)
        print "Please check that correct reference files are supplied"
        sys.exit(1)
elif args.bam1_reference == None or args.bam2_reference == None:
    # if only one of these are supplied
    print "When overriding REFERENCE settings using \
--bam1-reference/--bam2-reference, \
BOTH need to be specified separately."
    sys.exit(1)

elif args.bam1_reference == None or args.bam2_reference == None:
    print "both --bam1-reference and --bam2-reference must be specified when using these options"
    sys.exit(1)

else:
    # both have been supplied
    bam1_ref = args.bam1_reference
    bam2_ref = args.bam2_reference




    # check reference file and index
    for ref in [bam1_ref, bam2_ref]:
        if os.access(ref, os.R_OK) == False:
            print "Specified reference file ('%s') is either not present or \
readable" % ref
            sys.exit(1)
        ref_idx = ref + ".fai"
        if os.access(ref_idx, os.R_OK) == False:
            print "Reference fasta file ('%s') needs to be indexed by \
samtools" % ref
            sys.exit(1)


    # need to know whether these references have chr or not




    # ----------
    # compare chromosomes between reference and bam files

    bam1_haschr = False
    # get bam1 chromosomes
    bam_in = HTSeq.BAM_Reader(args.bam1)
    bam_header = bam_in.get_header_dict()["SQ"]
    for headline in bam_header:
        chr_name = headline["SN"]
        if "chr" in chr_name:
            bam1_haschr = True

    bam2_haschr = False
    # get bam2 chromosomes
    bam_in = HTSeq.BAM_Reader(args.bam2)
    bam_header = bam_in.get_header_dict()["SQ"]
    for headline in bam_header:
        chr_name = headline["SN"]
        if "chr" in chr_name:
            bam2_haschr = True


#-------------------------------------------
# Batch operations
# args.do_not_cache, args.recalculate

BATCH_WRITE_CACHE = True
BATCH_USE_CACHED  = True

if args.recalculate:
    BATCH_USE_CACHED = False
if args.do_not_cache:
    BATCH_WRITE_CACHE = False

if args.verbose:
    print "Using cached:  ", BATCH_USE_CACHED
    print "Writing cache: ", BATCH_WRITE_CACHE

if args.cache_dir != None:
    CACHE_DIR = args.cache_dir

if BATCH_WRITE_CACHE or BATCH_USE_CACHED:
    if CACHE_DIR == "":
        print """%s
No CACHE_DIR specified in configuration or at command line.
Cached operations requires this to work.
""" % CONFIG_ERROR
        sys.exit(1)

    # test if CACHE_DIR is writable
    try:
        test_cache_file = os.path.join(CACHE_DIR, "cache_test")
        test_cache = open( test_cache_file, "w"   )
        test_cache.write("")
        test_cache.close()
        os.remove(test_cache_file)
    except IOError, e:
        print """%s
Unable to write to specified cache directory ('%s')
Python error msg: %s
""" % (CONFIG_ERROR, CACHE_DIR, e)
        sys.exit(1)
CACHE_DIR = os.path.abspath(CACHE_DIR)


#===============================================================================
# Finished configuration and arguments checking and loading
#===============================================================================



#===============================================================================
# Variant calling
#===============================================================================
if args.verbose:
    print "\n=======================================\nCalling variants"

#-------------------------------------------------------------------------------
# first look for cached data if using BATCH_USE_CACHED is True
# cached file is named as the md5sum of:
# - the BAM path
# - number of variants compared
# - depth
# - bam file timestamp
# - bam header

m1 = md5()
bam1_path = os.path.abspath(args.bam1)
bam1_mtime = str(os.path.getmtime(bam1_path))
m1.update(bam1_path)
m1.update(str(NUMBER_OF_SNPS))
m1.update(str(DP_THRESH))
m1.update(bam1_mtime)
sam_cmd = [SAMTOOLS, "view", "-H", bam1_path]
sam_proc = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
for line in sam_proc.stdout:
    m1.update(line.strip("\n"))

m2 = md5()
bam2_path = os.path.abspath(args.bam2)
bam2_mtime = str(os.path.getmtime(bam2_path))
m2.update(bam2_path)
m2.update(str(NUMBER_OF_SNPS))
m2.update(str(DP_THRESH))
m2.update(bam2_mtime)
sam_cmd = [SAMTOOLS, "view", "-H", bam2_path]
sam_proc = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
for line in sam_proc.stdout:
    m2.update(line.strip("\n"))

bam1_cache_path = os.path.join(CACHE_DIR, m1.hexdigest())
bam2_cache_path = os.path.join(CACHE_DIR, m2.hexdigest())
bam1_is_cached = os.access(bam1_cache_path, os.R_OK)
bam2_is_cached = os.access(bam2_cache_path, os.R_OK)
if BATCH_USE_CACHED == False:
    bam1_is_cached = False
    bam2_is_cached = False

#-------------------------------------------
# SNPs file
# VCF_FILE, FILTER_VCF, NUMBER_OF_SNPS

# default 1KG VCF file doesn't use chr
# create intervals file
if args.verbose:
    print "Creating intervals file"

# f_itv = os.path.join(SCRATCH_DIR, "%s.target.intervals" % vcf_random_str)
f_itv = os.path.join(SCRATCH_DIR, "target.intervals")
temp_files.append(f_itv)
# convert variant VCF file to GATK intervals
convert_vcf_to_intervals(VCF_FILE, f_itv, 0, NUMBER_OF_SNPS, args.caller,
                         FILTER_VCF)

# Creating a version with 'chr'
f_itv_chr = f_itv + "_chr"
fin = open(f_itv, "r")
fout = open(f_itv_chr, "w")
for line in fin:
    fout.write("chr" + line)
fout.close()
temp_files.append(f_itv_chr)

# determining which one to use for each bam file
if bam1_haschr:
    bam1_itv = f_itv_chr
else:
    bam1_itv = f_itv

if bam2_haschr:
    bam2_itv = f_itv_chr
else:
    bam2_itv = f_itv

#-------------------------------------------------------------------------------
vcf1 = os.path.join(SCRATCH_DIR, "bam1.vcf")
vcf2 = os.path.join(SCRATCH_DIR, "bam2.vcf")

# pileup files, for VarScan
pup1 = os.path.join(SCRATCH_DIR, "bam1.pileup")
pup2 = os.path.join(SCRATCH_DIR, "bam2.pileup")

# lists
bam_list           = [args.bam1, args.bam2]
vcf_list           = [vcf1, vcf2]
interval_files_list = [bam1_itv, bam2_itv]
pup_list           = [pup1, pup2]
ref_list           = [bam1_ref, bam2_ref]
haschr_list        = [bam1_haschr, bam2_haschr]
cached_list        = [bam1_is_cached, bam2_is_cached]
cache_path_list    = [bam1_cache_path, bam2_cache_path]

# add temp files
temp_files += [vcf1, vcf2]
temp_files.append(vcf1+".idx")
temp_files.append(vcf2+".idx")
temp_files += pup_list

# Variant calling is not done in a single step,
# even though this is possible for some callers, because:
# 1. If the sample names are the same, this causes problems for GATK
# 2. They may have been mapped to different reference files
for i in [0,1]:
    if BATCH_USE_CACHED:
        if cached_list[i]:
            continue
    in_bam = bam_list[i]
    out_vcf = vcf_list[i]
    ref = ref_list[i]
    has_chr = haschr_list[i]
    interval_file = interval_files_list[i]

    if args.verbose:
        print "input bam: \t%s" % in_bam
        print "output vcf:\t%s" % out_vcf

    if args.caller == "gatk":        # GATK calling
        varcall_cmd = [JAVA, "-jar", "-Xmx%dg" % GATK_MEM,
                       "-XX:ParallelGCThreads=1", GATK, "-T",
                       "UnifiedGenotyper", "-R", ref, "-I", in_bam,
                       "-o", out_vcf]
        varcall_cmd += ["--output_mode", "EMIT_ALL_SITES", "-nt", str(GATK_NT),
                        "-L", interval_file]
        if args.verbose:
            print "GATK variant-calling command:\n(space in path not escaped \
here, but should be fine in actual call command)"
            print " ".join(varcall_cmd)
        varcall_proc = subprocess.Popen(varcall_cmd, stdout=subprocess.PIPE,
                                        stderr=STDERR_)

        varcall_proc.communicate()
    elif args.caller == "freebayes":    # Freebayes calling
        fout = open(out_vcf, "w")
        if FAST_FREEBAYES:
            varcall_cmd = [FREEBAYES, "--fasta-reference", ref, "--targets",
                           interval_file, "--no-indels", "--min-coverage",
                           str(DP_THRESH)]
            varcall_cmd += ["--report-all-haplotype-alleles",
                            "--report-monomorphic", in_bam]
            if args.verbose:
                print "Freebayes variant-calling command:"
                print " ".join(varcall_cmd)
            varcall_proc = subprocess.Popen(varcall_cmd, stdout=subprocess.PIPE,
                                            stderr=STDERR_)
            for line in varcall_proc.stdout:
                fout.write(line)
            fout.close()
        else:
            write_header = True
            itv_list = []
            fin = open(interval_file, "r")
            for line in fin:
                bits = line.strip("\n").split("\t")
                bits[1] = int(bits[1])
                bits[2] = int(bits[2]) + 1
                region_str = "%s:%d-%d" % (bits[0], bits[1], bits[2])

                varcall_cmd = [FREEBAYES, "--fasta-reference", ref, "--region",
                               region_str, "--no-indels", "--min-coverage",
                               str(DP_THRESH)]
                varcall_cmd += ["--report-all-haplotype-alleles",
                                "--report-monomorphic", in_bam]
                if args.verbose:
                    print "Freebayes variant-calling command:"
                    print " ".join(varcall_cmd)
                varcall_proc = subprocess.Popen(varcall_cmd,
                                                stdout=subprocess.PIPE,
                                                stderr=STDERR_)
                for line in varcall_proc.stdout:
                    if line.startswith("#"):
                        if write_header:
                            fout.write(line)
                    else:
                        if args.verbose:
                            print line.strip("\n")
                        fout.write(line)
                write_header = False
            fout.close()
    elif args.caller == "varscan":   # Varscan
        pup_file = pup_list[i]
        fout = open(pup_file, "w")
        # First need to do a pileup
        # not using the -l option in samtools mpileup because it doesn't seem to
        #   use indexed random BAM access, so ends up being much slower
        itv_list = []
        fin = open(interval_file, "r")
        for line in fin:
            region_str = line.strip("\n")
            sam_cmd = [SAMTOOLS, "mpileup", "-r", region_str, "-B", "-f",
                       ref, in_bam]
            if args.verbose:
                print "pileup command:\n%s" % " ".join(sam_cmd)
            sam_proc = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE,
                                        stderr=STDERR_)
            for line in sam_proc.stdout:
                bits = line.strip("\n").split("\t")
                if int(bits[3]) < args.dp_threshold:
                    continue
                fout.write(line)

        # calling with varscan
        varscan_cmd = [JAVA, "-XX:ParallelGCThreads=1",
                       "-Xmx%dg" % VARSCAN_MEM, "-jar", VARSCAN,
                       "mpileup2snp", pup_file,
                       "--output-vcf", "--min-coverage", str(DP_THRESH) ]
        if args.verbose:
            print "Varscan command:\n%s" % " ".join(varscan_cmd)
        fout = open(out_vcf, "w")
        varscan_proc = subprocess.Popen(varscan_cmd, stdout=subprocess.PIPE,
                                        stderr=STDERR_)
        for line in varscan_proc.stdout:
            fout.write(line)
        fout.close()

    # need to re-sort VCF file if genome has 'chr'
    # as the bam file chr order is not the same as the reference
    if has_chr:
        old_vcf = out_vcf + ".bak"
        new_vcf = out_vcf
        os.rename(out_vcf, old_vcf)
        ref_idx = ref + ".fai"
        sort_vcf_by_chrom_order(old_vcf, new_vcf, ref_idx)

if args.verbose:
    print """

Variant-calling finished
=============================
"""


#===============================================================================
# Finished variant calling
#===============================================================================

#-------------------------------------------------------------------------------







#===============================================================================
# Comparing variant data
#===============================================================================

if args.verbose:
    print """
=============================
Comparing VCF output

"""

#-------------------------------------------------------------------------------
# extract relevant VCF data into TSV files
tsv1 = os.path.join(SCRATCH_DIR, "bam1.tsv")
tsv2 = os.path.join(SCRATCH_DIR, "bam2.tsv")
tsv_list = [tsv1, tsv2]
temp_files += tsv_list

if args.verbose:
    print "Converting VCF to table"



for i in [0,1]:
    if BATCH_USE_CACHED:
        if cached_list[i]:
            continue
    in_vcf  = vcf_list[i]
    out_tsv = tsv_list[i]
    in_bam  = bam_list[i]
    has_chr = haschr_list[i]
    ref     = ref_list[i]

    if args.verbose:
        print """
input bam:  %s
reference:  %s
""" % (in_bam, ref)

# Using GATK to convert VCF to TSV
#     if args.caller == "gatk" or args.caller == "varscan":
#         gatk_cmd = [JAVA, "-jar", "-Xmx2g", "-XX:ParallelGCThreads=1",
#                     GATK, "-T", "VariantsToTable", "-R", ref, "-V", in_vcf,
#                     "-F", "CHROM", "-F", "POS"]
#         gatk_cmd += ["-F", "REF", "-F", "ALT", "-F", "QUAL", "-GF", "DP",
#                      "-GF", "AD", "-GF", "GT", "-o", out_tsv]
#     elif args.caller == "freebayes":
#         gatk_cmd = [JAVA, "-jar", "-Xmx2g", "-XX:ParallelGCThreads=1",
#                     GATK, "-T", "VariantsToTable", "-R", ref, "-V", in_vcf,
#                     "-F", "CHROM", "-F", "POS"]
#         gatk_cmd += ["-F", "REF", "-F", "ALT", "-F", "QUAL", "-GF", "DP",
#                      "-GF", "AO", "-GF", "GT", "-o", out_tsv]
#     if args.verbose:
#         print "GATK VariantsToTable command:\n%s" % " ".join(gatk_cmd)
#     gatk_proc = subprocess.Popen(gatk_cmd, stdout=subprocess.PIPE,
#                                  stderr=STDERR_)
#     gatk_proc.communicate()

    VCFtoTSV(in_vcf, out_tsv, args.caller)


if args.verbose:
    print "removing chr from tsv files"

for tsv in tsv_list:
    if BATCH_USE_CACHED:
        if cached_list[i]:
            continue
    f_temp = os.path.join(SCRATCH_DIR, "temp_file")
    fin = open(tsv, "r")
    fout = open(f_temp, "w")
    for line in fin:
        fout.write(line.lstrip("chr"))
    fout.close()
    os.rename(f_temp, tsv)
if args.verbose:
    print "Finished VCF conversion"

#-------------------------------------------------------------------------------
# Cache the tsv files

if BATCH_USE_CACHED:
    if bam1_is_cached:
        tsv1 = bam1_cache_path
    if bam2_is_cached:
        tsv2 = bam2_cache_path

if BATCH_WRITE_CACHE:
    for i in [0,1]:
        if cached_list[i] == False:
            try:
                if args.verbose:
                    print "copying BAM%d variant calling data to %s" % ( (i+1), cache_path_list[i])
                shutil.copyfile(tsv_list[i], cache_path_list[i])
            except IOError, e:
                print "Unable to write cache results for BAM%d" % (i+1)
                print "Python error msg:", e
                sys.exit(1)

#-------------------------------------------------------------------------------
# apply depth filter
if args.verbose:
    print "Variant comparison"
bam1_var = os.path.join(SCRATCH_DIR, "bam1.variants")
bam2_var = os.path.join(SCRATCH_DIR, "bam2.variants")
var_list = {}
common_vars = {}
temp_files += [bam1_var, bam2_var]

#-------------------------------------------------------------------------------
# first get list of passed variants in bam1
fin = open(tsv1, "r")
for line in fin:
    if line.startswith("CHROM\t"):
        continue
    bits = line.strip("\n").split("\t")

    if bits[5] == "NA":
        continue
    elif int(bits[5]) < DP_THRESH:
        continue
    else:
        # add variants to list
        var_list["\t".join(bits[:4])] = 1

#-------------------------------------------------------------------------------
# then parse second tsv file to get list of variants that passed in both bams
fin = open(tsv2, "r")
for line in fin:
    if line.startswith("CHROM\t"):
        continue
    bits = line.strip("\n").split("\t")

    var_ = "\t".join(bits[:4])
    if var_ in var_list:
        var_list[var_] = 2

#-------------------------------------------------------------------------------
# write out bam1 variants
fout = open(bam1_var, "w")
fin = open(tsv1, "r")
for line in fin:
    if line.startswith("CHROM"):
        continue
    bits = line.strip("\n").split("\t")
    out_line = "%s\t%s\t%s\t%s\t%s\n" % (bits[0], bits[1], bits[2],
                                         bits[3], bits[7])
    var_ = "\t".join(bits[:4])
    if var_ in var_list:
        if var_list[var_] == 2:
            fout.write(out_line)
fout.close()

#-------------------------------------------------------------------------------
# write out bam2 variants
fout = open(bam2_var, "w")
fin = open(tsv2, "r")
for line in fin:
    if line.startswith("CHROM"):
        continue
    bits = line.strip("\n").split("\t")
    out_line = "%s\t%s\t%s\t%s\t%s\n" % (bits[0], bits[1], bits[2], bits[3],
                                         bits[7])
    var_ = "\t".join(bits[:4])
    if var_ in var_list:
        if var_list[var_] == 2:
            fout.write(out_line)
fout.close()

#-------------------------------------------------------------------------------
# However, bam1_var and bam2_var may still be sorted differently
# so sort to the same way
for fvar in [bam1_var, bam2_var]:
#    f_sorted = os.path.join(SCRATCH_DIR,
#                            "%s.sorted_variants_file" % vcf_random_str)
    f_sorted = os.path.join(SCRATCH_DIR, "sorted_variants_file")
    sort_cmd = "sort -k1n -k2n '%s' > '%s' " % (fvar, f_sorted)
    if args.verbose:
        print sort_cmd
    sort_proc = subprocess.Popen([sort_cmd], shell=True, stdout=subprocess.PIPE,
                                 stderr=STDERR_)
    sort_proc.communicate()
    shutil.copy(f_sorted, fvar)

#-------------------------------------------------------------------------------
# Results variables

comm_het_ct = 0
comm_hom_ct = 0
diff_hom_ct = 0
diff_het_ct = 0
diff_het_hom_ct = 0   # het in 1, hom in 2
diff_hom_het_ct = 0   # hom in 1, het in 2
diff_1sub2_ct = 0   # bam1 GT is subset of bam2 GT  (bam1: hom, bam2: het)
diff_2sub1_ct = 0   # bam2 GT is subset of bam1 GT
ct_common = 0
ct_bam1 = 0
ct_bam2 = 0
ct_diff = 0
bam1_gt = {}
bam2_gt = {}
pos_list = []

#-------------------------------------------------------------------------------
# This is the actual comparison
fin = open(bam1_var, "r")
for line in fin:
    bits = line.strip("\n").split("\t")
    pos_ = "_".join(bits[:2])
    geno = bits[4]
    pos_list.append(pos_)
    bam1_gt[pos_] = geno

fin = open(bam2_var, "r")
for line in fin:
#    print line
    bits = line.strip("\n").split("\t")
    pos_ = "_".join(bits[:2])
    geno = bits[4]
    bam2_gt[pos_] = geno

for pos_ in pos_list:
    gt1 = bam1_gt[pos_]
    gt2 = bam2_gt[pos_]

    # if genotypes are the same
    if same_gt(gt1, gt2):
        ct_common += 1
        if is_hom(gt1):
            comm_hom_ct += 1
        else:
            comm_het_ct += 1
    else:
        ct_diff += 1

        # both are hom and different
        if is_hom(gt1) and is_hom(gt2):
            diff_hom_ct += 1
        # both are het and different
        elif is_hom(gt1) == False and is_hom(gt2) == False:
            diff_het_ct += 1
        # one is hom, one is het, test for subset
        elif is_hom(gt1):
            if is_subset(gt1, gt2):
#                print "%s is a subset of %s" % (gt1, gt2)
                diff_1sub2_ct += 1
            else:
                diff_hom_het_ct += 1
        elif is_hom(gt2):
            if is_subset(gt2, gt1):
#                print "%s is a subset of %s" % (gt2, gt1)
                diff_2sub1_ct += 1
            else:
                diff_het_hom_ct += 1
#-------------------------------------------------------------------------------







#===============================================================================
# write reports
#===============================================================================

if args.verbose:
    print "Writing output report"
total_compared = ct_common + ct_diff
frac_common = float(ct_common)/total_compared





# SAME PATIENT
# JUDGE_THRESHOLD = 0.95
# frac_common >= threshold & total_compared >= 50

# POSSIBLE RNA?
# check ratio of 1het-2sub and 1sub-2het

# frac_common <= 0.70 & total_compared >= 50





# determine whether they are from the same patient
judgement = "DIFFERENT SOURCES"
short_judgement = "Diff"
if frac_common >= JUDGE_THRESHOLD:
    judgement = "SAME SOURCE"
    short_judgement = "Same"

# but RNA-seq data...
# 1sub2 or 2sub should account for most of the differences in ct_diff

# p_1sub2 = pvalue(total_compared, diff_1sub2_ct, total_compared, ct_diff)
# p_2sub1 = pvalue(total_compared, diff_2sub1_ct, total_compared, ct_diff)
else:
    if float(diff_1sub2_ct)/ct_diff >= RNA_THRESHOLD:
        judgement = "Probably same sample. BAM1 is likely RNA-seq data of the same sample"
        short_judgement = "Same (BAM1=RNA)"
    if float(diff_2sub1_ct)/ct_diff >= RNA_THRESHOLD:
        judgement = "Probably same sample. BAM2 is likely RNA-seq data of the same sample"
        short_judgement = "Same (BAM2=RNA)"







# STANDARD FORMAT

# assume that the most space required is 4 digits,
# so pad numeric string to 6 spaces
diff_hom = ("%d" % diff_hom_ct).rjust(5)
diff_het = ("%d" % diff_het_ct).rjust(5)
diff_het_hom = ("%d" % diff_het_hom_ct).rjust(5)
diff_hom_het = ("%d" % diff_hom_het_ct).rjust(5)
diff_1sub2 = ("%d" % diff_1sub2_ct).rjust(5)
diff_2sub1 = ("%d" % diff_2sub1_ct).rjust(5)
std_report_str = """bam1:\t%s
bam2:\t%s
DP threshold: %d
________________________________________

Positions with same genotype:     %d
     breakdown:    hom: %d
                   het: %d
________________________________________

Positions with diff genotype:     %d
     breakdown:
                       BAM 1
               | het  | hom  | subset
        -------+------+------+-------
         het   |%s |%s |%s |
        -------+------+------+-------
BAM 2    hom   |%s |%s |   -  |
        -------+------+------+-------
         subset|%s |   -  |   -  |
________________________________________

Total sites compared: %d
Fraction of common: %f (%d/%d)
CONCLUSION: %s
"""  % (bam1_path, bam2_path, DP_THRESH,
        ct_common, comm_hom_ct, comm_het_ct,
        ct_diff, diff_het, diff_hom_het, diff_1sub2, diff_het_hom, diff_hom, diff_2sub1,
        total_compared, frac_common, ct_common, total_compared, judgement)

# SHORT FORMAT
short_report_str = """# BAM1\t BAM2\t DP_thresh\t FracCommon\t Same\t Same_hom\t Same_het\t Different\t 1het-2het\t 1het-2hom\t 1het-2sub\t 1hom-2het\t 1hom-2hom\t 1sub-2het\t Conclusion
%s\t%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s""" % (bam1_path,
       bam2_path, DP_THRESH, frac_common, ct_common, comm_hom_ct, comm_het_ct,
       ct_diff, diff_het_ct, diff_het_hom_ct, diff_2sub1_ct, diff_hom_het_ct,
       diff_hom_ct, diff_1sub2_ct, short_judgement)

# HTML FORMAT
if args.html:
    html_bam1_cached = "recalculated"
    if BATCH_USE_CACHED and bam1_is_cached:
        html_bam1_cached = "cached"
    html_bam2_cached = "recalculated"
    if BATCH_USE_CACHED and bam2_is_cached:
        html_bam2_cached = "cached"
    html_namespace = {"BAM1"        : bam1_path,
                      "BAM2"        : bam2_path,
                      "DP_THRESH"   : DP_THRESH,
                      "caller"      : args.caller,
                      "variants"    : VCF_FILE,
                      "bam1_cached" : html_bam1_cached,
                      "bam2_cached" : html_bam2_cached,

                      "total"          : total_compared,
                      "frac_common"    : frac_common*100,   # report % here
                      "ct_common"      : ct_common,
                      "comm_hom_ct"    : comm_hom_ct,
                      "comm_het_ct"    : comm_het_ct,

                      "ct_diff"        : ct_diff,
                      "diff_het_ct"    : diff_het_ct,
                      "diff_het_hom_ct": diff_het_hom_ct,
                      "diff_2sub1_ct"  : diff_2sub1_ct,
                      "diff_hom_het_ct": diff_hom_het_ct,
                      "diff_hom_ct"    : diff_hom_ct,
                      "diff_1sub2_ct"  : diff_1sub2_ct,
                      "judgement"      : judgement
                     }
    template_file = os.path.join(SCRIPT_DIR, "bam_matcher_html_template")
    template = Template(file=template_file, searchList=html_namespace)
    fout = open(html_path, "w")
    fout.write(str(template))
    fout.close()

# Writing output
fout = open(REPORT_PATH, "w")
if args.short_output:
    fout.write(short_report_str)
    print short_report_str
else:
    fout.write(std_report_str)
    print std_report_str
fout.close()

#===============================================================================
# house keeping
if args.debug == False:
    for f in os.listdir(SCRATCH_DIR):
        fpath = os.path.join(SCRATCH_DIR, f)
        os.remove(fpath)
    os.rmdir(SCRATCH_DIR)
else:
    print """
##DEBUG##
Temporary files were written to: %s
""" % SCRATCH_DIR

#===============================================================================
