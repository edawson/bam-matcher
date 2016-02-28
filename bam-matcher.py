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
import os
import subprocess
import sys
import vcf
import random
import string
import ConfigParser
import shutil
from Cheetah.Template import Template
from hashlib import md5
from bammatcher_methods import *

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
                             help="Default reference fasta file. Needs to be \
                             indexed with samtools faidx. Overrides config settings.")
    parser_grp6.add_argument("--ref-alternate", "-R2",  required=False,
                             help="Alternate reference fasta file. Needs to be \
                             indexed with samtools faidx. Overrides config settings.")
    parser_grp6.add_argument("--bam1-reference", "-B1R", required=False,
                             help="Reference fasta file for BAM1. Requires \
                             --bam2-reference/-B2R, overrides all other settings \
                             (i.e. will ignore default or alternate reference if specified).")
    parser_grp6.add_argument("--bam2-reference", "-B2R", required=False,
                             help="Reference fasta file for BAM2. Requires \
                             --bam1-reference/-B1R, overrides other settings \
                             (i.e. will ignore default or alternate reference if specified).")
    parser_grp6.add_argument("--chromosome-map", "-M", required=False,
                             help="Required when using alternate reference. \
                             Run BAM-matcher with --about-alternate-ref for more details.")
    parser_grp6.add_argument("--about-alternate-ref", "-A", required=False,
                             help="Print information about using --alternate-ref and --chromosome-map")

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

#===============================================================================
# Parsing configuration file and command line options
#-------------------------------------------------------------------------------
# Before calling handle_args(), need to check if --generate-config was called
GENERATE_CONFIG_SYMS = ["--generate-config", "--generate_config", "-G"]
ALTERNATE_REF_SYMS = ["--about-alternate-ref", "--about_alternate_ref", "-A"]
for idx, arg in enumerate(sys.argv):

    # generate config template
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
            print "%s\nThe specified path ('%s') for config template exists already." % (FILE_ERROR, config_template_output)
            print "Write to another file."
            sys.exit(1)
        fout = open(config_template_output, "w")
        fout.write(CONFIG_TEMPLATE_STR)
        fout.close()
        exit()

    # about alternate genome reference
    if arg in ALTERNATE_REF_SYMS:
        print ABOUT_ALTERNATE_REF_MSG
        exit()

# okay to parse the arguments now
args = handle_args()
#-------------------------------------------------------------------------------
# Some random-related stuff, this is for temp files
random.seed()
random_str = ""
for _ in range(12):
    random_str += random.choice(string.ascii_uppercase + string.digits)
temp_files = []
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

#-------------------------------------------------------------------------------
# Read in configuration file
config_file = ""
if args.config == None:
    config_file = os.path.join(SCRIPT_DIR, "bam-matcher.conf")
else:
    config_file = args.config

# Test if the config file exists
if check_file_read(config_file, "config", CONFIG_ERROR) == False: exit(1)

# Load config and check if all necessary bits are available
config = ConfigParser.ConfigParser()
try:
    config.read(config_file)
except ConfigParser.Error as e:
    print "%s\nUnspecified configuration file error. Please check configuration file.\n\nPython error msg:\n%s" % (CONFIG_ERROR, e)
    exit(1)

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
        exit(1)
#-------------------------------------------------------------------------------
# setting variables using the config file
GATK           = fetch_config_value(config, "VariantCallers", "GATK")
FREEBAYES      = fetch_config_value(config, "VariantCallers", "freebayes")
SAMTOOLS       = fetch_config_value(config, "VariantCallers", "samtools")
VARSCAN        = fetch_config_value(config, "VariantCallers", "varscan")
JAVA           = fetch_config_value(config, "VariantCallers", "java")
DP_THRESH      = fetch_config_value(config, "ScriptOptions", "DP_threshold")
NUMBER_OF_SNPS = fetch_config_value(config, "ScriptOptions", "number_of_SNPs")
FAST_FREEBAYES = fetch_config_value(config, "ScriptOptions", "fast_freebayes")
VCF_FILE       = fetch_config_value(config, "ScriptOptions", "VCF_file")
REFERENCE      = fetch_config_value(config, "GenomeReference", "REFERENCE")
REF_ALTERNATE  = fetch_config_value(config, "GenomeReference", "REF_ALTERNATE")
CHROM_MAP      = fetch_config_value(config, "GenomeReference", "CHROM_MAP")
GATK_MEM       = fetch_config_value(config, "VariantCallerParameters", "GATK_MEM")
GATK_NT        = fetch_config_value(config, "VariantCallerParameters", "GATK_nt")
VARSCAN_MEM    = fetch_config_value(config, "VariantCallerParameters", "VARSCAN_MEM")
CACHE_DIR      = fetch_config_value(config, "BatchOperations", "CACHE_DIR")

BATCH_RECALCULATE = False
BATCH_USE_CACHED  = True
BATCH_WRITE_CACHE = True
JUDGE_THRESHOLD   = 0.95
RNA_THRESHOLD     = 0.9

# don't write anything to standard output if not verbose
STDERR_ = open("/dev/null", "w")
if args.verbose:
    STDERR_ = None

#-------------------------------------------------------------------------------
# Input and output files
if args.verbose:
    print """
================================================================================
CHECKING INPUT AND OUTPUT
================================================================================
"""

# check input bams are readable and indexed
for bam_file in [args.bam1, args.bam2]:
    if check_file_read(bam_file, "BAM", FILE_ERROR) == False: exit(1)

    # check bam files are indexed:
    bam_index2 = bam_file.rstrip(".bam") + ".bai"
    bam_index1 = bam_file + ".bai"
    if (os.access(bam_index1, os.R_OK) == False and
        os.access(bam_index2, os.R_OK) == False):
        print "%s\nInput BAM file (%s) is either missing index or the index \
file is not readable." % (FILE_ERROR, bam_file)
        exit(1)

# ----------------------------
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

# ----------------------------
# check output directory is writable
REPORT_DIR = os.path.dirname(REPORT_PATH)
if REPORT_PATH != "/dev/null":
    if not check_file_write(REPORT_DIR, "output report directory", FILE_ERROR): exit(1)

# ----------------------------
# HTML path
html_path = ""
if args.html:
    html_path = REPORT_PATH + ".html"

# ----------------------------
# scratch dir
SCRATCH_DIR = "/tmp/%s" % random_str
if args.scratch_dir == None:
    if args.verbose:
        print "Making scratch directory: %s" % SCRATCH_DIR
    os.mkdir(SCRATCH_DIR)
else:
    SCRATCH_DIR = args.scratch_dir
    # if specified scratch exists, check it's writeable
    if os.path.isdir(SCRATCH_DIR):
        if check_file_write(SCRATCH_DIR, "scratch directory", FILE_ERROR) == False:
            exit(1)
    # if not, make it
    else:
        if args.verbose:
            print "Creating scratch directory: %s" % SCRATCH_DIR
        try:
            os.mkdir(SCRATCH_DIR)
        except OSError as e:
            print """%s
Unable to create specified scratch directory: %s
Python error: %s
""" % (FILE_ERROR, SCRATCH_DIR, e)
        except Exception as e:
            print """%s
Unknown error while trying to create scratch directory (%s)
Python error message: %s
""" % (FILE_ERROR, SCRATCH_DIR, e)

if args.verbose:
    print """

Input BAM files:
BAM1:           %s
BAM2:           %s

Output and scratch:
scratch dir:    %s
Output report:  %s

Input and output seem okay

""" % (args.bam1, args.bam2, SCRATCH_DIR, REPORT_PATH)




#===============================================================================
# Checking and validating settings and parameters
if args.verbose:
    print """
================================================================================
CHECKING SETTINGS AND PARAMETERS
================================================================================
"""

#-------------------------------------------
# Variants file

# is it overridden by args?
if args.vcf: VCF_FILE = args.vcf

# is it specified?
if VCF_FILE == "":
    print """%s
No variants file (VCF) has been specified.
Use --vcf/-V at command line or VCF_FILE in the configuration file.
""" % CONFIG_ERROR
    exit(1)

# is it readable?
VCF_FILE = os.path.abspath(VCF_FILE)
if check_file_read(VCF_FILE, "variants VCF", CONFIG_ERROR) == False: exit(1)

#-------------------------------------------
# DP threshold

# get from command line
if args.dp_threshold != None:
    DP_THRESH = args.dp_threshold
# get from config file
else:
    if DP_THRESH == "":
        print """%s
DP_threshold value was not specified in the config file or arguments.
Default value (15) will be used instead.
Setting DP_threshold = 15
""" % WARNING_MSG
        DP_THRESH = 15
    else:
        try:
            DP_THRESH = int(DP_THRESH)
        except ValueError, e:
            print """%s
DP_threshold value ('%s') in config file is not a valid integer.
""" % (CONFIG_ERROR, DP_THRESH)
            exit(1)

#-------------------------------------------
# Number of SNPs
# get from command line
if args.number_of_snps != None:
    NUMBER_OF_SNPS = args.number_of_snps
# get from config file
else:
    if NUMBER_OF_SNPS == "":
        NUMBER_OF_SNPS = 0
    else:
        try:
            NUMBER_OF_SNPS = int(NUMBER_OF_SNPS)
        except ValueError, e:
            print """%s
number_of_SNPs value ('%s') in config file is not a valid integer.
""" % (CONFIG_ERROR, NUMBER_OF_SNPS)
            sys.exit(1)

if NUMBER_OF_SNPS <= 200 and NUMBER_OF_SNPS > 0:
    print """%s
Using fewer than 200 SNPs is not recommended, may not be sufficient to
correctly discriminate between samples.
""" % WARNING_MSG

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
            print """%s
fast_freebayes was not set in the configuration file.
Default value (False) will be used instead.

Setting fast_freebayes = False
""" % WARNING_MSG
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

#-------------------------------------------
# GATK parameters
if args.caller == "gatk":
    # GATK_MEM
    if GATK_MEM == "":
        print """%s
GATK_MEM was not specified in the configuration file.
Default value (4G) will be used instead.

Setting GATK_MEM = 4G
""" % WARNING_MSG
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
        print """ %s GATK_nt was not specified in the configuration file.
Default value (1) will be used instead.

Setting GATK_nt = 1
""" % WARNING_MSG
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
    # get from config file
    if VARSCAN_MEM == "":
        print """
VARSCAN_MEM was not specified in the configuration file.
Default value (4GB) will be used instead.

Setting VARSCAN_MEM = 4GB
""" % WARNING_MSG
        VARSCAN_MEM = 4
    else:
        try:
            VARSCAN_MEM = int(VARSCAN_MEM)
        except ValueError, e:
            print """%s
VARSCAN_MEM value ('%s') in the config file is not a valid integer.
""" % (CONFIG_ERROR, VARSCAN_MEM)
            exit(1)


#===============================================================================
# New reference matching

bam1_ref = ""
bam2_ref = ""

# No references are specified anywhere
if REFERENCE == "" and args.reference == None and args.bam1_reference == None and args.bam2_reference == None:
    print """%s
No genome reference has been specified anywhere.
Need to do this in either the configuration file or at run time (--reference/-R).

ALTERNATE_REF (in config) or --alternate_ref/-A should only be used if there is
already a default genome reference speficied.
""" % CONFIG_ERROR
    exit(1)

# overriding config REFERENCE
if args.reference != None:
    REFERENCE = os.path.abspath(args.reference)
    print """%s
--reference/-R argument overrides config setting (REFERENCE).
Default reference = %s
""" % (WARNING_MSG, REFERENCE)

# overriding config REF_ALTERNATE
if args.ref_alternate != None:
    REF_ALTERNATE = os.path.abspath(args.ref_alternate)
    print """%s
--ref-alternate/-R2 argument overrides config setting (REF_ALTERNATE).
Alternate reference = %s
""" % (WARNING_MSG, REF_ALTERNATE)

# overriding config CHROM_MAP
if args.chromosome_map != None:
    CHROM_MAP = args.chromosome_map
    print """%s
--chromosome-map/-M argument overrides config setting (CHROM_MAP).
Chromosome map = %s
""" % (WARNING_MSG, CHROM_MAP)

# check that if bam1_reference or bam2_reference are both supplied
if (args.bam1_reference == None and args.bam2_reference != None) or (args.bam1_reference != None and args.bam2_reference == None):
    print """%s
When using --bam1_reference(-B1R) and --bam2_reference(-B2R), both need to be specified.
""" % (ARGUMENT_ERROR)
    exit(1)

# check that if ref_alternate is used, then chromosome_map is also supplied
if REF_ALTERNATE:
    if not CHROM_MAP:
        print """%s
When using an alternate genome reference (--ref-alternate), chromosome map
(--chromosome-map/-M or CHROM_MAP in config) must also be supplied.
For more details, run BAM-matcher with --about-alternate-ref/-A.
""" % (CONFIG_ERROR)
        exit(1)

# ---------------------------------------------------------
GT_DATA_CONVERT_TO = ""

using_default_reference = True
using_chrom_map = False

# When only 1 reference is available
if args.bam1_reference == None and args.bam2_reference == None and REF_ALTERNATE == "":
    # just use REFERENCE
    bam1_ref = REFERENCE
    bam2_ref = REFERENCE

    # check REFERENCE
    if check_file_read(REFERENCE, "default reference", FILE_ERROR) == False: exit(1)
    if not check_fasta_index(REFERENCE): exit(1)

    GT_DATA_CONVERT_TO = "REF"

# When using bam1_reference and bam2_reference
elif args.bam1_reference != None and args.bam2_reference != None:
    bam1_ref = args.bam1_reference
    bam2_ref = args.bam2_reference

    # check the reference files
    if not check_file_read(bam1_ref, "BAM1 reference", FILE_ERROR): exit(1)
    if not check_file_read(bam2_ref, "BAM2 reference", FILE_ERROR): exit(1)
    if not check_fasta_index(bam1_ref): exit(1)
    if not check_fasta_index(bam2_ref): exit(1)

    GT_DATA_CONVERT_TO = "B1REF"

    # if B1R != B2R, then need a chromosome map too
    if bam1_ref != bam2_ref:
        using_chrom_map = True
        # no CHROM_MAP
        if CHROM_MAP == "":
            print """%s
When --bam1-reference and --bam2-reference are supplied but are not the same,
a chromosome map (--CHROM-MAP/-M) is required.
""" % (CONFIG_ERROR)
            exit(1)

    using_default_reference = False

# When need to match BAM to correct reference file
else:
    # get BAM chroms
    bam1_chroms = get_chrom_names_from_BAM(args.bam1)
    bam2_chroms = get_chrom_names_from_BAM(args.bam2)

    # expect to have REFERENCE and ALTERNATE_REF
    # check ref files
    if not check_file_read(REFERENCE, "default genome reference FASTA", FILE_ERROR): exit(1)
    if not check_file_read(REF_ALTERNATE, "alternate genome reference FASTA", FILE_ERROR): exit(1)
    if not check_fasta_index(REFERENCE): exit(1)
    if not check_fasta_index(REF_ALTERNATE) : exit(1)

    # get ref chroms
    REF_CHROMS = get_chrom_names_from_REF(REFERENCE)
    ALT_CHROMS = get_chrom_names_from_REF(REF_ALTERNATE)

    # check chromsome_map
    if not check_file_read(CHROM_MAP, "chromosome map", FILE_ERROR): exit(1)

    # compare chromosomes
    MAP_REF_CHROMS, MAP_ALT_CHROMS, MAP_DEF2ALT, MAP_ALT2DEF = get_chrom_data_from_map(CHROM_MAP)
    n_chroms_expected = len(MAP_REF_CHROMS)


    # expect all REF_CHROMS to be in ref_chroms
    ref_chroms_diff = set(MAP_REF_CHROMS).difference(set(REF_CHROMS))
    alt_chroms_diff = set(MAP_ALT_CHROMS).difference(set(ALT_CHROMS))

    chrom_diffs = [ref_chroms_diff, alt_chroms_diff]
    for rid, reftype in enumerate(["default", "alternate"]):
        chrdiff = chrom_diffs[rid]
        if len(chrdiff) > 0: # - n_chroms_expected < 0:
            print """%s
Number of matching chromosomes in the %s reference genome (%d) is fewer than expected from the chromosome map (%d).\n
Missing chromosome: %r\n
Check that the correct genome reference files and chromosome map are used.
""" % (CONFIG_ERROR, reftype, n_chroms_expected-len(chrdiff), n_chroms_expected, ", ".join(chrdiff))
            exit(1)

    if args.verbose:
        print "Matching BAM files and genome references"
    bam1_REF_diff = set(MAP_REF_CHROMS).difference(set(bam1_chroms))
    bam2_REF_diff = set(MAP_REF_CHROMS).difference(set(bam2_chroms))
    bam1_ALT_diff = set(MAP_ALT_CHROMS).difference(set(bam1_chroms))
    bam2_ALT_diff = set(MAP_ALT_CHROMS).difference(set(bam2_chroms))

    if len(bam1_REF_diff) == 0:
        bam1_ref = REFERENCE
        if args.verbose:
            print "BAM1 (%s) is matched to default genome reference (%s)" % (args.bam1, REFERENCE)
    elif len(bam1_ALT_diff) == 0:
        bam1_ref = REF_ALTERNATE
        if args.verbose:
            print "BAM1 (%s) is matched to alternate genome reference (%s)" % (args.bam2, REF_ALTERNATE)
    else:
        print """%s
Cannot match BAM1 to a genome reference
BAM1 (%s) is missing:
- %d chromosomes against default reference (%s). Missing:
- %d chromosomes against alternate reference (%s). Missing:
""" % (CONFIG_ERROR, args.bam1, len(bam1_REF_diff), REFERENCE, len(bam1_ALT_diff), REF_ALTERNATE)
        exit(1)

    if len(bam2_REF_diff) == 0:
        bam2_ref = REFERENCE
        if args.verbose:
            print "BAM2 (%s) is matched to default genome reference (%s)" % (args.bam1, REFERENCE)
    elif len(bam2_ALT_diff) == 0:
        bam2_ref = REF_ALTERNATE
        if args.verbose:
            print "BAM2 (%s) is matched to alternate genome reference (%s)" % (args.bam2, REF_ALTERNATE)
    else:
        print """%s
Cannot match BAM2 to a genome reference
BAM2 (%s) is missing:
- %d chromosomes against default reference (%s). Missing:
- %d chromosomes against alternate reference (%s). Missing:
""" % (CONFIG_ERROR, args.bam2, len(bam2_REF_diff), REFERENCE, len(bam2_ALT_diff), REF_ALTERNATE)
        exit(1)

    GT_DATA_CONVERT_TO = "REF"

    # use chromosome map if either BAM1 or BAM2 are not using default REFERENCE
    if bam1_ref != REFERENCE or bam2_ref != REFERENCE:
        using_chrom_map = True

#===============================================================================
# Batch operations
# args.do_not_cache, args.recalculate

BATCH_WRITE_CACHE = True
BATCH_USE_CACHED  = True

if args.recalculate:
    BATCH_USE_CACHED = False
if args.do_not_cache:
    BATCH_WRITE_CACHE = False
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


#------------------------------------------
# Print config settings

if args.verbose:
    print """\n\nCONFIG SETTINGS
VCF file:          %s
DP_threshold:      %d
number_of_SNPs:    %d (if 0, all variants in VCF file will be used)

Caller:            %s""" % (VCF_FILE, DP_THRESH, NUMBER_OF_SNPS, args.caller)

    if args.caller == "freebayes":
        print "fast_freebayes:   ", FAST_FREEBAYES

    print """
default genome reference:   %s
alternate genome reference: %s
BAM1 matched to:            %s
BAM2 matched to:            %s""" % (REFERENCE, REF_ALTERNATE, bam1_ref, bam2_ref)

    if using_chrom_map:
        print """
using chromosome map:       %r
chromosome map:             %s
""" % (using_chrom_map, CHROM_MAP)

    print """
use cached wherever possible:     %r
write cache data for new samples: %r
""" % (BATCH_USE_CACHED, BATCH_WRITE_CACHE)

#===============================================================================
# Finished configuration and arguments checking and loading
#===============================================================================





#===============================================================================
# Variant calling
#===============================================================================

if args.verbose:
    print """
================================================================================
GENOTYPE CALLING
================================================================================
"""

#-------------------------------------------------------------------------------
# first look for cached data if using BATCH_USE_CACHED is True
# cached file is named as the md5sum of:
# - the BAM path
# - number of variants compared
# - depth
# - bam file timestamp
# - bam header
# - VCF file path (list of variants to check)

m1 = md5()
bam1_path = os.path.abspath(args.bam1)
bam1_mtime = str(os.path.getmtime(bam1_path))
m1.update(bam1_path)
m1.update(str(NUMBER_OF_SNPS))
m1.update(str(DP_THRESH))
m1.update(bam1_mtime)
m1.update(VCF_FILE)
m1.update(bam1_ref)
if GT_DATA_CONVERT_TO == "REF":
    m1.update(REFERENCE)
elif GT_DATA_CONVERT_TO == "B1R":
    m1.update(bam1_ref)
for line in get_bam_header(bam1_path):
    m1.update(line)

m2 = md5()
bam2_path = os.path.abspath(args.bam2)
bam2_mtime = str(os.path.getmtime(bam2_path))
m2.update(bam2_path)
m2.update(str(NUMBER_OF_SNPS))
m2.update(str(DP_THRESH))
m2.update(bam2_mtime)
m2.update(VCF_FILE)
m2.update(bam2_ref)
if GT_DATA_CONVERT_TO == "REF":
    m2.update(REFERENCE)
elif GT_DATA_CONVERT_TO == "B1R":
    m2.update(bam1_ref)
for line in get_bam_header(bam2_path):
    m2.update(line)

bam1_cache_path = os.path.join(CACHE_DIR, m1.hexdigest())
bam2_cache_path = os.path.join(CACHE_DIR, m2.hexdigest())
bam1_is_cached = os.access(bam1_cache_path, os.R_OK)
bam2_is_cached = os.access(bam2_cache_path, os.R_OK)
if BATCH_USE_CACHED == False:
    bam1_is_cached = False
    bam2_is_cached = False
print "BAM1 is cached:", bam1_is_cached
print "BAM2 is cached:", bam2_is_cached

# Only check callers if not using cached data
# Test caller binaries
if bam1_is_cached == False or bam2_is_cached==False:
    if args.verbose:
        print "\n---------------\nChecking caller\n---------------"
    if args.caller == "freebayes" and not JAVA:
        print "%s\nJava command was not specified.\nDo this in the configuration file" % CONFIG_ERROR
        sys.exit(1)
    caller_check_log = os.path.join(SCRATCH_DIR, "caller_check.log")
    if args.caller == "gatk":
        check_caller(args.caller, GATK,      JAVA, args.verbose, logfile=caller_check_log)
    elif args.caller == "freebayes":
        check_caller(args.caller, FREEBAYES, JAVA, args.verbose, logfile=caller_check_log)
    elif args.caller == "varscan":
        check_caller(args.caller, VARSCAN,   JAVA, args.verbose, logfile=caller_check_log, SAMTL=SAMTOOLS)
    if args.verbose:
        print "Caller settings seem okay.\n"
else:
    if args.verbose:
        print """
---------------
Checking caller
---------------
Using cached data for both BAM files, so don't need to test caller.
"""
#-------------------------------------------------------------------------------
# generating intervals file for variant calling - only required if not using cached data
# VCF_FILE, NUMBER_OF_SNPS
if bam1_is_cached == False or bam2_is_cached == False:
    if args.verbose:
        print "Creating intervals file"
    bam1_itv = os.path.join(SCRATCH_DIR, "bam1.intervals")
    bam2_itv = os.path.join(SCRATCH_DIR, "bam2.intervals")
    temp_files.append(bam1_itv)
    temp_files.append(bam2_itv)

    # if using default reference (or alternate reference),
    # assume that the VCF file is same as default
    if using_default_reference:
        # if using chrom_map, assume that at least one BAM is not mapped to REFERENCE
        if using_chrom_map:
            _, _, MAP_DEF2ALT, _ = get_chrom_data_from_map(CHROM_MAP)
            if bam1_ref == REFERENCE:
                convert_vcf_to_intervals(VCF_FILE, bam1_itv, 0, NUMBER_OF_SNPS, args.caller)
            else:
                convert_vcf_to_intervals(VCF_FILE, bam1_itv, 0, NUMBER_OF_SNPS, args.caller, cmap=MAP_DEF2ALT)
            if bam2_ref == REFERENCE:
                convert_vcf_to_intervals(VCF_FILE, bam2_itv, 0, NUMBER_OF_SNPS, args.caller)
            else:
                convert_vcf_to_intervals(VCF_FILE, bam2_itv, 0, NUMBER_OF_SNPS, args.caller, cmap=MAP_DEF2ALT)
        # else, assume that both are mapped to REFERENCE
        else:
            convert_vcf_to_intervals(VCF_FILE, bam1_itv, 0, NUMBER_OF_SNPS, args.caller)
            convert_vcf_to_intervals(VCF_FILE, bam2_itv, 0, NUMBER_OF_SNPS, args.caller)
    else:
        # ===============================================
        # ===============================================
        # ===============================================
        # STILL NEED TO DO THE CASE for -B1R/-B2R
        # ===============================================
        # ===============================================
        # ===============================================
        pass
    interval_files_list = [bam1_itv, bam2_itv]

#-------------------------------------------------------------------------------
# Caller output files
if args.verbose:
    print "\n----------------\nCalling variants\n----------------"
vcf1 = os.path.join(SCRATCH_DIR, "bam1.vcf")
vcf2 = os.path.join(SCRATCH_DIR, "bam2.vcf")
# pileup files, for VarScan
pup1 = os.path.join(SCRATCH_DIR, "bam1.pileup")
pup2 = os.path.join(SCRATCH_DIR, "bam2.pileup")
# lists
bam_list           = [args.bam1, args.bam2]
vcf_list           = [vcf1, vcf2]
pup_list           = [pup1, pup2]
ref_list           = [bam1_ref, bam2_ref]
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
    if args.verbose:
        print "\nGenotype calling for BAM %d" % (i+1)
    if BATCH_USE_CACHED:
        if cached_list[i]:
            if args.verbose:
                print "BAM %d has cached genotype data" % (i+1)
            continue
    in_bam = bam_list[i]
    out_vcf = vcf_list[i]
    ref = ref_list[i]
    interval_file = interval_files_list[i]
    if args.verbose:
        print "input bam: \t%s\noutput vcf:\t%s" % (in_bam, out_vcf)
    # set up caller log file
    caller_log_file = os.path.join(SCRATCH_DIR, "caller%d.log" % i)
    temp_files.append(caller_log_file)
    caller_log = open(caller_log_file, "w")
    # ----------------------------------------------------------
    # Genotype calling with GATK
    if args.caller == "gatk":
        varcall_cmd = [JAVA, "-jar", "-Xmx%dg" % GATK_MEM,
                       "-XX:ParallelGCThreads=1", GATK, "-T",
                       "UnifiedGenotyper", "-R", ref, "-I", in_bam,
                       "-o", out_vcf]
        varcall_cmd += ["--output_mode", "EMIT_ALL_SITES", "-nt", str(GATK_NT),
                        "-L", interval_file]
        if args.verbose:
            print "\nGATK variant-calling command (space in path not escaped here, but should be fine in actual call command):\n"
            print " ".join(varcall_cmd) + "\n"
        varcall_proc = subprocess.Popen(varcall_cmd, stdout=subprocess.PIPE, stderr=caller_log)
        varcall_proc.communicate()
        varcall_proc_returncode = varcall_proc.returncode
        caller_log.close()

        if args.verbose:
            fin = open(caller_log_file, "r")
            for line in fin:
                print line.strip()

        # check calling was successful
        if varcall_proc_returncode != 0:
            print_caller_failure_message(" ".join(varcall_cmd), caller_log_file)
            exit(1)
        else:
            if args.verbose:
                print "\nVariant calling successful."

    # ----------------------------------------------------------
    # Genotype calling with Freebayes
    elif args.caller == "freebayes":
        fout = open(out_vcf, "w")
        # -----------------------
        # fast-Freebayes, single intervals file
        if FAST_FREEBAYES:
            varcall_cmd = [FREEBAYES, "--fasta-reference", ref, "--targets",
                           interval_file, "--no-indels", "--min-coverage",
                           str(DP_THRESH)]
            varcall_cmd += ["--report-all-haplotype-alleles", "--report-monomorphic", in_bam]
            if args.verbose:
                print "Freebayes variant-calling command:\n" + " ".join(varcall_cmd)

            varcall_proc = subprocess.Popen(varcall_cmd, stdout=fout, stderr=caller_log)
            varcall_proc.communicate()
            varcall_proc_returncode = varcall_proc.returncode
            fout.close()
            caller_log.close()

            if args.verbose:
                fin = open(caller_log_file, "r")
                for line in fin:
                    print line.strip()

            # check calling was successful
            if varcall_proc_returncode != 0:
                print_caller_failure_message(" ".join(varcall_cmd), caller_log_file)
                exit(1)
            else:
                if args.verbose:
                    print "\nVariant calling successful.\n\n"
        else:
        # -----------------------
        # slow Freebayes,calling each site separately
            write_header = True
            caller_log.write("""## SLOW FREEBAYES CALLING
## Calling each variant position separately, as some versions of Freebayes sometimes fail with --targets
""")
            fin = open(interval_file, "r")
            for line in fin:
                bits = line.strip("\n").split("\t")
                bits[1] = int(bits[1])
                bits[2] = int(bits[2]) + 1
                region_str = "%s:%d-%d" % (bits[0], bits[1], bits[2])
                varcall_cmd = [FREEBAYES, "--fasta-reference", ref, "--region",
                               region_str, "--no-indels", "--min-coverage",
                               str(DP_THRESH)]
                varcall_cmd += ["--report-all-haplotype-alleles", "--report-monomorphic", in_bam]
                caller_log.write("FREEBAYES COMMAND:\n" + " ".join(varcall_cmd))
                if args.verbose:
                    print "Freebayes variant-calling command:\n" + " ".join(varcall_cmd)
                varcall_proc = subprocess.Popen(varcall_cmd, stdout=subprocess.PIPE, stderr=caller_log)
                varcall_proc.communicate()

                # write output, this is different to fastfreebayes,
                # because we need to strip the header lines after the first call
                have_header = False
                for line in varcall_proc.stdout:
                    if line.startswith("#"):
                        have_header = True
                        if write_header:
                            fout.write(line)
                    else:
                        if args.verbose:
                            print line.strip("\n")
                        fout.write(line)
                # only switch write_header to False after confirming that header lines have been written once
                if have_header:
                    write_header = False

                # if a call failed
                if varcall_proc.returncode != 0:
                    fout.close()
                    caller_log.close()
                    print_caller_failure_message(" ".join(varcall_cmd), caller_log_file)
                    exit(1)

            fout.close()
            caller_log.close()
            if args.verbose:
                print "\nVariant calling successful.\n\n"

    # ----------------------------------------------------------
    # Genotype calling with VarScan2
    elif args.caller == "varscan":

        # Generate pileup with samtools
        if args.verbose:
            "Generate pileup using samtools.\n"
        pup_file = pup_list[i]
        pup_out = open(pup_file, "w")
        pileup_log_file = os.path.join(SCRATCH_DIR, "samtools_mpileup.log")
        pileup_log = open(pileup_log_file, "w")
        temp_files.append(pileup_log_file)

        # First need to do a pileup
        # not using the -l option in samtools mpileup because it doesn't seem to
        #   use indexed random BAM access, so ends up being much slower
        fin = open(interval_file, "r")
        for line in fin:
            region_str = line.strip("\n")
            sam_cmd = [SAMTOOLS, "mpileup", "-r", region_str, "-B", "-f", ref, in_bam]
            if args.verbose:
                print "pileup command:\n%s" % " ".join(sam_cmd)
            sam_proc = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE, stderr=pileup_log)

            # write samtools mpileup output
            for line in sam_proc.stdout:
                if args.verbose:
                    print line
                bits = line.strip("\n").split("\t")
                if int(bits[3]) < args.dp_threshold:
                    continue
                pup_out.write(line)
            sam_proc.communicate()

            # catch errors
            if sam_proc.returncode != 0:
                pup_out.close()
                pileup_log.close()
                print_caller_failure_message(" ".join(sam_cmd), caller_log_file)
                exit(1)
        pup_out.close()
        pileup_log.close()

        # calling with varscan
        fout = open(out_vcf, "w")
        varcall_cmd = [JAVA, "-XX:ParallelGCThreads=1",
                       "-Xmx%dg" % VARSCAN_MEM, "-jar", VARSCAN,
                       "mpileup2snp", pup_file,
                       "--output-vcf", "--min-coverage", str(DP_THRESH) ]
        if args.verbose:
            print "Varscan command:\n%s" % " ".join(varcall_cmd)
        varcall_proc = subprocess.Popen(varcall_cmd, stdout=fout, stderr=caller_log)
        varcall_proc.communicate()

        fout.close()
        caller_log.close()
        if args.verbose:
            fin = open(caller_log_file, "r")
            for line in fin:
                print line.strip()

        # check calling was successful
        if varcall_proc.returncode != 0:
            fout.close()
            caller_log.close()
            print_caller_failure_message(" ".join(varcall_cmd), caller_log_file)
            exit(1)
        else:
            if args.verbose:
                print "\nVariant calling successful.\n\n"
    # ----------------------------------------------------

if args.verbose:
    print """

Variant-calling finished
"""

#===============================================================================
# Finished variant calling
#===============================================================================






#===============================================================================
# Comparing variant data
#===============================================================================

tsv1 = os.path.join(SCRATCH_DIR, "bam1.tsv")
tsv2 = os.path.join(SCRATCH_DIR, "bam2.tsv")
tsv_list = [tsv1, tsv2]
temp_files += tsv_list

if args.verbose:
    print """
================================================================================
GENOTYPE DATA COMPARISON
================================================================================
"""

#-------------------------------------------------------------------------------
# 1. convert and reorder VCF file to REFERENCE
if args.verbose:
    print "Converting VCF to table"

# convert and re-order
# only need to do this if alternate reference is used
if using_chrom_map:
    _, _, _, MAP_ALT2DEF = get_chrom_data_from_map(CHROM_MAP)

    for i in [0, 1]:
        # don't bother if using cached
        if cached_list[i] == True:
            continue

        # only necessary if not using REFREENCE
        if ref_list[i] == REFERENCE:
            continue

        # converting chrom names
        ftemp = os.path.join(SCRATCH_DIR, "temp_file")
        temp_files.append(ftemp)
        fout = open(ftemp, "w")
        fin = open(vcf_list[i])
        for line in fin:
            # ignore contig information, as we are changing reference
            if line.startswith("##contig="):
                continue
            elif line.startswith("##reference=file://"):
                fout.write("##reference=file://" + REFERENCE + "\n")
            elif line.startswith("#"):
                fout.write(line)
            else:
                bits = line.strip().split("\t")
                bits[0] = MAP_ALT2DEF[bits[0]]
                fout.write("\t".join(bits) + "\n")
        fout.close()

        shutil.copyfile(vcf_list[i], vcf_list[i]+".original")
        temp_files.append(vcf_list[i]+".original")
        # re-sorting VCF file
        sort_vcf_by_chrom_order(ftemp, vcf_list[i], REFERENCE+".fai")



#-------------------------------------------------------------------------------
# extract relevant VCF data into TSV files
for i in [0,1]:
    # if cached
    if BATCH_USE_CACHED:
        if cached_list[i]:
            if args.verbose:
                print "BAM %d has cached genotype data" % (i+1)
            continue
    # otherwise, convert
    in_vcf  = vcf_list[i]
    out_tsv = tsv_list[i]
    in_bam  = bam_list[i]
    ref     = ref_list[i]

    if args.verbose:
        print """
input bam:  %s
reference:  %s
""" % (in_bam, ref)
    VCFtoTSV(in_vcf, out_tsv, args.caller)

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
    out_line = "%s\t%s\t%s\t%s\t%s\n" % (bits[0], bits[1], bits[2], bits[3], bits[7])
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
variants:\t%s
depth threshold: %d
________________________________________

Positions with same genotype:   %d
     breakdown:    hom: %d
                   het: %d
________________________________________

Positions with diff genotype:   %d
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
"""  % (bam1_path, bam2_path, VCF_FILE, DP_THRESH,
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
