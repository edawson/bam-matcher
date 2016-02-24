#!/usr/bin/env python
'''
Created on 23/02/2016

@author: Paul Wang (ppswang@gmail.com)

Generate example data set from CML samples

- need to combine from multiple samples as a way to annonymise patient ID

'''


import os
import HTSeq
import vcf
import sys

# ==============================
def change_RG(aligned_read, new_RG):
    bits = aligned_read.get_sam_line().split("\t")
    for i in range(len(bits)):
        if bits[i].startswith("RG:"):
            bits[i] = "RG:Z:%s" % new_RG
            break
    return "\t".join(bits)


# ==============================





BAMDIR = "/data/sacgf/molpath/data/aligned/CML_WES/bam_files"
BAM_FILES = {}
bam_idx = 0
for f in os.listdir(BAMDIR):
    if f.endswith(".bam"):
        BAM_FILES[bam_idx] = f
        bam_idx += 1

# number of variants to get:
VAR_N = 300

VCF_FILE = "1KG_1500_exon_variants_noX.vcf"
invcf = vcf.Reader(filename=VCF_FILE)

SAM1 = "sample1.sam"
SAM2 = "sample2.sam"

HEADER = """@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:1	LN:249250621
@SQ	SN:2	LN:243199373
@SQ	SN:3	LN:198022430
@SQ	SN:4	LN:191154276
@SQ	SN:5	LN:180915260
@SQ	SN:6	LN:171115067
@SQ	SN:7	LN:159138663
@SQ	SN:8	LN:146364022
@SQ	SN:9	LN:141213431
@SQ	SN:10	LN:135534747
@SQ	SN:11	LN:135006516
@SQ	SN:12	LN:133851895
@SQ	SN:13	LN:115169878
@SQ	SN:14	LN:107349540
@SQ	SN:15	LN:102531392
@SQ	SN:16	LN:90354753
@SQ	SN:17	LN:81195210
@SQ	SN:18	LN:78077248
@SQ	SN:19	LN:59128983
@SQ	SN:20	LN:63025520
@SQ	SN:21	LN:48129895
@SQ	SN:22	LN:51304566
@SQ	SN:X	LN:155270560
@SQ	SN:Y	LN:59373566
@SQ	SN:MT	LN:16569
@SQ	SN:GL000207.1	LN:4262
@SQ	SN:GL000226.1	LN:15008
@SQ	SN:GL000229.1	LN:19913
@SQ	SN:GL000231.1	LN:27386
@SQ	SN:GL000210.1	LN:27682
@SQ	SN:GL000239.1	LN:33824
@SQ	SN:GL000235.1	LN:34474
@SQ	SN:GL000201.1	LN:36148
@SQ	SN:GL000247.1	LN:36422
@SQ	SN:GL000245.1	LN:36651
@SQ	SN:GL000197.1	LN:37175
@SQ	SN:GL000203.1	LN:37498
@SQ	SN:GL000246.1	LN:38154
@SQ	SN:GL000249.1	LN:38502
@SQ	SN:GL000196.1	LN:38914
@SQ	SN:GL000248.1	LN:39786
@SQ	SN:GL000244.1	LN:39929
@SQ	SN:GL000238.1	LN:39939
@SQ	SN:GL000202.1	LN:40103
@SQ	SN:GL000234.1	LN:40531
@SQ	SN:GL000232.1	LN:40652
@SQ	SN:GL000206.1	LN:41001
@SQ	SN:GL000240.1	LN:41933
@SQ	SN:GL000236.1	LN:41934
@SQ	SN:GL000241.1	LN:42152
@SQ	SN:GL000243.1	LN:43341
@SQ	SN:GL000242.1	LN:43523
@SQ	SN:GL000230.1	LN:43691
@SQ	SN:GL000237.1	LN:45867
@SQ	SN:GL000233.1	LN:45941
@SQ	SN:GL000204.1	LN:81310
@SQ	SN:GL000198.1	LN:90085
@SQ	SN:GL000208.1	LN:92689
@SQ	SN:GL000191.1	LN:106433
@SQ	SN:GL000227.1	LN:128374
@SQ	SN:GL000228.1	LN:129120
@SQ	SN:GL000214.1	LN:137718
@SQ	SN:GL000221.1	LN:155397
@SQ	SN:GL000209.1	LN:159169
@SQ	SN:GL000218.1	LN:161147
@SQ	SN:GL000220.1	LN:161802
@SQ	SN:GL000213.1	LN:164239
@SQ	SN:GL000211.1	LN:166566
@SQ	SN:GL000199.1	LN:169874
@SQ	SN:GL000217.1	LN:172149
@SQ	SN:GL000216.1	LN:172294
@SQ	SN:GL000215.1	LN:172545
@SQ	SN:GL000205.1	LN:174588
@SQ	SN:GL000219.1	LN:179198
@SQ	SN:GL000224.1	LN:179693
@SQ	SN:GL000223.1	LN:180455
@SQ	SN:GL000195.1	LN:182896
@SQ	SN:GL000212.1	LN:186858
@SQ	SN:GL000222.1	LN:186861
@SQ	SN:GL000200.1	LN:187035
@SQ	SN:GL000193.1	LN:189789
@SQ	SN:GL000194.1	LN:191469
@SQ	SN:GL000225.1	LN:211173
@SQ	SN:GL000192.1	LN:547496
@SQ	SN:NC_007605	LN:171823
"""

# write header
sam1_out = open(SAM1, "w")
sam2_out = open(SAM2, "w")

sam1_out.write(HEADER)
sam1_out.write("@RG\tID:sample1\tSM:sample1\tPL:Illumina\n")
sam2_out.write(HEADER)
sam2_out.write("@RG\tID:sample2\tSM:sample2\tPL:Illumina\n")

sample_idx = 0
sample_N = len(BAM_FILES)
var_ct = 0
for var in invcf:
    # write SAM1

    print var.CHROM, var.POS, var.REF, var.ALT

    SAM1_done = False
    SAM2_done = False

    inbam1 = HTSeq.BAM_Reader(os.path.join(BAMDIR, BAM_FILES[sample_idx]))
    sample_idx += 1
    sample_idx = sample_idx % sample_N
    inbam2 = HTSeq.BAM_Reader(os.path.join(BAMDIR, BAM_FILES[sample_idx]))
    sample_idx += 1
    sample_idx = sample_idx % sample_N

    SAM1_reads = []
    SAM2_reads = []

    for read in inbam1.fetch(region="%s:%d-%d" % (var.CHROM, var.POS, var.POS)):
        if read.pcr_or_optical_duplicate:
            continue
        if read.proper_pair == False:
            continue
        SAM1_reads.append(read)

    for read in inbam2.fetch(region="%s:%d-%d" % (var.CHROM, var.POS, var.POS)):
        if read.pcr_or_optical_duplicate:
            continue
        if read.proper_pair == False:
            continue
        SAM2_reads.append(read)

    # don't write anything if neither samples are sufficiently covered
    if len(SAM1_reads) < 10 or len(SAM2_reads) < 10:
        continue

    if len(SAM1_reads) > 100 or len(SAM2_reads) > 100:
        continue

    print "sample 1 writing %d reads" % len(SAM1_reads)
    for ct, read in enumerate(SAM1_reads):
        # print "%d/%d" % (ct+1, len(SAM1_reads))
        sys.stdout.write("\b\b\b\b\b\b\b\b\b\b\b\b%d/%d" % (ct+1, len(SAM1_reads)))
        sys.stdout.flush()

        # need to replace read group
        sam1_out.write(change_RG(read, "sample1") + "\n")
        # bits = read.get_sam_line().split("\t")
        # for i in range(len(bits)):
        #     if bits[i].startswith("RG:"):
        #         bits[i] = "RG:Z:sample1"
        # sam1_out.write("\t".join(bits) + "\n")

        # get paired mate
        mate_pos = read.mate_start
        mate_found = False
        for read2 in inbam1.fetch(region="%s:%d-%d" % (mate_pos.chrom, mate_pos.pos+1, mate_pos.pos+1)):
            if read2.read.name == read.read.name:
                mate_found = True
                sam1_out.write(change_RG(read2, "sample1") + "\n")
                break
        if mate_found == False:
            print "ERROR: Cannot find mate read"
            exit(1)

    print "\b\b\b\b\b\b\b\b\b\b\b\bsample 2 writing %d reads" % len(SAM2_reads)
    for ct, read in enumerate(SAM2_reads):
        # print "%d/%d" % (ct+1, len(SAM2_reads))
        sys.stdout.write("\b\b\b\b\b\b\b\b\b\b\b\b%d/%d" % (ct+1, len(SAM2_reads)))
        sys.stdout.flush()

        # need to replace read group
        sam2_out.write(change_RG(read, "sample2") + "\n")

        # get paired mate
        mate_pos = read.mate_start
        mate_found = False
        for read2 in inbam2.fetch(region="%s:%d-%d" % (mate_pos.chrom, mate_pos.pos+1, mate_pos.pos+1)):
            if read2.read.name == read.read.name:
                mate_found = True
                sam2_out.write(change_RG(read2, "sample2") + "\n")
                break
        if mate_found == False:
            print "ERROR: Cannot find mate read"
            exit(1)

    var_ct += 1
    print "\b\b\b\b\b\b\b\b\b\b\b\bwrote %d sites"  % var_ct
    if var_ct >= VAR_N:
        break

sam1_out.close()
sam2_out.close()


os.system("samtools view -bS sample1.sam | samtools sort - > sample1.bam")
os.system("samtools index sample1.bam")
os.system("samtools view -bS sample2.sam | samtools sort - > sample2.bam")
os.system("samtools index sample2.bam")
