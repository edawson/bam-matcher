#!/bin/bash
## This script takes in the required arguments for BAM-Matcher.conf, fills out the config,
## and passes it to bam-matcher

bmdir=$4

firstbam=$1
secondbam=$2
reportname=$(basename ${firstbam} .bam).$(basename $secondbam .bam).bam-matcher-report.txt
ref=$3
cache=$4
vcf=/app/bam-matcher/1kg.exome.highAF.7550.vcf

# Fill in the cache 
sed -i "s/_CACHETARGET_/${cache}/g" ${bmdir}/bam-matcher.default.conf
# Fill in the reference
sed -i "s/_FASTATARGET_/${ref}/g" ${bmdir}/bam-matcher.default.conf
# Fill in the VCF
sed -i "s/_VCFTARGET_/${vcf}/g" ${bmdir}/bam-matcher.default.conf

cp ${bmdir}/bam-matcher.default.conf ${bmdir}/bam-matcher.conf


${bmdir} bam-matcher.py --bam1 ${firstbam} --bam2 ${secondbam} -o ${reportname}
