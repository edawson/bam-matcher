#!/bin/bash
## This script takes in the required arguments for BAM-Matcher.conf, fills out the config,
## and passes it to bam-matcher


firstbam=$1
secondbam=$2
reportname=$(basename ${firstbam} .bam).$(basename $secondbam .bam).bam-matcher-report.txt
ref=$3
bmdir=$4
cache=$5
vcf=${bmdir}1kg.exome.highAF.7550.vcf

# Fill in the cache 
# pass @ as the sed delimiter (thanks https://stackoverflow.com/questions/9366816/sed-fails-with-unknown-option-to-s-error)
sed -i "s@_CACHETARGET_@${cache}@g" ${bmdir}/bam-matcher.default.conf
# Fill in the reference
sed -i "s@_FASTATARGET_@${ref}@g" ${bmdir}/bam-matcher.default.conf
# Fill in the VCF
sed -i "s@_VCFTARGET_@${vcf}@g" ${bmdir}/bam-matcher.default.conf

cp ${bmdir}/bam-matcher.default.conf ${bmdir}/bam-matcher.conf


${bmdir}/bam-matcher.py --bam1 ${firstbam} --bam2 ${secondbam} -o ${reportname}
