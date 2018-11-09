#!/bin/bash

chrom=$1
start=$2
family=$3

./src/ExamineDenovo \
    --strvcf /storage/mgymrek/ssc-imputation/filtered_vcfs/hipstr.chr${chrom}.allfilters.vcf.gz \
    --fam /home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam \
    --locus ${chrom}:${start} \
    --family ${family}
