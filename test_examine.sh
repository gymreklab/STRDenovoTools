#!/bin/bash

chrom=$1
start=$2
family=$3

./src/ExamineDenovo \
    --strvcf /storage/s1saini/hipstr_genomewide/chr${chrom}/hipstr_calls_${chrom}.vcf.gz \
    --fam /home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam \
    --locus ${chrom}:${start} \
    --family ${family}
