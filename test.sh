#!/bin/bash

chrom=10

./src/STRDenovoTools \
    --strvcf /storage/s1saini/hipstr_rerun/chr${chrom}/hipstr.chr${chrom}.with.1kg.vcf.gz \
    --fam /home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam \
    --mutation-models /storage/mgymrek/ssc-denovos/mutea-results/predicted_str_mutrates_GRCh37.bed \
    --default-prior -5.0 \
    --max-num-alleles 25 \
    --require-all-children \
    --require-num-children 2 \
    --min-coverage 10 \
    --min-score 0.8 \
    --posterior-threshold 0.9 \
    --out test \
    --combine-alleles-by-length \
    --min-span-cov 10 \
    --min-supp-reads 2 \
    --region ${chrom}:383896-383896 \
    --round-alleles \
    --include-invariant \
    --output-all-loci
#    --region ${chrom}:125907115-125907115 \
# 20:52477355-52477366 \
#    --region 22:17631413-17631426
#    --region 22:39350442-39350443 \
#    --family 13924

#    --region 21:32734150-32734151
#     --period 3 \

#    --family 14658 \
#    --region 21:14606632-14606633 \
