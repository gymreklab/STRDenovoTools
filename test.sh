#!/bin/bash

./src/STRDenovoTools \
    --strvcf /storage/s1saini/hipstr_genomewide/chr22/hipstr_calls_22.vcf.gz \
    --fam /home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam \
    --max-num-alleles 25 \
    --require-all-children \
    --require-num-children 2 \
    --min-coverage 10 \
    --min-score 0.9 \
    --posterior-threshold 0.9 \
    --out test \
    --combine-alleles-by-length \
    --min-span-cov 10 \
    --min-supp-reads 2 
#    --region 22:17631413-17631426
#    --region 22:39350442-39350443 \
#    --family 13924

#    --region 21:32734150-32734151
#     --period 3 \

#    --family 14658 \
#    --region 21:14606632-14606633 \
