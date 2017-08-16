#!/bin/bash

./src/STRDenovoTools \
    --strvcf /storage/s1saini/hipstr_genomewide/chr21/hipstr_calls_21.vcf.gz \
    --fam /home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam \
    --max-num-alleles 25 \
    --require-all-children \
    --require-num-children 2 \
    --min-coverage 10 \
    --min-score 0.9 \
    --posterior-threshold 0.9 \
    --out test \
    --combine-alleles-by-length

#     --region 21:32734150-32734151 \
#     --period 3 \
