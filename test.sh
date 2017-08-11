#!/bin/bash

./src/STRDenovoTools \
    --strvcf /storage/s1saini/hipstr_genomewide/chr21/hipstr_calls_21.vcf.gz \
    --fam /home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam \
    --max-num-alleles 25 \
    --posterior-threshold 0.9 \
    --out test

#     --region 21:32734150-32734151 \
