#!/bin/bash

./src/ExamineDenovo \
    --strvcf /storage/s1saini/hipstr_genomewide/chr21/hipstr_calls_21.vcf.gz \
    --fam /home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam \
    --family 13889 \
    --locus 19:529942
