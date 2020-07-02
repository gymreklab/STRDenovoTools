#!/bin/bash

### test naive for expansions
./src/CookieMonSTR \
    --chrX \
    --region chrX:0-34026969 \
    --strvcf /storage/ileena/pg/gangstr_chrX/pg_quad_NA12879_82_hg19_chrX.sorted.vcf.gz \
    --fam family.ped \
    --max-num-alleles 100 \
    --round-alleles --include-invariant \
    --output-all-loci \
    --gangstr \
    --require-all-children \
    --naive \
    --naive-expansions-frr 3,8 \
    --min-num-encl-child 3 --max-perc-encl-parent 0.05 --min-encl-match 0.9 --min-total-encl 10 \
    --out test_naive_results --verbose --debug

### test model based
# ./src/CookieMonSTR \
#     --chrX \
#     --region chrX:0-34026969 \
#     --strvcf /storage/ileena/pg/gangstr_chrX/pg_quad_NA12879_82_hg19_chrX.sorted.vcf.gz \
#     --fam family.ped \
#     --max-num-alleles 100 \
#     --round-alleles --include-invariant \
#     --output-all-loci \
#     --gangstr \
#     --require-all-children \
#     --default-prior -3 \
#     --posterior-threshold 0.5 \
#     --min-num-encl-child 3 --max-perc-encl-parent 0.05 --min-encl-match 0.9 --min-total-encl 10 \
#     --out test_model_results --verbose --debug
