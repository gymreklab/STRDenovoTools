#!/bin/bash

VCF=/storage/ileena/ssc-gangstr-denovos/vcf/merged/phase2/phase2_1.filtered.vcf.gz
FAM=/storage/ileena/denovos4/metadata/ssc_4phases_ids.ped

./src/CookieMonSTR \
    --strvcf ${VCF} \
    --fam ${FAM} --family 14686 \
    --include-invariant \
    --output-all-loci \
    --gangstr \
    --require-all-children \
    --out test \
    --naive --min-num-encl-child 2 --max-num-encl-parent 0
