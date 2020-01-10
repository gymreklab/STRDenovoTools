#!/bin/bash

VCF=test.vcf.gz #/storage/ileena/ssc-gangstr-denovos/vcf/merged/phase2/phase2_1.filtered.vcf.gz
FAM=/storage/ileena/denovos4/metadata/ssc_4phases_ids.ped

./src/CookieMonSTR \
    --strvcf ${VCF} \
    --fam ${FAM} --family 14686 \
    --include-invariant \
    --output-all-loci \
    --gangstr \
    --require-all-children \
    --out debug \
    --naive --min-num-encl-child 3 --max-perc-encl-parent 0.05 --min-encl-match 0.9 --min-total-encl 10 --debug --region chr1:145928571-145928571
