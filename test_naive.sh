#!/bin/bash

#VCF=test.vcf.gz #/storage/ileena/ssc-gangstr-denovos/vcf/merged/phase2/phase2_1.filtered.vcf.gz
#FAM=test.fam #/storage/ileena/denovos4/metadata/ssc_4phases_ids.ped
#VCF=/storage/ileena/ssc-gangstr-denovos/vcf/merged/phase2/phase2_1.filtered.vcf.gz
#FAM=/storage/ileena/denovos4/metadata/ssc_4phases_ids.ped
#VCF=test.vcf.gz
#FAM=test.fam
 #--family 14686 \


FAM=/storage/ileena/denovos4/metadata/ssc_4phases_ids.ped

# Naive method
for chrom in $(seq 1 22)
do
#chrom=22
#    VCF=/storage/ileena/ssc-gangstr-denovos/vcf/merged/phase2/phase2_${chrom}.filtered.PASS_only.vcf.gz 
    #VCF=/storage/ileena/ssc-gangstr-denovos/vcf/phase2/11963_${chrom}.sorted.vcf.gz
    VCF=/storage/ileena/ssc-gangstr-denovos/vcf/merged/phase2/phase2_${chrom}.filtered.vcf.gz
    ./src/CookieMonSTR \
	--strvcf ${VCF} \
	--fam ${FAM}  \
	--include-invariant \
	--output-all-loci \
	--gangstr \
	--require-all-children \
	--out test2-naive-chr${chrom} \
	--naive --min-num-encl-child 3 --max-perc-encl-parent 0.05 --min-encl-match 0.9 --min-total-encl 10 --naive-expansions-frr 3,8 #--debug --region chr1:6733191-6733191
done

# Model based method
#./src/CookieMonSTR \
#    --strvcf ${VCF} \
#    --fam ${FAM} \
#    --include-invariant \
#    --output-all-loci \
#    --gangstr \
#    --require-all-children \
#    --out test-model \
#    --min-num-encl-child 3 --max-perc-encl-parent 0.05 --min-encl-match 0.9 --min-total-encl 10 --combine-alleles #--debug --region chr1:231069412-231069412
