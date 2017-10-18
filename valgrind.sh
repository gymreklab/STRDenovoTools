#!/bin/bash

source ~/workspace/ssc-imputation/denovos2/params.sh

chrom=20
valgrind --leak-check=full \
    /home/mgymrek/workspace/STRDenovoTools/src/STRDenovoTools \
    --strvcf /storage/s1saini/hipstr_rerun/chr${chrom}/hipstr.chr${chrom}.with.1kg.vcf.gz \
    --fam /home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam \
    --max-num-alleles ${MAXALLELES} \
    --require-all-children \
    --require-num-children 2 \
    --min-coverage ${MINCOV} \
    --min-score ${MINSCORE} \
    --min-span-coverage ${MINSPANCOV} \
    --min-supp-reads ${MINSUPPREADS} \
    --posterior-threshold ${PTHRESH} \
    --out test \
    --region 20:52467355-52477366 \
    --combine-alleles-by-length 
