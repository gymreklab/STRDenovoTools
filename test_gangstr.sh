#!/bin/bash

GANGSTRFILE=/storage/resources/datasets/PlatinumGenomesDbGaP/gangstr/pathogenic_gangstr_012819_filtered.vcf.gz
MUTMODEL=/storage/mgymrek/ssc-denovos/mutea-results/predicted_str_mutrates_GRCh37.bed

MINCOV=10
MINSCORE=0.9
MINSPANCOV=10
MINSUPPREADS=2
PTHRESH=0.8
MAXALLELES=100

./src/STRDenovoTools \
    --strvcf ${GANGSTRFILE} \
    --fam /storage/ileena/pg/denovos/metadata/fam/1 \
    --max-num-alleles ${MAXALLELES} \
    --min-coverage ${MINCOV} \
    --min-score ${MINSCORE} \
    --posterior-threshold ${PTHRESH} \
    --round-alleles --include-invariant \
    --output-all-loci \
    --debug \
    --out test \
    --gangstr
