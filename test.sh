#!/bin/bash

# VCF=/storage/ileena/ssc-gangstr-denovos/vcf/phase2/11963_22.filtered.vcf.gz
VCF=/storage/ileena/tmp/test_vcf.vcf.gz
FAM=/storage/ileena/denovos4/metadata/ssc_4phases_ids.ped
MUT=/storage/ileena/denovos4/denovos_mean_priors/mean_str_mutrates_hg38_ver16.bed
OUT=/storage/ileena/tmp/test_denovos

MAXALLELES=100
MINCOV=10
MINSCORE=0.9
MINSPANCOV=10
MINSUPPREADS=2
PTHRESH=0.8

/home/ileena/STRDenovoTools/src/CookieMonSTR \
    --strvcf ${VCF} \
    --fam ${FAM} \
    --max-num-alleles ${MAXALLELES} \
    --posterior-threshold ${PTHRESH} \
    --round-alleles --include-invariant \
    --output-all-loci \
    --mutation-models ${MUT} \
    --default-beta 0 --default-central 0 \
    --gangstr \
    --require-all-children \
    --out ${OUT}
echo "Results: "${OUT}
grep 28161724 /storage/ileena/tmp/test_denovos.all_mutations.tab

#### Expected result:
#### chrom	pos	period	prior	family	child	phenotype	posterior	newallele	mutsize	inparents	poocase	isnew	case_count	ctrl_count	unk_count	child_gt	mat_gt	pat_gt
#### chr22	28161724	5	1.92576e-05	11963	SSC03646	1	0.800138	-5	-5	0	4	1	0	0	0	32764	0	32607
#### chrom	pos	end	period	num_alleles_bylength	num_alleles_byseq	het_by_length	het_by_seq	total_children	total_mutations	total_mutation_rate	affected_children	affected_mutations	affected_new_mutations	affected_mutation_rate	unaffected_children	unaffected_mutations	unaffected_new_mutationunaffected_mutation_rate	p-value	children_with_mutations
#### chr22	28161724	28161753	5	3	3	0.595041	0.595041	2	1	0.5	1	0	0	0	1	1	1	1	1	11963:1:-5:0,0,0
