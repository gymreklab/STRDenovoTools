# MonSTR: de novo tandem repeat mutations discovery.

#### Authors: * Ileena Mitra <ilmitra@ucsd.edu> * Melissa Gymrek <mgymrek@ucsd.edu> <br>  
License: MIT

Toolkit for calling and analyzing de novo tandem repeat mutations.
This repository is under construction but is planned for official release soon.

The core code for calling de novos is modified from that originally written by Thomas Willems in the [HipSTR](https://github.com/tfwillems/HipSTR) repository.

## Docker
The Dockerfile in this container sets up a Docker with GangSTR and dependencies installed. You can pull this image from https://hub.docker.com/r/ileena/cookiemonstr.

## Install

You can install MonSTR using the following commands:

```
git clone https://github.com/gymreklab/STRDenovoTools  
cd STRDenovoTools
mkdir build
cd build
cmake ..
make
sudo make install
```

## Usage
```
MonSTR --help
```

Demo example:
```
MonSTR \
    --strvcf test.vcf.gz \
    --fam family.ped \
    --max-num-alleles 100 \
    --round-alleles --include-invariant \
    --output-all-loci \
    --gangstr \
    --require-all-children \
    --naive \
    --naive-expansions-frr 3,8 \
    --min-num-encl-child 3 --max-perc-encl-parent 0.05 --min-encl-match 0.9 --min-total-encl 10 \
    --out test_results.txt
```

The file used for --fam should be in pedigree file format (see more information here: https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).

The output file (--out) will be a tab delimited text file with the following columns:
chromosome,	start position,	TR repeat unit,	prior used for model, family ID, child ID, phenotype, posterior mutation probability, de novo allele observed, mutation size, observed in parents, parent of orgin, is rare allele, number of times observed in cases, in controls,	in unknown phenotype, child's genotype, mother's genotype, father's genotype, number of enclosing reads in child, in mother, in father, enclosing reads in parents, allele type from mother, allele type from father.

## Cite 
If you use this method for your research, please cite:
Mitra, et al. 2020 https://www.biorxiv.org/content/10.1101/2020.03.04.974170v1
