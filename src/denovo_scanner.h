/*
Copyright (C) 2017 Melissa Gymrek <mgymrek@ucsd.edu>

This file is part of STRDenovoTools.

STRDenovoTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

STRDenovoTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with STRDenovoTools.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TRIO_DENOVO_SCANNER_H_
#define TRIO_DENOVO_SCANNER_H_

#include <assert.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "src/mutation_priors.h"
#include "src/options.h"
#include "src/pedigree.h"
#include "src/vcf_reader.h"

using namespace std;

class DenovoResult {
 public:
  DenovoResult(const std::string& family_id,
	       const std::string& mother_id,
	       const std::string& father_id,
	       const std::string& child_id,
	       const int& phenotype,
	       const double& total_ll_no_mutation,
	       const double& total_ll_one_denovo,
	       const double& log10prior);
  virtual ~DenovoResult();

  const int& get_phenotype() const {return phenotype_;}
  const double& get_posterior() const {return posterior_;}
  const double& get_prior() const {return log10_prior_mutation_;}
  const std::string& get_child_id() const {return child_id_;}
  const std::string& get_mother_id() const {return mother_id_;}
  const std::string& get_father_id() const {return father_id_;}
  const std::string& get_family_id() const {return family_id_;}

 private:
  void CalculatePosterior();
  std::string family_id_;
  std::string mother_id_;
  std::string father_id_;
  std::string child_id_;
  int phenotype_;
  double total_ll_no_mutation_;
  double total_ll_one_denovo_;
  double posterior_;
  double log10_prior_mutation_;
};

class TrioDenovoScanner {
 private:
  static std::string START_KEY, END_KEY, PERIOD_KEY;
  PedigreeSet pedigree_set_;
  Options options_;
  ofstream locus_summary_;
  ofstream all_mutations_file_;

 public:
  TrioDenovoScanner(const PedigreeSet& pedigree_set,
		    const Options& options)
    : pedigree_set_(pedigree_set), options_(options),
    locus_summary_(options.outprefix + ".locus_summary.tab"),
    all_mutations_file_(options.outprefix + ".all_mutations.tab") {
    locus_summary_ << "chrom\tpos\tend\tperiod\tref_allele_size\t"
		   << "num_alleles_bylength\tnum_alleles_byseq\thet_by_length\thet_by_seq\t"
		   << "total_children\ttotal_mutations\ttotal_mutation_rate\t"
		   << "affected_children\taffected_mutations\taffected_new_mutations\taffected_mutation_rate\t"
		   << "unaffected_children\tunaffected_mutations\tunaffected_new_mutations\tunaffected_mutation_rate\t"
		   << "p-value\tchildren_with_mutations\n";
    all_mutations_file_ << "chrom\tpos\tperiod\tprior\tfamily\tchild\tphenotype\t"
       << "posterior\tnewallele\tmutsize\tinparents\tpoocase\tisnew\tcase_count\t"
       << "ctrl_count\tunk_count\tchild_gt\tmat_gt\tpat_gt\n";
  }
  virtual ~TrioDenovoScanner();

  void scan(VCF::VCFReader& strvcf, MutationPriors& priors);
  void naive_scan(VCF::VCFReader& strvcf);
  bool check_mutation(const VCF::Variant& str_variant,
		      const int32_t& mother_ind,
		      const int32_t& father_ind,
		      const int32_t& child_ind);
  void summarize_results(std::vector<DenovoResult>& dnr,
			 VCF::Variant& str_variant);
  /*
    poocase: describes inheritance pattern
    -1: unknown
    1: Mendelian (no denovo)
    2, 21: New allele from father
       2: allele a in mother only, allele b not in father
       21: allele b in mother only, allele a not in father
    3, 31: new allele from mother
       3: allele a in father only, allele b not in mother
       31: allele b in father only, allele a not in mother
    4, 41: unclear (should this happen?)
       4: allele a not in either
       41: allele b not in either
    5: none of the above (shouldn't happen)
   */
  void GetMutationInfo(const VCF::Variant& variant, const std::string& mother_id,
		       const std::string& father_id, const std::string& child_id,
		       std::string* new_allele, std::string* new_allele_raw, std::string* mut_size,
		       bool* new_allele_in_parents, int* poocase,
           std::string* child_gt, std::string* mat_gt, std::string* pat_gt);

};


#endif
