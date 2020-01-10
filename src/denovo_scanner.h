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
	       const int32_t mother_ind,
	       const int32_t father_ind,
	       const int32_t child_ind,
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
  const std::string& get_child_gt() const {return child_gt_;}
  const std::string& get_mat_gt() const {return mat_gt_;}
  const std::string& get_pat_gt() const {return pat_gt_;}
  const int get_new_allele() const {return new_allele_;}
  const int get_mut_size() const {return mut_size_;}
  const int get_poocase() const {return poocase_;}
  const bool get_new_allele_in_parents() const {return new_allele_in_parents_;}
  const int get_encl_reads_child() const {return encl_reads_child_;}
  const int get_encl_reads_parent() const {return encl_reads_parent_;}
  const int get_encl_reads_mother() const {return encl_reads_mother_;}
  const int get_encl_reads_father() const {return encl_reads_father_;}

  const void zero_posterior() {posterior_ = 0;}
  const void set_child_gt(const std::string& child_gt) {child_gt_ = child_gt;}
  const void set_mat_gt(const std::string& mat_gt) {mat_gt_ = mat_gt;}
  const void set_pat_gt(const std::string& pat_gt) {pat_gt_ = pat_gt;}

  void GetMutationInfo(const Options& options, const VCF::Variant& variant,
		       bool* filter_mutation);
  // 0=unknown, 1=shorter, 2=longer
  void TestTransmission(int* long_mother, int* long_father);
  int GetTrans(const int& child, const int& p1, const int& p2);

 private:
  static std::string PERIOD_KEY;
  void CalculatePosterior();
  std::string family_id_;
  std::string mother_id_;
  std::string father_id_;
  std::string child_id_;
  int32_t mother_ind_;
  int32_t father_ind_;
  int32_t child_ind_;
  int phenotype_;
  double total_ll_no_mutation_;
  double total_ll_one_denovo_;
  double posterior_;
  double log10_prior_mutation_;

  // Genotype info
  std::string child_gt_;
  std::string mat_gt_;
  std::string pat_gt_;
  int child_gt_a_, child_gt_b_;
  int mat_gt_a_, mat_gt_b_;
  int pat_gt_a_, pat_gt_b_;
  void GetRepcn(const VCF::Variant& variant, const int32_t& sample_ind,
		int* repcn_a, int* repcn_b);
  int GetMutSize(const int& new_allele, const int& a1, const int& a2);
  int GetMutSize(const int& new_allele, const int& a1, const int& a2, const int& a3, const int& a4);
  void GetEnclosing(const std::string& enclstring, int& new_allele,
		    const int32_t& repcn_a, const int32_t& repcn_b,
		    int* encl_newallele, int* encl_total, int* encl_match);
		    

  // New allele info and POO
  int new_allele_ = 0;
  int mut_size_ = 0;
  int poocase_ = 0;
  bool new_allele_in_parents_ = false;

  // Enclosing read info
  int encl_reads_child_ = 0; // num encl reads matching new allele in child
  int encl_reads_parent_ = 0 ; // num in parent the allele was inherited from
  int encl_reads_mother_ = 0;
  int encl_reads_father_ = 0;
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
			<< "ctrl_count\tunk_count\tchild_gt\tmat_gt\tpat_gt\t"
			<< "encl_child\tencl_mother\tencl_father\tencl_parent"
			<< "long_mother\tlong_father\n";
  }
  virtual ~TrioDenovoScanner();

  void scan(VCF::VCFReader& strvcf, MutationPriors& priors);
  void naive_scan(VCF::VCFReader& strvcf);
  void summarize_results(std::vector<DenovoResult>& dnr,
			 VCF::Variant& str_variant);
};


#endif
