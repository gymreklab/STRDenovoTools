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

#include <stdlib.h>

#include <cfloat>
#include <sstream>
#include <vector>

#include "src/denovo_allele_priors.h"
#include "src/denovo_scanner.h"
#include "src/mathops.h"
#include "src/mutation_model.h"
#include "src/vcf_input.h"

std::string TrioDenovoScanner::START_KEY   = "START";
std::string TrioDenovoScanner::END_KEY     = "END";

TrioDenovoScanner::~TrioDenovoScanner() {}

void TrioDenovoScanner::scan(VCF::VCFReader& strvcf) {
  VCF::Variant str_variant;
  while (strvcf.get_next_variant(str_variant)) {
    // Initial checks
    int num_alleles = str_variant.num_alleles();
    if (num_alleles <= 1) {
      continue;
    }
    if (str_variant.num_samples() == str_variant.num_missing()) {
      continue;
    }
    if (num_alleles > options_.max_num_alleles) {
      continue;
    }
    
    // Get locus info
    int32_t start; str_variant.get_INFO_value_single_int(START_KEY, start);
    int32_t end; str_variant.get_INFO_value_single_int(END_KEY, end);
    std::stringstream ss;
    ss << "Processing STR region " << str_variant.get_chromosome() << ":" << start << "-" << end
       << " with " << num_alleles << " alleles";
    PrintMessageDieOnError(ss.str(), M_PROGRESS);

    // Set up
    UnphasedGL unphased_gls(str_variant);
    MutationModel mut_model(str_variant);
    DiploidGenotypePrior* dip_gt_priors;
    std::vector<NuclearFamily> families = pedigree_set_.get_families();
    if (options_.use_pop_priors) {
      dip_gt_priors = new PopulationGenotypePrior(str_variant, families);
    } else {
      dip_gt_priors = new UniformGenotypePrior(str_variant, families);
    }
    const double LOG_ONE_FOURTH = -log10(4);
    const double LOG_TWO        = log10(2);

    // Scan each family
    for (auto family_iter = families.begin(); family_iter != families.end(); family_iter++) {
      bool scan_for_denovo = unphased_gls.has_sample(family_iter->get_mother()) && unphased_gls.has_sample(family_iter->get_father());
      // Scan each child in the family
      for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); ++child_iter) {
	if (!scan_for_denovo || !unphased_gls.has_sample(*child_iter)) {
	  continue;
	}
	// To accelerate computations, we will ignore configurations that make a neglible contribution (< 0.01%) to the total LL
	// For mutational scenarios, we aggregate 1/4*A^2*(A+1)^2*4*2*A values. Therefore, to ignore a configuration with LL=X:
	// X*A^3*(A+1)^2*2 < TOTAL/10000;
	// logX < log(TOTAL) - log(10000*A^3*(A+1)^2*2) = log(TOTAL) - [log(10000) + 3log(A) + 2log(A+1) + log(2)];
	double MIN_CONTRIBUTION   = 4 + 3*log10(num_alleles) + 2*log(num_alleles+1) + LOG_TWO;
	double ll_no_mutation_max = -DBL_MAX/2, ll_no_mutation_total = 0.0;
	double ll_one_denovo_max  = -DBL_MAX/2, ll_one_denovo_total  = 0.0;
	double ll_one_other_max   = -DBL_MAX/2, ll_one_other_total   = 0.0;
	int mother_gl_index       = unphased_gls.get_sample_index(family_iter->get_mother());
	int father_gl_index       = unphased_gls.get_sample_index(family_iter->get_father());
	int child_gl_index        = unphased_gls.get_sample_index(*child_iter);

	// Iterate over all maternal genotypes
	for (int mat_i = 0; mat_i < num_alleles; mat_i++){
	  for (int mat_j = 0; mat_j <= mat_i; mat_j++){
	    double mat_ll = dip_gt_priors->log_unphased_genotype_prior(mat_j, mat_i, family_iter->get_mother()) + unphased_gls.get_gl(mother_gl_index, mat_j, mat_i);

	    // Iterate over all paternal genotypes
	    for (int pat_i = 0; pat_i < num_alleles; pat_i++){
	      for (int pat_j = 0; pat_j <= pat_i; pat_j++){
		double pat_ll    = dip_gt_priors->log_unphased_genotype_prior(pat_j, pat_i, family_iter->get_father()) + unphased_gls.get_gl(father_gl_index, pat_j, pat_i);
		double config_ll = mat_ll + pat_ll + LOG_ONE_FOURTH;

		// Iterate over all 4 possible inheritance patterns for the child
		for (int mat_index = 0; mat_index < 2; ++mat_index){
		  int mat_allele = (mat_index == 0 ? mat_i : mat_j);
		  for (int pat_index = 0; pat_index < 2; ++pat_index){
		    int pat_allele = (pat_index == 0 ? pat_i : pat_j);

		    double no_mutation_config_ll = config_ll + unphased_gls.get_gl(child_gl_index, std::min(mat_allele, pat_allele), std::max(mat_allele, pat_allele));
		    update_streaming_log_sum_exp(no_mutation_config_ll, ll_no_mutation_max, ll_no_mutation_total);

		    // All putative mutations to the maternal allele
		    double max_ll_mat_mut = config_ll + unphased_gls.get_max_gl_allele_fixed(child_gl_index, pat_allele) + mut_model.max_log_prior_mutation(mat_allele);
		    if (max_ll_mat_mut > std::min(ll_one_denovo_max, ll_one_other_max)-MIN_CONTRIBUTION){
		      for (int mut_allele = 0; mut_allele < num_alleles; mut_allele++){
			if (mut_allele == mat_allele)
			  continue;
			double prob = config_ll + unphased_gls.get_gl(child_gl_index, std::min(mut_allele, pat_allele), std::max(mut_allele, pat_allele))
			  + mut_model.log_prior_mutation(mat_allele, mut_allele);
			if (mut_allele != mat_i && mut_allele != mat_j && mut_allele != pat_i && mut_allele != pat_j)
			  update_streaming_log_sum_exp(prob, ll_one_denovo_max, ll_one_denovo_total);
			else
			  update_streaming_log_sum_exp(prob, ll_one_other_max, ll_one_other_total);
		      }
		    }

		    // All putative mutations to the paternal allele
		    double max_ll_pat_mut = config_ll + unphased_gls.get_max_gl_allele_fixed(child_gl_index, mat_allele) + mut_model.max_log_prior_mutation(pat_allele);
		    if (max_ll_pat_mut > std::min(ll_one_denovo_max, ll_one_other_max)-MIN_CONTRIBUTION){
		      for (int mut_allele = 0; mut_allele < num_alleles; mut_allele++){
			if (mut_allele == pat_allele)
			  continue;
			double prob = config_ll + unphased_gls.get_gl(child_gl_index, std::min(mat_allele, mut_allele), std::max(mat_allele, mut_allele))
			  + mut_model.log_prior_mutation(pat_allele, mut_allele);
			if (mut_allele != mat_i && mut_allele != mat_j && mut_allele != pat_i && mut_allele != pat_j)
			  update_streaming_log_sum_exp(prob, ll_one_denovo_max, ll_one_denovo_total);
			else
			  update_streaming_log_sum_exp(prob, ll_one_other_max, ll_one_other_total);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}

	// Compute total LL for each scenario and add it to the VCF
	double total_ll_no_mutation = finish_streaming_log_sum_exp(ll_no_mutation_max, ll_no_mutation_total);
	double total_ll_one_denovo  = finish_streaming_log_sum_exp(ll_one_denovo_max,  ll_one_denovo_total);
	double total_ll_one_other   = finish_streaming_log_sum_exp(ll_one_other_max,   ll_one_other_total);
      }
    }
  }
}
