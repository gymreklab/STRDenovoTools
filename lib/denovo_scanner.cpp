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

#include "denovo_allele_priors.h"
#include "denovo_scanner.h"
#include "locus_inspector.h"
#include "mathops.h"
#include "mutation_model.h"
#include "vcf_input.h"

#include "htslib/kfunc.h"

std::string TrioDenovoScanner::START_KEY   = "START";
std::string TrioDenovoScanner::END_KEY     = "END";
std::string TrioDenovoScanner::PERIOD_KEY  = "PERIOD";
std::string DenovoResult::PERIOD_KEY  = "PERIOD";

TrioDenovoScanner::~TrioDenovoScanner() {
  locus_summary_.close();
  all_mutations_file_.close();
}

bool TrioDenovoScanner::GetFollowsMI(const int& gt_mother_a, const int& gt_mother_b,
				     const int& gt_father_a, const int& gt_father_b,
				     const int& gt_child_a, const int& gt_child_b,
				     bool is_chrx, const int& child_sex) {
  bool follows_MI = false;
  if (is_chrx && child_sex == SEX_MALE) {
    if (gt_child_a != gt_child_b) {
      return follows_MI; // ThIs shouldn't really happen, should be homozygous
    }
    if ((gt_child_a == gt_mother_a || gt_child_a == gt_mother_b)) {
      follows_MI = true;
    }
  } else {
    if ((gt_child_a == gt_mother_a || gt_child_a == gt_mother_b) && (gt_child_b == gt_father_a || gt_child_b == gt_father_b)) {
      follows_MI = true;
    } else if ((gt_child_a == gt_father_a || gt_child_a == gt_father_b) && (gt_child_b == gt_mother_a || gt_child_b == gt_mother_b)) {
      follows_MI = true;
    }
  }
  return follows_MI;
}

void TrioDenovoScanner::naive_scan(VCF::VCFReader& strvcf) {
  if (options_.debug) PrintMessageDieOnError("Scanning for de novos using naive method...", M_PROGRESS);
  VCF::Variant str_variant;
  std::vector<NuclearFamily> families = pedigree_set_.get_families();
  std::vector<DenovoResult> denovo_results;
  while (strvcf.get_next_variant(&str_variant, options_.round_alleles)) {
    // Clear leftovers
    denovo_results.clear();
    // Get locus info
    if (options_.debug) PrintMessageDieOnError("Getting locus info...", M_PROGRESS);
    int32_t start;
    if (options_.gangstr) {
      start = str_variant.get_position();
    } else {
      str_variant.get_INFO_value_single_int(START_KEY, start);
    }
    int32_t end; str_variant.get_INFO_value_single_int(END_KEY, end);
    int32_t period; str_variant.get_INFO_value_single_int(PERIOD_KEY, period);
    if (options_.period != 0 && period != options_.period) {
      str_variant.destroy_record();
      continue;
    }
    std::stringstream ss;
    ss << "Processing STR region " << str_variant.get_chromosome() << ":" << start << "-" << end;
    PrintMessageDieOnError(ss.str(), M_PROGRESS);
    if (options_.debug)  PrintMessageDieOnError("Checking if we have any samples...", M_PROGRESS);
    if (str_variant.num_samples() == str_variant.num_missing()) {
      std::stringstream ss;
      ss << "Found " << str_variant.num_samples() << " missing: " << str_variant.num_missing() << endl;
      if (options_.debug) PrintMessageDieOnError(ss.str(), M_PROGRESS);
      str_variant.destroy_record();
      continue;
    }
    // Scan each family
    if (options_.debug) PrintMessageDieOnError("Scanning families...", M_PROGRESS);
    for (auto family_iter = families.begin(); family_iter != families.end(); family_iter++) {
      if (!options_.family.empty() and family_iter->get_family_id() != options_.family) {
	continue;
      }
      if (str_variant.sample_call_missing(family_iter->get_mother()) ||
	  str_variant.sample_call_missing(family_iter->get_father())) {
        continue; // No point if there are no calls for parents
      }
      // First check we have children we need
      bool scan_for_denovo = true;
      int num_children_nonmissing = 0;
      for (auto child_iter = family_iter->get_children().begin();
	   child_iter != family_iter->get_children().end(); child_iter++) {
	if (options_.require_all_children) {
	  scan_for_denovo = scan_for_denovo && !str_variant.sample_call_missing(*child_iter);
	}
	if (!str_variant.sample_call_missing(*child_iter)) {
	  num_children_nonmissing += 1;
	}
      }
      ss.str("");
      ss << "Found " << num_children_nonmissing << " nonmissing children . Scanning for de novo=" << scan_for_denovo;
      if (options_.debug) PrintMessageDieOnError(ss.str(), M_PROGRESS);
      scan_for_denovo = scan_for_denovo && (num_children_nonmissing > 0);

      // Get parent gentoype info
      int gt_mother_a, gt_mother_b, gt_father_a, gt_father_b, gt_child_a, gt_child_b;
      str_variant.get_genotype(family_iter->get_mother(), gt_mother_a, gt_mother_b);
      str_variant.get_genotype(family_iter->get_father(), gt_father_a, gt_father_b);

      // Scan each child in the family
      for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); child_iter++) {
	if (options_.debug) PrintMessageDieOnError("Scanning child...", M_PROGRESS);
	if (!scan_for_denovo || str_variant.sample_call_missing(*child_iter)) {
	  continue;
	}
	// Get child genotype info
	str_variant.get_genotype(*child_iter, gt_child_a, gt_child_b);

	// Test for Mendelian inheritance
	bool follows_MI = GetFollowsMI(gt_mother_a, gt_mother_b,
				       gt_father_a, gt_father_b,
				       gt_child_a, gt_child_b,
				       options_.chrX,
				       family_iter->get_child_sex(*child_iter));

	// Add to results
	if (follows_MI) { // Set dummy values to force posterior to 0
	  DenovoResult dnr(family_iter->get_family_id(), family_iter->get_mother(), family_iter->get_father(), (*child_iter),
			   family_iter->get_child_phenotype(*child_iter),
         family_iter->get_child_sex(*child_iter),
			   strvcf.get_sample_index(family_iter->get_mother()),
			   strvcf.get_sample_index(family_iter->get_father()),
			   strvcf.get_sample_index(*child_iter),
			   0, 0, -8);
	  denovo_results.push_back(dnr);
	} else { // Set dummy values to force posterior to 1
	  DenovoResult dnr(family_iter->get_family_id(), family_iter->get_mother(), family_iter->get_father(), (*child_iter),
			   family_iter->get_child_phenotype(*child_iter),
         family_iter->get_child_sex(*child_iter),
			   strvcf.get_sample_index(family_iter->get_mother()),
			   strvcf.get_sample_index(family_iter->get_father()),
			   strvcf.get_sample_index(*child_iter),
			   -10000, 0, -8);
	  denovo_results.push_back(dnr);
	}
      }
    }
    summarize_results(denovo_results, str_variant);
  }
}

void TrioDenovoScanner::scan(VCF::VCFReader& strvcf,
			     MutationPriors& priors) {
  if (options_.debug) PrintMessageDieOnError("Scanning for de novos using model based method...", M_PROGRESS);
  VCF::Variant str_variant;
  std::vector<NuclearFamily> families = pedigree_set_.get_families();
  std::vector<DenovoResult> denovo_results;
  while (strvcf.get_next_variant(&str_variant, options_.round_alleles)) {
    // Clear leftovers
    denovo_results.clear();
    // Get locus info
    bool dummy_models = false; // Use if only 1 allele but still including
    int num_alleles = str_variant.num_alleles();
    if (options_.combine_alleles) {
      num_alleles = str_variant.num_alleles_by_length(options_.round_alleles);
    }
    if (options_.gangstr) {
      num_alleles = str_variant.num_gangstr_alleles();
    }
    if (num_alleles <= 1) {
      if (options_.include_invariant && num_alleles == 1) {
	dummy_models = true;
      } else {
	str_variant.destroy_record();
	continue;
      }
    }
    if (options_.debug) PrintMessageDieOnError("Getting locus info...", M_PROGRESS);
    int32_t start;
    if (options_.gangstr) {
      start = str_variant.get_position();
    } else {
      str_variant.get_INFO_value_single_int(START_KEY, start);
    }
    int32_t end; str_variant.get_INFO_value_single_int(END_KEY, end);
    int32_t period; str_variant.get_INFO_value_single_int(PERIOD_KEY, period);
    if (options_.period != 0 && period != options_.period) {
      str_variant.destroy_record();
      continue;
    }
    std::stringstream ss;
    ss << "Processing STR region " << str_variant.get_chromosome() << ":" << start << "-" << end
       << " with " << num_alleles << " alleles.";
    if (options_.combine_alleles) {
      ss << " Total allele seqs: " << str_variant.num_alleles();
    }
    PrintMessageDieOnError(ss.str(), M_PROGRESS);
    if (options_.debug)  PrintMessageDieOnError("Checking if we have any samples...", M_PROGRESS);
    if (str_variant.num_samples() == str_variant.num_missing()) {
      std::stringstream ss;
      ss << "Found " << str_variant.num_samples() << " missing: " << str_variant.num_missing() << endl;
      if (options_.debug) PrintMessageDieOnError(ss.str(), M_PROGRESS);
      str_variant.destroy_record();
      continue;
    }
    if (options_.debug) PrintMessageDieOnError("Checking if we have too many alleles...", M_PROGRESS);
    if (num_alleles > options_.max_num_alleles) {
      str_variant.destroy_record();
      continue;
    }

    // Set up
    if (options_.debug) PrintMessageDieOnError("Set up GLs...", M_PROGRESS);
    GL* unphased_gls;
    UnphasedLengthGL unphased_length_gls(str_variant, options_, dummy_models);
    UnphasedGL unphased_seq_gls(str_variant, options_, dummy_models);
    GangSTRGL gangstr_gls(str_variant, options_, dummy_models);
    MutationModel mut_model(str_variant, priors, options_, dummy_models);
    DiploidGenotypePrior* dip_gt_priors;
    if (options_.combine_alleles) {
      if (options_.debug) PrintMessageDieOnError("Set up length GLs...", M_PROGRESS);
      unphased_gls = &unphased_length_gls;
      if (options_.use_pop_priors) {
	dip_gt_priors = new PopulationGenotypeLengthPrior(str_variant, families, options_.round_alleles);
      } else {
	dip_gt_priors = new UniformGenotypeLengthPrior(str_variant, families, options_.round_alleles);
      }
      if (options_.debug) PrintMessageDieOnError("Done setting up priors...", M_PROGRESS);
    } else if (options_.gangstr) {
      unphased_gls = &gangstr_gls;
      dip_gt_priors = new UniformGenotypePrior(str_variant, families, true); // TODO support priors for gangstr
    } else {
      unphased_gls = &unphased_seq_gls;
      if (options_.use_pop_priors) {
	dip_gt_priors = new PopulationGenotypePrior(str_variant, families);
      } else {
	dip_gt_priors = new UniformGenotypePrior(str_variant, families, false);
      }
    }
    double prior_rate = priors.GetPrior(str_variant.get_chromosome(), start);
    const double LOG_ONE_FOURTH = -log10(4);
    const double LOG_TWO        = log10(2);

    // Scan each family
    if (options_.debug) PrintMessageDieOnError("Scanning families...", M_PROGRESS);
    for (auto family_iter = families.begin(); family_iter != families.end(); family_iter++) {
      if (!options_.family.empty() and family_iter->get_family_id() != options_.family) {
	continue;
      }
      bool scan_for_denovo = true;
      if (dummy_models) {
	scan_for_denovo = !str_variant.sample_call_missing(family_iter->get_mother()) && !str_variant.sample_call_missing(family_iter->get_father());
      } else {
	scan_for_denovo = unphased_gls->has_sample(family_iter->get_mother()) && unphased_gls->has_sample(family_iter->get_father());
      }
      if (options_.debug) PrintMessageDieOnError("Check children...", M_PROGRESS);
      int num_children_nonmissing = 0;

      ss.str("");
      ss << "Before checking children, Scan for de novo=" << scan_for_denovo << " "
	 << "Mother has data=" << unphased_gls->has_sample(family_iter->get_mother()) << " " << family_iter->get_mother() << " "
	 << "Father has data=" << unphased_gls->has_sample(family_iter->get_father()) << " " << family_iter->get_father();
      if (options_.debug) PrintMessageDieOnError(ss.str(), M_PROGRESS);
      for (auto child_iter = family_iter->get_children().begin();
	   child_iter != family_iter->get_children().end(); child_iter++) {
	if (dummy_models) {
	  if (options_.require_all_children) {
	    scan_for_denovo = scan_for_denovo && !str_variant.sample_call_missing(*child_iter);
	  }
	  if (!str_variant.sample_call_missing(*child_iter)) {
	    num_children_nonmissing += 1;
	  }
	} else {
	  if (options_.require_all_children) {
	    scan_for_denovo = scan_for_denovo && unphased_gls->has_sample(*child_iter);
	  }
	  if (unphased_gls->has_sample(*child_iter)) {
	    num_children_nonmissing += 1;
	  }
	}
      }
      ss.str("");
      ss << "Found " << num_children_nonmissing << " nonmissing children . Scanning for de novo=" << scan_for_denovo;
      if (options_.debug) PrintMessageDieOnError(ss.str(), M_PROGRESS);
      scan_for_denovo = scan_for_denovo && (num_children_nonmissing > 0);

      // Scan each child in the family
      for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); child_iter++) {
	double total_ll_no_mutation, total_ll_one_denovo;
	if (options_.debug) PrintMessageDieOnError("Scanning child...", M_PROGRESS);
	if (dummy_models) {
	  if (options_.debug) PrintMessageDieOnError("Inside dummy...", M_PROGRESS);
	  if (!scan_for_denovo) {
	    continue;
	  }
	  total_ll_no_mutation = 0;
	  total_ll_one_denovo = -DBL_MAX/2;
	} else {
	  if (options_.debug) PrintMessageDieOnError("Checking child...", M_PROGRESS);
	  if (!scan_for_denovo || !unphased_gls->has_sample(*child_iter)) {
	    continue;
	  }

	  // To accelerate computations, we will ignore configurations that make a neglible contribution (< 0.01%) to the total LL
      	  // For mutational scenarios, we aggregate 1/4*A^2*(A+1)^2*4*2*A values. Therefore, to ignore a configuration with LL=X:
      	  // X*A^3*(A+1)^2*2 < TOTAL/10000;
      	  // logX < log(TOTAL) - log(10000*A^3*(A+1)^2*2) = log(TOTAL) - [log(10000) + 3log(A) + 2log(A+1) + log(2)];
      	  double MIN_CONTRIBUTION   = 4 + 3*log10(num_alleles) + 2*log(num_alleles+1) + LOG_TWO;
      	  double ll_no_mutation_max = -DBL_MAX/2, ll_no_mutation_total = 0.0;
      	  double ll_one_denovo_max  = -DBL_MAX/2, ll_one_denovo_total  = 0.0;
      	  if (options_.debug) PrintMessageDieOnError("Get GL indices...", M_PROGRESS);
      	  int mother_gl_index       = unphased_gls->get_sample_index(family_iter->get_mother());
      	  int father_gl_index       = unphased_gls->get_sample_index(family_iter->get_father());
      	  int child_gl_index        = unphased_gls->get_sample_index(*child_iter);
      	  // Iterate over all maternal genotypes
      	  if (options_.debug) PrintMessageDieOnError("Compute likelihoods...", M_PROGRESS);
      	  for (int mat_i = 0; mat_i < num_alleles; mat_i++){
      	    for (int mat_j = 0; mat_j <= mat_i; mat_j++){
      	      double mat_ll = dip_gt_priors->log_unphased_genotype_prior(mat_j, mat_i, family_iter->get_mother()) + unphased_gls->get_gl(mother_gl_index, mat_j, mat_i);

	      // Iterate over all paternal genotypes
	      for (int pat_i = 0; pat_i < num_alleles; pat_i++){
		for (int pat_j = 0; pat_j <= pat_i; pat_j++){

		  if (options_.chrX && (pat_i != pat_j)) { continue;} // ignore case of father heterozygous for chrX

		  double pat_ll    = dip_gt_priors->log_unphased_genotype_prior(pat_j, pat_i, family_iter->get_father()) + unphased_gls->get_gl(father_gl_index, pat_j, pat_i);
		  double config_ll = mat_ll + pat_ll + LOG_ONE_FOURTH;

		  // Compute total LL for each scenario
		  // Iterate over all 4 possible inheritance patterns for the child
		  for (int mat_index = 0; mat_index < 2; ++mat_index) {
		    int mat_allele = (mat_index == 0 ? mat_i : mat_j);
		    for (int pat_index = 0; pat_index < 2; ++pat_index) {
		      int pat_allele = (pat_index == 0 ? pat_i : pat_j);

		      double no_mutation_config_ll;
		      if (options_.chrX && family_iter->get_child_sex(*child_iter)==SEX_MALE) {
			no_mutation_config_ll = config_ll + unphased_gls->get_gl(child_gl_index, mat_allele, mat_allele); // homozygous for mom allele
		      } else {
			no_mutation_config_ll = config_ll + unphased_gls->get_gl(child_gl_index, std::min(mat_allele, pat_allele), std::max(mat_allele, pat_allele));
		      }
		      update_streaming_log_sum_exp(no_mutation_config_ll, ll_no_mutation_max, ll_no_mutation_total);

		      // All putative mutations to the maternal allele
		      double max_ll_mat_mut = config_ll + unphased_gls->get_max_gl_allele_fixed(child_gl_index, pat_allele) +
			mut_model.max_log_prior_mutation(str_variant.GetSizeFromLengthAllele(mat_allele));
		      if (max_ll_mat_mut > ll_one_denovo_max-MIN_CONTRIBUTION) {
			for (int mut_allele = 0; mut_allele < num_alleles; mut_allele++) {
			  if (mut_allele == mat_allele)
			    continue;
			  double prob;
			  if (!(options_.chrX && family_iter->get_child_sex(*child_iter)==SEX_MALE)) {
			    prob = config_ll + unphased_gls->get_gl(child_gl_index, std::min(mut_allele, pat_allele),
									   std::max(mut_allele, pat_allele))
			      + mut_model.log_prior_mutation(str_variant.GetSizeFromLengthAllele(mat_allele),
							     str_variant.GetSizeFromLengthAllele(mut_allele));
			  } else {
			    prob = config_ll + unphased_gls->get_gl(child_gl_index, mut_allele, mut_allele)
			      + mut_model.log_prior_mutation(str_variant.GetSizeFromLengthAllele(mat_allele),
							     str_variant.GetSizeFromLengthAllele(mut_allele));
			  }
			  update_streaming_log_sum_exp(prob, ll_one_denovo_max, ll_one_denovo_total);
			}
		      }

		      // All putative mutations to the paternal allele - skip if male chrX
		      if (!(options_.chrX && family_iter->get_child_sex(*child_iter)==SEX_MALE)) {
			double max_ll_pat_mut = config_ll + unphased_gls->get_max_gl_allele_fixed(child_gl_index, mat_allele) +
			  mut_model.max_log_prior_mutation(str_variant.GetSizeFromLengthAllele(pat_allele));
			if (max_ll_pat_mut > ll_one_denovo_max - MIN_CONTRIBUTION){
			  for (int mut_allele = 0; mut_allele < num_alleles; mut_allele++){
			    if (mut_allele == pat_allele)
			      continue;
			    double prob = config_ll + unphased_gls->get_gl(child_gl_index, std::min(mat_allele, mut_allele), std::max(mat_allele, mut_allele))
			      + mut_model.log_prior_mutation(str_variant.GetSizeFromLengthAllele(pat_allele),
							     str_variant.GetSizeFromLengthAllele(mut_allele));
			    update_streaming_log_sum_exp(prob, ll_one_denovo_max, ll_one_denovo_total);
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  total_ll_no_mutation = finish_streaming_log_sum_exp(ll_no_mutation_max, ll_no_mutation_total);
	  total_ll_one_denovo  = finish_streaming_log_sum_exp(ll_one_denovo_max,  ll_one_denovo_total);
	}
	// Add to results
	DenovoResult dnr(family_iter->get_family_id(),
       family_iter->get_mother(), family_iter->get_father(), (*child_iter),
			 family_iter->get_child_phenotype(*child_iter),
       family_iter->get_child_sex(*child_iter),
			 strvcf.get_sample_index(family_iter->get_mother()),
			 strvcf.get_sample_index(family_iter->get_father()),
			 strvcf.get_sample_index(*child_iter),
			 total_ll_no_mutation, total_ll_one_denovo, prior_rate);
	if (options_.debug) PrintMessageDieOnError("Adding dnr result", M_PROGRESS);
	denovo_results.push_back(dnr);
      }
    }
    summarize_results(denovo_results, str_variant);
    if (options_.debug) PrintMessageDieOnError("Done summarizing...", M_PROGRESS);
    delete dip_gt_priors;
    denovo_results.clear();

    // Destroy vcf_record_ from variant
    str_variant.destroy_record();
    if (options_.debug) PrintMessageDieOnError("Destroy...", M_PROGRESS);
  }
}

void DenovoResult::GetRepcn(const VCF::Variant& variant, const int32_t& sample_ind,
			    int* repcn_a, int* repcn_b) {
  if (variant.sample_call_missing(sample_ind)) {
    // Return -1, -1 if the variant not called. This shouldn't happen
    // since we check for no call earlier
    PrintMessageDieOnError("Encountered unexpected missing genotype in GetRepcn", M_WARNING);
    *repcn_a = -1;
    *repcn_b = -1;
    return;
  }
  int32_t period; variant.get_INFO_value_single_int(PERIOD_KEY, period);
  int gt_a, gt_b;
  variant.get_genotype(sample_ind, gt_a, gt_b);
  *repcn_a = ((int)variant.get_allele(gt_a).size())/period;
  *repcn_b = ((int)variant.get_allele(gt_b).size())/period;
}

int DenovoResult::GetMaxFlankAllele(const std::string flnkstring) {
  if (flnkstring == "NULL") {
    return 0;
  }
  int maxallele = 0;
  std::vector<std::string> flnkreads_items, items;
  split_by_delim(flnkstring, '|', flnkreads_items);
  for (auto item_iter = flnkreads_items.begin(); item_iter != flnkreads_items.end(); item_iter++) {
    items.clear();
    split_by_delim(*item_iter, ',' ,items);
    if (items.size() != 2) {
      PrintMessageDieOnError("Invalid flankstring encountered " + flnkstring, M_ERROR);
    }
    int a = stoi(items[0]);
    int count = stoi(items[1]);
    if (a > maxallele) {
      maxallele = a;
    }
  }
  return maxallele;
}

int DenovoResult::GetFlankLargerThan(const std::string flnkstring, const int& allele) {
  if (flnkstring == "NULL") {
    return 0;
  }
  int totalcount = 0;
  std::vector<std::string> flnkreads_items, items;
  split_by_delim(flnkstring, '|', flnkreads_items);
  for (auto item_iter = flnkreads_items.begin(); item_iter != flnkreads_items.end(); item_iter++) {
    items.clear();
    split_by_delim(*item_iter, ',' ,items);
    if (items.size() != 2) {
      PrintMessageDieOnError("Invalid flankstring encountered " + flnkstring, M_ERROR);
    }
    int a = stoi(items[0]);
    int count = stoi(items[1]);
    if (a > allele) {
      totalcount += count;
    }
  }
  return totalcount;
}

int DenovoResult::GetTrans(const int& child, const int& p1, const int& p2) {
  if (child == p1) {
    if (p1 == p2) {
      return 0;
    } else if (p1<p2) {
      return 1;
    } else {
      return 2;
    }
  }
  if (child == p2) {
    if (p1 == p2) {
      return 0;
    } else if (p2<p1) {
      return 1;
    } else {
      return 2;
    }
  }
  return 2;
}

// Must be run after GetMutationInfo
// 0=unknown, 1=shorter, 2=longer
void DenovoResult::TestTransmission(int* long_mother, int* long_father) {
  if ((child_gt_a_ == mat_gt_a_ || child_gt_a_ == mat_gt_b_) &&
      (child_gt_b_ == pat_gt_a_ || child_gt_b_ == pat_gt_b_)) {
    // a came from mother, b from father
    *long_mother = GetTrans(child_gt_a_, mat_gt_a_, mat_gt_b_);
    *long_father = GetTrans(child_gt_b_, pat_gt_a_, pat_gt_b_);
    return;
  }
  if ((child_gt_a_ == pat_gt_a_ || child_gt_a_ == pat_gt_b_) &&
      (child_gt_b_ == mat_gt_a_ || child_gt_b_ == mat_gt_b_)) {
    // b came from mother, a from father
    *long_mother = GetTrans(child_gt_b_, mat_gt_a_, mat_gt_b_);
    *long_father = GetTrans(child_gt_a_, pat_gt_a_, pat_gt_b_);
    return;
  }
  // Else, we aren't sure
  *long_mother = 0;
  *long_father = 0;
}

int DenovoResult::GetMutSize(const int& new_allele, const int& a1, const int& a2) {
  int diff1 = new_allele - a1;
  int diff2 = new_allele - a2;
  if (abs(diff1) < abs(diff2)) {
    return diff1;
  } else {
    return diff2;
  }
}

int DenovoResult::GetMutSize(const int& new_allele, const int& a1, const int& a2, const int& a3, const int& a4) {
  int diff1 = new_allele - a1;
  int diff2 = new_allele - a2;
  int diff3 = new_allele - a3;
  int diff4 = new_allele - a4;
  int min_size = std::min(std::min(abs(diff1), abs(diff2)), std::min(abs(diff3), abs(diff4)));
  if (min_size == abs(diff1)) {
    return diff1;
  } else if (min_size == abs(diff2)) {
    return diff2;
  } else if (min_size == abs(diff2)) {
    return diff3;
  } else {
    return diff4;
  }
}

void DenovoResult::GetFRR(const std::string& rcstring, int* frrcount) {
  std::vector<std::string> items;
  split_by_delim(rcstring, ',', items);
  if (items.size() != 4) {
    PrintMessageDieOnError("Invalid RC string encountered " + rcstring, M_ERROR);
  }
  *frrcount = stoi(items[2]);
}

void DenovoResult::GetEnclosing(const std::string& enclstring, int& new_allele,
				const int32_t& repcn_a, const int32_t& repcn_b,
				int* encl_newallele, int* encl_total, int* encl_match) {
  *encl_newallele = 0;
  *encl_total = 0;
  *encl_match = 0;
  if (enclstring == "NULL") {
    return;
  }
  std::vector<std::string> enclreads_items, items;
  split_by_delim(enclstring, '|', enclreads_items);
  for (auto item_iter = enclreads_items.begin(); item_iter != enclreads_items.end(); item_iter++) {
    split_by_delim(*item_iter, ',', items);
    if (items.size() != 2) {
      PrintMessageDieOnError("Invalid enclreads " + enclstring, M_ERROR);
    }
    int allele = stoi(items[0]);
    int count = stoi(items[1]);
    *encl_total += count;
    if (allele == new_allele) {
      *encl_newallele = count;
    }
    if (allele == repcn_a || allele == repcn_b) {
      *encl_match += count;
    }
    items.clear();
  }
}


/*
  Infer parent of origin and mutation size info

  Sets:
   new_allele_
   mut_size_
   poocase_
   new_allele_in_parents_

  Returns false if Mendelian, true otherwise

      poocase: describes inheritance pattern
    0: unknown
    1: Mendelian (no denovo)
    2: New allele from father
    3: New allele from mother
    4: Unclear

*/
bool DenovoResult::GetPOOMutationInfo(const bool& chrX) {
  if (! vcfinfo_set_) {
    PrintMessageDieOnError("Can't run GetPOOMutationInfo before setting VCF info", M_ERROR);
  }

  new_allele_ = 0;
  mut_size_ = 0;
  poocase_ = 0; // 1=Mendelian 2=new allele from father, 3=new allele from mother, 4=not in anyone
  new_allele_in_parents_ = false;

  if (chrX && child_sex_ == SEX_MALE) {
    // Case 1: Mendelian
    if ((child_gt_a_ == mat_gt_a_ || child_gt_a_  == mat_gt_b_)) {
      poocase_ = 1;
      new_allele_in_parents_ = true;
      new_allele_ = 0;
      mut_size_ = 0;
      return false;
    } else { // Case 3: New allele (only check a) from mother
      new_allele_ = child_gt_a_;
      poocase_ = 3;
      mut_size_ = GetMutSize(new_allele_, mat_gt_a_, mat_gt_b_);
    }
    return true;
  }
  else {
    // Case 1: Mendelian
    if ((child_gt_a_ == mat_gt_a_ || child_gt_a_ == mat_gt_b_) &&
        (child_gt_b_ == pat_gt_a_ || child_gt_b_ == pat_gt_b_)) {
      poocase_ = 1;
      new_allele_in_parents_ = true;
      new_allele_ = 0;
      mut_size_ = 0;
      return false;
    }
    if ((child_gt_a_ == pat_gt_a_  || child_gt_a_ == pat_gt_b_ ) &&
        (child_gt_b_ == mat_gt_a_  || child_gt_b_ == mat_gt_b_)) {
      poocase_ = 1;
      new_allele_in_parents_ = true;
      new_allele_ = 0;
      mut_size_ = 0;
      return false;
    }
    // Case 2: new allele from father.
    // Allele a in mother only, allele b not in father
    if ((child_gt_a_ == mat_gt_a_ || child_gt_a_ == mat_gt_b_) &&
        (child_gt_a_ != pat_gt_a_ && child_gt_a_ != pat_gt_b_) &&
        (child_gt_b_ != pat_gt_a_  && child_gt_b_ != pat_gt_b_)) {
      new_allele_ = child_gt_b_;
      poocase_ = 2;
      mut_size_ = GetMutSize(new_allele_, pat_gt_a_, pat_gt_b_);
    }
    // Allele b in mother only, allele a not in father
    if ((child_gt_b_  == mat_gt_a_ || child_gt_b_  == mat_gt_b_) &&
        (child_gt_b_  != pat_gt_a_ && child_gt_b_  != pat_gt_b_) &&
        (child_gt_a_ != pat_gt_a_ && child_gt_a_ != pat_gt_b_)) {
      new_allele_ = child_gt_a_;
      poocase_ = 2;
      mut_size_ = GetMutSize(new_allele_, pat_gt_a_, pat_gt_b_);
    }
    // Case 3: New allele from mother
    // Allele a in father only, allele b not in mother
    if ((child_gt_a_ == pat_gt_a_ || child_gt_a_ == pat_gt_b_) &&
        (child_gt_a_ != mat_gt_a_ && child_gt_a_ != mat_gt_b_) &&
        (child_gt_b_ != mat_gt_a_ && child_gt_b_ != mat_gt_b_)) {
      new_allele_ = child_gt_b_;
      poocase_ = 3;
      mut_size_ = GetMutSize(new_allele_, mat_gt_a_, mat_gt_b_);
    }
    // Allele b in father only, allele a not in mother
    if ((child_gt_b_ == pat_gt_a_ || child_gt_b_ == pat_gt_b_) &&
        (child_gt_b_ != mat_gt_a_ && child_gt_b_ != mat_gt_b_) &&
        (child_gt_a_ != mat_gt_a_ && child_gt_a_ != mat_gt_b_)) {
      new_allele_ = child_gt_a_;
      poocase_ = 3;
      mut_size_ = GetMutSize(new_allele_, mat_gt_a_, mat_gt_b_);
    }
    // Case 4: new allele not in either parent at all and we haven't figured it out yet
    if (poocase_ == 0 && child_gt_a_ != mat_gt_a_ && child_gt_a_ != mat_gt_b_ &&
        child_gt_a_ != pat_gt_a_ && child_gt_a_ != pat_gt_b_) {
      new_allele_ = child_gt_a_;
      poocase_ = 4;
      mut_size_ = GetMutSize(new_allele_, mat_gt_a_, mat_gt_b_, pat_gt_a_, pat_gt_b_);
    }
    if (poocase_ == 0 && child_gt_b_ != mat_gt_a_ && child_gt_b_ != mat_gt_b_ &&
        child_gt_b_ != pat_gt_a_ && child_gt_b_ != pat_gt_b_) {
      new_allele_ = child_gt_b_;
      poocase_ = 4;
      mut_size_ = GetMutSize(new_allele_, mat_gt_a_, mat_gt_b_, pat_gt_a_, pat_gt_b_);
    }

    if (new_allele_ == mat_gt_a_ || new_allele_ == mat_gt_b_ ||
        new_allele_ == pat_gt_a_ || new_allele_ == pat_gt_b_) {
      new_allele_in_parents_ = true;
    }
    return true;
  }
}

/*
  Set info extracted from the VCF file
*/
void DenovoResult::SetVCFInfo(const VCF::Variant& variant) {
  // *** First set genotypes *** //
  int repcn_mother_a, repcn_mother_b, repcn_father_a, repcn_father_b, repcn_child_a, repcn_child_b;
  GetRepcn(variant, mother_ind_, &repcn_mother_a, &repcn_mother_b);
  GetRepcn(variant, father_ind_, &repcn_father_a, &repcn_father_b);
  GetRepcn(variant, child_ind_, &repcn_child_a, &repcn_child_b);
  child_gt_ = std::to_string(repcn_child_a) + "," + std::to_string(repcn_child_b);
  mat_gt_ = std::to_string(repcn_mother_a) + "," + std::to_string(repcn_mother_b);
  pat_gt_ = std::to_string(repcn_father_a) + "," + std::to_string(repcn_father_b);
  child_gt_a_ = repcn_child_a; child_gt_b_ = repcn_child_b;
  mat_gt_a_ = repcn_mother_a; mat_gt_b_ = repcn_mother_b;
  pat_gt_a_ = repcn_father_a; pat_gt_b_ = repcn_father_b;

  vcfinfo_set_ = true;
}

/*
  Check read filters

  Return true if passing, else false
*/
bool DenovoResult::CheckReadFilters(const Options& options, const VCF::Variant& variant) {
  if (!vcfinfo_set_) {
    SetVCFInfo(variant);
  }
  // *** Get enclosing read info to set filter *** //
  std::vector<std::string> enclreads;
  variant.get_FORMAT_value_single_string(ENCLREADS_KEY, enclreads);
  int encl_child, total_encl_child, child_encl_match;
  int encl_mother, total_encl_mother, mother_encl_match;
  int encl_father, total_encl_father, father_encl_match;
  float encl_reads_perc_parent = 0;;

  GetEnclosing(enclreads[child_ind_], new_allele_, child_gt_a_, child_gt_b_, &encl_child, &total_encl_child, &child_encl_match);
  GetEnclosing(enclreads[mother_ind_], new_allele_, mat_gt_a_, mat_gt_b_, &encl_mother, &total_encl_mother, &mother_encl_match);
  GetEnclosing(enclreads[father_ind_], new_allele_, pat_gt_a_, pat_gt_b_, &encl_father, &total_encl_father, &father_encl_match);

  // Set in dnr_iter
  encl_reads_child_ = encl_child;
  encl_reads_mother_ = encl_mother;
  encl_reads_father_ = encl_father;

  if (poocase_ == 2) {
    encl_reads_parent_ = encl_reads_father_;
    encl_reads_perc_parent = (float)encl_father/(float)total_encl_father;
  } else if (poocase_ == 3) {
    encl_reads_parent_ = encl_reads_mother_;
    encl_reads_perc_parent = (float)encl_mother/(float)total_encl_mother;
  } else {
    if (encl_reads_father_ > encl_reads_mother_) {
      encl_reads_parent_ = encl_reads_father_;
      encl_reads_perc_parent = (float)encl_father/(float)total_encl_father;
    } else {
      encl_reads_parent_ = encl_reads_mother_;
      encl_reads_perc_parent = (float)encl_mother/(float)total_encl_mother;
    }
  }

  // Total enclosing
  if ((total_encl_child < options.min_total_encl) ||
      (total_encl_mother < options.min_total_encl) ||
      (total_encl_father < options.min_total_encl)) {
    if (options.debug) {cerr << "reject based on total enclosing reads" << endl;}
    return false;
  }
  // Messiness
  if (((float)child_encl_match/(float)total_encl_child < options.min_encl_match) ||
      ((float)mother_encl_match/(float)total_encl_mother < options.min_encl_match) ||
      ((float)father_encl_match/(float)total_encl_father < options.min_encl_match)) {
    if (options.debug) {cerr << "reject based on messy enclosing reads" << endl;}
    return false;
  }
  // Enclosing matching new allele
  if (encl_child < options.min_num_encl_child) {
    if (options.debug) {cerr << "reject based on child encl: " << encl_child << endl;}
    return false;
  }
  if (encl_reads_parent_ > options.max_num_encl_parent) {
    if (options.debug) {cerr << "reject based on num encl in parents" << endl;}
    return false;
  }
  if (encl_reads_perc_parent > options.max_perc_encl_parent) {
    if (options.debug) {cerr << "reject based on perc encl in parents" << endl;}
    return false;
  }
  return true;
}

/*
  Apply naive expansion detection
  Return false if no expansion or filtered
*/
bool DenovoResult::NaiveExpansionDetection(const Options& options, const VCF::Variant& variant) {
  // Extract enclosing
  std::vector<std::string> enclreads;
  variant.get_FORMAT_value_single_string(ENCLREADS_KEY, enclreads);
  int encl_child, total_encl_child, child_encl_match;
  int encl_mother, total_encl_mother, mother_encl_match;
  int encl_father, total_encl_father, father_encl_match;
  float encl_reads_perc_parent = 0;;

  GetEnclosing(enclreads[child_ind_], new_allele_, child_gt_a_, child_gt_b_, &encl_child, &total_encl_child, &child_encl_match);
  GetEnclosing(enclreads[mother_ind_], new_allele_, mat_gt_a_, mat_gt_b_, &encl_mother, &total_encl_mother, &mother_encl_match);
  GetEnclosing(enclreads[father_ind_], new_allele_, pat_gt_a_, pat_gt_b_, &encl_father, &total_encl_father, &father_encl_match);

  // Extract number of FRRs from RC field and FLNKREADS field
  std::vector<std::string> readcounts, flnkreads;
  variant.get_FORMAT_value_single_string(RC_KEY, readcounts);
  variant.get_FORMAT_value_single_string(FLNKREADS_KEY, flnkreads);
  int frr_child, frr_mother, frr_father;
  GetFRR(readcounts[child_ind_], &frr_child);
  GetFRR(readcounts[mother_ind_], &frr_mother);
  GetFRR(readcounts[father_ind_], &frr_father);

  // First, if we see no enclosing reads at all that is a bad sign and we should just quit.
  // Note this might not be always desirable
  if (enclreads[child_ind_] == "NULL" && enclreads[mother_ind_] == "NULL" && enclreads[father_ind_] == "NULL") {
    return false;
  }
  // Test if child has many FRRs with none in parents
  int max_parent_frr = 0; // TODO make an option?
  if (frr_mother <= max_parent_frr && frr_father <= max_parent_frr && frr_child >= options.min_exp_frr) {
    posterior_ = -1; // TODO for now, so we can easily find
    return true;
  }

  // *** Apply naive expansion detection - second looking at flanks *** //
  // First get max parent allele size either in flanks or enclosing
  int max_parent_allele = std::max(std::max(GetMaxFlankAllele(enclreads[mother_ind_]), GetMaxFlankAllele(enclreads[father_ind_])),
           std::max(GetMaxFlankAllele(flnkreads[mother_ind_]), GetMaxFlankAllele(flnkreads[father_ind_])));
  // Then get num child flank reads > max parent allele size
  int num_large_child_flank = GetFlankLargerThan(flnkreads[child_ind_], max_parent_allele);
  if (num_large_child_flank >= options.min_exp_flnk) {
    posterior_ = -1; // TODO for now, so we can easily find
    return true;
  }
  return false;
}

/*
  Get a lot of metadata about the mutation

  1. Infer parent of origin/mutationsize/new allele
  2. Check if mutation should be filtered
  3. Apply naive expansion detection

 */
void DenovoResult::GetMutationInfo(const Options& options, const VCF::Variant& variant,
				   bool* filter_mutation) {
  // Set up
  *filter_mutation = false;
  SetVCFInfo(variant);

  // *** Figure out new allele and POO *** //
  if (!GetPOOMutationInfo(options.chrX)) {
    return; // Mendelian
  }


  // *** Check if we should filter *** //
  if ((child_gt_a_ == new_allele_) && (child_gt_b_ == new_allele_) & options.filter_hom) {
      *filter_mutation = true;
      return;
  }
  if (!CheckReadFilters(options, variant)) {
    *filter_mutation = true;
    return;
  }

  // If mutation passed already, we're good
  if (!filter_mutation) {return;}

  // *** Apply naive expansion detection - first looking at FRR *** //
  // If not, and we didn't return yet, we're doing expansion detection.
  if (!options.naive_expansion_detection) {return;}
  if (!NaiveExpansionDetection(options, variant)) {
    *filter_mutation = true;
    return;
  }
}

void TrioDenovoScanner::summarize_results(std::vector<DenovoResult>& dnr,
					  VCF::Variant& str_variant) {
  LocusInspector locus_inspector;
  if (dnr.empty()) {
    stringstream ss;
    ss << " Skipping - no called children";
    if (options_.debug) PrintMessageDieOnError(ss.str(), M_WARNING);
    return;
  }
  if (options_.debug) PrintMessageDieOnError("Get locus info...", M_PROGRESS);
  // Get locus sinfo
  int32_t start;
  if (options_.gangstr) {
    start = str_variant.get_position();
  } else {
    str_variant.get_INFO_value_single_int(START_KEY, start);
  }
  int32_t end; str_variant.get_INFO_value_single_int(END_KEY, end);
  int32_t period; str_variant.get_INFO_value_single_int(PERIOD_KEY, period);
  int ref_allele_size = (int)str_variant.get_allele(0).size();
  if (options_.debug) PrintMessageDieOnError("Get mutation info...", M_PROGRESS);
  // Get mutation info
  int total_children = 0;
  int total_unaffected = 0;
  int total_affected = 0;
  int num_mutations = 0;
  int num_mutations_unaffected = 0;
  int num_mutations_affected = 0;
  int num_new_affected = 0;
  int num_new_unaffected = 0;
  std::vector<std::string>children_with_mutations;
  // Transmission info
  std::vector<int> unaff_tdt_counts_mat(3, 0);
  std::vector<int> aff_tdt_counts_mat(3, 0);
  std::vector<int> unaff_tdt_counts_pat(3, 0);
  std::vector<int> aff_tdt_counts_pat(3, 0);
  for (auto dnr_iter = dnr.begin(); dnr_iter != dnr.end(); dnr_iter++) {
    total_children++;
    bool filter_mutation = false;
    dnr_iter->GetMutationInfo(options_, str_variant, &filter_mutation);
    if (filter_mutation) {
      dnr_iter->zero_posterior();
    }
    int count_control = 0;
    int count_case = 0;
    int count_unknown = 0;
    bool is_new = false;
    int long_mother = 0;
    int long_father = 0;
    dnr_iter->TestTransmission(&long_mother, &long_father);
    if (dnr_iter->get_poocase() > 1) {
      // Compare to rest of samples
      int new_allele_gb = dnr_iter->get_new_allele()*period - ref_allele_size;
      locus_inspector.GetAlleleCountByPhenotype(str_variant, pedigree_set_.get_families(),
						new_allele_gb, &count_control, &count_case, &count_unknown,
						options_.combine_alleles);
      is_new = (count_unknown == 0);
    } else {
      // Get transmission info. Don't bother if there was a mutation
      dnr_iter->TestTransmission(&long_mother, &long_father);
      if (dnr_iter->get_phenotype() == PT_CASE) {
	aff_tdt_counts_mat[long_mother] += 1;
	aff_tdt_counts_pat[long_father] += 1;
      } else {
	unaff_tdt_counts_mat[long_mother] += 1;
	unaff_tdt_counts_pat[long_father] += 1;
      }
    }
    if (dnr_iter->get_posterior() > options_.posterior_threshold || options_.outputall) {
      all_mutations_file_ << str_variant.get_chromosome() << "\t" << start << "\t"
			  << period  << "\t" << pow(10, dnr_iter->get_prior()) << "\t"
			  << dnr_iter->get_family_id() << "\t" << dnr_iter->get_child_id() << "\t"
			  << dnr_iter->get_phenotype() << "\t" << dnr_iter->get_posterior() << "\t"
			  << dnr_iter->get_new_allele() << "\t" << dnr_iter->get_mut_size() << "\t"
			  << dnr_iter->get_new_allele_in_parents() << "\t" << dnr_iter->get_poocase() << "\t"
			  << is_new << "\t" << count_case << "\t" << count_control << "\t" << count_unknown << "\t"
			  << dnr_iter->get_child_gt() << "\t" << dnr_iter->get_mat_gt() << "\t" << dnr_iter->get_pat_gt() << "\t"
			  << dnr_iter->get_encl_reads_child() << "\t" << dnr_iter->get_encl_reads_mother() << "\t"
			  << dnr_iter->get_encl_reads_father() << "\t" << dnr_iter->get_encl_reads_parent() << "\t"
			  << long_mother << "\t" << long_father
			  << "\n";
      all_mutations_file_.flush();
    }
    if (dnr_iter->get_posterior() > options_.posterior_threshold) {
      num_mutations++;
      children_with_mutations.push_back(dnr_iter->get_family_id() + ":" +
					std::to_string(dnr_iter->get_phenotype()) + ":" + std::to_string(dnr_iter->get_new_allele()) + ":" +
					std::to_string(count_control)+","+std::to_string(count_case) + ","+std::to_string(count_unknown));
      if (dnr_iter->get_phenotype() == PT_CONTROL) {
	num_mutations_unaffected++;
	if (is_new) {
	  num_new_unaffected++;
	}
      }
      if (dnr_iter->get_phenotype() == PT_CASE) {
	num_mutations_affected++;
	if (is_new) {
	  num_new_affected++;
	}
      }
    }
    if (dnr_iter->get_phenotype() == PT_CONTROL) {
      total_unaffected++;
    }
    if (dnr_iter->get_phenotype() == PT_CASE) {
      total_affected++;
    }
  }
  double total_mutation_rate = (double)(num_mutations)/(double)total_children;
  double affected_mutation_rate = (double)(num_mutations_affected)/(double)total_affected;
  double unaffected_mutation_rate = (double)(num_mutations_unaffected)/(double)total_unaffected;
  // Perform fisher exact test
  double fisher_left_p, fisher_right_p, fisher_twosided_p;
  int n11 = num_mutations_affected;
  int n12 = total_affected - num_mutations_affected;
  int n21 = num_mutations_unaffected;
  int n22 = total_unaffected - num_mutations_unaffected;
  kt_fisher_exact(n11, n12, n21, n22, &fisher_left_p, &fisher_right_p, &fisher_twosided_p);
  double p_value = fisher_right_p;
  // Output status
  if (num_mutations > 0) {
    std::stringstream ss;
    ss << "   Analyzed " << total_children << " children\n"
       << "      Total mutation rate: " << total_mutation_rate
       << " (" << num_mutations << "/" << total_children << ")\n"
       << "      Affected mutation rate: " << affected_mutation_rate
       << " (" << num_mutations_affected << "/" << total_affected << ")\n"
       << "      Unaffected mutation rate: " << unaffected_mutation_rate
       << " (" << num_mutations_unaffected << "/" << total_unaffected << ")\n"
       << " P-value: " << p_value;
    PrintMessageDieOnError(ss.str(), M_PROGRESS);
  }
  // Set transmission info (counts for 0=unknown/1=short/2=long mother:father)
  std::string aff_tdt = std::to_string(aff_tdt_counts_mat[0]) + "/" +
    std::to_string(aff_tdt_counts_mat[1]) + "/" +
    std::to_string(aff_tdt_counts_mat[2]) + ":" +
    std::to_string(aff_tdt_counts_pat[0]) + "/" +
    std::to_string(aff_tdt_counts_pat[1]) + "/" +
    std::to_string(aff_tdt_counts_pat[2]);
  std::string unaff_tdt = std::to_string(unaff_tdt_counts_mat[0]) + "/" +
    std::to_string(unaff_tdt_counts_mat[1]) + "/" +
    std::to_string(unaff_tdt_counts_mat[2]) + ":" +
    std::to_string(unaff_tdt_counts_pat[0]) + "/" +
    std::to_string(unaff_tdt_counts_pat[1]) + "/" +
    std::to_string(unaff_tdt_counts_pat[2]);
  // Write output file
  std::string children_with_mutations_string = ".";
  if (children_with_mutations.size() > 0) {
    join(&children_with_mutations_string, children_with_mutations, ",");
  }
  locus_summary_ << str_variant.get_chromosome() << "\t"
		 << start << "\t" << end << "\t" << period << "\t" << ref_allele_size << "\t"
		 << str_variant.num_alleles_by_length(options_.round_alleles) << "\t" << str_variant.num_alleles()  << "\t"
		 << str_variant.heterozygosity_by_length(options_.round_alleles) << "\t" << str_variant.heterozygosity() << "\t"
		 << total_children << "\t" << num_mutations << "\t" << total_mutation_rate << "\t"
		 << total_affected << "\t" << num_mutations_affected << "\t" << num_new_affected << "\t" << affected_mutation_rate << "\t"
		 << total_unaffected << "\t" << num_mutations_unaffected << "\t" << num_new_unaffected << "\t" << unaffected_mutation_rate << "\t"
		 << p_value << "\t" << children_with_mutations_string << "\t"
		 << aff_tdt << "\t" << unaff_tdt << "\n";
  locus_summary_.flush();
}

DenovoResult::DenovoResult(const std::string& family_id,
			   const std::string& mother_id,
			   const std::string& father_id,
			   const std::string& child_id,
			   const int& phenotype,
         const int& child_sex,
			   const int32_t mother_ind,
			   const int32_t father_ind,
			   const int32_t child_ind,
			   const double& total_ll_no_mutation,
			   const double& total_ll_one_denovo,
			   const double& log10prior) {
  family_id_ = family_id;
  mother_id_ = mother_id;
  father_id_ = father_id;
  child_id_ = child_id;
  mother_ind_ = mother_ind;
  father_ind_ = father_ind;
  child_ind_ = child_ind;
  phenotype_ = phenotype;
  child_sex_ = child_sex;
  total_ll_no_mutation_ = total_ll_no_mutation;
  total_ll_one_denovo_ = total_ll_one_denovo;
  log10_prior_mutation_ = log10prior;
  CalculatePosterior();
}

void DenovoResult::CalculatePosterior() {
  double log10_prior_nomut = log10(1-pow(10, log10_prior_mutation_));
  // Note, all likelihoods in log10, but fast_log_sum_exp assumes log_e
  // denom is p(data|nomut)*p(nomut) + p(data|mut)*p(mut)
  // log10 p(data|nomut)*p(nomut) = total_ll_no_mutation_+ log10_prior_nomut
  // log10 p(data|mut)*p(mut) = total_ll_one_denovo_+ log10_prior_nomut
  double denom = fast_log_sum_exp((total_ll_one_denovo_  + log10_prior_mutation_)*log(10),
				                          (total_ll_no_mutation_ + log10_prior_nomut)*log(10)) *log10(exp(1)); // Converts denom to log10
  double posterior = pow(10, total_ll_one_denovo_+log10_prior_mutation_ - denom);
  posterior_ = posterior;
}

DenovoResult::~DenovoResult() {}
