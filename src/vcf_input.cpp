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

#include <algorithm>
#include <assert.h>
#include <cfloat>
#include <math.h>

#include "src/mathops.h"
#include "src/region.h"
#include "src/vcf_input.h"

const std::string GT_KEY          = "GT";
const std::string UNPHASED_GL_KEY = "GL";
const std::string PHASED_GL_KEY   = "PHASEDGL";
const std::string COVERAGE_KEY    = "DP";
const std::string SCORE_KEY       = "Q";
std::string START_INFO_TAG        = "START";
std::string STOP_INFO_TAG         = "END";

// Because HipSTR extends putative STR regions if there are nearby indels, the STR coordinates in the VCF may
// not exactly match the original reference region coordinates. As a result, when looking for a particular STR region,
// we look for entries a window around the locus. The size of this window is controlled by this parameter
const int32_t pad = 50;

bool read_vcf_alleles(VCF::VCFReader* ref_vcf, const Region& region, std::vector<std::string>& alleles, int32_t& pos){
  assert(alleles.size() == 0 && ref_vcf != NULL);
  int32_t pad_start = (region.start() < pad ? 0 : region.start()-pad);
  if (!ref_vcf->set_region(region.chrom(), pad_start, region.stop()+pad)){
    // Retry setting region if chr is in chromosome name
    if (region.chrom().size() <= 3 || region.chrom().substr(0, 3).compare("chr") != 0
	|| !ref_vcf->set_region(region.chrom().substr(3), pad_start, region.stop()+pad)){
      pos     = -1;
      return false;
    }
  }
   
  // Extract STR and ensure the coordinates match
  VCF::Variant variant;
  while (ref_vcf->get_next_variant(variant)){
    // Skip variants without the appropriate INFO fields (as they're not STRs)
    if (!variant.has_info_field(START_INFO_TAG) || !variant.has_info_field(STOP_INFO_TAG))
      continue;

    int32_t str_start, str_stop;
    variant.get_INFO_value_single_int(START_INFO_TAG, str_start);
    variant.get_INFO_value_single_int(STOP_INFO_TAG, str_stop);
    if (str_start == region.start()+1 && str_stop == region.stop()){
      pos = variant.get_position()-1;
      alleles.insert(alleles.end(), variant.get_alleles().begin(), variant.get_alleles().end());
      return true;
    }
    if (variant.get_position() > region.start()+pad)
      break;
  }

  pos = -1;
  return false;
}

bool UnphasedGL::build(const VCF::Variant& variant){
  std::vector< std::vector<float> > values;
  std::vector<int32_t> coverage_values;
  std::vector<float> score_values;
  variant.get_FORMAT_value_multiple_floats(UNPHASED_GL_KEY, values);
  variant.get_FORMAT_value_single_int(COVERAGE_KEY, coverage_values);
  variant.get_FORMAT_value_single_float(SCORE_KEY, score_values);
  num_samples_         = 0;
  num_alleles_         = variant.num_alleles();
  int vcf_sample_index = 0;

  const std::vector<std::string>& samples = variant.get_samples();
  for (auto sample_iter = samples.begin(); sample_iter != samples.end(); ++sample_iter, ++vcf_sample_index){
    if (variant.sample_call_missing(vcf_sample_index)) {
      continue;
    }
    // Apply filters
    if (coverage_values[vcf_sample_index] < options_.min_coverage ||
	score_values[vcf_sample_index] < options_.min_score) {
      continue;
    }
    unphased_gls_.push_back(values[vcf_sample_index]);
    sample_indices_[*sample_iter] = num_samples_++;

    std::vector<float> max_allele_gl(num_alleles_, -DBL_MAX/2);
    int gl_index = 0;
    for (int i = 0; i < num_alleles_; ++i){
      for (int j = 0; j <= i; ++j, ++gl_index){
	max_allele_gl[i] = std::max(max_allele_gl[i], unphased_gls_.back()[gl_index]);
	max_allele_gl[j] = std::max(max_allele_gl[j], unphased_gls_.back()[gl_index]);
      }
    }
    max_gls_.push_back(max_allele_gl);
  }

  return true;
}

void UnphasedLengthGL::convert_gl_to_length(const std::vector<float>& gl_vals,
					    const VCF::Variant& variant,
					    std::vector<float>* gl_by_length) {
  // Initialize
  gl_by_length->clear();
  int numgls = (num_alleles_-1)*(num_alleles_)/2 + num_alleles_;
  for (int i = 0; i < numgls; i++) {
    gl_by_length->push_back(-DBL_MAX/2);
  }
  // Iterate through GLs, add to appropriate length based GL
  for (int a1 = 0; a1 < num_seq_alleles_; a1++) {
    for (int a2 = 0; a2 <= a1; a2++) {
      int a1_length = gt_to_allele_index_.at(a1);
      int a2_length = gt_to_allele_index_.at(a2);
      int oldindex = a1*(a1+1)/2+a2;
      int newindex = a1_length*(a1_length+1)/2 + a2_length;
      float addval = gl_vals[oldindex];
      float currentval = (*gl_by_length)[newindex];
      (*gl_by_length)[newindex] = log_sum_exp(currentval, addval);
    }
  }
}

bool UnphasedLengthGL::build(const VCF::Variant& variant){
  std::vector< std::vector<float> > values;
  std::vector<int32_t> coverage_values;
  std::vector<float> score_values;
  variant.get_FORMAT_value_multiple_floats(UNPHASED_GL_KEY, values);
  variant.get_FORMAT_value_single_int(COVERAGE_KEY, coverage_values);
  variant.get_FORMAT_value_single_float(SCORE_KEY, score_values);
  num_samples_         = 0;
  num_alleles_         = variant.num_alleles_by_length();
  num_seq_alleles_     = variant.num_alleles();
  int vcf_sample_index = 0;

  // Build map of GT->allele length. Assume alleles ordered by length
  // Except reference allele, which is always first
  const std::vector<std::string> alleles = variant.get_alleles();
  int ref_allele_size = (int) alleles.front().size();
  allele_sizes_.resize(num_alleles_, 0);
  allele_sizes_[0] = ref_allele_size;
  gt_to_allele_index_[0] = 0;
  int allele_index = 1;
  int prev_index = 0;
  int allele_size = (int)alleles[1].size();
  for (int i = 1; i < alleles.size(); i++) {
    int len = (int)alleles[i].size();
    assert(len >= allele_size);
    if (len > allele_size) {
      if (len == ref_allele_size) { // If this alleles is ref length
	allele_size = len;
	prev_index = allele_index;
	allele_index = 0;
      } else if (allele_size == ref_allele_size) { // If previous allele was ref length
	allele_size = len;
	allele_index = prev_index + 1;
      } else {
	allele_size = len;
	allele_index++;
      }
    }
    gt_to_allele_index_[i] = allele_index;
    allele_sizes_[allele_index] = allele_size;
  }
  const std::vector<std::string>& samples = variant.get_samples();
  for (auto sample_iter = samples.begin(); sample_iter != samples.end(); ++sample_iter, ++vcf_sample_index){
    if (variant.sample_call_missing(vcf_sample_index)) {
      continue;
    }
    // Apply filters
    if (coverage_values[vcf_sample_index] < options_.min_coverage ||
	score_values[vcf_sample_index] < options_.min_score) {
      continue;
    }
    std::vector<float> gl_by_length;
    convert_gl_to_length(values[vcf_sample_index], variant, &gl_by_length);
    unphased_gls_.push_back(gl_by_length);
    sample_indices_[*sample_iter] = num_samples_++;

    std::vector<float> max_allele_gl(num_alleles_, -DBL_MAX/2);
    int gl_index = 0;
    for (int i = 0; i < num_alleles_; ++i){
      for (int j = 0; j <= i; ++j, ++gl_index){
	max_allele_gl[i] = std::max(max_allele_gl[i], unphased_gls_.back()[gl_index]);
	max_allele_gl[j] = std::max(max_allele_gl[j], unphased_gls_.back()[gl_index]);
      }
    }
    max_gls_.push_back(max_allele_gl);
  }

  return true;
}

bool PhasedGL::build(const VCF::Variant& variant){
  std::vector< std::vector<float> > values;
  variant.get_FORMAT_value_multiple_floats(PHASED_GL_KEY, values);
  num_samples_         = 0;
  num_alleles_         = variant.num_alleles();
  int vcf_sample_index = 0;

  const std::vector<std::string>& samples = variant.get_samples();
  for (auto sample_iter = samples.begin(); sample_iter != samples.end(); ++sample_iter, ++vcf_sample_index){
    if (variant.sample_call_missing(vcf_sample_index))
      continue;
    phased_gls_.push_back(values[vcf_sample_index]);
    sample_indices_[*sample_iter] = num_samples_++;

    std::vector<float> max_gl_one(num_alleles_, -DBL_MAX/2), max_gl_two(num_alleles_, -DBL_MAX/2);
    int gl_index = 0;
    for (int i = 0; i < num_alleles_; ++i){
      for (int j = 0; j < num_alleles_; ++j, ++gl_index){
	max_gl_one[i] = std::max(max_gl_one[i], phased_gls_.back()[gl_index]);
	max_gl_two[j] = std::max(max_gl_two[j], phased_gls_.back()[gl_index]);
      }
    }
    max_gls_one_.push_back(max_gl_one);
    max_gls_two_.push_back(max_gl_two);
  }

  return true;
}
