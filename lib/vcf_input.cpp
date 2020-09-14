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

#include "mathops.h"
#include "region.h"
#include "vcf_input.h"

const std::string GT_KEY          = "GT";
const std::string UNPHASED_GL_KEY = "GL";
const std::string GANGSTR_GL_KEY  = "GGL";
const std::string PHASED_GL_KEY   = "PHASEDGL";
const std::string COVERAGE_KEY    = "DP";
const std::string SCORE_KEY       = "Q";
const std::string MALLREADS_KEY   = "MALLREADS";
const std::string GB_KEY          = "GB";
std::string START_INFO_TAG        = "START";
std::string STOP_INFO_TAG         = "END";
const std::string REPCN_KEY       = "REPCN";
const std::string ENCLREADS_KEY   = "ENCLREADS";
const std::string FLNKREADS_KEY   = "FLNKREADS";
const std::string RC_KEY          = "RC";

// Because HipSTR extends putative STR regions if there are nearby indels, the STR coordinates in the VCF may
// not exactly match the original reference region coordinates. As a result, when looking for a particular STR region,
// we look for entries a window around the locus. The size of this window is controlled by this parameter
const int32_t pad = 50;

/*
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
*/

bool GangSTRGL::build(const VCF::Variant& variant) {
  std::vector< std::vector<float> > values;
  std::vector<int32_t> coverage_values;
  std::vector<float> score_values;
  variant.get_FORMAT_value_single_int(COVERAGE_KEY, coverage_values);
  variant.get_FORMAT_value_single_float(SCORE_KEY, score_values);
  // Get GLs - they are stored as a string 
  std::vector<std::string> glstring;
  std::vector<std::string> glvals;
  variant.get_FORMAT_value_single_string(GANGSTR_GL_KEY, glstring);
  for (std::vector<std::string>::iterator it=glstring.begin(); it<glstring.end(); it++) {
    glvals.clear();
    std::vector<float> sample_gl_vals;
    split_by_delim(*it, ',', glvals);
    for (std::vector<std::string>::iterator it2=glvals.begin(); it2<glvals.end(); it2++) {
      sample_gl_vals.push_back(atof((*it2).c_str()));
    }
    values.push_back(sample_gl_vals);
  }

  // More setup
  num_samples_         = 0;
  num_alleles_         = variant.num_gangstr_alleles();
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

bool UnphasedGL::build(const VCF::Variant& variant){
  std::vector< std::vector<float> > values;
  std::vector<int32_t> coverage_values;
  std::vector<float> score_values;
  std::vector<std::string> mallreads_values;
  std::vector<std::string> gb_values;
  variant.get_FORMAT_value_multiple_floats(UNPHASED_GL_KEY, values);
  variant.get_FORMAT_value_single_int(COVERAGE_KEY, coverage_values);
  variant.get_FORMAT_value_single_float(SCORE_KEY, score_values);
  variant.get_FORMAT_value_single_string(MALLREADS_KEY, mallreads_values);
  variant.get_FORMAT_value_single_string(GB_KEY, gb_values);
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
	score_values[vcf_sample_index] < options_.min_score || 
	GetSpanCount(mallreads_values[vcf_sample_index]) < options_.min_span_cov ||
	GetMinAlleleCount(mallreads_values[vcf_sample_index], gb_values[vcf_sample_index]) < options_.min_supp_reads) {
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
					    const std::string& sample,
					    std::vector<float>* gl_by_length) {
  // Initialize
  gl_by_length->clear();
  int numgls = (num_alleles_-1)*(num_alleles_)/2 + num_alleles_;
  for (int i = 0; i < numgls; i++) {
    gl_by_length->push_back(-DBL_MAX/2);
  }
  // Iterate through GLs, add to appropriate length based GL
  //  std::cerr << "debugging GL - " << sample << std::endl;
  for (int a1 = 0; a1 < num_seq_alleles_; a1++) {
    for (int a2 = 0; a2 <= a1; a2++) {
      int a1_length = variant.GetLengthIndexFromGT(a1);
      int a2_length = variant.GetLengthIndexFromGT(a2);
      if (a2_length > a1_length) { // Always require a1_length is max
	a1_length = variant.GetLengthIndexFromGT(a2);
	a2_length = variant.GetLengthIndexFromGT(a1);
      }
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
  std::vector<std::string> mallreads_values;
  std::vector<std::string> gb_values;
  variant.get_FORMAT_value_multiple_floats(UNPHASED_GL_KEY, values);
  variant.get_FORMAT_value_single_int(COVERAGE_KEY, coverage_values);
  variant.get_FORMAT_value_single_float(SCORE_KEY, score_values);
  variant.get_FORMAT_value_single_string(MALLREADS_KEY, mallreads_values);
  variant.get_FORMAT_value_single_string(GB_KEY, gb_values);
  num_samples_         = 0;
  num_alleles_         = variant.num_alleles_by_length(options_.round_alleles);
  num_seq_alleles_     = variant.num_alleles();
  int vcf_sample_index = 0;
  const std::vector<std::string>& samples = variant.get_samples();

  for (auto sample_iter = samples.begin(); sample_iter != samples.end(); ++sample_iter, ++vcf_sample_index){
    if (variant.sample_call_missing(vcf_sample_index)) {
      continue;
    }
    // Apply filters
    if (coverage_values[vcf_sample_index] < options_.min_coverage ||
	score_values[vcf_sample_index] < options_.min_score || 
	GetSpanCount(mallreads_values[vcf_sample_index]) < options_.min_span_cov ||
	GetMinAlleleCount(mallreads_values[vcf_sample_index], gb_values[vcf_sample_index]) < options_.min_supp_reads) {
      continue;
    }
    std::vector<float> gl_by_length;
    convert_gl_to_length(values[vcf_sample_index], variant, (*sample_iter), &gl_by_length);
    unphased_gls_.push_back(gl_by_length);
    sample_indices_[*sample_iter] = num_samples_;

    std::vector<float> max_allele_gl(num_alleles_, -DBL_MAX/2);
    int gl_index = 0;
    for (int i = 0; i < num_alleles_; ++i){
      for (int j = 0; j <= i; ++j, ++gl_index){
	max_allele_gl[i] = std::max(max_allele_gl[i], unphased_gls_.back()[gl_index]);
	max_allele_gl[j] = std::max(max_allele_gl[j], unphased_gls_.back()[gl_index]);
      }
    }
    max_gls_.push_back(max_allele_gl);
    num_samples_++;
  }

  return true;
}

UnphasedLengthGL::~UnphasedLengthGL() {}

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

int GL::GetSpanCount(const std::string& mallreads) {
  int totalreads = 0;
  std::vector<std::string> items;
  split_by_delim(mallreads, ';', items);
  for (auto iter = items.begin(); iter != items.end(); iter++) {
    std::vector<std::string> ainfo;
    split_by_delim((*iter), '|', ainfo);
    if (ainfo.size() != 2) {
      return 0;
    }
    totalreads += atoi(ainfo[1].c_str());
  }
  return totalreads;
}

int GL::GetMinAlleleCount(const std::string& mallreads, const std::string& gbstring) {
  std::vector<std::string> alleles;
  std::vector<std::string> items;
  split_by_delim(gbstring, '|', alleles);
  if (alleles.size() != 2) {
    return 0;
  }
  split_by_delim(mallreads, ';', items);
  int suppreads = INT_MAX;
  bool found_a0 = false;
  bool found_a1 = false;
  for (auto iter = items.begin(); iter != items.end(); iter++) {
    std::vector<std::string> ainfo;
    split_by_delim((*iter), '|', ainfo);
    if (ainfo.size() != 2) {
      return 0;
    }
    int cov = atoi(ainfo[1].c_str());
    if (ainfo[0] == alleles[0]) {
      found_a0 = true;
    }
    if (ainfo[0] == alleles[1]) {
      found_a1 = true;
    }
    if (ainfo[0] == alleles[0] || ainfo[0] == alleles[1]) {
      if (cov < suppreads) {
	suppreads = cov;
      }
    }
  }
  if (!found_a0 or !found_a1) {
    suppreads = 0;
  }
  return suppreads;
}
