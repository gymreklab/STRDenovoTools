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

#ifndef VCF_INPUT_H_
#define VCF_INPUT_H_

#include <assert.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "src/options.h"
#include "src/region.h"
#include "src/vcf_reader.h"

extern const std::string GENOTYPE_KEY;
extern const std::string UNPHASED_GL_KEY;
extern const std::string PHASED_GL_KEY;
extern const std::string COVERAGE_KEY;
extern const std::string SCORE_KEY;
extern const std::string MALLREADS_KEY;
extern const std::string GB_KEY;
extern std::string START_INFO_TAG;
extern std::string STOP_INFO_TAG;
extern const int32_t pad;

bool read_vcf_alleles(VCF::VCFReader* ref_vcf, const Region& region, std::vector<std::string>& alleles, int32_t& pos);

class GL {
 protected:
  int num_alleles_;
  int num_seq_alleles_;
  int num_samples_;
  std::map<std::string, int> sample_indices_;

 public:
  GL(){
    num_alleles_ = 0;
    num_samples_ = 0;
  }

  bool has_sample(const std::string& sample) const {
    return sample_indices_.find(sample) != sample_indices_.end();
  }

  int get_sample_index(const std::string& sample) const {
    auto sample_iter = sample_indices_.find(sample);
    return (sample_iter == sample_indices_.end() ? -1 : sample_iter->second);
  }

  // Get number of reads spanning the STR
  int GetSpanCount(const std::string& mallreads);
  // Get min number of reads supporting a called allele
  int GetMinAlleleCount(const std::string& mallreads, const std::string& gbstring);

  virtual float get_gl(int sample_index, int gt_a, int gt_b) const = 0;
  virtual float get_max_gl_allele_fixed(int sample_index, int gt_a) const = 0;
};

class UnphasedGL : public GL {
 private:
  std::vector< std::vector<float> > unphased_gls_;
  std::vector< std::vector<float> > max_gls_;

  bool build(const VCF::Variant& variant);
  Options options_;

 public:
  explicit UnphasedGL(const VCF::Variant& variant,
		      const Options& options): options_(options) {
    if (!variant.has_format_field(UNPHASED_GL_KEY))
      PrintMessageDieOnError("Required FORMAT field " + UNPHASED_GL_KEY + " not present in VCF", M_ERROR);
    if (!build(variant))
      PrintMessageDieOnError("Failed to construct UnphasedGL instance from VCF record", M_ERROR);
  }

  float get_gl(int sample_index, int min_gt, int max_gt) const {
    assert(min_gt <= max_gt);
    return unphased_gls_[sample_index][max_gt*(max_gt+1)/2 + min_gt];
  }

  /*
   * For the relevant sample, returns the maximum unphased GL of all genotypes
   * that contain GT_A as an allele
   */
  float get_max_gl_allele_fixed(int sample_index, int gt_a) const {
    return max_gls_[sample_index][gt_a];
  }
};

class UnphasedLengthGL : public GL {
 private:
  std::vector< std::vector<float> > unphased_gls_;
  std::vector< std::vector<float> > max_gls_;

  bool build(const VCF::Variant& variant);
  Options options_;

 public:
  explicit UnphasedLengthGL(const VCF::Variant& variant,
			    const Options& options): options_(options) {
    if (!variant.has_format_field(UNPHASED_GL_KEY))
      PrintMessageDieOnError("Required FORMAT field " + UNPHASED_GL_KEY + " not present in VCF", M_ERROR);
    if (!build(variant))
      PrintMessageDieOnError("Failed to construct UnphasedGL instance from VCF record", M_ERROR);
  }
  ~UnphasedLengthGL();

  float get_gl(int sample_index, int min_gt, int max_gt) const {
    assert(min_gt <= max_gt);
    return unphased_gls_[sample_index][max_gt*(max_gt+1)/2 + min_gt];
  }

  /*
    Convert to new GL field, combining alleles of the same length
   */
  void convert_gl_to_length(const std::vector<float>& gl_vals,
			    const VCF::Variant& variant,
			    const std::string& sample,
			    std::vector<float>* gl_by_length);

  /*
   * For the relevant sample, returns the maximum unphased GL of all genotypes
   * that contain GT_A as an allele
   */
  float get_max_gl_allele_fixed(int sample_index, int gt_a) const {
    return max_gls_[sample_index][gt_a];
  }
};

class PhasedGL : public GL {
 private:
  std::vector< std::vector<float> > phased_gls_;
  std::vector< std::vector<float> > max_gls_one_;
  std::vector< std::vector<float> > max_gls_two_;

  bool build(const VCF::Variant& variant);

 public:
  explicit PhasedGL(const VCF::Variant& variant){
    if (!variant.has_format_field(PHASED_GL_KEY))
      PrintMessageDieOnError("Required FORMAT field " + PHASED_GL_KEY + " not present in VCF", M_ERROR);
    if (!build(variant))
      PrintMessageDieOnError("Failed to construct PhasedGL instance from VCF record", M_ERROR);
  }

  float get_gl(int sample_index, int gt_a, int gt_b) const {
    return phased_gls_[sample_index][gt_a*num_alleles_ + gt_b];
  }

  /*
   * For the relevant sample, returns the maximum phased GL of all genotypes
   * of the form GT_A | X, where X is any valid allele
   */
  float get_max_gl_allele_one_fixed(int sample_index, int gt_a) const {
    return max_gls_one_[sample_index][gt_a];
  }

  /*
   * For the relevant sample, returns the maximum phased GL of all genotypes
   * of the form X | GT_A, where X is any valid allele
   */
  float get_max_gl_allele_two_fixed(int sample_index, int gt_b) const {
    return max_gls_two_[sample_index][gt_b];
  }
};

#endif
