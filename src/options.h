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

#ifndef SRC_OPTIONS_H__
#define SRC_OPTIONS_H__

#include <vector>
#include <string>

#include <stdint.h>

class Options {
 public:
  Options();
  virtual ~Options();

  // Input/output paths
  std::string strvcf;
  std::string famfile;
  std::string outprefix;
  std::string priors_file;

  // Filtering samples
  int require_num_children;

  // Filtering calls
  bool require_all_children;
  int32_t min_coverage;
  double min_score;
  int32_t min_span_cov;
  int32_t min_supp_reads;

  // Filtering loci
  std::string region;
  int max_num_alleles;
  int period;

  // Options for denovo calling
  bool use_pop_priors;
  bool combine_alleles;
  double posterior_threshold;

  // Other options
  bool verbose;

  // Options for output
  bool outputall;
  
  // Options for ExamineDenovo
  std::string family;
  std::string locus;

  // Mutation model
  double default_prior;
  double default_beta;
  double default_pgeom;
  int default_central;
};

#endif  // SRC_OPTIONS_H__
