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

#include "src/options.h"

using namespace std;

Options::Options() {
  strvcf = "";
  famfile = "";
  outprefix = "";
  priors_file = "";
  gangstr = false;
  outputall = false;
  combine_alleles = false;
  round_alleles = false;
  use_pop_priors = false;
  verbose = false;
  require_all_children = false;
  require_num_children = 0;

  max_num_alleles = 25;
  region = "";
  posterior_threshold = 0.9;
  include_invariant = false;
  period = 0;
  min_coverage = 0;
  min_score = 0.0;
  min_span_cov = 0;
  min_supp_reads = 0;

  family = "";
  locus = "";

  default_prior = -5.0;
  default_beta = 0.0;
  default_pgeom = 1.0;
  default_central = 0;

  naive = false;
  min_num_encl_child = 0;
  max_num_encl_parent = 10000;
  max_perc_encl_parent = 1.0;
  min_encl_match = 0.0;
  min_total_encl = 0;
  filter_hom = false;
  naive_expansion_detection = false;
  min_exp_frr = 1;
  min_exp_flnk = 10;
  chrX = false;

  debug = false;
}

Options::~Options() {}

