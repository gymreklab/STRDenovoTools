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

#ifndef MUTATION_MODEL_H_
#define MUTATION_MODEL_H_

#include <math.h>

#include "src/mutation_priors.h"
#include "src/options.h"
#include "src/vcf_reader.h"

class MutationModel {
public:
MutationModel(const VCF::Variant& str_variant, MutationPriors& priors, const Options& options);
~MutationModel();
  /*
   * Log10-likelihood of mutating from the parental to the child allele,
   * given that a mutation occurred
   */
  double log_prior_mutation(const int& parental_allele, const int& child_allele);

  double max_log_prior_mutation(const int& parental_allele);

private:
   bool combine_alleles;
   bool round_alleles;
   double log_mut_prior_;
   double beta;
   double geomp;
   int central_allele;
   int period;
   int ref_allele_size;
};

#endif
