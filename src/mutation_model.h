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

#include "src/vcf_reader.h"

class MutationModel {
  double log_mut_prior_;

 public:
  explicit MutationModel(const VCF::Variant& str_variant){
    assert(str_variant.num_alleles() > 1);
    
    // The allele on each haplotype can mutate to N-1 alleles, so assuming a 
    // uniform prior each mutation has a prior of 1/(2*(N-1))
    log_mut_prior_ = -log10(2) - log10(str_variant.num_alleles()-1);
  }

  /*
   * Log10-likelihood of mutating from the parental to the child allele,
   * given that a mutation occurred
   */
  double log_prior_mutation(int parental_allele, int child_allele) const {
    return log_mut_prior_;
  }

  double max_log_prior_mutation(int parental_allele) const {
    return log_mut_prior_;
  }
};

#endif
