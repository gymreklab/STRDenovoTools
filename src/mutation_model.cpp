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

#include <cfloat>

#include "math.h"

#include "src/mutation_model.h"

MutationModel::MutationModel(const VCF::Variant& str_variant,
			     MutationPriors& priors,
			     const Options& options, const bool& dummy_models) {
  combine_alleles = options.combine_alleles;
  round_alleles = options.round_alleles;
  if (options.combine_alleles) {
    if (!dummy_models) {
      assert(str_variant.num_alleles_by_length(options.round_alleles) > 1);
      log_mut_prior_ = -log10(2) - log10(str_variant.num_alleles_by_length(options.round_alleles)-1);
    }
    std::string chrom = str_variant.get_chromosome();
    int32_t start;
    str_variant.get_INFO_value_single_int("START", start);
    beta = priors.GetBeta(chrom, start);
    geomp = priors.GetGeomp(chrom, start);
    central_allele = priors.GetCentralAllele(chrom, start);
    str_variant.get_INFO_value_single_int("PERIOD", period);
    ref_allele_size = str_variant.GetSizeFromLengthAllele(0);
  } else {
    if (! dummy_models) {
      assert(str_variant.num_alleles() > 1);
    
      // The allele on each haplotype can mutate to N-1 alleles, so assuming a 
      // uniform prior each mutation has a prior of 1/(2*(N-1))
      log_mut_prior_ = -log10(2) - log10(str_variant.num_alleles()-1);
    }
  } 
}

/*
  Use mutation model to get probability of a mutation step size
  parental and child alleles are in bp relative to reference allele
  central_allele is in num. repeat units compared to reference
 */
double MutationModel::log_prior_mutation(const int& parental_allele, const int& child_allele) {
  if (!combine_alleles) {
    return log_mut_prior_;
  }
  int parental_allele_centered = (parental_allele)/period - central_allele;
  int child_allele_centered = (child_allele)/period - central_allele;
  double up_prob = (1-beta*geomp*parental_allele_centered)/2;
  if (up_prob > 1) {
    up_prob = 0.9999; // boundary state
    PrintMessageDieOnError("Encountered boundary allele " + std::to_string(parental_allele_centered), M_WARNING);
  }
  if (up_prob < 0) {
    up_prob = 0.0001; // boundary state
    PrintMessageDieOnError("Encountered boundary allele " + std::to_string(parental_allele_centered), M_WARNING);
  }
  int k = child_allele_centered - parental_allele_centered;
  double logp;
  if (k == 0 && !round_alleles) {
    // Happens for non-unit mutations that round to same unit size
    k = (child_allele-parental_allele > 0) ? 1 : -1;
  }
  if (k > 0) {
    logp = log10(up_prob*geomp*pow(1-geomp, k-1));
  } else if (k < 0) {
    logp = log10((1-up_prob)*geomp*pow(1-geomp, -1*k-1));
  } else {
    PrintMessageDieOnError("Encountered mutation of length 0", M_ERROR);
  }
  if (isnan(logp)) {
    //    std::cerr << k << " " << up_prob << " " << beta << " " << geomp << " " << parental_allele_centered <<  std::endl;
    PrintMessageDieOnError("Encountered nan mutation prob", M_ERROR);
  }
  return logp;
}

double MutationModel::max_log_prior_mutation(const int& parental_allele) {
  if (!combine_alleles) {
    return log_mut_prior_;
  }
  // Max prob is always one unit away. Max direction depends on allele size relative to center
  return max(log_prior_mutation(parental_allele, parental_allele-period),
	     log_prior_mutation(parental_allele, parental_allele+period));
}

MutationModel::~MutationModel() {}
