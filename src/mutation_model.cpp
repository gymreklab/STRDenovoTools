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

#include "src/mutation_model.h"

MutationModel::MutationModel(const VCF::Variant& str_variant,
			     MutationPriors& priors,
			     const Options& options) {
  combine_alleles = options.combine_alleles;
  if (options.combine_alleles) {
    assert(str_variant.num_alleles_by_length() > 1);
    log_mut_prior_ = -log10(2) - log10(str_variant.num_alleles_by_length()-1);
    std::string chrom = str_variant.get_chromosome();
    int32_t start;
    str_variant.get_INFO_value_single_int("START", start);
    beta = priors.GetBeta(chrom, start);
    geomp = priors.GetGeomp(chrom, start);
    central_allele = priors.GetCentralAllele(chrom, start);
    str_variant.get_INFO_value_single_int("PERIOD", period);
    ref_allele_size = str_variant.GetSizeFromLengthAllele(0);
  } else {
    assert(str_variant.num_alleles() > 1);
    
    // The allele on each haplotype can mutate to N-1 alleles, so assuming a 
    // uniform prior each mutation has a prior of 1/(2*(N-1))
    log_mut_prior_ = -log10(2) - log10(str_variant.num_alleles()-1);
  } 
}

/*
  Use mutation model to get probability of a mutation step size
  parental and child alleles are in total bp
  central_allele is in num. repeat units compared to reference
 */
double MutationModel::log_prior_mutation(const int& parental_allele, const int& child_allele) {
  if (!combine_alleles) {
    return log_mut_prior_;
  }
  int parental_allele_centered = (parental_allele - ref_allele_size)/period - central_allele;
  int child_allele_centered = (child_allele - ref_allele_size)/period - central_allele;
  double up_prob = (1-beta*geomp*parental_allele_centered)/2;
  int k = child_allele_centered - parental_allele_centered;
  //  std::cerr << parental_allele << " " << child_allele << " " << period << " " << parental_allele_centered << " " << child_allele_centered << std::endl;
  if (k > 0) {
    return log10(up_prob*geomp*pow(1-geomp, k-1));
  } else if (k < 0) {
    return log10((1-up_prob)*geomp*pow(1-geomp, -1*k-1));
  } else {
    // If k=0, set mutation prob to 0
    // This happens for non-unit mutations. e.g. period 5 but get mutation of 1bp
    // TODO how to handle these?
    return -DBL_MAX/2;
    //PrintMessageDieOnError("Encountered mutation of length 0", M_WARNING);
  }
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
