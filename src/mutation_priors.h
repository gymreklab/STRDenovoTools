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

#ifndef MUTATION_PRIORS_H_
#define MUTATION_PRIORS_H_

#include <map>
#include <string>

using namespace std;

struct LocusParams {
  double log10mutrate;
  double beta;
  double geomp;
  int central_allele;
};

class MutationPriors {
 public:
  MutationPriors(const double& _default_prior,
		 const double& _default_beta,
		 const double& _default_pgeom,
		 const double& _default_central,
		 const std::string& _priors_file);
  virtual ~MutationPriors();

  // Get mutation params for the locus
  double GetPrior(const std::string& chrom, const int32_t& start);
  double GetBeta(const std::string& chrom, const int32_t& start);
  double GetGeomp(const std::string& chrom, const int32_t& start);
  int GetCentralAllele(const std::string& chrom, const int32_t& start);

 private:
  std::map< std::pair<std::string, int32_t>, LocusParams> priors_map;
  double default_prior;
  double default_beta;
  double default_geomp;
  double default_central;

  // Load priors from file
  void LoadPriors(const std::string& priors_file);
};

#endif
