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

#include <fstream>
#include <sstream>

#include "src/common.h"
#include "src/mutation_priors.h"

MutationPriors::MutationPriors(const double& _default_prior,
			       const double& _default_beta,
			       const double& _default_pgeom,
			       const double& _default_central,
			       const std::string& _priors_file) {
  default_prior = _default_prior;
  default_beta = _default_beta;
  default_geomp = _default_pgeom;
  default_central = _default_central;
  if (!_priors_file.empty()) {
    LoadPriors(_priors_file);
  }
}

double MutationPriors::GetPrior(const std::string& chrom, const int32_t& start) {
  std::pair<std::string, int32_t> key(chrom, start);
  if (priors_map.find(key) == priors_map.end()) {
    return default_prior;
  } else {
    return priors_map[key].log10mutrate;
  }
}

double MutationPriors::GetBeta(const std::string& chrom, const int32_t& start) {
  std::pair<std::string, int32_t> key(chrom, start);
  if (priors_map.find(key) == priors_map.end()) {
    return default_beta;
  } else {
    return priors_map[key].beta;
  }
}

double MutationPriors::GetGeomp(const std::string& chrom, const int32_t& start) {
  std::pair<std::string, int32_t> key(chrom, start);
  if (priors_map.find(key) == priors_map.end()) {
    return default_geomp;
  } else {
    return priors_map[key].geomp;
  }
}

int MutationPriors::GetCentralAllele(const std::string& chrom, const int32_t& start) {
  std::pair<std::string, int32_t> key(chrom, start);
  if (priors_map.find(key) == priors_map.end()) {
    return default_central;
  } else {
    return priors_map[key].central_allele;
  }
}

void MutationPriors::LoadPriors(const std::string& priors_file) {
  std::ifstream input(priors_file.c_str());
  if (!input.is_open()) {
    PrintMessageDieOnError("Failed to open priors file", M_ERROR);
  }
  std::string line;
  while (std::getline(input, line)) {
    std::istringstream iss(line);
    std::string chrom;
    int32_t start, stop;
    int period, central_allele;
    double rate, beta, geomp;
    if (!(iss >> chrom >> start >> stop >> period >> rate >> beta >> geomp >> central_allele)) {
      PrintMessageDieOnError("Improperly formatted priors file. Should have chrom, start, stop, period, rate, beta, geomp, central_allele (tab delimited).\n Bad line: " + line, M_ERROR);
    }
    if (start < 1)      PrintMessageDieOnError("Improperly formatted region file. \n Region has a START < 1, but START must be >= 1\n Bad line: " + line, M_ERROR);
    if (stop <= start)  PrintMessageDieOnError("Improperly formatted region file. \n Region has a STOP <= START. Bad line: " + line, M_ERROR);
    if (period < 1)     PrintMessageDieOnError("Improperly formatted region file. \n Region has a PERIOD < 1. Bad line: " + line, M_ERROR);
    if (period > 99)     PrintMessageDieOnError("Improperly formatted region file. \n Region has a PERIOD > 99. Bad line: " + line, M_ERROR);
    std::pair <std::string, int32_t> locus(chrom, start);
    LocusParams lp;
    lp.log10mutrate = rate;
    lp.beta = beta;
    lp.geomp = geomp;
    lp.central_allele = central_allele;
    priors_map[locus] = lp;
  }
}

MutationPriors::~MutationPriors() {}
