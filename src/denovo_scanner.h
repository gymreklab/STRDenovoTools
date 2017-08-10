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

#ifndef TRIO_DENOVO_SCANNER_H_
#define TRIO_DENOVO_SCANNER_H_

#include <assert.h>

#include <iostream>
#include <vector>
#include <string>

#include "src/options.h"
#include "src/pedigree.h"
#include "src/vcf_reader.h"

class DenovoResult {
 public:
  DenovoResult(const std::string& family_id,
	       const std::string& child_id,
	       const int& phenotype,
	       const double& total_ll_no_mutation,
	       const double& total_ll_one_denovo);
  virtual ~DenovoResult();

  const int& get_phenotype() const {return phenotype_;}
  const double& get_posterior() const {return posterior_;}
  
 private:
  void CalculatePosterior();
  std::string family_id_;
  std::string child_id_;
  int phenotype_;
  double total_ll_no_mutation_;
  double total_ll_one_denovo_;
  double posterior_;
};

class TrioDenovoScanner {
 public:
  TrioDenovoScanner(const PedigreeSet& pedigree_set,
		    const Options& options)
    : pedigree_set_(pedigree_set), options_(options) {}
  virtual ~TrioDenovoScanner();

  void scan(VCF::VCFReader& strvcf);
  void summarize_results(std::vector<DenovoResult>& dnr);
 private:
  static std::string START_KEY, END_KEY;
  PedigreeSet pedigree_set_;
  Options options_;
};


#endif
