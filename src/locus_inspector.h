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

#ifndef SRC_LOCUS_INSPECTOR_H__
#define SRC_LOCUS_INSPECTOR_H__

#include <string>

#include "src/pedigree.h"
#include "src/vcf_reader.h"

using namespace std;

class LocusInspector {
 public:
  LocusInspector();
  virtual ~LocusInspector();

  void Inspect(VCF::Variant& str_variant, const std::string& sample_id, const int& sample_index,
	       std::vector<NuclearFamily> families,
	       const std::string& label, const int& status, const bool& combine_alleles);
  void GetAlleleCountByPhenotype(VCF::Variant& variant, std::vector<NuclearFamily> families,
				 int allele,
				 int* count_control, int* count_case, int* count_unknown, const bool& combine_alleles);
 private:
  int GetAlleleCount(VCF::Variant& str_variant, int allele, const bool& combine_alleles);
};

#endif  // SRC_LOCUS_INSPECTOR_H__
