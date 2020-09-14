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

#include <sstream>

#include "common.h"
#include "locus_inspector.h"

LocusInspector::LocusInspector() {}

void LocusInspector::GetAlleleCountByPhenotype(VCF::Variant& variant, std::vector<NuclearFamily> families,
					       int allele,
					       int* count_control, int* count_case, int* count_unknown,
					       const bool& combine_alleles) {
  *count_control = 0;
  *count_case = 0;
  *count_unknown = 0;
  for (auto niter = families.begin(); niter != families.end(); niter++) {
    int m_gt_a, m_gt_b, f_gt_a, f_gt_b;
    std::string mother = (*niter).get_mother();
    std::string father = (*niter).get_father();
    if (!variant.sample_call_missing(mother)) {
      if (combine_alleles) {
	variant.get_length_genotype(mother, m_gt_a, m_gt_b);
      } else {
	variant.get_genotype(mother, m_gt_a, m_gt_b);
      }
      if (m_gt_a == allele) {
	(*count_unknown)++;
      }
      if (m_gt_b == allele) {
	(*count_unknown)++;
      }
    }
    if (!variant.sample_call_missing(father)) {
      if (combine_alleles) {
	variant.get_length_genotype(father, f_gt_a, f_gt_b);
      } else {
	variant.get_genotype(father, f_gt_a, f_gt_b);
      }
      if (f_gt_a == allele) {
	(*count_unknown)++;
      }
      if (f_gt_b == allele) {
	(*count_unknown)++;
      }
    }
    for (size_t i = 0; i < (*niter).num_children(); i++) {
      int c_gt_a, c_gt_b;
      std:string child = (*niter).get_children()[i];
      if (variant.sample_call_missing(child)) {
	continue;
      }
      int pt = (*niter).GetChildrenStatus()[i];
      if (combine_alleles) {
	variant.get_length_genotype(child, c_gt_a, c_gt_b);
      } else {
	variant.get_genotype(child, c_gt_a, c_gt_b);
      }
      if (pt == PT_CASE) {
	if (c_gt_a == allele) {
	  (*count_case)++;
	}
	if (c_gt_b == allele) {
	  (*count_case)++;
	}
      }
      else if (pt == PT_CONTROL) {
	if (c_gt_a == allele) {
	  (*count_control)++;
	}
	if (c_gt_b == allele) {
	  (*count_control)++;
	}
      } else {
	if (c_gt_a == allele) {
	  (*count_unknown)++;
	}
	if (c_gt_b == allele) {
	  (*count_unknown)++;
	}
      }
    }
  }
}

int LocusInspector::GetAlleleCount(VCF::Variant& str_variant, int allele, const bool& combine_alleles) {
  int count = 0;
  for (int i = 0; i < str_variant.num_samples(); i++) {
    int gt_a, gt_b;
    if (!str_variant.sample_call_missing(i)) {
      if (combine_alleles ) {
	str_variant.get_length_genotype(i, gt_a, gt_b);
      } else {
	str_variant.get_genotype(i, gt_a, gt_b);
      }
      if (gt_a == allele) {
	count++;
      }
      if (gt_b == allele) {
	count++;
      }
    }
  }
  return count;
}

void LocusInspector::Inspect(VCF::Variant& str_variant, const std::string& sample_id,
			     const int& sample_index,
			     std::vector<NuclearFamily> families,
			     const std::string& label, const int& status, const bool& combine_alleles) {
  stringstream ss;
  ss << "******* " << label << " (" << sample_id << ") " << "Status=" << status << " ******\n";
  if (str_variant.sample_call_missing(sample_index)) {
    ss << "No call\n";
    PrintMessageDieOnError(ss.str(), M_PROGRESS);
    return;
  }
  // Get Allele (GB and sequence)
  std::vector<std::string> gb_vals;
  std::string allele1, allele2;
  int gt_1, gt_2;
  str_variant.get_genotype(sample_index, gt_1, gt_2);
  str_variant.get_FORMAT_value_single_string("GB", gb_vals);
  allele1 = str_variant.get_alleles()[gt_1];
  allele2 = str_variant.get_alleles()[gt_2];
  int allele1_count = GetAlleleCount(str_variant, gt_1, combine_alleles);
  int allele2_count = GetAlleleCount(str_variant, gt_2, combine_alleles);
  int allele1_count_control, allele1_count_case, allele1_count_unknown;
  int allele2_count_control, allele2_count_case, allele2_count_unknown;
  GetAlleleCountByPhenotype(str_variant, families, gt_1,
			    &allele1_count_control, &allele1_count_case, &allele1_count_unknown, combine_alleles);
  GetAlleleCountByPhenotype(str_variant, families, gt_2,
			    &allele2_count_control, &allele2_count_case, &allele2_count_unknown, combine_alleles);
  ss << "GB: " << gb_vals[sample_index] << "\t (GT: " << gt_1 << "|" << gt_2 << ")\t" 
     << "Count: " << allele1_count << ", " << allele2_count << "\n"
     << "Count by phenotype (ctrl,case,unknown): "
     << allele1_count_control << "," << allele1_count_case << "," << allele1_count_unknown << "\t"
     << allele2_count_control << "," << allele2_count_case << "," << allele2_count_unknown << "\n"
     << allele1 << "\t" << allele2 << "\n";
  // Get quality and othe rinfo (Q, DP, ALLREADS, MALLREADS)
  std::vector<std::string> allreads_vals;
  std::vector<std::string> mallreads_vals;
  std::vector<int32_t> depth_vals;
  std::vector<float> score_vals;
  str_variant.get_FORMAT_value_single_string("ALLREADS", allreads_vals);
  str_variant.get_FORMAT_value_single_string("MALLREADS", mallreads_vals);
  str_variant.get_FORMAT_value_single_int("DP", depth_vals);
  str_variant.get_FORMAT_value_single_float("Q", score_vals);
  ss << "\tDP: " << depth_vals[sample_index] << "\t"
     << "\tQ: " << score_vals[sample_index] << "\n"
     << "\tALLREADS: " << allreads_vals[sample_index] << "\n"
     << "\tMALLREADS: " << mallreads_vals[sample_index] << "\n";
  PrintMessageDieOnError(ss.str(), M_PROGRESS);
}
LocusInspector::~LocusInspector() {}
