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

using namespace std;

#include <getopt.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

#include "src/common.h"
#include "src/denovo_scanner.h"
#include "src/options.h"
#include "src/pedigree.h"
#include "src/vcf_reader.h"

bool file_exists(const std::string& path){
  return (access(path.c_str(), F_OK) != -1);
}

void show_help() {
  std::stringstream help_msg;
  help_msg << "\nSTRDenovoTools [OPTIONS]"
	   << " --str-vcf <STR VCF file>"
	   << " --fam <pedigree file>"
	   << " --out <outprefix>"
	   << "\n\nOptions:\n"
	   << "********* Mutation model ***********************\n"
	   << "********* Filtering calls **********************\n"
	   << "--min-coverage <INT>       Discard calls with less than this much coverage\n"
	   << "--min-span-coverage <INT>  Discard calls with less than this many spanning reads\n"
	   << "--min-score <FLOAT>        Discard calls with less than this score\n"
	   << "--min-supp-reads <INT>     Discard calls with an allele supported by less than this many reads\n"
	   << "--require-all-children     Discard loci in family where not all children have calls\n" 
	   << "********* Filtering samples ********************\n"
	   << "--family <STR>             Restrict to analyzing this family\n"
	   << "--require-num-children <INT> Require family to have this many children.\n"
	   << "********* Filtering loci ***********************\n"
	   << "--region <STR>             Restrict to loci in this region (chrom:start-end). \n"
	   << "--period <INT>             Restrict to loci with this motif length.\n"
	   << "--max-num-alleles <INT>    Filter loci with more than this many alleles. \n"
	   << "********* Parameters for de novo calling *******\n"
	   << "--combine-alleles-by-length Collapse alleles of the same length to one. \n"
	   << "--use-pop-priors            Get genotype priors from population. \n"
	   << "--posterior-threshold <FLOAT>  Cutoff to call something de novo. \n"
	   << "********* Other options ************************\n"
	   << "-h,--help      display this help screen\n"
	   << "-v,--verbose   print out useful progress messages\n"
	   << "--version      print out the version of this software\n\n";
  cerr << help_msg.str();
  exit(1);
}

void parse_commandline_options(int argc, char* argv[], Options* options) {
  enum LONG_OPTIONS {
    OPT_COMBINEALLELES,
    OPT_FAM,
    OPT_FAMILY,
    OPT_HELP,
    OPT_MAXNUMALLELES,
    OPT_MINCOVERAGE,
    OPT_MINSCORE,
    OPT_MINSPANCOV,
    OPT_MINSUPPREADS,
    OPT_OUT,
    OPT_PERIOD,
    OPT_POSTERIORTHRESHOLD,
    OPT_REGION,
    OPT_REQUIREALLCHILDREN,
    OPT_REQUIRENUMCHILDREN,
    OPT_STRVCF,
    OPT_USEPOPPRIORS,
    OPT_VERBOSE,
    OPT_VERSION,
  };
  static struct option long_options[] = {
    {"combine-alleles-by-length", 0, 0, OPT_COMBINEALLELES},
    {"fam", 1, 0, OPT_FAM},
    {"family", 1, 0, OPT_FAMILY},
    {"help", 0, 0, OPT_HELP},
    {"max-num-alleles", 1, 0, OPT_MAXNUMALLELES},
    {"min-coverage", 1, 0, OPT_MINCOVERAGE},
    {"min-score", 1, 0, OPT_MINSCORE},
    {"min-span-coverage", 1, 0, OPT_MINSPANCOV},
    {"min-supp-reads", 1, 0, OPT_MINSUPPREADS},
    {"out", 1, 0, OPT_OUT},
    {"period", 1, 0, OPT_PERIOD},
    {"posterior-threshold", 1, 0, OPT_POSTERIORTHRESHOLD},
    {"region", 1, 0, OPT_REGION},
    {"require-all-children", 0, 0, OPT_REQUIREALLCHILDREN},
    {"require-num-children", 1, 0, OPT_REQUIRENUMCHILDREN},
    {"strvcf", 1, 0, OPT_STRVCF},
    {"use-pop-priors", 0, 0, OPT_USEPOPPRIORS},
    {"verbose", 0, 0, OPT_VERBOSE},
    {"version", 0, 0, OPT_VERSION},
    {NULL, no_argument, NULL, 0},
  };
  int ch;
  int option_index = 0;
  ch = getopt_long(argc, argv, "hv?",
                   long_options, &option_index);
  while (ch != -1) {
    switch (ch) {
    case OPT_COMBINEALLELES:
      options->combine_alleles++;
      break;
    case OPT_FAM:
      options->famfile = optarg;
      break;
    case OPT_FAMILY:
      options->family = optarg;
      break;
    case OPT_HELP:
    case 'h':
      show_help();
    case OPT_MAXNUMALLELES:
      options->max_num_alleles = atoi(optarg);
      break;
    case OPT_MINCOVERAGE:
      options->min_coverage = atoi(optarg);
      break;
    case OPT_MINSCORE:
      options->min_score = atof(optarg);
      break;
    case OPT_MINSPANCOV:
      options->min_span_cov = atoi(optarg);
      break;
    case OPT_MINSUPPREADS:
      options->min_supp_reads = atoi(optarg);
      break;
    case OPT_OUT:
      options->outprefix = optarg;
      break;
    case OPT_PERIOD:
      options->period = atoi(optarg);
      break;
    case OPT_POSTERIORTHRESHOLD:
      options->posterior_threshold = atof(optarg);
      break;
    case OPT_REGION:
      options->region = optarg;
      break;
    case OPT_REQUIREALLCHILDREN:
      options->require_all_children++;
      break;
    case OPT_REQUIRENUMCHILDREN:
      options->require_num_children = atoi(optarg);
      break;
    case OPT_STRVCF:
      options->strvcf = optarg;
      break;
    case OPT_USEPOPPRIORS:
      options->use_pop_priors++;
      break;
    case OPT_VERBOSE:
    case 'v':
      options->verbose++;
      break;
    case OPT_VERSION:
      cerr << _GIT_VERSION << endl;
      exit(0);
    case '?':
      show_help();
    default:
      show_help();
    };
    ch = getopt_long(argc, argv, "hv?",
		     long_options, &option_index);
  };
  // Leftover arguments are errors
  if (optind < argc) {
    PrintMessageDieOnError("Unnecessary leftover arguments", M_ERROR);
  }
  // Perform other error checking
  if (options->strvcf.empty()) {
    PrintMessageDieOnError("No --strvcf specified", M_ERROR);
  }
  if (options->famfile.empty()) {
    PrintMessageDieOnError("No --fam option specified", M_ERROR);
  }
  if (options->outprefix.empty()) {
    PrintMessageDieOnError("No --out option specified", M_ERROR);
  }
}

int main(int argc, char* argv[]) {
  // Set up
  double total_time = clock();
  Options options;
  parse_commandline_options(argc, argv, &options);

  // Load STR VCF file
  if (!file_exists(options.strvcf)) {
    PrintMessageDieOnError("STR vcf file " + options.strvcf + " does not exist.", M_ERROR);
  }
  if (!file_exists(options.strvcf + ".tbi")) {
    PrintMessageDieOnError("No .tbi index found for " + options.strvcf, M_ERROR);
  }
  VCF::VCFReader strvcf(options.strvcf);
  std::set<std::string> str_samples(strvcf.get_samples().begin(), strvcf.get_samples().end());
  if (!options.region.empty()) {
    strvcf.set_region(options.region);
  }

  // Extract nuclear families for samples with data
  if (!file_exists(options.famfile)) {
    PrintMessageDieOnError("FAM file " + options.famfile + " does not exist.", M_ERROR);
  }
  PedigreeSet pedigree_set;
  if (!pedigree_set.ExtractFamilies(options.famfile, str_samples, options.require_num_children)) {
    PrintMessageDieOnError("Error extracting families from the pedigree.", M_ERROR);
  }
  pedigree_set.PrintStatus();

  // Run denovo analysis, based on trios
  TrioDenovoScanner denovo_scanner(pedigree_set, options);
  denovo_scanner.scan(strvcf);

  // Tear down
  total_time = (clock() - total_time)/CLOCKS_PER_SEC;
  stringstream ss;
  ss << "Total run time = " << total_time << " sec";
  PrintMessageDieOnError(ss.str(), M_PROGRESS);
  return 0;
}
