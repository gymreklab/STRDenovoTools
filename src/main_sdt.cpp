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

#include "monstr/common.h"
#include "monstr/denovo_scanner.h"
#include "monstr/mutation_priors.h"
#include "monstr/options.h"
#include "monstr/pedigree.h"
#include "monstr/vcf_reader.h"
#include "MonSTRConfig.h"

bool file_exists(const std::string& path){
  return (access(path.c_str(), F_OK) != -1);
}

void show_help() {
  std::stringstream help_msg;
  help_msg << "\nSTRDenovoTools [OPTIONS]"
	   << "--str-vcf <STR VCF file>"
	   << "--fam <pedigree file>"
	   << "--out <outprefix>"
	   << "\n\nOptions:\n"
	   << "--gangstr                 Indicates input VCF is from GangSTR\n"
	   << "********* Naive mutation calling ***************\n"
	   << "--naive                    Use naive mutation calling based on hard calls.\n"
	   << "                           Only works with GangSTR input.\n"
	   << "--min-num-encl-child       Require this many enclosing reads supporting de novo allele.\n" 
	   << "                           Only works with GangSTR input.\n"
	   << "--max-num-encl-parent      Discard if more than this many enclosing reads \n"
	   << "                           support de novo allele in parent. Only works with GangSTR input.\n"
	   << "--max-perc-encl-parent     Discard if more than this percentage enclosing reads \n"
	   << "                           support de novo allele in parent. Only works with GangSTR input.\n"
	   << "--min-encl-match           Discard if fewer than this percentage of enclosing reads match the call.\n"
	   << "                           Only works with GangSTR input.\n"
	   << "--min-total-encl           Require this many total enclosing reads in each sample\n"
	   << "                           Only works with GangSTR input.\n"
	   << "--filter-hom               Filter calls where child is homozygous for new allele\n"
	   << "--naive-expansions-frr <int1,int2> Use naive method to detect expansions.\n"
	   << "                           Look for <int1> FRRs in child and none in parent. If not:\n"
	   << "                           Look for <int2> flanks in child greater than largest parent allele\n"
	   << "                           Only works with GangSTR input.\n"
	   << "********* Mutation model ***********************\n"
	   << "--default-prior <FLOAT>    Default log10 mutation rate to use as prior\n"
	   << "--default-beta <FLOAT>     Default value to use for length constraint\n"
	   << "--default-geomp <FLOAT>    Default value to use for step size parameter\n"
	   << "--default-central <INT>    Default central allele\n"
	   << "--mutation-models <STR>    Use per-locus mutation parameters as priors\n"
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
	   << "--round-alleles             Round allele lengths to nearest repeat unit. \n"
	   << "--use-pop-priors            Get genotype priors from population. \n"
	   << "--posterior-threshold <FLOAT>  Cutoff to call something de novo. \n"
	   << "--include-invariant         Output info for loci even if no alt alleles present.\n"
	   << "--chrX                      Call denovos on chrX.\n"
	   << "********* Filtering loci ***********************\n"
	   << "--output-all-loci           Output all loci to mutations file regardless of score.\n"
	   << "********* Other options ************************\n"
	   << "-h,--help      display this help screen\n"
	   << "-v,--verbose   print out useful progress messages\n"
	   << "--debug        print out lots of debugging messages\n"
	   << "--version      print out the version of this software\n\n";
  cerr << help_msg.str();
  exit(1);
}

void parse_commandline_options(int argc, char* argv[], Options* options) {
  enum LONG_OPTIONS {
    OPT_CHRX,
    OPT_NAIVEEXPANSIONFRR,
    OPT_NAIVE,
    OPT_MINNUMENCLCHILD,
    OPT_MAXNUMENCLPARENT,
    OPT_MAXPERCENCLPARENT,
    OPT_MINENCLMATCH,
    OPT_MINTOTALENCL,
    OPT_FILTERHOM,
    OPT_COMBINEALLELES,
    OPT_DEBUG,
    OPT_DEFAULTPRIOR,
    OPT_DEFAULTBETA,
    OPT_DEFAULTPGEOM,
    OPT_DEFAULTCENTRAL,
    OPT_FAM,
    OPT_FAMILY,
    OPT_GANGSTR,
    OPT_HELP,
    OPT_INCLUDEINVARIANT,
    OPT_MAXNUMALLELES,
    OPT_MINCOVERAGE,
    OPT_MINSCORE,
    OPT_MINSPANCOV,
    OPT_MINSUPPREADS,
    OPT_OUT,
    OPT_OUTPUTALL,
    OPT_PERIOD,
    OPT_POSTERIORTHRESHOLD,
    OPT_PRIORS,
    OPT_REGION,
    OPT_REQUIREALLCHILDREN,
    OPT_REQUIRENUMCHILDREN,
    OPT_ROUNDALLELES,
    OPT_STRVCF,
    OPT_USEPOPPRIORS,
    OPT_VERBOSE,
    OPT_VERSION,
  };
  static struct option long_options[] = {
    {"chrX", 0, 0, OPT_CHRX},
    {"naive-expansions-frr", 1, 0, OPT_NAIVEEXPANSIONFRR},
    {"naive", 0, 0, OPT_NAIVE},
    {"min-num-encl-child", 1, 0, OPT_MINNUMENCLCHILD},
    {"max-num-encl-parent", 1, 0, OPT_MAXNUMENCLPARENT},
    {"max-perc-encl-parent", 1, 0, OPT_MAXPERCENCLPARENT},
    {"min-encl-match", 1, 0, OPT_MINENCLMATCH},
    {"min-total-encl", 1, 0, OPT_MINTOTALENCL},
    {"filter-hom", 0, 0, OPT_FILTERHOM},
    {"combine-alleles-by-length", 0, 0, OPT_COMBINEALLELES},
    {"debug", 0, 0, OPT_DEBUG},
    {"default-prior", 1, 0, OPT_DEFAULTPRIOR},
    {"default-beta", 1, 0, OPT_DEFAULTBETA},
    {"default-pgeom", 1, 0, OPT_DEFAULTPGEOM},
    {"default-central", 1, 0, OPT_DEFAULTCENTRAL},
    {"fam", 1, 0, OPT_FAM},
    {"family", 1, 0, OPT_FAMILY},
    {"gangstr", 0, 0, OPT_GANGSTR},
    {"help", 0, 0, OPT_HELP},
    {"include-invariant", 0, 0, OPT_INCLUDEINVARIANT},
    {"max-num-alleles", 1, 0, OPT_MAXNUMALLELES},
    {"min-coverage", 1, 0, OPT_MINCOVERAGE},
    {"min-score", 1, 0, OPT_MINSCORE},
    {"min-span-coverage", 1, 0, OPT_MINSPANCOV},
    {"min-supp-reads", 1, 0, OPT_MINSUPPREADS},
    {"mutation-models", 1, 0, OPT_PRIORS},
    {"out", 1, 0, OPT_OUT},
    {"output-all-loci", 0, 0, OPT_OUTPUTALL},
    {"period", 1, 0, OPT_PERIOD},
    {"posterior-threshold", 1, 0, OPT_POSTERIORTHRESHOLD},
    {"region", 1, 0, OPT_REGION},
    {"require-all-children", 0, 0, OPT_REQUIREALLCHILDREN},
    {"require-num-children", 1, 0, OPT_REQUIRENUMCHILDREN},
    {"round-alleles", 0, 0, OPT_ROUNDALLELES},
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
  std::vector<std::string> eitems;
  while (ch != -1) {
    switch (ch) {
    case OPT_CHRX:
      options->chrX = true;
      break;
    case OPT_NAIVEEXPANSIONFRR:
      options->naive_expansion_detection = true;
      eitems.clear();
      split_by_delim(optarg, ',', eitems);
      options->min_exp_frr = stoi(eitems[0]);
      options->min_exp_flnk = stoi(eitems[1]);
      break;
    case OPT_NAIVE:
      options->naive = true;
      break;
    case OPT_MINNUMENCLCHILD:
      options->min_num_encl_child = atoi(optarg);
      break;
    case OPT_MINENCLMATCH:
      options->min_encl_match = atof(optarg);
      break;
    case OPT_MAXNUMENCLPARENT:
      options->max_num_encl_parent = atoi(optarg);
      break;
    case OPT_MAXPERCENCLPARENT:
      options->max_perc_encl_parent = atof(optarg);
      break;
    case OPT_MINTOTALENCL:
      options->min_total_encl = atoi(optarg);
      break;
    case OPT_FILTERHOM:
      options->filter_hom = true;
      break;
    case OPT_COMBINEALLELES:
      options->combine_alleles = true;
      break;
    case OPT_DEBUG:
      options->debug = true;
      break;
    case OPT_DEFAULTPRIOR:
      options->default_prior = atof(optarg);
      break;
    case OPT_DEFAULTBETA:
      options->default_beta = atof(optarg);
      break;
    case OPT_DEFAULTPGEOM:
      options->default_pgeom = atof(optarg);
      break;
    case OPT_DEFAULTCENTRAL:
      options->default_central = atoi(optarg);
      break;
    case OPT_FAM:
      options->famfile = optarg;
      break;
    case OPT_FAMILY:
      options->family = optarg;
      break;
    case OPT_GANGSTR:
      options->gangstr = true;
      break;
    case OPT_HELP:
    case 'h':
      show_help();
    case OPT_INCLUDEINVARIANT:
      options->include_invariant = true;
      break;
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
    case OPT_OUTPUTALL:
      options->outputall = true;
      break;
    case OPT_PERIOD:
      options->period = atoi(optarg);
      break;
    case OPT_POSTERIORTHRESHOLD:
      options->posterior_threshold = atof(optarg);
      break;
    case OPT_PRIORS:
      options->priors_file = optarg;
      break;
    case OPT_REGION:
      options->region = optarg;
      break;
    case OPT_REQUIREALLCHILDREN:
      options->require_all_children = true;
      break;
    case OPT_REQUIRENUMCHILDREN:
      options->require_num_children = atoi(optarg);
      break;
    case OPT_ROUNDALLELES:
      options->round_alleles = true;
      break;
    case OPT_STRVCF:
      options->strvcf = optarg;
      break;
    case OPT_USEPOPPRIORS:
      options->use_pop_priors = true;
      break;
    case OPT_VERBOSE:
    case 'v':
      options->verbose = true;
      break;
    case OPT_VERSION:
      cout << argv[0] << " Version " << MonSTR_VERSION_MAJOR << "."
	   << MonSTR_VERSION_MINOR << endl;
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

  // Check some options
  if (options.gangstr) {
    if (options.combine_alleles) {
      PrintMessageDieOnError("combine-alleles ignored for GangSTR files", M_WARNING);
      options.combine_alleles = false;
    }
  }

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
  PrintMessageDieOnError("Running de novo analysis...", M_PROGRESS);
  TrioDenovoScanner denovo_scanner(pedigree_set, options);
  if (options.naive) {
    denovo_scanner.naive_scan(strvcf);
  } else {
    // Set up priors
    PrintMessageDieOnError("Opening priors file...", M_PROGRESS);
    MutationPriors priors(options.default_prior, options.default_beta,
			  options.default_pgeom, options.default_central,
			  options.priors_file);
    denovo_scanner.scan(strvcf, priors);
  }

  // Tear down
  total_time = (clock() - total_time)/CLOCKS_PER_SEC;
  stringstream ss;
  ss << "Total run time = " << total_time << " sec";
  PrintMessageDieOnError(ss.str(), M_PROGRESS);
  return 0;
}
