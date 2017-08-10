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

#include <iostream>
#include <sstream>

#include "src/common.h"
#include "src/options.h"

void show_help() {
  std::stringstream help_msg;
  help_msg << "\nSTRDenovoTools [OPTIONS]"
	   << " --str-vcf <STR VCF file>"
	   << " --fam <pedigree file>"
	   << " --out <outprefix>"
	   << "\n\nOptions:\n"
	   << "********* Mutation model ***********************\n"
	   << "********* Filtering calls **********************\n"
	   << "********* Filtering samples ********************\n"
	   << "********* Filtering loci ***********************\n"
	   << "********* Parameters for de novo calling *******\n"
	   << "--combine-alleles-by-length Collapse alleles of the same length to one. \n"
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
    OPT_HELP,
    OPT_OUT,
    OPT_STRVCF,
    OPT_VERBOSE,
    OPT_VERSION,
  };
  static struct option long_options[] = {
    {"combine-alleles-by-length", 0, 0, OPT_COMBINEALLELES},
    {"fam", 1, 0, OPT_FAM},
    {"help", 0, 0, OPT_HELP},
    {"out", 1, 0, OPT_OUT},
    {"strvcf", 1, 0, OPT_STRVCF},
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
    case OPT_HELP:
    case 'h':
      show_help();
    case OPT_OUT:
      options->outprefix = optarg;
      break;
    case OPT_STRVCF:
      options->strvcf = optarg;
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
  Options options;
  parse_commandline_options(argc, argv, &options);
  return 0;
}
