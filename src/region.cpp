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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "src/region.h"

void readRegions(const std::string& input_file, uint32_t max_regions, const std::string& chrom_limit, std::vector<Region>& regions, std::ostream& logger){
  logger << "Reading region file " << input_file << std::endl;
  std::ifstream input(input_file.c_str());
  if (!input.is_open())
    PrintMessageDieOnError("Failed to open region file", M_ERROR);

  regions.clear();
  std::string line;
  int32_t num_regions = 0;
  while (std::getline(input, line) && regions.size() < max_regions){
    num_regions++;
    std::istringstream iss(line);
    std::string chrom, name;
    int32_t start, stop;
    int period;
    double ref_copy;
    if (!(iss >> chrom >> start >> stop >> period >> ref_copy))
      PrintMessageDieOnError("Improperly formatted region file. \nRequired format is tab-delimited columns CHROM START STOP PERIOD NCOPIES\n Bad line: " + line, M_ERROR);
    if (start < 1)      PrintMessageDieOnError("Improperly formatted region file. \n Region has a START < 1, but START must be >= 1\n Bad line: " + line, M_ERROR);
    if (stop <= start)  PrintMessageDieOnError("Improperly formatted region file. \n Region has a STOP <= START. Bad line: " + line, M_ERROR);
    if (period < 1)     PrintMessageDieOnError("Improperly formatted region file. \n Region has a PERIOD < 1. Bad line: " + line, M_ERROR);
    if (period > 99)     PrintMessageDieOnError("Improperly formatted region file. \n Region has a PERIOD > 99. Bad line: " + line, M_ERROR);

    if (!chrom_limit.empty() && chrom.compare(chrom_limit) != 0)
      continue;
    if (iss >> name)
      regions.push_back(Region(chrom, start-1, stop, period, name));
    else
      regions.push_back(Region(chrom, start-1, stop, period));
  }
  input.close();
  logger << "Region file contains " << num_regions << " regions";
  if (!chrom_limit.empty())
    logger << ", of which " << regions.size() << " were located on the requested chromosome";
  logger << "\n" << std::endl;

  if (!chrom_limit.empty() && regions.empty())
    PrintMessageDieOnError("Region file " + input_file + " did not contain any regions on the requested chromosome: " + chrom_limit, M_ERROR);
}

void orderRegions(std::vector<Region>& regions){
  std::sort(regions.begin(), regions.end());
}

void orderRegions(std::vector<Region>& input_regions, std::vector< std::vector<Region> >& output_regions, std::map<std::string, int>& chrom_order){
  output_regions.clear();
  chrom_order.clear();
  int chrom_count = 0;
  for (auto iter = input_regions.begin(); iter != input_regions.end(); iter++){
    if (chrom_order.find(iter->chrom()) == chrom_order.end()){
      chrom_order[iter->chrom()] = chrom_count++;
      output_regions.push_back(std::vector<Region>());
    }
    output_regions[chrom_order[iter->chrom()]].push_back(*iter);
  }
  for (unsigned int i = 0; i < output_regions.size(); i++)
    std::sort(output_regions[i].begin(), output_regions[i].end());
}
