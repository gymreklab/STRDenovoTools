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

#ifndef REGION_H_
#define REGION_H_

#include <assert.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "src/common.h"

class Region{
 private:
  std::string chrom_, name_;
  int32_t start_, stop_;
  int period_;
 public:
  Region(const std::string& chrom, int32_t start, int32_t stop, int period)
    : chrom_(chrom){
    assert(stop > start);
    start_ = start; stop_ = stop; period_ = period;
  }

  Region(const std::string& chrom, int32_t start, int32_t stop, int period, const std::string& name)
    : chrom_(chrom), name_(name){
    assert(stop > start);
    start_ = start; stop_ = stop; period_ = period;
  }

  const std::string& chrom() const { return chrom_;  }
  const std::string& name()  const { return name_;   }
  int32_t start()            const { return start_;  }
  int32_t  stop()            const { return stop_;   }
  int     period()           const { return period_; }
  Region*   copy()           const { return new Region(chrom_, start_, stop_, period_, name_); }

  void set_start(int32_t start){ start_ = start; }
  void set_stop(int32_t stop)  { stop_  = stop;  }

  std::string str() const {
    std::stringstream ss;
    ss << chrom_ << ":" << start_ << "-" << stop_;
    return ss.str();
  }

  bool operator<(const Region &r) const {
    if (chrom_.compare(r.chrom()) != 0)
      return chrom_.compare(r.chrom()) < 0;
    if (start_ != r.start())
      return start_ < r.start();
    if (stop_ != r.stop())
      return stop_ < r.stop();
    return false;
  }

};

void readRegions(const std::string& input_file, uint32_t max_regions, const std::string& chrom_limit, std::vector<Region>& regions, std::ostream& logger);

void orderRegions(std::vector<Region>& regions);

void orderRegions(std::vector<Region>& input_regions, 
		  std::vector< std::vector<Region> >& output_regions, 
		  std::map<std::string, int>& chrom_order);

class RegionGroup {
  std::vector<Region> regions_;
  std::string chrom_;
  int32_t start_;
  int32_t stop_;

 public:
  explicit RegionGroup(const Region& region){
    regions_.push_back(region);
    chrom_ = region.chrom();
    start_ = region.start();
    stop_  = region.stop();
  }

  const std::vector<Region>& regions() const {
    return regions_;
  }

  const std::string& chrom() const { return chrom_;          }
  int32_t  start()           const { return start_;          }
  int32_t  stop()            const { return stop_;           }
  int      num_regions()     const { return regions_.size(); }

  RegionGroup* copy() const {
    RegionGroup* clone = new RegionGroup(regions_[0]);
    for (int i = 1; i < regions_.size(); i++)
      clone->add_region(regions_[i]);
    return clone;
  }

  void add_region(const Region& region){
    if (region.chrom().compare(chrom_) != 0)
      PrintMessageDieOnError("RegionGroup can only consist of regions on a single chromosome", M_ERROR);
    start_ = std::min(start_, region.start());
    stop_  = std::max(stop_,  region.stop());
    regions_.push_back(region);
    std::sort(regions_.begin(), regions_.end());
  }
};

#endif
