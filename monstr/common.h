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

#ifndef SRC_COMMON_H__
#define SRC_COMMON_H__

#include <string>
#include <vector>

// Print msg, exit if error
enum MSGTYPE {
  M_ERROR = 0,
  M_WARNING = 1,
  M_DEBUG = 2,
  M_PROGRESS = 3
};
void PrintMessageDieOnError(const std::string& msg,
                            MSGTYPE msgtype);

void join(std::string* target_string,
	  const std::vector<std::string>& items,
	  const std::string& delim);

void split_by_delim(const std::string &s, char delim, 
		    std::vector<std::string>& substrings);

#endif  // SRC_COMMON_H__
