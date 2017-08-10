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

#include <err.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>

#include "src/common.h"

using namespace std;

void PrintMessageDieOnError(const string& msg, MSGTYPE msgtype) {
  string typestring = "";
  switch (msgtype) {
  case M_ERROR:
    typestring = "ERROR: ";
    break;
  case M_WARNING:
    typestring = "WARNING: ";
    break;
  case M_PROGRESS:
    typestring = "ProgressMeter: ";
    break;
  case M_DEBUG:
    typestring = "DEBUG: ";
    break;
  default:
    errx(1, "Invalid message type. This should never happen");
  }
  stringstream ss;
  ss  << "[STRDenovoTools"
      << "-" << _GIT_VERSION << "] " << typestring << msg << endl;
  cerr << ss.str();

  if (msgtype == M_ERROR) {
    exit(1);
  }
}
