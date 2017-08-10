# Copyright (C) 2017 Melissa Gymrek <mgymrek@ucsd.edu>
#
# This file is part of STRDenovoTools.
#
# STRDenovoTools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STRDenovoTools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with STRDenovoTools.  If not, see <http://www.gnu.org/licenses/>.

HTS_LIBS = -Ilib/htslib -lhts
NLOP_LIBS = -lnlopt -lm

AM_CPPFLAGS = -I$(top_srcdir) $(GSL_CFLAGS)
AM_LDFLAGS = $(GSL_LIBS) $(NLOP_LIBS)

bin_PROGRAMS = STRDenovoTools
noinst_LIBRARIES = libsdt.a

libsdt_a_CPPFLAGS = $(AM_CPPFLAGS)
libsdt_a_SOURCES = \
	common.h common.cpp \
	options.h options.cpp

STRDenovoTools_LDFLAGS = $(AM_LDFLAGS) $(HTS_LIBS)
STRDenovoTools_SOURCES = main_sdt.cpp
STRDenovoTools_CPPFLAGS = $(AM_CPPFLAGS)
STRDenovoTools_CFLAGS = $(AM_CFLAGS)
STRDenovoTools_CXXFLAGS = $(AM_CXXFLAGS)
STRDenovoTools_LDADD = libsdt.a $(GSL_LIBS)

