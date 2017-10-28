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

bin_PROGRAMS = STRDenovoTools ExamineDenovo
noinst_LIBRARIES = libsdt.a

libsdt_a_CPPFLAGS = $(AM_CPPFLAGS)
libsdt_a_SOURCES = \
	common.h common.cpp \
	denovo_allele_priors.h denovo_allele_priors.cpp \
	denovo_scanner.h denovo_scanner.cpp \
	fastonebigheader.h \
	locus_inspector.h locus_inspector.cpp \
	mathops.h mathops.cpp \
	mutation_model.h mutation_model.cpp \
	mutation_priors.h mutation_priors.cpp \
	options.h options.cpp \
	pedigree.h pedigree.cpp \
	region.h region.cpp \
	vcf_input.h vcf_input.cpp \
	vcf_reader.h vcf_reader.cpp

STRDenovoTools_LDFLAGS = $(AM_LDFLAGS) $(HTS_LIBS)
STRDenovoTools_SOURCES = main_sdt.cpp
STRDenovoTools_CPPFLAGS = $(AM_CPPFLAGS)
STRDenovoTools_CFLAGS = $(AM_CFLAGS)
STRDenovoTools_CXXFLAGS = $(AM_CXXFLAGS)
STRDenovoTools_LDADD = libsdt.a $(GSL_LIBS)

ExamineDenovo_LDFLAGS = $(AM_LDFLAGS) $(HTS_LIBS)
ExamineDenovo_SOURCES = main_examine.cpp
ExamineDenovo_CPPFLAGS = $(AM_CPPFLAGS)
ExamineDenovo_CFLAGS = $(AM_CFLAGS)
ExamineDenovo_CXXFLAGS = $(AM_CXXFLAGS)
ExamineDenovo_LDADD = libsdt.a $(GSL_LIBS)

