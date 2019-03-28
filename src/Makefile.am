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

AM_CPPFLAGS = -I$(top_srcdir)/src/ -I$(top_srcdir) $(GSL_CFLAGS) $(HTSLIB_CFLAGS)
AM_LDFLAGS = $(GSL_LIBS) $(HTSLIB_LIBS)

bin_PROGRAMS = CookieMonSTR

CookieMonSTR_SOURCES = \
	common.h common.cpp \
	denovo_allele_priors.h denovo_allele_priors.cpp \
	denovo_scanner.h denovo_scanner.cpp \
	fastonebigheader.h \
	locus_inspector.h locus_inspector.cpp \
	main_sdt.cpp \
	mathops.h mathops.cpp \
	mutation_model.h mutation_model.cpp \
	mutation_priors.h mutation_priors.cpp \
	options.h options.cpp \
	pedigree.h pedigree.cpp \
	region.h region.cpp \
	vcf_input.h vcf_input.cpp \
	vcf_reader.h vcf_reader.cpp

CookieMonSTR_CPPFLAGS = $(AM_CPPFLAGS)
CookieMonSTR_CFLAGS = $(AM_CFLAGS)
CookieMonSTR_CXXFLAGS = $(AM_CXXFLAGS)
CookieMonSTR_LDFLAGS = $(AM_LDFLAGS) $(LT_LDFLAGS)
CookieMonSTR_LDADD = $(AM_LDFLAGS) $(LT_LDFLAGS)
