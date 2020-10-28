/*
Copyright (C) 2020 Melissa Gymrek <mgymrek@ucsd.edu>

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

#include "lib/denovo_scanner.h"

#include "gtest/gtest.h"

#include <string>

// TODO fill in tests for these functions
// DenovoResult::GetRepcn
// DenovoResult::GetTrans
// DenovoResult::TestTransmission
// DenovoResult::GetMutSize
// DenovoResult::GetFRR
// DenovoResult::GetEnclosing
// DenovoResult::GetMutationInfo
// DenovoResult::CalculatePosterior
// TrioDenovoScanner::naive_scan
// TrioDenovoScanner::scan

TEST(GetFollowsMI, GetFollowsMI) {
	// Test TrioDenovoScanner::GetFollowsMI(mother_a, mother_b, father_a, father_b, child_a, child_b, is_chrx, child_sex)

	// Set up TrioDenovoScanner dummy object
	Options options;
	PedigreeSet pset;
	TrioDenovoScanner tds(pset, options);

	// Autosomal
	ASSERT_TRUE(tds.GetFollowsMI(5, 5, 5, 5, 5, 5, false, SEX_MALE));
	ASSERT_TRUE(tds.GetFollowsMI(5, 6, 7, 8, 5, 7, false, SEX_MALE));
	ASSERT_TRUE(tds.GetFollowsMI(5, 6, 7, 8, 6, 7, false, SEX_MALE));
	ASSERT_TRUE(tds.GetFollowsMI(5, 6, 7, 8, 6, 8, false, SEX_MALE));
	ASSERT_FALSE(tds.GetFollowsMI(5, 5, 5, 5, 6, 5, false, SEX_MALE));
	ASSERT_FALSE(tds.GetFollowsMI(5, 5, 5, 5, 5, 6, false, SEX_MALE));
	ASSERT_FALSE(tds.GetFollowsMI(5, 5, 5, 5, 6, 6, false, SEX_MALE));
	ASSERT_FALSE(tds.GetFollowsMI(7, 7, 5, 5, 5, 5, false, SEX_MALE));
	ASSERT_FALSE(tds.GetFollowsMI(5, 5, 7, 7, 5, 5, false, SEX_MALE));

	// ChrX - male
	ASSERT_TRUE(tds.GetFollowsMI(5, 5, 6, 6, 5, 5, true, SEX_MALE));
	ASSERT_TRUE(tds.GetFollowsMI(5, 7, 6, 6, 7, 7, true, SEX_MALE));
	ASSERT_FALSE(tds.GetFollowsMI(5, 7, 6, 6, 6, 7, true, SEX_MALE));
	ASSERT_FALSE(tds.GetFollowsMI(5, 7, 6, 6, 6, 6, true, SEX_MALE));

	// ChrX - female
	ASSERT_FALSE(tds.GetFollowsMI(5, 5, 6, 6, 5, 5, true, SEX_FEMALE));
	ASSERT_FALSE(tds.GetFollowsMI(5, 7, 6, 6, 7, 7, true, SEX_FEMALE));
	ASSERT_TRUE(tds.GetFollowsMI(5, 7, 6, 6, 6, 7, true, SEX_FEMALE));
	ASSERT_TRUE(tds.GetFollowsMI(5, 7, 6, 6, 5, 6, true, SEX_FEMALE));
}

TEST(GetMaxFlankAllele, GetMaxFlankAllele) {
	// Test DenovoResult::GetMaxFlankAllele(const std::string flnkstring)
	// flnkstring is NULL, or allele1,reads1|allele2,reads2...

	// Set up dummy DenovoResult
	DenovoResult dnr("family", "mother", "father", "chid", 1, SEX_FEMALE, 0, 0, 0, 0, 0, -8);

	std::string flnkstring = "NULL";
	ASSERT_EQ(dnr.GetMaxFlankAllele(flnkstring), 0);

	flnkstring = "10,6";
	ASSERT_EQ(dnr.GetMaxFlankAllele(flnkstring), 10);

	flnkstring = "10,6|12,2";
	ASSERT_EQ(dnr.GetMaxFlankAllele(flnkstring), 12);

	flnkstring = "25,100|10,6|12,2";
	ASSERT_EQ(dnr.GetMaxFlankAllele(flnkstring), 25);

	flnkstring = "dummystring";
	ASSERT_DEATH(dnr.GetMaxFlankAllele(flnkstring), "");
}

TEST(GetFlankLargerThan, GetFlankLargerThan) {
	// Test DenovoResult::GetFlankLargerThan(const std::string flnkstring, const int& allele)
	// flnkstring is NULL, or allele1,reads1|allele2,reads2...
	// Count how many flanks support a count > allele

	DenovoResult dnr("family", "mother", "father", "chid", 1, SEX_FEMALE, 0, 0, 0, 0, 0, -8);

	std::string flnkstring = "NULL";
	ASSERT_EQ(dnr.GetFlankLargerThan(flnkstring, 0), 0);

	flnkstring = "10,6";
	ASSERT_EQ(dnr.GetFlankLargerThan(flnkstring, 8), 6);
	ASSERT_EQ(dnr.GetFlankLargerThan(flnkstring, 10), 0);
	ASSERT_EQ(dnr.GetFlankLargerThan(flnkstring, 12), 0);

	flnkstring = "25,100|10,6|12,2";
	ASSERT_EQ(dnr.GetFlankLargerThan(flnkstring, 8), 108);
	ASSERT_EQ(dnr.GetFlankLargerThan(flnkstring, 10), 102);
	ASSERT_EQ(dnr.GetFlankLargerThan(flnkstring, 12), 100);
	ASSERT_EQ(dnr.GetFlankLargerThan(flnkstring, 30), 0);

	flnkstring = "dummystring";
	ASSERT_DEATH(dnr.GetFlankLargerThan(flnkstring, 0), "");
}