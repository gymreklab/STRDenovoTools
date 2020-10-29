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
#include "lib/vcf_reader.h"
#include "MonSTRConfig.h"

#include "gtest/gtest.h"

#include <math.h>
#include <sstream>
#include <string>


// ********* tests below don't rely on reading in test files ***** //
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
	ASSERT_DEATH(dnr.GetMaxFlankAllele(flnkstring), "Invalid flankstring encountered .*");
}

TEST(GetFlankLargerThan, GetFlankLargerThan) {
	// Test DenovoResult::GetFlankLargerThan(const std::string flnkstring, const int& allele)
	// flnkstring is NULL, or allele1,reads1|allele2,reads2...
	// Count how many flanks support a count > allele

	// Set up dummy DenovoResult
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
	ASSERT_DEATH(dnr.GetFlankLargerThan(flnkstring, 0), "Invalid flankstring encountered .*");
}


TEST(TestTransmission, TestTransmission) {
	// Test DenovoResult::TestTransmission(int* long_mother, int* long_father)
	// 0=unknown, 1=shorter, 2=longer

	// Set up dummy DenovoResult
	DenovoResult dnr("family", "mother", "father", "chid", 1, SEX_FEMALE, 0, 0, 0, 0, 0, -8);

	// Set up result variables
	int long_mother, long_father;

	// Parents homozygous - unknown
	dnr.set_child_gt(10, 12);
	dnr.set_mat_gt(10, 10);
	dnr.set_pat_gt(12, 12);
	dnr.TestTransmission(&long_mother, &long_father);
	ASSERT_EQ(long_mother, 0);
	ASSERT_EQ(long_father, 0);

	// Mother unknown, father long
	dnr.set_child_gt(10, 12);
	dnr.set_mat_gt(10, 10);
	dnr.set_pat_gt(10, 12);
	dnr.TestTransmission(&long_mother, &long_father);
	ASSERT_EQ(long_mother, 0);
	ASSERT_EQ(long_father, 2);

	// Mother unknown, father short
	dnr.set_child_gt(10, 10);
	dnr.set_mat_gt(10, 10);
	dnr.set_pat_gt(10, 12);
	dnr.TestTransmission(&long_mother, &long_father);
	ASSERT_EQ(long_mother, 0);
	ASSERT_EQ(long_father, 1);

	// Mother long, father unknown
	dnr.set_child_gt(12, 15);
	dnr.set_mat_gt(10, 15);
	dnr.set_pat_gt(12, 12);
	dnr.TestTransmission(&long_mother, &long_father);
	ASSERT_EQ(long_mother, 2);
	ASSERT_EQ(long_father, 0);

	// Mother short, father unknown
	dnr.set_child_gt(12, 10);
	dnr.set_mat_gt(10, 15);
	dnr.set_pat_gt(12, 12);
	dnr.TestTransmission(&long_mother, &long_father);
	ASSERT_EQ(long_mother, 1);
	ASSERT_EQ(long_father, 0);

	// Mother long, father long
	dnr.set_child_gt(12, 15);
	dnr.set_mat_gt(10, 15);
	dnr.set_pat_gt(10, 12);
	dnr.TestTransmission(&long_mother, &long_father);
	ASSERT_EQ(long_mother, 2);
	ASSERT_EQ(long_father, 2);

	// Mother long, father short
	dnr.set_child_gt(15, 10);
	dnr.set_mat_gt(10, 15);
	dnr.set_pat_gt(10, 12);
	dnr.TestTransmission(&long_mother, &long_father);
	ASSERT_EQ(long_mother, 2);
	ASSERT_EQ(long_father, 1);

	// Mother short, father long
	dnr.set_child_gt(10, 12);
	dnr.set_mat_gt(10, 15);
	dnr.set_pat_gt(10, 12);
	dnr.TestTransmission(&long_mother, &long_father);
	ASSERT_EQ(long_mother, 1);
	ASSERT_EQ(long_father, 2);

	// Mother short, father short
	dnr.set_child_gt(10, 10);
	dnr.set_mat_gt(10, 15);
	dnr.set_pat_gt(10, 12);
	dnr.TestTransmission(&long_mother, &long_father);
	ASSERT_EQ(long_mother, 1);
	ASSERT_EQ(long_father, 1);

	// Not Mendelian, unknown
	dnr.set_child_gt(10, 21);
	dnr.set_mat_gt(10, 15);
	dnr.set_pat_gt(10, 12);
	dnr.TestTransmission(&long_mother, &long_father);
	ASSERT_EQ(long_mother, 0);
	ASSERT_EQ(long_father, 0);
}

TEST(TestMutSize, TestMutSize) {
	// Set up dummy DenovoResult
	DenovoResult dnr("family", "mother", "father", "chid", 1, SEX_FEMALE, 0, 0, 0, 0, 0, -8);

	// Test int GetMutSize(const int& new_allele, const int& a1, const int& a2, const int& a3, const int& a4); (autosome)
	ASSERT_EQ(dnr.GetMutSize(10, 10, 10, 10, 15), 0);
	ASSERT_EQ(dnr.GetMutSize(15, 10, 10, 10, 15), 0);
	ASSERT_EQ(dnr.GetMutSize(11, 10, 10, 10, 15), 1);
	ASSERT_EQ(dnr.GetMutSize(12, 10, 10, 10, 15), 2);
	ASSERT_EQ(dnr.GetMutSize(13, 10, 10, 10, 15), -2);
	ASSERT_EQ(dnr.GetMutSize(14, 10, 10, 10, 15), -1);
	ASSERT_EQ(dnr.GetMutSize(5, 10, 10, 10, 15), -5);
	ASSERT_EQ(abs(dnr.GetMutSize(11, 10, 10, 12, 12)), 1);


	// Test int DenovoResult::GetMutSize(const int& new_allele, const int& a1, const int& a2) (chrX)
	ASSERT_EQ(dnr.GetMutSize(10, 10, 15), 0);
	ASSERT_EQ(dnr.GetMutSize(12, 10, 15), 2);
	ASSERT_EQ(dnr.GetMutSize(8, 10, 15), -2);
	ASSERT_EQ(dnr.GetMutSize(16, 10, 15), 1);
}

TEST(TestGetFRR, TestGetFRR) {
	// Test DenovoResult::GetFRR(const std::string& rcstring, int* frrcount)
	// RC string has "(enclosing, spanning, FRR, flanking)"

	// Set up dummy DenovoResult
	DenovoResult dnr("family", "mother", "father", "chid", 1, SEX_FEMALE, 0, 0, 0, 0, 0, -8);
	int frrcount;

	std::string rcstring = "10,20,30,40";
	dnr.GetFRR(rcstring, &frrcount);
	ASSERT_EQ(frrcount, 30);

	rcstring = "10,20,0,10";
	dnr.GetFRR(rcstring, &frrcount);
	ASSERT_EQ(frrcount, 0);

	// Try an illegal rcstring
	rcstring = "";
	ASSERT_DEATH(dnr.GetFRR(rcstring, &frrcount), "Invalid RC string encountered .*");
	rcstring = "0,0,";
	ASSERT_DEATH(dnr.GetFRR(rcstring, &frrcount), "Invalid RC string encountered .*");
}

TEST(TestGetEnclosing, TestGetEnclosing) {
	// Test DenovoResult::GetEnclosing(const std::string& enclstring, int& new_allele,
	//			const int32_t& repcn_a, const int32_t& repcn_b,
	//			int* encl_newallele, int* encl_total, int* encl_match) 

	// Set up dummy DenovoResult
	DenovoResult dnr("family", "mother", "father", "chid", 1, SEX_FEMALE, 0, 0, 0, 0, 0, -8);
	int encl_newallele, encl_total, encl_match;

	std::string enclstring = "4,5|10,6|11,2";
	int repcn_a = 4;
	int repcn_b = 10;
	int new_allele=10;
	dnr.GetEnclosing(enclstring, new_allele, repcn_a, repcn_b,
		&encl_newallele, &encl_total, &encl_match);
	ASSERT_EQ(encl_newallele, 6);
	ASSERT_EQ(encl_total, 13);
	ASSERT_EQ(encl_match, 11);

	enclstring = "NULL";
	dnr.GetEnclosing(enclstring, new_allele, repcn_a, repcn_b,
		&encl_newallele, &encl_total, &encl_match);
	ASSERT_EQ(encl_newallele, 0);
	ASSERT_EQ(encl_total, 0);
	ASSERT_EQ(encl_match, 0);

	enclstring = "4,5|10,6|11,2";
	new_allele = 13;
	dnr.GetEnclosing(enclstring, new_allele, repcn_a, repcn_b,
		&encl_newallele, &encl_total, &encl_match);
	ASSERT_EQ(encl_newallele, 0);
	ASSERT_EQ(encl_total, 13);
	ASSERT_EQ(encl_match, 11);

	enclstring = "dummystring";
	ASSERT_DEATH(dnr.GetEnclosing(enclstring, new_allele, repcn_a, repcn_b,
		&encl_newallele, &encl_total, &encl_match), "Invalid enclreads .*");
}

TEST(TestCalculatePosterior, TestCalculatePosterior) {
	// Test DenovoResult::CalculatePosterior() 

	// Dummy case that forces posterior to 0
	double total_ll_no_mutation = 0;
	double total_ll_one_denovo = 0;
	double log10prior = -8;
	DenovoResult dnr("family", "mother", "father", "chid", 1, SEX_FEMALE, 0, 0, 0, total_ll_no_mutation, total_ll_one_denovo, log10prior);
	dnr.CalculatePosterior();
	ASSERT_NEAR(dnr.get_posterior(), 0, 0.00001);

	// Dummy case that forces posterior to 1
	total_ll_no_mutation = -10000;
	total_ll_one_denovo = 0;
	log10prior = -8;
	DenovoResult dnr2("family", "mother", "father", "chid", 1, SEX_FEMALE, 0, 0, 0, total_ll_no_mutation, total_ll_one_denovo, log10prior);
	dnr2.CalculatePosterior();
	ASSERT_NEAR(dnr2.get_posterior(), 1, 0.00001);

	// zero posterior
	dnr2.zero_posterior();
	ASSERT_EQ(dnr2.get_posterior(), 0);

	// More complicated example
	// posterior = P(mut|data) = P(data|mut)*P(mut)/denom
	// denom = P(data|mut)*P(mut) + P(data|nomut)*P(nomut)
	// Using P(mut)=0.1, P(nomut)=0.9, P(data|mut)=0.20, P(data|nomut)=0.10
	//    num = 0.20*0.1 = 0.02 
	//    denom = 0.20*0.1+0.10*0.9 = 0.11
	//    post = 0.02/0.11 = 0.18181818
	total_ll_no_mutation = log10(0.10);
	total_ll_one_denovo = log10(0.20);
	log10prior = log10(0.10);
	DenovoResult dnr3("family", "mother", "father", "chid", 1, SEX_FEMALE, 0, 0, 0, total_ll_no_mutation, total_ll_one_denovo, log10prior);
	dnr3.CalculatePosterior();
	ASSERT_NEAR(dnr3.get_posterior(), 0.1818, 0.01);
}

// ********* DenovoScannerTest tests need test files ***** //

class DenovoScannerTest : public ::testing::Test {
    protected:
    	void SetUp() override {
    		std::stringstream ss;
			ss << SOURCE_DIR << "/TestData/Q3.filtered.vcf.gz";
			vcffile = ss.str();

			std::stringstream ssped;
			ssped << SOURCE_DIR << "/TestData/TestQuadFam.ped";
			pedfile = ss.str();
    	}
    std::string vcffile;
    std::string pedfile;
};

TEST_F(DenovoScannerTest, GetRepcn) {
	// DenovoResult::GetRepcn(const VCF::Variant& variant, const int32_t& sample_ind,
	//			    int* repcn_a, int* repcn_b)

	// Set up dummy DenovoResult
	DenovoResult dnr("family", "mother", "father", "chid", 1, SEX_FEMALE, 0, 0, 0, 0, 0, -8);

	// Set up VCF
	VCF::VCFReader strvcf(vcffile);
	VCF::Variant str_variant;
	int repcn_a, repcn_b;

	// First variant is nocall in all samples
	strvcf.get_next_variant(&str_variant, true);
	dnr.GetRepcn(str_variant, 0, &repcn_a, &repcn_b);
	ASSERT_EQ(repcn_a, -1);
	ASSERT_EQ(repcn_b, -1);

	// Second variant is period 1. 
	// All samples should have 18,18
	strvcf.get_next_variant(&str_variant, true);
	dnr.GetRepcn(str_variant, 0, &repcn_a, &repcn_b);
	ASSERT_EQ(repcn_a, 18);
	ASSERT_EQ(repcn_b, 18);

	dnr.GetRepcn(str_variant, 1, &repcn_a, &repcn_b);
	ASSERT_EQ(repcn_a, 18);
	ASSERT_EQ(repcn_b, 18);

	dnr.GetRepcn(str_variant, 2, &repcn_a, &repcn_b);
	ASSERT_EQ(repcn_a, 18);
	ASSERT_EQ(repcn_b, 18);

	dnr.GetRepcn(str_variant, 3, &repcn_a, &repcn_b);
	ASSERT_EQ(repcn_a, 18);
	ASSERT_EQ(repcn_b, 18);

	// Third variant is period 4
	// All samples should have 5,5
	strvcf.get_next_variant(&str_variant, true);
	dnr.GetRepcn(str_variant, 0, &repcn_a, &repcn_b);
	ASSERT_EQ(repcn_a, 5);
	ASSERT_EQ(repcn_b, 5);

	dnr.GetRepcn(str_variant, 1, &repcn_a, &repcn_b);
	ASSERT_EQ(repcn_a, 5);
	ASSERT_EQ(repcn_b, 5);

	dnr.GetRepcn(str_variant, 2, &repcn_a, &repcn_b);
	ASSERT_EQ(repcn_a, 5);
	ASSERT_EQ(repcn_b, 5);

	dnr.GetRepcn(str_variant, 3, &repcn_a, &repcn_b);
	ASSERT_EQ(repcn_a, 5);
	ASSERT_EQ(repcn_b, 5);
}

TEST_F(DenovoScannerTest, GetMutationInfo) {
	// TODO - test DenovoResult::GetMutationInfo

	// Set up VCF
	VCF::VCFReader strvcf(vcffile);
}

TEST_F(DenovoScannerTest, NaiveScan) {
	// TODO - test TrioDenovoScanner::naive_scan

	// Set up VCF
	VCF::VCFReader strvcf(vcffile);
}

TEST_F(DenovoScannerTest, Scan) {
	// TODO - test TrioDenovoScanner::scan

	// Set up VCF
	VCF::VCFReader strvcf(vcffile);
}
