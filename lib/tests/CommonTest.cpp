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

#include "lib/common.h"

#include "gtest/gtest.h"

#include <string>
#include <vector>

TEST(CommonTests, CommonTests_Join) {
	std::vector<std::string> join_items;
	std::string joined_string;
	join_items.push_back("1");
	join_items.push_back("2");
	join(&joined_string, join_items, ",");
	ASSERT_EQ(joined_string, "1,2");

	join_items.clear();
	join_items.push_back("");
	join_items.push_back("0");
	join(&joined_string, join_items, ":");
	ASSERT_EQ(joined_string, ":0");
}

TEST(CommonTests, CommonTests_SplitByDelim) {
	std::string mystring = "1,2,3,4";
	std::vector<std::string> split_items;
	split_by_delim(mystring, ',', split_items);
	ASSERT_EQ(split_items.size(), 4);
	ASSERT_EQ(split_items[0], "1");
	ASSERT_EQ(split_items[1], "2");
	ASSERT_EQ(split_items[2], "3");
	ASSERT_EQ(split_items[3], "4");
}