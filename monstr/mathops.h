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

#ifndef MATH_OPS_H_
#define MATH_OPS_H_

#include <math.h>
#include <cstdint>
#include <algorithm>
#include <vector>

extern const double LOG_ONE_HALF;
extern const double TOLERANCE;
extern const double LOG_E_BASE_10;

void precompute_integer_logs();

double int_log(int val);

double sum(const double* begin, const double* end);

double sum(const std::vector<double>& vals);

int sum (const std::vector<bool>& vals);

double log_sum_exp(const double* begin, const double* end);

double log_sum_exp(double log_v1, double log_v2);

double log_sum_exp(double log_v1, double log_v2, double log_v3);

double log_sum_exp(const std::vector<double>& log_vals);

double expected_value(const double* log_probs, const std::vector<int>& vals);

double expected_value(const std::vector<double>& log_likelihoods, const std::vector<int>& vals);

void update_streaming_log_sum_exp(double log_val, double& max_val, double& total);

double finish_streaming_log_sum_exp(double max_val, double total);

void fast_update_streaming_log_sum_exp(double log_val, double& max_val, double& total);

double fast_finish_streaming_log_sum_exp(double max_val, double total);

// To accelerate logsumexp, ignore values if they're 1/1000th or less than the maximum value
const double LOG_THRESH = log(0.001);

double fast_log_sum_exp(double log_v1, double log_v2);
double fast_log_sum_exp(const std::vector<double>& log_vals);

#endif
