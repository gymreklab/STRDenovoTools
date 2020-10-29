/*
This file is resued from Thomas Willems's HipSTR: https://github.com/tfwillems/HipSTR
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
