#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "package.h"

void factorial(std::size_t, mpz_t);

double log_sum(double, double);

double log_sum(Eigen::VectorXd&);

void assert_valid_probability_data(Eigen::MatrixXd&, std::string msg);

void assert_valid_probability_data(double, std::string msg);

void assert_equal_value(double, double, std::string msg);

void extract_k_fold_indices(Rcpp::List&,
  std::vector<std::vector<std::vector<std::size_t>>>&);

void extract_xgboost_parameters(Rcpp::List&,
  std::vector<std::vector<std::string>>&,
  std::vector<std::vector<std::string>>&);

#endif
