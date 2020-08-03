#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "package.h"

void factorial(std::size_t, mpz_t);

void set_seed(double);

void log_matrix(Eigen::MatrixXd&);

void log_1m_matrix(Eigen::MatrixXd&);

double log_sum(double, double);

double log_sum(Eigen::VectorXd&);

double mean_value(Eigen::VectorXd&);

double standard_error_value(Eigen::VectorXd&);

double variance_value(Eigen::VectorXd&);

void assert_valid_probability_data(Eigen::MatrixXd&, std::string msg);

void assert_valid_probability_data(double, std::string msg);

template<typename T>
inline void assert_equal_value(T x, T y, std::string msg) {
  if (std::abs(x - y) > 1.0e-15)
    Rcpp::stop(msg);
  return;
}

template<typename T>
inline void assert_gt_value(T x, T y, std::string msg) {
  if (!(x > y))
    Rcpp::stop(msg);
  return;
}

void extract_k_fold_indices(Rcpp::List&,
  std::vector<std::vector<std::vector<std::size_t>>>&);

void extract_list_of_list_of_indices(
  Rcpp::List&, std::vector<std::vector<std::size_t>>&);

void extract_xgboost_parameters(Rcpp::List&,
  std::vector<std::vector<std::string>>&,
  std::vector<std::vector<std::string>>&);

double convolve_binomial(Eigen::VectorXd&, int);

#endif
