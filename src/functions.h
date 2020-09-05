#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "package.h"

void factorial(std::size_t, mpz_t);

void create_reverse_lookup_id(
  std::vector<bool>&, std::vector<std::size_t>&);

void calculate_survey_tss(
  Eigen::MatrixXd&, std::vector<std::size_t>&,
  Eigen::VectorXd&, Eigen::VectorXd&, Eigen::MatrixXd&);

void set_seed(double);

void log_matrix(Eigen::MatrixXd&);

void log_1m_matrix(Eigen::MatrixXd&);

double log_sum(double, double);

double log_substract(double, double);

double log_sum(Eigen::VectorXd&);

double mean_value(Eigen::VectorXd&);

double standard_error_value(Eigen::VectorXd&);

double variance_value(Eigen::VectorXd&);

void assert_valid_probability_data(Eigen::MatrixXd&, std::string msg);

void assert_valid_probability_data(double, std::string msg);

void extract_k_fold_y_data_from_indices(
  std::vector<std::vector<std::vector<std::size_t>>>&,
  std::vector<std::size_t>&,
  std::vector<std::vector<Eigen::VectorXf>>&);

void extract_k_fold_train_w_data_from_indices(
  Eigen::MatrixXd&,
  std::vector<std::vector<std::vector<std::size_t>>>&,
  std::vector<std::size_t>&,
  std::vector<std::vector<Eigen::VectorXf>>&);

void extract_k_fold_test_w_data_from_indices(
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  std::vector<std::vector<std::vector<std::size_t>>>&,
  std::vector<std::size_t>&,
  std::vector<std::vector<Eigen::VectorXf>>&);

void extract_k_fold_x_data_from_indices(
  MatrixXfRM&,
  std::vector<std::vector<std::vector<std::size_t>>>&,
  std::vector<std::size_t>&,
  std::vector<std::vector<MatrixXfRM>>&);

void initialize_DMatrixHandle(
  Eigen::VectorXf&, Eigen::VectorXf&, MatrixXfRM&, DMatrixHandle&);

void initialize_DMatrixHandle(
  Eigen::VectorXf&, MatrixXfRM&, DMatrixHandle&);

void initialize_DMatrixHandle(
  MatrixXfRM&, DMatrixHandle&);

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
