#include "functions.h"

void factorial(std::size_t x, mpz_t out) {
  if (x == 0 || x == 1) {
    // return simple answers quickly
    mpz_set_ui(out, x);
    return;
  } else {
    // run calculations for larger numbers
    mpz_set_ui(out, 1);
    for (std::size_t i = x; i > 0; --i)
      mpz_mul_ui(out, out, i);
    return;
  }
}

void log_matrix(Eigen::MatrixXd &x) {
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wenum-compare"
  x.array() = x.array().log();
  #pragma GCC diagnostic pop
  return;
}

void log_1m_matrix(Eigen::MatrixXd &x) {
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wenum-compare"
  x.array() = (1.0 - x.array()).array().log();
  #pragma GCC diagnostic pop
  return;
}

double log_sum(double u, double v) {
  // https://statmodeling.stat.columbia.edu/2016/06/11/log-sum-of-exponentials
  const double m = std::max(u, v);
  return m + std::log(std::exp(u - m) + std::exp(v - m));
}

double log_sum(Eigen::VectorXd &x) {
  double m = x.maxCoeff();
  return m + std::log((x.array() - m).exp().sum());
}

double mean_value(Eigen::VectorXd &x) {
  return x.mean();
}

double standard_error_value(Eigen::VectorXd &x) {
  return std::sqrt(variance_value(x) / static_cast<double>(x.size()));
}

double variance_value(Eigen::VectorXd &x) {
  return (x.array() - x.mean()).array().pow(2).sum() /
         static_cast<double>(x.size() - 1);
}

// [[Rcpp::export]]
double rcpp_standard_error_value(Eigen::VectorXd &x) {
  return standard_error_value(x);
}

// [[Rcpp::export]]
double rcpp_log_sum(Eigen::VectorXd &x) {
  return log_sum(x);
}

void assert_valid_probability_data(Eigen::MatrixXd &x, std::string msg) {
  if ((x.maxCoeff() > (1.0 + 1.0e-15)) | (x.minCoeff() < (0.0 - 1.0e-15)))
    Rcpp::stop(msg);
  return;
}

void assert_valid_probability_data(double x, std::string msg) {
  if ((x > (1.0 + 1.0e-15)) | (x < (0.0 - 1.0e-15)))
    Rcpp::stop(msg);
  return;
}

void extract_k_fold_indices(Rcpp::List &x,
  std::vector<std::vector<std::vector<std::size_t>>> &out) {
  const std::size_t n1 = x.size();
  std::size_t n2;
  std::size_t n3;
  Rcpp::IntegerVector idx;
  Rcpp::List curr_sub_x;
  out.resize(n1);
  for (std::size_t i = 0; i < n1; ++i) {
    curr_sub_x = Rcpp::as<Rcpp::List>(x[i]);
    n2 = curr_sub_x.size();
    out[i].resize(n2);
    for (std::size_t j = 0; j < n2; ++j) {
      idx = Rcpp::as<Rcpp::IntegerVector>(curr_sub_x[j]);
      n3 = idx.size();
      out[i][j].resize(n3);
      for (std::size_t k = 0; k < n3; ++k) {
        out[i][j][k] = idx[k] - 1;
      }
    }
  }
  return;
}

void extract_xgboost_parameters(Rcpp::List &x,
  std::vector<std::vector<std::string>> &out_names,
  std::vector<std::vector<std::string>> &out_values) {
  const std::size_t n1 = x.size();
  Rcpp::CharacterVector curr_sub_names;
  Rcpp::List curr_sub_list;
  std::size_t n2;
  out_names.resize(n1);
  out_values.resize(n1);
  for (std::size_t i = 0; i < n1; ++i) {
    curr_sub_list = Rcpp::as<Rcpp::List>(x[i]);
    n2 = curr_sub_list.size();
    curr_sub_names = curr_sub_list.names();
    out_names[i].reserve(n2);
    out_values[i].reserve(n2);
    for (std::size_t p = 0; p < n2; ++p) {
      out_names[i].push_back(Rcpp::as<std::string>(curr_sub_names[p]));
      out_values[i].push_back(Rcpp::as<std::string>(curr_sub_list[p]));
    }
  }
  return;
}

/* Copyright (c) 2015 by Xgboost Contributors */
// This file contains the customization implementations of R module
// to change behavior of libxgboost
namespace xgboost {
namespace common {

// redirect the nath functions.
bool CheckNAN(double v) {
  return ISNAN(v);
}
#if !defined(XGBOOST_USE_CUDA)
double LogGamma(double v) {
  return R::lgammafn(v);
}
#endif  // !defined(XGBOOST_USE_CUDA)
// customize random engine.
void CustomGlobalRandomEngine::seed(CustomGlobalRandomEngine::result_type val) {
  // ignore the seed
}

// use R's PRNG to replace
CustomGlobalRandomEngine::result_type
CustomGlobalRandomEngine::operator()() {
  return static_cast<result_type>(
      std::floor(unif_rand() * CustomGlobalRandomEngine::max()));
}
}  // namespace common
}  // namespace xgboost
