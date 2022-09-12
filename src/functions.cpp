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

void create_reverse_lookup_id(
  std::vector<bool> &in, std::vector<std::size_t> &out) {
  /// e.g. if in = [0, 1, 0, 1]
  ///         out = [?, 0, ?, 1]
  ///  note ? = initial number in out, defaults to 1e-5 to deliberately crash
  ///           if there are coding errors
  const std::size_t n = in.size();
  std::size_t counter = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (in[i]) {
      out[i] = counter;
      ++counter;
    }
  }
  return;
}

void calculate_survey_tss(
  Eigen::MatrixXd &nij, std::vector<std::size_t> &survey_features_idx,
  Eigen::VectorXd &survey_sensitivity, Eigen::VectorXd &survey_specificity,
  Eigen::MatrixXd &out) {
  // declare constants
  const std::size_t n_f_survey = survey_features_idx.size();
  const std::size_t n_pu = nij.cols();
  // declare temporary loop variables
  double curr_sens, curr_spec;
  // calculate TSS for survey data at each planning unit for each
  // species that we are surveying
  // note that these TSS account for multiple surveys per site,
  // assuming that each survey is independent
  for (std::size_t i = 0; i < n_f_survey; ++i) {
    for (std::size_t j = 0; j < n_pu; ++j) {
      // calculate sensitivity
      curr_sens = 1.0;
      for (std::size_t r = 0; r < nij(survey_features_idx[i], j); ++r)
        curr_sens *= (1.0 - survey_sensitivity[survey_features_idx[i]]);
      curr_sens = 1.0 - curr_sens;
      // calculate specificity
      curr_spec = 1.0;
      for (std::size_t r = 0; r < nij(survey_features_idx[i], j); ++r)
        curr_spec *= (1.0 - survey_specificity[survey_features_idx[i]]);
      curr_spec = 1.0 - curr_spec;
      // calculate TSS
      out(i, j) = curr_sens + curr_spec - 1.0;
    }
  }
  // return void
  return;
}

void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

void log_matrix(Eigen::MatrixXd &x) {
  x.array() = x.array().log();
  return;
}

void log_1m_matrix(Eigen::MatrixXd &x) {
  x.array() = (1.0 - x.array()).array().log();
  return;
}

double log_sum(double u, double v) {
  // https://statmodeling.stat.columbia.edu/2016/06/11/log-sum-of-exponentials
  const double m = std::max(u, v);
  return m + std::log(std::exp(u - m) + std::exp(v - m));
}

double log_subtract(double u, double v) {
 // https://stackoverflow.com/a/778273/3483791
  return u + std::log1p(-std::exp(v - u));
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
  if ((x.maxCoeff() > (1.0 + 1.0e-15)) || (x.minCoeff() < (0.0 - 1.0e-15)))
    Rcpp::stop(msg);
  return;
}

void assert_valid_probability_data(double x, std::string msg) {
  if ((x > (1.0 + 1.0e-15)) || (x < (0.0 - 1.0e-15)))
    Rcpp::stop(msg);
  return;
}

void extract_list_of_list_of_indices(
  Rcpp::List &x, std::vector<std::vector<std::size_t>> &out) {
  const std::size_t n1 = x.size();
  std::size_t n2;
  Rcpp::IntegerVector idx;
  out.resize(n1);
  for (std::size_t i = 0; i < n1; ++i) {
    idx = Rcpp::as<Rcpp::IntegerVector>(x[i]);
    n2 = idx.size();
    out[i].resize(n2);
    for (std::size_t k = 0; k < n2; ++k) {
        out[i][k] = idx[k] - 1;
    }
  }
  return;
}
