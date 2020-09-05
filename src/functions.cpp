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

double log_substract(double u, double v) {
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

double convolve_binomial(Eigen::VectorXd& x, int threshold) {
  // initialize
  std::size_t n = x.size() + 1;
  std::size_t n_copy = x.size();
  Eigen::VectorXd z(n);
  z.setZero();
  z[0] = 1.0;
  Eigen::VectorXd new_z(n);
  new_z[0] = 0.0;
  // perform convolution
  for (std::size_t i = 0; i < x.size(); ++i) {
    std::copy(z.data(), z.data() + n_copy, new_z.data() + 1);
    z = ((1.0 - x[i]) * z.array()) + (x[i] * new_z.array());
  }
  // calculate probability
  return (1.0 - z.segment(0, threshold).sum());
}


void extract_k_fold_y_data_from_indices(
  std::vector<std::vector<std::vector<std::size_t>>> &idx,
  std::vector<std::size_t> &survey_features_idx,
  std::vector<std::vector<Eigen::VectorXf>> &out) {
  // preallocate output
  std::size_t n_pu_in_species_fold;
  out.resize(survey_features_idx.size());
  for (std::size_t i = 0; i < survey_features_idx.size(); ++i) {
    out[i].resize(idx[survey_features_idx[i]].size());
    for (std::size_t j = 0; j < idx[survey_features_idx[i]].size(); ++j) {
      out[i][j].resize(idx[survey_features_idx[i]][j].size() * 2);
    }
  }
  // populate data
  for (std::size_t i = 0; i < survey_features_idx.size(); ++i) {
    for (std::size_t j = 0; j < idx[survey_features_idx[i]].size(); ++j) {
      n_pu_in_species_fold = idx[survey_features_idx[i]][j].size();
      for (std::size_t k = 0; k < n_pu_in_species_fold; ++k) {
        out[i][j][k] = 1.0;
        out[i][j][k + n_pu_in_species_fold] = 0.0;
      }
    }
  }
  // return void
  return;
}

void extract_k_fold_train_w_data_from_indices(
  Eigen::MatrixXd &x,
  std::vector<std::vector<std::vector<std::size_t>>> &idx,
  std::vector<std::size_t> &survey_features_idx,
  std::vector<std::vector<Eigen::VectorXf>> &out) {
  // preallocate output
  std::size_t n_pu_in_species_fold;
  out.resize(survey_features_idx.size());
  for (std::size_t i = 0; i < survey_features_idx.size(); ++i) {
    out[i].resize(idx[survey_features_idx[i]].size());
    for (std::size_t j = 0; j < idx[survey_features_idx[i]].size(); ++j) {
      out[i][j].resize(idx[survey_features_idx[i]][j].size() * 2);
    }
  }
  // populate data
  for (std::size_t i = 0; i < survey_features_idx.size(); ++i) {
    for (std::size_t j = 0; j < idx[survey_features_idx[i]].size(); ++j) {
      n_pu_in_species_fold = idx[survey_features_idx[i]][j].size();
      for (std::size_t k = 0; k < n_pu_in_species_fold; ++k) {
        out[i][j][k] = static_cast<float>(
          100.0 * x(survey_features_idx[i], idx[survey_features_idx[i]][j][k]));
        out[i][j][k + n_pu_in_species_fold] = static_cast<float>(
          100.0 * (1.0 - x(survey_features_idx[i],
                           idx[survey_features_idx[i]][j][k])));

      }

    }
  }
  // return void
  return;
}

void extract_k_fold_test_w_data_from_indices(
  Eigen::MatrixXd &dij,
  Eigen::MatrixXd &nij,
  std::vector<std::vector<std::vector<std::size_t>>> &idx,
  std::vector<std::size_t> &survey_features_idx,
  std::vector<std::vector<Eigen::VectorXf>> &out) {
  // preallocate output
  std::size_t n_pu_in_species_fold;
  out.resize(survey_features_idx.size());
  for (std::size_t i = 0; i < survey_features_idx.size(); ++i) {
    out[i].resize(idx[survey_features_idx[i]].size());
    for (std::size_t j = 0; j < idx[survey_features_idx[i]].size(); ++j) {
      out[i][j].resize(idx[survey_features_idx[i]][j].size() * 2);
    }
  }
  // populate data
  for (std::size_t i = 0; i < survey_features_idx.size(); ++i) {
    for (std::size_t j = 0; j < idx[survey_features_idx[i]].size(); ++j) {
      n_pu_in_species_fold = idx[survey_features_idx[i]][j].size();
      // copy data
      for (std::size_t k = 0; k < n_pu_in_species_fold; ++k) {
        out[i][j][k] = static_cast<float>(
          dij(survey_features_idx[i], idx[survey_features_idx[i]][j][k]));
        out[i][j][k + n_pu_in_species_fold] = static_cast<float>(
          1.0 - dij(survey_features_idx[i], idx[survey_features_idx[i]][j][k]));
      }
    }
  }
  // return void
  return;
}

void extract_k_fold_x_data_from_indices(
  MatrixXfRM &x,
  std::vector<std::vector<std::vector<std::size_t>>> &idx,
  std::vector<std::size_t> &survey_features_idx,
  std::vector<std::vector<MatrixXfRM>> &out) {
  // preallocate output
  std::size_t n_pu_in_species_fold;
  out.resize(survey_features_idx.size());
  for (std::size_t i = 0; i < survey_features_idx.size(); ++i) {
    out[i].resize(idx[survey_features_idx[i]].size() * 2);
  }
  // populate data
  for (std::size_t i = 0; i < survey_features_idx.size(); ++i) {
    for (std::size_t j = 0; j < idx[survey_features_idx[i]].size(); ++j) {
      n_pu_in_species_fold = idx[survey_features_idx[i]][j].size();
      out[i][j].resize(n_pu_in_species_fold * 2, x.cols());
      for (std::size_t k = 0; k < n_pu_in_species_fold; ++k) {
        out[i][j].row(k).array() =
          x.row(idx[survey_features_idx[i]][j][k]).array();
        out[i][j].row(k + n_pu_in_species_fold).array() =
          x.row(idx[survey_features_idx[i]][j][k]).array();
      }
    }
  }
  // return void
  return;
}


void initialize_DMatrixHandle(
  Eigen::VectorXf &y, Eigen::VectorXf &w, MatrixXfRM &x, DMatrixHandle &h) {
  XGDMatrixCreateFromMat((float *) x.data(), x.rows(), x.cols(), -99999.0f, &h);
  XGDMatrixSetFloatInfo(h, "label", y.data(), y.size());
  XGDMatrixSetFloatInfo(h, "weight", w.data(), w.size());
  return;
}

void initialize_DMatrixHandle(
  Eigen::VectorXf &w, MatrixXfRM &x, DMatrixHandle &h) {
  XGDMatrixCreateFromMat((float *) x.data(), x.rows(), x.cols(), -99999.0f, &h);
  XGDMatrixSetFloatInfo(h, "weight", w.data(), w.size());
  return;
}

void initialize_DMatrixHandle(
  MatrixXfRM &x, DMatrixHandle &h) {
  XGDMatrixCreateFromMat((float *) x.data(), x.rows(), x.cols(), -99999.0f, &h);
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
