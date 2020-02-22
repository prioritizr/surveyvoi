#include "rcpp_states.h"

void nth_state(mpz_t n, Eigen::MatrixXd &matrix) {
  // if n is equal to zero then simply set all values to zero and exit
  if (mpz_cmp_ui(n, 0) == 0) {
    matrix.setZero();
    return;
  }
  // initialization
  std::size_t n_cells = matrix.size();
  mpz_t mask, bitwise_and, n2;
  mpz_init(mask);
  mpz_init(bitwise_and);
  mpz_init(n2);
  mpz_set(n2, n);
  mpz_set_ui(mask, 1U << (n_cells - 1));
  // store n'th state by converting the integer to binary representation,
  // see https://stackoverflow.com/a/31578829/3483791
  auto j = matrix.data();
  for (std::size_t i = 0; i < n_cells; i++) {
    mpz_and(bitwise_and, n2, mask);
    *j = (mpz_cmp_ui(bitwise_and, 0) > 0) ? 1.0 : 0.0;
    mpz_mul_2exp(n2, n2, 1);
    ++j;
  }
  // clean-up
  mpz_clear(mask);
  mpz_clear(bitwise_and);
  mpz_clear(n2);
  // return void
  return;
}

void nth_state(
  mpz_t n,
  Eigen::MatrixXd &matrix,
  std::vector<std::size_t> &idx) {
  // if n is equal to zero then simply set all values to zero and exit
  if (mpz_cmp_ui(n, 0) == 0) {
    for (auto itr = idx.cbegin(); itr != idx.cend(); ++itr)
      matrix(*itr) = 0.0;
    return;
  }
  // initialization
  std::size_t n_cells = idx.size();
  mpz_t mask, bitwise_and, n2;
  mpz_init(mask);
  mpz_init(bitwise_and);
  mpz_init(n2);
  mpz_set(n2, n);
  mpz_set_ui(mask, 1U << (n_cells - 1));
  // store n'th state by converting the integer to binary representation,
  // see https://stackoverflow.com/a/31578829/3483791
  auto j = idx.cbegin();
  for (std::size_t i = 0; i < n_cells; i++) {
    // original c++ implementation
    /// matrix(*j) = (n & mask) ? 1.0 : 0.0;
    /// n <<= 1;
    /// ++j;

    // gmp implementation for large integers
    mpz_and(bitwise_and, n2, mask);
    matrix(*j) = (mpz_cmp_ui(bitwise_and, 0) > 0) ? 1.0 : 0.0;
    mpz_mul_2exp(n2, n2, 1);
    ++j;
  }
  // clean-up
  mpz_clear(mask);
  mpz_clear(bitwise_and);
  mpz_clear(n2);
  // return void
  return;
}

void n_states(std::size_t x, mpz_t out) {
  // init
  mpz_set_ui(out, 1);
  mpz_t v1, v2, v3;
  mpz_init(v1);
  mpz_init(v2);
  mpz_init(v3);
  // main
  for (std::size_t i = 1; i < x; ++i) {
    // out += (factorial(x) / (factorial(i) * factorial(x - i)));
    factorial(x, v1);
    factorial(i, v2);
    factorial(x - i, v3);
    mpz_mul(v2, v2, v3);
    mpz_divexact(v1, v1, v2);
    mpz_add(out, out, v1);
  }
  // clean-up
  mpz_clear(v1);
  mpz_clear(v2);
  mpz_clear(v3);
  // return void
  return;
}

std::size_t n_states(std::size_t x) {
  mpz_t tmp;
  mpz_init(tmp);
  n_states(x, tmp);
  std::size_t out = mpz_get_ui(tmp);
  mpz_clear(tmp);
  return out;
}

void which_state(
  Eigen::MatrixXd &matrix,
  std::vector<std::size_t> &idx,
  mpz_t out) {
  // init
  mpz_t tmp;
  mpz_init(tmp);
  mpz_set_ui(out, 0);
  // main
  for (auto itr = idx.cbegin(); itr != idx.cend(); ++itr) {
    // out = out * 2 + matrix(*itr);
    mpz_mul_ui(tmp, out, 2);
    mpz_add_ui(out, tmp, static_cast<std::size_t>(matrix(*itr)));
  }
  // clean-up
  mpz_clear(tmp);
  // return void
  return;
}

std::size_t which_state(
  Eigen::MatrixXd &matrix,
  std::vector<std::size_t> &idx) {
  mpz_t tmp;
  mpz_init(tmp);
  which_state(matrix, idx, tmp);
  std::size_t out = mpz_get_ui(tmp);
  mpz_clear(tmp);
  return out;
}

void which_feature_state(
  Eigen::MatrixXd &oij,
  std::vector<std::size_t> &features,
  std::vector<std::size_t> &pu_survey_solution_idx,
  std::vector<std::size_t> &out) {
  // init
  const std::size_t n_f_total = oij.rows();
  const std::size_t n_f = features.size();
  const std::size_t n_pu_survey = pu_survey_solution_idx.size();
  // find out which model coefficients should be used for making predictions,
  // in other words, find out what outcome has been generated for the feature.
  // this involves finding the cell numbers of rij that correspond to
  // the planning units that have been selected for surveying
  std::vector<std::size_t> curr_feature_pu_survey_solution_idx(n_pu_survey);
  for (std::size_t i = 0; i < n_f; ++i) {
    for (std::size_t j = 0; j < n_pu_survey; ++j) {
      curr_feature_pu_survey_solution_idx[j] =
        (n_f_total * pu_survey_solution_idx[j]) + features[i];
      out[i] = which_state(oij, curr_feature_pu_survey_solution_idx);
    }
  }
  // return void
  return;
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_nth_state_sparse(
  std::size_t n, Eigen::MatrixXd matrix, std::vector<std::size_t> idx) {
  // initialization
  for(auto& i : idx)
    i -= 1;
  mpz_t m;
  mpz_init(m);
  mpz_set_ui(m, n);
  // main processing
  nth_state(m, matrix, idx);
  // clean-up
  mpz_clear(m);
  return matrix;
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_nth_state(std::size_t n, Eigen::MatrixXd matrix) {
  // initialization
  mpz_t m;
  mpz_init(m);
  mpz_set_ui(m, n);
  // main processing
  nth_state(m, matrix);
  // clean-up
  mpz_clear(m);
  return matrix;
}

// [[Rcpp::export]]
std::size_t rcpp_n_states(std::size_t n) {
  return n_states(n);
}

// [[Rcpp::export]]
std::size_t rcpp_which_state(
  Eigen::MatrixXd &matrix, std::vector<std::size_t> &idx) {
  return which_state(matrix, idx);
}
