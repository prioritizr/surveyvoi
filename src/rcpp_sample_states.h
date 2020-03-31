#ifndef SAMPLE_STATES_H
#define SAMPLE_STATES_H

#include "package.h"
#include "rcpp_states.h"
#include "functions.h"

void sample_n_weighted_states_with_replacement(
  std::size_t, Eigen::MatrixXd&, std::vector<mpz_class>&);

void sample_n_uniform_states_with_replacement(
  std::size_t, Eigen::MatrixXd&, std::vector<mpz_class>&);

void sample_n_uniform_states_without_replacement(
  std::size_t, Eigen::MatrixXd&, std::vector<mpz_class>&);

#endif
