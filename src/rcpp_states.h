#ifndef STATES_H
#define STATES_H

#include "package.h"
#include "functions.h"

void n_states(std::size_t, mpz_t);

void nth_state(mpz_t, Eigen::MatrixXd&);

void nth_state_sparse(mpz_t, Eigen::MatrixXd&, std::vector<std::size_t>&);

void which_state(Eigen::MatrixXd&, mpz_t);

std::size_t which_state(Eigen::MatrixXd&);

void which_state_sparse(Eigen::MatrixXd&, std::vector<std::size_t>&, mpz_t);

std::size_t which_state_sparse(Eigen::MatrixXd&, std::vector<std::size_t>&);

void sample_k_weighted_nth_states(
  std::size_t, Eigen::MatrixXd&, std::vector<mpz_t>&);

void sample_k_uniform_nth_states(
  std::size_t, Eigen::MatrixXd&, std::vector<mpz_t>&);

void sample_k_uniform_no_replacement_nth_states(
  std::size_t, Eigen::MatrixXd&, std::vector<mpz_t>&);

void which_feature_state(
  Eigen::MatrixXd&, std::vector<std::size_t>&, std::vector<std::size_t>&,
  std::vector<std::size_t>&);

#endif
