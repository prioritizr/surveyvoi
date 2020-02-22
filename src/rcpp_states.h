#ifndef STATES_H
#define STATES_H

#include "package.h"
#include "functions.h"

void nth_state(mpz_t, Eigen::MatrixXd&, std::vector<std::size_t>&);

void nth_state(mpz_t, Eigen::MatrixXd&);

std::size_t n_states(std::size_t);

void n_states(std::size_t, mpz_t);

void which_state(Eigen::MatrixXd&, std::vector<std::size_t>&, mpz_t);

std::size_t which_state(Eigen::MatrixXd&, std::vector<std::size_t>&);

void which_feature_state(
  Eigen::MatrixXd&, std::vector<std::size_t>&, std::vector<std::size_t>&,
  std::vector<std::size_t>&);

#endif
