#pragma once
#ifndef UPDATE_POSTERIOR_PROBABILITY_MATRIX_H
#define UPDATE_POSTERIOR_PROBABILITY_MATRIX_H

#include "package.h"
#include "rcpp_probability.h"

void update_posterior_probability_matrix(
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  std::vector<std::size_t>&,
  Eigen::VectorXd&,
  Eigen::VectorXd&,
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  Eigen::MatrixXd&);

#endif
