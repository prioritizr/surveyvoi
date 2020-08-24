#pragma once
#ifndef POSTERIOR_PROBABILITY_MATRIX_H
#define POSTERIOR_PROBABILITY_MATRIX_H

#include "package.h"
#include "functions.h"
#include "rcpp_probability.h"

void initialize_posterior_probability_matrix(
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  std::vector<bool>&,
  std::vector<bool>&,
  std::vector<std::size_t>&,
  Eigen::VectorXd&,
  Eigen::VectorXd&,
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  Eigen::MatrixXd&);

void update_model_posterior_probabilities(
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  std::vector<bool>&,
  std::vector<bool>&,
  std::vector<std::size_t>&,
  Eigen::VectorXd&,
  Eigen::VectorXd&,
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  Eigen::MatrixXd&);

#endif
