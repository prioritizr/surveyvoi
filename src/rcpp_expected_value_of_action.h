#pragma once
#ifndef EXPECTED_VALUE_OF_ACTION_H
#define EXPECTED_VALUE_OF_ACTION_H

#include "package.h"
#include "functions.h"

double expected_value_of_action(
  std::vector<bool>&, Eigen::MatrixXd&, Eigen::VectorXi&);

double approx_expected_value_of_action(
  std::vector<std::vector<double>>&,
  Rcpp::IntegerVector&,
  Rcpp::IntegerVector&,
  std::vector<std::size_t>&
);

double exact_expected_value_of_action(
  Eigen::MatrixXd&, Rcpp::IntegerVector&);

double log_proxy_expected_value_of_action(Eigen::MatrixXd&);

#endif
