#ifndef EXPECTED_VALUE_OF_ACTION_H
#define EXPECTED_VALUE_OF_ACTION_H

#include "package.h"
#include "functions.h"

double expected_value_of_action(
  std::vector<bool>&, Eigen::MatrixXd&, Eigen::VectorXi&);

double approx_expected_value_of_action(
  Eigen::MatrixXd&, Rcpp::IntegerVector&);

#endif
