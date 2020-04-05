#ifndef EXPECTED_VALUE_OF_ACTION_H
#define EXPECTED_VALUE_OF_ACTION_H

#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_conservation_benefit.h"

double expected_value_of_action(
  std::vector<bool>&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&);

double expected_value_of_action(
  std::vector<bool>&, Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::VectorXd&,
  Eigen::VectorXd&);

#endif
