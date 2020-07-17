#ifndef APPROX_EXPECTED_VALUE_OF_ACTION_H
#define APPROX_EXPECTED_VALUE_OF_ACTION_H

#include "package.h"
#include "functions.h"
#include "rcpp_conservation_value.h"
#include "rcpp_states.h"
#include "rcpp_sample_states.h"
#include "rcpp_probability.h"

double approx_expected_value_of_action(
  std::vector<bool>&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&, std::size_t);

#endif
