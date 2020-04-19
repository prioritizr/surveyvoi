#ifndef APPROX_EXPECTED_VALUE_OF_ACTION_H
#define APPROX_EXPECTED_VALUE_OF_ACTION_H

#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_sample_states.h"
#include "rcpp_probability.h"
#include "rcpp_conservation_benefit.h"

double approx_expected_value_of_action(
  std::vector<bool>&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::VectorXd&,
  Eigen::VectorXd&, Eigen::VectorXd&, std::vector<mpz_class>&);

void approx_expected_value_of_action_values(
  std::vector<bool>&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::VectorXd&,
  Eigen::VectorXd&, Eigen::VectorXd&, std::vector<mpz_class>&,
  std::vector<double>&, std::vector<double>&, std::vector<double>&, double);

#endif
