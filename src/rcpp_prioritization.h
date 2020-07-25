#ifndef PRIORITIZATION_H
#define PRIORITIZATION_H

#include "package.h"
#include "rcpp_expected_value_of_action.h"

void prioritization(
  Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&,
  Eigen::VectorXi&, double, std::vector<bool>&);

#endif
