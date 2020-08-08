#ifndef HEURISTIC_PRIORITIZATION_H
#define HEURISTIC_PRIORITIZATION_H

#include "package.h"
#include "rcpp_expected_value_of_action.h"

void greedy_heuristic_prioritization(
  Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&,
  Eigen::VectorXi&, double, std::vector<bool>&);

#endif
