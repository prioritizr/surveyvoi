#ifndef HEURISTIC_PRIORITIZATION_H
#define HEURISTIC_PRIORITIZATION_H

#include "package.h"
#include "rcpp_expected_value_of_action.h"
#include "rcpp_expected_value_of_action.h"

void stingy_heuristic_prioritization(
  Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&,
  Eigen::VectorXi&, double, std::vector<bool>&);

void greedy_heuristic_prioritization(
  Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&,
  Eigen::VectorXi&, double, std::vector<bool>&);

#endif
