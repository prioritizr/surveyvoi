#ifndef PROB_H
#define PROB_H

#include "package.h"
#include "functions.h"

void total_probability_of_positive_result(
  Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::MatrixXd&);

void total_probability_of_negative_result(
  Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::MatrixXd&);

void total_probability_of_positive_model_result(
  Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::MatrixXd&);

void total_probability_of_negative_model_result(
  Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::MatrixXd&);

double log_probability_of_outcome(
  Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&,
  std::vector<std::size_t>&);

double log_probability_of_state(
  Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&);

#endif
