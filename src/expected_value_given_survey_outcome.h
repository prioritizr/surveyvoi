#pragma once
#ifndef EXPECTED_VALUE_OF_ACTION_GIVEN_SURVEY_OUTCOME_H
#define EXPECTED_VALUE_OF_ACTION_GIVEN_SURVEY_OUTCOME_H

#include "package.h"
#include "functions.h"
#include "rcpp_sample_states.h"
#include "rcpp_probability.h"
#include "rcpp_posterior_probability_matrix.h"
#include "rcpp_expected_value_of_action.h"

double exact_expected_value_given_survey_outcome_with_model_estimates(
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  std::vector<bool>&,
  std::vector<std::size_t>&,
  std::vector<std::size_t>&,
  Eigen::VectorXd&,
  Eigen::VectorXd&,
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  Eigen::VectorXi&);

double approx_expected_value_given_survey_outcome_with_model_estimates(
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  std::vector<bool>&,
  std::vector<std::size_t>&,
  std::vector<std::size_t>&,
  Eigen::VectorXd&,
  Eigen::VectorXd&,
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  Eigen::VectorXi&,
  std::size_t,
  std::string&,
  std::size_t);

#endif
