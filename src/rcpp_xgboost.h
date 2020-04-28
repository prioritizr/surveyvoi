#ifndef XGBOOST_MODEL_H
#define XGBOOST_MODEL_H

#include "package.h"
#include "functions.h"

void fit_xgboost_models_and_assess_performance(
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  MatrixXfRM &,
  std::vector<std::size_t>&,
  std::vector<mpz_class>&,
  std::vector<std::vector<std::string>>&,
  std::vector<std::vector<std::string>>&,
  std::vector<std::size_t>&,
  std::vector<std::vector<std::vector<std::size_t>>>&,
  std::vector<std::vector<std::vector<std::size_t>>>&,
  model_beta_map&,
  model_performance_map&,
  std::vector<std::vector<BoosterHandle> *>&,
  Eigen::VectorXd&,
  Eigen::VectorXd&);

void predict_xgboost_model(
  BoosterHandle&,
  DMatrixHandle&,
  Eigen::VectorXf&);

void xgboost_model_sensitivity_and_specificity(
  Eigen::VectorXf&,
  Eigen::VectorXf&,
  DMatrixHandle&,
  BoosterHandle&,
  double&,
  double&);

#endif
