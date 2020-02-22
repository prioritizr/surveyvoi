#ifndef XGBOOST_MODEL_H
#define XGBOOST_MODEL_H

#include "package.h"
#include "functions.h"

void fit_xgboost_models_and_assess_performance(
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  MatrixXfRM&,
  std::vector<std::size_t>&,
  std::size_t,
  std::vector<std::vector<std::string>>&,
  std::vector<std::vector<std::string>>&,
  std::vector<std::size_t>&,
  std::vector<std::vector<std::vector<std::size_t>>> &xgb_train_folds,
  std::vector<std::vector<std::vector<std::size_t>>> &xgb_test_folds,
  Eigen::Array<std::vector<BoosterHandle>, Eigen::Dynamic, Eigen::Dynamic>&,
  Eigen::MatrixXd&,
  Eigen::MatrixXd&);

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
