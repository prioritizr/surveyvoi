#ifndef XGBOOST_MODEL_H
#define XGBOOST_MODEL_H

#include "package.h"
#include "functions.h"

void fit_xgboost_models_and_assess_performance(
  Eigen::MatrixXd&,
  Eigen::MatrixXd&,
  std::vector<std::size_t>&,
  std::vector<mpz_class>&,
  std::vector<std::string>&,
  MatrixXs&,
  std::vector<std::size_t>&,
  std::vector<std::size_t>&,
  std::vector<std::vector<MatrixXfRM>>&,
  std::vector<std::vector<Eigen::VectorXf>>&,
  std::vector<std::vector<Eigen::VectorXf>>&,
  std::vector<std::vector<MatrixXfRM>>&,
  std::vector<std::vector<Eigen::VectorXf>>&,
  std::vector<std::vector<Eigen::VectorXf>>&,
  std::vector<MatrixXfRM>&,
  model_yhat_map&,
  model_performance_map&,
  Eigen::VectorXd&,
  Eigen::VectorXd&);

double xgboost_model_tss(
  Eigen::VectorXf&,
  Eigen::VectorXf&,
  DMatrixHandle&,
  int,
  BoosterHandle&);

void predict_xgboost_model(
  int,
  BoosterHandle&,
  DMatrixHandle&,
  Eigen::VectorXf&);

void xgboost_model_sensitivity_and_specificity(
  Eigen::VectorXf&,
  Eigen::VectorXf&,
  DMatrixHandle&,
  int,
  BoosterHandle&,
  double&,
  double&);

#endif
