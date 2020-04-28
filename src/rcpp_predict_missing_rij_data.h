#ifndef PREDICT_MISSING_RIJ_DATA_H
#define PREDICT_MISSING_RIJ_DATA_H

#include "package.h"
#include "rcpp_xgboost.h"
#include "rcpp_states.h"

void predict_missing_rij_data(
  Eigen::MatrixXd&,
  MatrixXfRM&,
  std::vector<std::size_t>&,
  std::vector<std::size_t>&,
  std::vector<std::vector<BoosterHandle> *> &);

#endif
