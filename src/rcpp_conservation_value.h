#ifndef CONSERVATION_BENEFIT_H
#define CONSERVATION_BENEFIT_H

#include "package.h"
#include "functions.h"

double conservation_value_amount(
  double, double, double, double, double);

double conservation_value_state(
  Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&,
  double);

#endif
