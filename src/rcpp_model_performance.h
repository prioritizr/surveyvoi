#pragma once
#ifndef MODEL_PERFORMANCE_H
#define MODEL_PERFORMANCE_H

#include "package.h"
// [[Rcpp::depends(nloptr)]]

void model_sensitivity_and_specificity(
  Eigen::VectorXf&, Eigen::VectorXf&, Eigen::VectorXf&,
  double, double, double&, double&);

void formula_sensitivity_and_specificity(
  double, double, double, double, double, double, double&, double&);

void maxlik_sensitivity_and_specificity(
  double, double, double, double, double, double, double&, double&);

double nll(unsigned n, const double*, double*, void*);

/* typedefs */
typedef struct {
    mpfr_t tp, fp, fn, tn, t1, t2, t3, t4;
    double dse, dsp;
} nll_data;

#endif
