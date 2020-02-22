#ifndef PACKAGE_H
#define PACKAGE_H

/* Set plugins */
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppNumerical)]]

/* Load header files */
// R includes package
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppNumerical.h>
#include <progress.hpp>
#include <progress_bar.hpp>
// GMP library for large numbers
#include <gmp.h>
#include <gmpxx.h>
// R internal functions
#include <Rinternals.h>
#include <R_ext/Random.h>
#include <Rmath.h>
// xgboost library
#include <dmlc/logging.h>
#include <dmlc/omp.h>
#include <xgboost/c_api.h>
#include "src/common/random.h"
// standard C++/C libraries
#include <vector>
#include <string>
#include <utility>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <stdio.h>
#include <stdarg.h>

/* Import namespaces */
using namespace Rcpp;
using namespace Eigen;
using namespace Numer;
using namespace dmlc;

/* typedefs */
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXfRM;

#endif
