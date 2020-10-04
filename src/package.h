#pragma once
#ifndef PACKAGE_H
#define PACKAGE_H

/* Set plugins */
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

/* Load header files */
// R includes package
#include <Rcpp.h>
#include <RcppEigen.h>
#include <PoissonBinomial.h>
// R internal functions
#include <Rinternals.h>
#include <R_ext/Random.h>
#include <Rmath.h>
// standard C++/C libraries
#include <vector>
#include <string>
#include <utility>
#include <cstring>
#include <cstdio>
#include <unordered_map>
#include <utility>
#include <sstream>
#include <stdio.h>
#include <stdarg.h>
// GMP library for large integers
#include <gmp.h>
#include <gmpxx.h>
// MPFR library for large and small floats (used in maximum likelihood calc)
#include <mpfr.h>

/* Import namespaces */
using namespace Rcpp;
using namespace Eigen;

/* typedefs */
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXfRM;

typedef Eigen::Matrix<std::string, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

#endif
