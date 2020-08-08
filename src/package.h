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
#include <PoissonBinomial.h>
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
#include <unordered_map>
#include <utility>
#include <sstream>
#include <stdio.h>
#include <stdarg.h>

/* Import namespaces */
using namespace Rcpp;
using namespace Eigen;
using namespace Numer;
using namespace dmlc;

/* hash functions */
typedef std::pair<std::size_t, mpz_class> model_key;

struct model_key_hash {
    unsigned long operator()(const model_key& key) const {
      std::string v = std::to_string(key.first) + "-" + key.second.get_str();
      return std::hash<std::string>{}(v);
    }
};

struct model_key_equal {
    bool operator()(const model_key& t1, const model_key& t2) const {
        return (t1.first == t2.first) && (cmp(t1.second, t2.second) == 0);
    }
};

/* typedefs */
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXfRM;

typedef Eigen::Matrix<BoosterHandle, Eigen::Dynamic, Eigen::Dynamic> MatrixXbh;

typedef Eigen::Matrix<std::string, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

typedef std::unordered_map<std::pair<std::size_t, mpz_class>, std::vector<BoosterHandle>*, model_key_hash, model_key_equal> model_beta_map;

typedef std::unordered_map<std::pair<std::size_t, mpz_class>, std::vector<BoosterHandle>*, model_key_hash, model_key_equal>::iterator model_beta_iterator;

typedef std::unordered_map<std::pair<std::size_t, mpz_class>, std::pair<double, double>, model_key_hash, model_key_equal> model_performance_map;

typedef std::unordered_map<std::pair<std::size_t, mpz_class>, Eigen::VectorXd, model_key_hash, model_key_equal> model_yhat_map;

typedef std::unordered_map<std::pair<std::size_t, mpz_class>, Eigen::VectorXd, model_key_hash, model_key_equal>::iterator model_yhat_iterator;

#endif
