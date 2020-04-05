#include "rcpp_conservation_benefit.h"

double conservation_benefit_amount(
  double x, double preweight, double postweight, double target) {
  bool met = x >= target;
  return (!met * (preweight * (x / target))) +
         (met * (preweight + (postweight * (x - target))));
}

double conservation_benefit_state(
  Eigen::MatrixXd &x, Eigen::VectorXd &preweight, Eigen::VectorXd &postweight,
  Eigen::VectorXd &target) {
  Eigen::VectorXd v = x.rowwise().sum();
  Eigen::VectorXd is_met = (v.array() >= target.array()).cast<double>();
  Eigen::VectorXd not_met_value =
    preweight.array() * (v.array() / target.array()).array();
  Eigen::VectorXd met_value =
    preweight.array() +
    (postweight.array() * (v.array() - target.array())).array();
  return ((1.0 - is_met.array()).array() * not_met_value.array()).sum() +
         (is_met.array() * met_value.array()).sum();
}

// [[Rcpp::export]]
double rcpp_conservation_benefit_state(
  Eigen::MatrixXd x, Eigen::VectorXd preweight, Eigen::VectorXd postweight,
  Eigen::VectorXd target) {
  return conservation_benefit_state(x, preweight, postweight, target);
}

// [[Rcpp::export]]
double rcpp_conservation_benefit_amount(
  double x, double preweight, double postweight, double target) {
  return conservation_benefit_amount(x, preweight, postweight, target);
}
