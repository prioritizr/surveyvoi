#include "rcpp_conservation_benefit.h"

double conservation_benefit_amount(
  double x, double preweight, double postweight, double target, double total) {
  bool met = x >= target;
  return (!met * (preweight * (x / target))) +
         (met * (preweight + (postweight * ((x - target) / (total - target)))));
}

double conservation_benefit_state(
  Eigen::MatrixXd &x, Eigen::VectorXd &preweight, Eigen::VectorXd &postweight,
  Eigen::VectorXd &target, double total) {
  Eigen::VectorXd v = x.rowwise().sum();
  Eigen::VectorXd is_met = (v.array() >= target.array()).cast<double>();
  Eigen::VectorXd not_met_value =
    preweight.array() * (v.array() / target.array()).array();
  Eigen::VectorXd met_value(is_met.size());
  for (std::size_t i = 0; i < met_value.size(); ++i) {
    if (std::abs(total - target[i]) > 1.0e-15) {
      met_value[i] = (v[i] - target[i]) / (total - target[i]);
    } else {
      met_value[i] = 1.0;
    }
  }
  met_value = preweight.array() + (postweight.array() * met_value.array());
  return ((1.0 - is_met.array()).array() * not_met_value.array()).sum() +
         (is_met.array() * met_value.array()).sum();
}

// [[Rcpp::export]]
double rcpp_conservation_benefit_state(
  Eigen::MatrixXd x, Eigen::VectorXd preweight, Eigen::VectorXd postweight,
  Eigen::VectorXd target, double total) {
  return conservation_benefit_state(x, preweight, postweight, target, total);
}

// [[Rcpp::export]]
double rcpp_conservation_benefit_amount(
  double x, double preweight, double postweight, double target, double total) {
  return conservation_benefit_amount(x, preweight, postweight, target, total);
}
