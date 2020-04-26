#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_sample_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"
#include "rcpp_posterior_probability_matrix.h"
#include "rcpp_predict_missing_rij_data.h"

// [[Rcpp::export]]
Rcpp::NumericVector
  rcpp_approx_expected_value_of_decision_given_survey_scheme_n_states(
  Eigen::MatrixXd rij,
  Eigen::MatrixXd pij,
  Eigen::MatrixXd wij,
  std::vector<bool> survey_features,
  Eigen::VectorXd survey_sensitivity,
  Eigen::VectorXd survey_specificity,
  std::vector<bool> pu_survey_solution,
  Eigen::VectorXd pu_survey_status,
  Eigen::VectorXd pu_survey_costs,
  Eigen::VectorXd pu_purchase_costs,
  Eigen::VectorXd pu_purchase_locked_in,
  Eigen::MatrixXf pu_env_data,
  Rcpp::List xgb_parameters,
  Rcpp::List xgb_train_folds,
  Rcpp::List xgb_test_folds,
  std::vector<std::size_t> n_xgb_nrounds,
  Eigen::VectorXd obj_fun_preweight,
  Eigen::VectorXd obj_fun_postweight,
  Eigen::VectorXd obj_fun_target,
  std::size_t n_approx_obj_fun_points,
  double total_budget,
  double optim_gap,
  std::size_t n_approx_replicates,
  std::size_t n_approx_outcomes_per_replicate,
  std::string method_approx_outcomes) {
  Rcpp::stop("TODO");
}
