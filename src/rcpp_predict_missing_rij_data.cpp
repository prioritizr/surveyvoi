#include "rcpp_predict_missing_rij_data.h"

void predict_missing_rij_data(
  Eigen::MatrixXd &rij,
  std::vector<std::size_t> &survey_features_idx,
  std::vector<mpz_class> &feature_outcome_idx,
  std::vector<std::size_t> &pu_model_prediction_idx,
  model_yhat_map &model_yhat) {

  // initialization
  const std::size_t n_f = survey_features_idx.size();
  const std::size_t n_pu_model_prediction = pu_model_prediction_idx.size();
  model_key curr_key;
  model_yhat_iterator curr_itr;
  Eigen::VectorXd curr_yhat;

  // main processing
  for (std::size_t i = 0; i < n_f; ++i) {
    // initialize unordered map key
    curr_key = std::make_pair(survey_features_idx[i], feature_outcome_idx[i]);
    // check if model has already been fit, and if so then skip model fitting
    curr_itr = model_yhat.find(curr_key);
    if (curr_itr == model_yhat.end())
      Rcpp::stop("trying to make predictions from non-existant model");
    // extract predictions
    Eigen::VectorXd &curr_yhat = curr_itr->second;
    for (std::size_t j = 0; j < n_pu_model_prediction; ++j) {
      rij(survey_features_idx[i], pu_model_prediction_idx[j]) = curr_yhat[j];
    }
  }

  // return void
  return;
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_predict_missing_rij_data(
  Eigen::MatrixXd rij,
  Eigen::MatrixXd wij,
  Eigen::MatrixXf pu_env_data_raw,
  std::vector<bool> survey_features,
  std::vector<std::size_t> pu_model_prediction_idx,
  Rcpp::List xgb_parameters,
  std::vector<std::size_t> n_xgb_nrounds,
  Rcpp::List xgb_train_folds,
  Rcpp::List xgb_test_folds) {

  // init
  MatrixXfRM pu_env_data = pu_env_data_raw;
  const std::size_t n_pu = rij.cols();
  const std::size_t n_f = rij.rows();
  const std::size_t n_vars = pu_env_data_raw.cols();
  const std::size_t n_pu_model_prediction = pu_model_prediction_idx.size();
  std::vector<std::size_t> survey_features_idx;
  survey_features_idx.reserve(n_f);
  for (std::size_t i = 0; i < n_f; ++i)
    if (survey_features[i])
      survey_features_idx.push_back(i);
  survey_features_idx.shrink_to_fit();
  const std::size_t n_f_survey = survey_features_idx.size();
  const std::size_t n_f_outcomes = 1;
  std::vector<mpz_class> feature_outcome_idx(n_f_survey, 0);

  // prepare objects for modelling
  model_yhat_map model_yhat;
  model_performance_map model_performance;
  Eigen::VectorXd curr_sensitivity(n_f_survey);
  Eigen::VectorXd curr_specificity(n_f_survey);

  // subset environmental data
  MatrixXfRM pu_predict_env_data(n_pu_model_prediction, n_vars);
  for (std::size_t i = 0; i < n_pu_model_prediction; ++i)
    pu_predict_env_data.row(i).array() =
      pu_env_data.row(pu_model_prediction_idx[i]).array();

  // extract xgboost parameters
  std::vector<std::vector<std::string>> xgb_parameter_names;
  std::vector<std::vector<std::string>> xgb_parameter_values;
  extract_xgboost_parameters(xgb_parameters, xgb_parameter_names,
                             xgb_parameter_values);

  // format xgboost fold indices
  std::vector<std::vector<std::vector<std::size_t>>>
    xgb_train_folds2;
  std::vector<std::vector<std::vector<std::size_t>>>
    xgb_test_folds2;
  extract_k_fold_indices(xgb_train_folds, xgb_train_folds2);
  extract_k_fold_indices(xgb_test_folds, xgb_test_folds2);

  // fit models
  fit_xgboost_models_and_assess_performance(
    rij, wij, pu_env_data, pu_predict_env_data,
    survey_features_idx, feature_outcome_idx,
    xgb_parameter_names, xgb_parameter_values, n_xgb_nrounds,
    xgb_train_folds2, xgb_test_folds2,
    model_yhat, model_performance,
    curr_sensitivity, curr_specificity);

  // predict missing values
  predict_missing_rij_data(
    rij, survey_features_idx, feature_outcome_idx, pu_model_prediction_idx,
    model_yhat);

  // exports
  return rij;
}
