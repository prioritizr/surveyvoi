#include "rcpp_predict_missing_rij_data.h"

void predict_missing_rij_data(
  Eigen::MatrixXd &pij,
  std::vector<std::size_t> &survey_features_idx,
  std::vector<mpz_class> &feature_outcome_idx,
  std::vector<std::vector<std::size_t>> &pu_model_prediction_idx,
  model_yhat_map &model_yhat) {

  // initialization
  const std::size_t n_f = survey_features_idx.size();
  std::size_t n_pu_predict;
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
    n_pu_predict = pu_model_prediction_idx[survey_features_idx[i]].size();
    Eigen::VectorXd &curr_yhat = curr_itr->second;
    for (std::size_t j = 0; j < n_pu_predict; ++j) {
      pij(survey_features_idx[i],
          pu_model_prediction_idx[survey_features_idx[i]][j]) = curr_yhat[j];
    }
  }

  // return void
  return;
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_predict_missing_rij_data(
  Eigen::MatrixXd dij,
  Eigen::MatrixXd nij,
  Eigen::MatrixXd pij,
  Eigen::MatrixXf pu_env_data_raw,
  std::vector<bool> survey_features,
  Eigen::VectorXd survey_sensitivity,
  Eigen::VectorXd survey_specificity,
  std::vector<std::string> xgb_parameter_names,
  Rcpp::CharacterMatrix xgb_parameter_values,
  std::vector<std::size_t> n_xgb_rounds,
  std::vector<std::size_t> n_xgb_early_stopping_rounds,
  Rcpp::List xgb_train_folds,
  Rcpp::List xgb_test_folds,
  Rcpp::List pu_model_prediction) {

  // init
  MatrixXfRM pu_env_data = pu_env_data_raw;
  const std::size_t n_f = survey_features.size();
  std::vector<std::size_t> survey_features_idx;
  survey_features_idx.reserve(n_f);
  for (std::size_t i = 0; i < n_f; ++i)
    if (survey_features[i])
      survey_features_idx.push_back(i);
  survey_features_idx.shrink_to_fit();
  const std::size_t n_f_survey = survey_features_idx.size();
  const std::size_t n_f_outcomes = 1;
  std::vector<mpz_class> feature_outcome_idx(n_f_survey, 0);

  // extract xgboost parameter values
  MatrixXs xgb_parameter_values2(
    xgb_parameter_values.rows(), xgb_parameter_values.cols());
  for (std::size_t i = 0; i != xgb_parameter_values2.size(); ++i)
    xgb_parameter_values2(i) = Rcpp::as<std::string>(xgb_parameter_values[i]);

  // prepare objects for modelling
  model_yhat_map model_yhat;
  model_performance_map model_performance;
  Eigen::VectorXd curr_sensitivity(n_f_survey);
  Eigen::VectorXd curr_specificity(n_f_survey);

  // format xgboost fold indices
  std::vector<std::vector<std::vector<std::size_t>>>
    xgb_train_folds2;
  std::vector<std::vector<std::vector<std::size_t>>>
    xgb_test_folds2;
  extract_k_fold_indices(xgb_train_folds, xgb_train_folds2);
  extract_k_fold_indices(xgb_test_folds, xgb_test_folds2);

  // prepare xgboost data structures for model training
  std::vector<std::vector<Eigen::VectorXf>> train_y;
  std::vector<std::vector<Eigen::VectorXf>> train_w;
  std::vector<std::vector<MatrixXfRM>> train_x;
  extract_k_fold_y_data_from_indices(
    xgb_train_folds2, survey_features_idx, train_y);
  extract_k_fold_train_w_data_from_indices(
    pij, xgb_train_folds2, survey_features_idx, train_w);
  extract_k_fold_x_data_from_indices(
    pu_env_data, xgb_train_folds2, survey_features_idx, train_x);

  // prepare xgboost data structures for model evaluation
  std::vector<std::vector<Eigen::VectorXf>> test_y;
  std::vector<std::vector<Eigen::VectorXf>> test_w;
  std::vector<std::vector<MatrixXfRM>> test_x;
  extract_k_fold_y_data_from_indices(
    xgb_test_folds2, survey_features_idx, test_y);
  extract_k_fold_test_w_data_from_indices(
    dij, nij, xgb_test_folds2, survey_features_idx, test_w);
  extract_k_fold_x_data_from_indices(
    pu_env_data, xgb_test_folds2, survey_features_idx, test_x);

  // format pu model prediction indices
  std::vector<std::vector<std::size_t>> pu_model_prediction_idx;
  extract_list_of_list_of_indices(pu_model_prediction, pu_model_prediction_idx);

  // prepare pu prediction data
  std::vector<MatrixXfRM> pu_predict_env_data(n_f_survey);
  std::size_t curr_n;
  std::size_t curr_row;
  const std::size_t n_vars = pu_env_data_raw.cols();
  for (std::size_t i = 0; i < n_f_survey; ++i) {
    /// prepare matrix
    curr_n = pu_model_prediction_idx[i].size();
    pu_predict_env_data[i].resize(curr_n, n_vars);
    /// store environmental values for feature needing predictions
    for (std::size_t j = 0; j < curr_n; ++j) {
      curr_row = pu_model_prediction_idx[survey_features_idx[i]][j];
      pu_predict_env_data[i].row(j) = pu_env_data.row(curr_row);
    }
  }

  // fit models
  fit_xgboost_models_and_assess_performance(
    survey_features_idx,
    survey_sensitivity, survey_specificity,
    feature_outcome_idx,
    xgb_parameter_names, xgb_parameter_values2,
    n_xgb_rounds, n_xgb_early_stopping_rounds,
    train_x, train_y, train_w, test_x, test_y, test_w, pu_predict_env_data,
    model_yhat, model_performance,
    curr_sensitivity, curr_specificity);

  // predict missing values
  predict_missing_rij_data(
    pij, survey_features_idx, feature_outcome_idx, pu_model_prediction_idx,
    model_yhat);

  // exports
  return pij;
}
