#include "rcpp_predict_missing_rij_data.h"

void predict_missing_rij_data(
  Eigen::MatrixXd &rij,
  MatrixXfRM &x,
  std::vector<std::size_t> &survey_features_idx,
  std::vector<std::size_t> &pu_model_prediction_idx,
  std::vector<std::vector<BoosterHandle> *> &models) {

  // initialization
  const std::size_t n_f = survey_features_idx.size();
  const std::size_t n_vars = x.cols();
  const std::size_t n_pu_model_prediction = pu_model_prediction_idx.size();
  std::size_t curr_n_folds;
  Eigen::MatrixXf predict_y;
  Eigen::VectorXf tmp_predict_y;

  // preliminary processing
  /// prepare matrix for predicting missing data using model
  MatrixXfRM predict_x(n_pu_model_prediction, n_vars);
  for (std::size_t i = 0; i < n_pu_model_prediction; ++i)
    predict_x.row(i).array() = x.row(pu_model_prediction_idx[i]).array();

  /// create xgboost matrix
  DMatrixHandle predict_x_handle;
  XGDMatrixCreateFromMat((float *) predict_x.data(), n_pu_model_prediction,
                         n_vars, -1.0f, &predict_x_handle);

  // main processing
  for (std::size_t i = 0; i < n_f; ++i) {
    /// determine number of folds for the i'th species
    curr_n_folds = models[i]->size();
    /// allocate memory for storing predictions
    predict_y.resize(n_pu_model_prediction, curr_n_folds);
    tmp_predict_y.resize(n_pu_model_prediction);
    /// generate predictions from each fold
    for (std::size_t k = 0; k < curr_n_folds; ++k) {
      predict_xgboost_model(models[i]->at(k), predict_x_handle, tmp_predict_y);
      predict_y.col(k) = tmp_predict_y.array();
    }
    /// calculate average predictions across folds
    for (std::size_t j = 0; j < n_pu_model_prediction; ++j) {
      rij(survey_features_idx[i], pu_model_prediction_idx[j]) =
        static_cast<double>(predict_y.row(j).sum()) /
        static_cast<double>(curr_n_folds);
    }
  }

  // free xgboost pointers
  XGDMatrixFree(predict_x_handle);

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
  // init
  MatrixXfRM pu_env_data = pu_env_data_raw;
  const std::size_t n_pu = rij.cols();
  const std::size_t n_f = rij.rows();
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
  model_beta_map model_beta;
  model_performance_map model_performance;
  std::vector<std::vector<BoosterHandle> *> curr_models(n_f_survey);
  Eigen::VectorXd curr_sensitivity(n_f_survey);
  Eigen::VectorXd curr_specificity(n_f_survey);

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
    rij, wij, pu_env_data,
    survey_features_idx, feature_outcome_idx,
    xgb_parameter_names, xgb_parameter_values, n_xgb_nrounds,
    xgb_train_folds2, xgb_test_folds2,
    model_beta, model_performance,
    curr_models, curr_sensitivity, curr_specificity);

  // predict missing values
  predict_missing_rij_data(
    rij, pu_env_data, survey_features_idx, pu_model_prediction_idx,
    curr_models);

  // clean-up
  std::size_t curr_n_folds;
  for (auto itr = model_beta.begin(); itr != model_beta.end(); ++itr) {
    curr_n_folds = itr->second->size();
    for (std::size_t k = 0; k < curr_n_folds; ++k)
      XGBoosterFree(itr->second->at(k));
    delete itr->second;
  }

  // exports
  return rij;
}
