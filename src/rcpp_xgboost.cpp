#include "rcpp_xgboost.h"

void fit_xgboost_models_and_assess_performance(
  Eigen::MatrixXd &rij,
  Eigen::MatrixXd &wij,
  MatrixXfRM &x,
  std::vector<MatrixXfRM> &predict_x,
  std::vector<std::size_t> &survey_features_idx,
  std::vector<mpz_class> &feature_outcome_idx,
  std::vector<std::vector<std::string>> &xgb_parameter_names,
  std::vector<std::vector<std::string>> &xgb_parameter_values,
  std::vector<std::size_t> &n_xgb_nrounds,
  std::vector<std::vector<std::vector<std::size_t>>> &xgb_train_folds,
  std::vector<std::vector<std::vector<std::size_t>>> &xgb_test_folds,
  model_yhat_map &model_yhat,
  model_performance_map &model_performance,
  Eigen::VectorXd &output_model_sensitivity,
  Eigen::VectorXd &output_model_specificity) {

  // initialization
  const std::size_t n_pu = rij.cols();
  const std::size_t n_f = survey_features_idx.size();
  const std::size_t n_vars = x.cols();

  // prepare looping variables
  /// class prototypes
  model_key curr_key;
  std::pair<double, double> curr_performance;
  Eigen::VectorXd curr_predictions;
  Eigen::VectorXf curr_fold_predictions;
  /// counters
  std::size_t curr_n_pu_predict;
  std::size_t curr_n_folds;
  std::size_t curr_n_k_train;
  std::size_t curr_n_k_test;
  /// training dta
  Eigen::VectorXf train_y;
  Eigen::VectorXf train_w;
  MatrixXfRM train_x;
  /// test data
  Eigen::VectorXf test_y;
  Eigen::VectorXf test_w;
  MatrixXfRM test_x;
  /// performance maetrics
  Eigen::VectorXd fold_sensitivity;
  Eigen::VectorXd fold_specificity;

  // fit models if needed
  std::size_t curr_n_xgb_parameters;
  for (std::size_t i = 0; i < n_f; ++i) {
    // initialize unordered map key
    curr_key = std::make_pair(survey_features_idx[i], feature_outcome_idx[i]);
    // check if model has already been fit, and if so then skip model fitting
    if (model_yhat.find(curr_key) != model_yhat.end())
      continue;
    // determine number of folds for the i'th species
    curr_n_folds = xgb_train_folds[i].size();
    // allocate memory for storing performance of each fold
    fold_sensitivity.resize(curr_n_folds);
    fold_specificity.resize(curr_n_folds);
    // allocate memory for i'th species predictions
    curr_n_pu_predict = predict_x[i].rows();
    curr_predictions.resize(curr_n_pu_predict);
    curr_fold_predictions.resize(curr_n_pu_predict);
    // reset predictions vector
    curr_predictions.setZero();
    for (std::size_t k = 0; k < curr_n_folds; ++k) {
      // initialize new model
      BoosterHandle curr_model;
      // determine number of training and testing observations
      // for i'th species, k'th fold
      curr_n_k_train = xgb_train_folds[i][k].size();
      curr_n_k_test = xgb_test_folds[i][k].size();
      // prepare train labels
      train_y.resize(curr_n_k_train);
      for (std::size_t j = 0; j < curr_n_k_train; ++j)
        train_y[j] = static_cast<float>(
          rij(survey_features_idx[i], xgb_train_folds[i][k][j]));
      // prepare test labels
      test_y.resize(curr_n_k_test);
      for (std::size_t j = 0; j < curr_n_k_test; ++j)
        test_y[j] = static_cast<float>(
          rij(survey_features_idx[i], xgb_test_folds[i][k][j]));
      // prepare train weights
      train_w.resize(curr_n_k_train);
      for (std::size_t j = 0; j < curr_n_k_train; ++j)
        train_w[j] = static_cast<float>(
          wij(survey_features_idx[i], xgb_train_folds[i][k][j]));
      // prepare test weights
      test_w.resize(curr_n_k_test);
      for (std::size_t j = 0; j < curr_n_k_test; ++j)
        test_w[j] = static_cast<float>(
          wij(survey_features_idx[i], xgb_test_folds[i][k][j]));
      // prepare training matrix
      train_x.resize(curr_n_k_train, n_vars);
      for (std::size_t j = 0; j < curr_n_k_train; ++j)
          train_x.row(j).array() = x.row(xgb_train_folds[i][k][j]).array();
      // prepare test matrix
      test_x.resize(curr_n_k_test, n_vars);
      for (std::size_t j = 0; j < curr_n_k_test; ++j)
          test_x.row(j).array() = x.row(xgb_test_folds[i][k][j]).array();
      // prepare xgboost data handles for model training
      DMatrixHandle train_x_handle[1];
      XGDMatrixCreateFromMat((float *) train_x.data(), curr_n_k_train,
                             n_vars, -1.0f, &train_x_handle[0]);
      XGDMatrixSetFloatInfo(train_x_handle[0], "label", train_y.data(),
                            curr_n_k_train);
      XGDMatrixSetFloatInfo(train_x_handle[0], "weight", train_w.data(),
                            curr_n_k_train);
      // create the booster and set parameters
      XGBoosterCreate(train_x_handle, 1, &curr_model);
      curr_n_xgb_parameters =
        xgb_parameter_names[survey_features_idx[i]].size();
      for (std::size_t p = 0; p < curr_n_xgb_parameters; ++p) {
        XGBoosterSetParam(
          curr_model,
          xgb_parameter_names[survey_features_idx[i]][p].c_str(),
          xgb_parameter_values[survey_features_idx[i]][p].c_str());
      }
      // train the model
      for (std::size_t iter = 0; iter < n_xgb_nrounds[survey_features_idx[i]];
           iter++)
        XGBoosterUpdateOneIter(curr_model, iter, train_x_handle[0]);
      // prepare xgboost data handles for model testing
      DMatrixHandle test_x_handle[1];
      XGDMatrixCreateFromMat((float *) test_x.data(), curr_n_k_test,
                             n_vars, -1.0f, &test_x_handle[0]);
      // store model performance values
      xgboost_model_sensitivity_and_specificity(
        test_y, test_w, test_x_handle[0], curr_model,
        fold_sensitivity[k], fold_specificity[k]);
      // prepare predictions matrix
      DMatrixHandle predict_x_handle[1];
      XGDMatrixCreateFromMat((float *) predict_x[i].data(), curr_n_pu_predict,
                             n_vars, -1.0f, &predict_x_handle[0]);
      // generate predictions
      predict_xgboost_model(curr_model, predict_x_handle[0],
                            curr_fold_predictions);
      curr_predictions.array() +=
        curr_fold_predictions.cast<double>().array();
      // clean up
      XGDMatrixFree(train_x_handle[0]);
      XGDMatrixFree(test_x_handle[0]);
      XGDMatrixFree(predict_x_handle[0]);
      XGBoosterFree(curr_model);
    }
    // calculate average model performance across the different folds
    curr_performance = std::make_pair(
      fold_sensitivity.sum() / static_cast<double>(curr_n_folds),
      fold_specificity.sum() / static_cast<double>(curr_n_folds));
    // store model performance
    model_performance[curr_key] = curr_performance;
    // sotre model predictions
    curr_predictions.array() /= static_cast<double>(curr_n_folds);
    model_yhat[curr_key] = curr_predictions;
  }

  // extract performance statistics for output
  for (std::size_t i = 0; i < n_f; ++i) {
    // create key
    curr_key = std::make_pair(survey_features_idx[i], feature_outcome_idx[i]);
    // extract performance values
    curr_performance = model_performance.find(curr_key)->second;
    output_model_sensitivity[i] = curr_performance.first;
    output_model_specificity[i] = curr_performance.second;
  }

  // return void;
  return;
}

void predict_xgboost_model(BoosterHandle &model, DMatrixHandle &x_handle,
                           Eigen::VectorXf &out) {
  bst_ulong out_len;
  const float *raw_out;
  XGBoosterPredict(model, x_handle, 0, 0, 0, &out_len, &raw_out);
  for (unsigned int i = 0; i < out_len; ++i)
    out(i) = raw_out[i];
  return;
}

void xgboost_model_sensitivity_and_specificity(
  Eigen::VectorXf &y, Eigen::VectorXf &w, DMatrixHandle &x_handle,
  BoosterHandle &model, double &sens, double &spec) {
  // initialization
  Eigen::VectorXf yhat(y.size());

  // generate model predictions
  predict_xgboost_model(model, x_handle, yhat);

  // generate confusion table
  double total_positive = static_cast<double>((y.array() * w.array()).sum());
  double total_negative = static_cast<double>(w.sum()) - total_positive;
  double true_positive = static_cast<double>(
    (w.array() *
    ((yhat.array() >= 0.5) && (y.array() >= 0.5)).cast<float>()).sum()
  );
  double true_negative = static_cast<double>(
    (w.array() *
     ((yhat.array() < 0.5) && (y.array() < 0.5)).cast<float>()).sum()
  );

  // calculate model sensitivity
  if (std::abs(total_positive) > 1.0e-10) {
    // if there are actual positives then calculate correclty
    // sens = true_positive / (true_positive + false_negative);
    sens = true_positive / total_positive;
  } else {
    // otherwise, if there are no positives then report if the predictions
    // got this right
    sens = static_cast<double>((yhat.array() >= 0.5).count() == 0);
  }

  // calculate model specificity
  if (std::abs(total_negative) > 1.0e-10) {
    // if there are actual negatives then calculate correclty
    // spec = true_negative / (true_negative + false_positive);
    spec = true_negative / total_negative;
  } else {
    // otherwise, if there are no negatives then report if the predictions
    // got this right
    spec = static_cast<double>((yhat.array() < 0.5).count() == 0);
  }

  // clamp values to (1e-10) and (1 - 1e-10) to avoid numerical issues
  // with probabilities that are exactly zero and one
  spec = std::max(spec, 1.0e-10);
  sens = std::max(sens, 1.0e-10);
  spec = std::min(spec, 1.0 - 1.0e-10);
  sens = std::min(sens, 1.0 - 1.0e-10);

  // return void
  return;
}

// [[Rcpp::export]]
Rcpp::List rcpp_fit_xgboost_models_and_assess_performance(
  Eigen::MatrixXd rij,
  Eigen::MatrixXd wij,
  Eigen::MatrixXf pu_env_data_raw,
  std::vector<bool> survey_features,
  Rcpp::List xgb_parameters,
  std::vector<std::size_t> n_xgb_nrounds,
  Rcpp::List xgb_train_folds,
  Rcpp::List xgb_test_folds) {

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
  model_yhat_map model_yhat;
  model_performance_map model_performance;
  Eigen::VectorXd curr_sensitivity(n_f_survey);
  Eigen::VectorXd curr_specificity(n_f_survey);

  // prepare pu prediction data
  std::vector<MatrixXfRM> pu_predict_env_data(n_f_survey, pu_env_data);

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

  // exports
   return Rcpp::List::create(
     Rcpp::Named("sens") = Rcpp::wrap(curr_sensitivity),
     Rcpp::Named("spec") = Rcpp::wrap(curr_specificity));
}

// [[Rcpp::export]]
Eigen::VectorXf rcpp_xgboost(
  Eigen::VectorXf y,
  Eigen::MatrixXf train_x_raw,
  Eigen::MatrixXf predict_x_raw,
  Rcpp::List xgb_parameters,
  std::size_t n_xgb_nrounds) {

  // init
  BoosterHandle booster_handle;
  const std::size_t nvars = train_x_raw.cols();
  const std::size_t nobs_train = train_x_raw.rows();
  const std::size_t nobs_predict = predict_x_raw.rows();
  const std::size_t n_xgb_parameters = xgb_parameters.size();
  Eigen::VectorXf predict_y(nobs_predict);

  // convert the data to row-major format
  MatrixXfRM train_x_rm = train_x_raw;
  MatrixXfRM predict_x_rm = predict_x_raw;

  // prepare training data
  DMatrixHandle train_x_handle[1];
  XGDMatrixCreateFromMat((float *) train_x_rm.data(), nobs_train,
                         nvars, -1.0f, &train_x_handle[0]);
  XGDMatrixSetFloatInfo(train_x_handle[0], "label", (float *) y.data(),
                        nobs_train);

  // prepare predict data
  DMatrixHandle predict_x_handle;
  XGDMatrixCreateFromMat((float *) predict_x_rm.data(), nobs_predict,
                         nvars, -1.0f, &predict_x_handle);

  // create booster handle
  XGBoosterCreate(train_x_handle, 1, &booster_handle);

  // set parameters
  std::vector<std::string> xgb_parameter_names(n_xgb_parameters);
  Rcpp::CharacterVector xgb_parameter_names_chr = xgb_parameters.names();
  std::vector<std::string> xgb_parameter_values(n_xgb_parameters);
  for (std::size_t i = 0; i < n_xgb_parameters; ++i) {
    xgb_parameter_names[i] = Rcpp::as<std::string>(xgb_parameter_names_chr[i]);
    xgb_parameter_values[i] = Rcpp::as<std::string>(xgb_parameters[i]);
  }
  for (std::size_t i = 0; i < n_xgb_parameters; ++i)
    XGBoosterSetParam(booster_handle, xgb_parameter_names[i].c_str(),
                      xgb_parameter_values[i].c_str());

  // train model
  for (std::size_t iter = 0; iter < n_xgb_nrounds; iter++)
    XGBoosterUpdateOneIter(booster_handle, iter, train_x_handle[0]);

  // predictions
  const float *raw_out;
  bst_ulong out_len;
  XGBoosterPredict(booster_handle, predict_x_handle, 0, 0, 0, &out_len,
                   &raw_out);
  for (std::size_t i = 0; i < nobs_predict; ++i)
    predict_y[i] = raw_out[i];

  // clean up
  XGDMatrixFree(train_x_handle[0]);
  XGDMatrixFree(predict_x_handle);
  XGBoosterFree(booster_handle);

  // return result
  return predict_y;
}
