#include "rcpp_xgboost.h"

void fit_xgboost_models_and_assess_performance(
  Eigen::MatrixXd &rij,
  Eigen::MatrixXd &wij,
  std::vector<std::size_t> &survey_features_idx,
  std::vector<mpz_class> &feature_outcome_idx,
  std::vector<std::string> &xgb_parameter_names,
  MatrixXs &xgb_parameter_values,
  std::vector<std::size_t> &n_xgb_rounds,
  std::vector<std::size_t> &n_xgb_early_stopping_rounds,
  std::vector<std::vector<MatrixXfRM>> &train_x,
  std::vector<std::vector<Eigen::VectorXf>> &train_y,
  std::vector<std::vector<Eigen::VectorXf>> &train_w,
  std::vector<std::vector<MatrixXfRM>> &test_x,
  std::vector<std::vector<Eigen::VectorXf>> &test_y,
  std::vector<std::vector<Eigen::VectorXf>> &test_w,
  std::vector<MatrixXfRM> &predict_x,
  model_yhat_map &model_yhat,
  model_performance_map &model_performance,
  Eigen::VectorXd &output_model_sensitivity,
  Eigen::VectorXd &output_model_specificity) {

  // initialization
  const std::size_t n_pu = rij.cols();
  const std::size_t n_f = survey_features_idx.size();
  const std::size_t n_xgb_parameters = xgb_parameter_names.size();
  const std::size_t n_tuning_parameter_combs = xgb_parameter_values.rows();
  // prepare looping variables
  /// class prototypes
  model_key curr_key;
  std::pair<double, double> curr_performance;
  Eigen::VectorXd curr_predictions;
  Eigen::VectorXf curr_fold_predictions;
  /// counters
  std::size_t curr_n_pu_predict;
  std::size_t curr_n_folds;
  std::size_t curr_n_max_rounds;
  std::size_t curr_n_early_rounds;
  /// performance metrics
  double curr_score, best_score, curr_det, spw;
  std::size_t curr_no_improvement_iterations;
  std::size_t curr_best_parameter_comb;
  Eigen::VectorXd models_avg_performance;
  Eigen::VectorXd fold_sensitivity;
  Eigen::VectorXd fold_specificity;
  Eigen::MatrixXi best_ntree_limit;
  Eigen::MatrixXd models_performance;
  /// string variable
  std::string spw2;
  std::string spw_name = "scale_pos_weight";
  // Main processing
  // fit models if needed
  for (std::size_t i = 0; i < n_f; ++i) {
    // initialize unordered map key
    curr_key = std::make_pair(survey_features_idx[i], feature_outcome_idx[i]);
    // check if model has already been fit, and if so then skip model fitting
    if (model_yhat.find(curr_key) != model_yhat.end())
      continue;
    // calculate constants for i'th species
    curr_n_folds = train_y[i].size();
    curr_n_max_rounds = n_xgb_rounds[survey_features_idx[i]];
    curr_n_early_rounds = n_xgb_early_stopping_rounds[survey_features_idx[i]];
    // initialize models for i'th species
    MatrixXbh models(n_tuning_parameter_combs, curr_n_folds);
    models_performance.resize(n_tuning_parameter_combs, curr_n_folds);
    models_performance.setZero();
    best_ntree_limit.resize(n_tuning_parameter_combs, curr_n_folds);
    best_ntree_limit.setZero();
    // iterate over each fold
    for (std::size_t k = 0; k < curr_n_folds; ++k) {
      // prepare training data for booster
      DMatrixHandle train_x_handle[1];
      initialize_DMatrixHandle(
        train_y[i][k], train_w[i][k], train_x[i][k], train_x_handle[0]);
      // prepare testing data for booster
      DMatrixHandle test_x_handle[1];
      initialize_DMatrixHandle(
        test_x[i][k], test_x_handle[0]);
      // calculate spw for the k'th fold
      curr_det = train_y[i][k].sum();
      spw = (static_cast<double>(train_y[i][k].size()) - curr_det) / curr_det;
      spw2 = std::to_string(spw);
      // iterate over each combination of tunning parameters
      for (std::size_t t = 0; t < n_tuning_parameter_combs; ++t) {
        // create the booster
        XGBoosterCreate(train_x_handle, 1, &models(t, k));
        // set parameters
        XGBoosterSetParam(models(t, k), spw_name.c_str(), spw2.c_str());
        for (std::size_t p = 0; p < n_xgb_parameters; ++p) {
          XGBoosterSetParam(
            models(t, k),
            xgb_parameter_names[p].c_str(),
            xgb_parameter_values(t, p).c_str());
        }
        // train the model
        curr_no_improvement_iterations = 0;
        best_score = -10.0;
        for (std::size_t iter = 0; iter < curr_n_max_rounds; iter++) {
          // update model
          XGBoosterUpdateOneIter(models(t, k), iter, train_x_handle[0]);
          // calculate performance of model
          curr_score = xgboost_model_tss(
            test_y[i][k], test_w[i][k], test_x_handle[0], 0, models(t, k));
          // check if the updated model has not improved
          if (curr_score <= best_score) {
            ++curr_no_improvement_iterations;
          } else {
            best_ntree_limit(t, k) = iter;
            best_score = curr_score;
            curr_no_improvement_iterations = 0;
          }
          // exit early if model has not improved for a set number of iterations
          if (curr_no_improvement_iterations >= curr_n_early_rounds) {
            break;
          }
        }
        // store the model performance
        models_performance(t, k) = best_score;
        // increment ntree limit because it uses base-1 indexing for some reason
        ++best_ntree_limit(t, k);
      }
      // clean up fold training and test matrices
      XGDMatrixFree(test_x_handle[0]);
      XGDMatrixFree(train_x_handle[0]);
    }
    // determine the best model
    models_avg_performance = models_performance.rowwise().sum();
    curr_best_parameter_comb =
      std::max_element(
        models_avg_performance.data(),
       models_avg_performance.data() + models_avg_performance.size()) -
      models_avg_performance.data();
    // allocate memory for calculating and storing performance of each fold
    fold_sensitivity.resize(curr_n_folds);
    fold_specificity.resize(curr_n_folds);
    DMatrixHandle best_test_x_handle[1];
    // calculate model sensitivity and specificity for best model
    for (std::size_t k = 0; k < curr_n_folds; ++k) {
      initialize_DMatrixHandle(
        test_x[i][k], best_test_x_handle[0]);
      xgboost_model_sensitivity_and_specificity(
        test_y[i][k], test_w[i][k], best_test_x_handle[0],
        best_ntree_limit(curr_best_parameter_comb, k),
        models(curr_best_parameter_comb, k),
        fold_sensitivity[k], fold_specificity[k]);
      XGDMatrixFree(best_test_x_handle[0]);
    }
    // calculate average model performance across the different folds
    curr_performance = std::make_pair(
      fold_sensitivity.sum() / static_cast<double>(curr_n_folds),
      fold_specificity.sum() / static_cast<double>(curr_n_folds));
    model_performance[curr_key] = curr_performance;
    // allocate memory for i'th species predictions
    curr_n_pu_predict = static_cast<std::size_t>(predict_x[i].rows());
    curr_predictions.resize(curr_n_pu_predict);
    curr_fold_predictions.resize(curr_n_pu_predict);
    DMatrixHandle predict_x_handle[1];
    initialize_DMatrixHandle(predict_x[i], predict_x_handle[0]);
    // generate predictions
    curr_predictions.setZero();
    for (std::size_t k = 0; k < curr_n_folds; ++k) {
      predict_xgboost_model(
        best_ntree_limit(curr_best_parameter_comb, k),
        models(curr_best_parameter_comb, k),
        predict_x_handle[0], curr_fold_predictions);
      curr_predictions.array() += curr_fold_predictions.cast<double>().array();
    }
    curr_predictions.array() /= static_cast<double>(curr_n_folds);
    // store predictions
    model_yhat[curr_key] = curr_predictions;
    // clean up
    XGDMatrixFree(predict_x_handle[0]);
    for (std::size_t k = 0; k < curr_n_folds; ++k)
      for (std::size_t t = 0; t < n_tuning_parameter_combs; ++t)
        XGBoosterFree(models(t, k));
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

void predict_xgboost_model(
  int ntree_limit, BoosterHandle &model, DMatrixHandle &x_handle,
  Eigen::VectorXf &out) {
  bst_ulong out_len;
  const float *raw_out;
  XGBoosterPredict(model, x_handle, 0, ntree_limit, 0, &out_len, &raw_out);
  for (unsigned int i = 0; i < out_len; ++i)
    out(i) = raw_out[i];
  return;
}

double xgboost_model_tss(
  Eigen::VectorXf &y, Eigen::VectorXf &w, DMatrixHandle &x_handle,
  int ntree_limit, BoosterHandle &model) {
  // initialization
  Eigen::VectorXf yhat(y.size());
  double sens, spec;

  // calculate model sensitivity and specificty
  xgboost_model_sensitivity_and_specificity(
    y, w, x_handle, ntree_limit, model, sens, spec);

  // return result
  return sens + spec - 1.0;
}

void xgboost_model_sensitivity_and_specificity(
  Eigen::VectorXf &y, Eigen::VectorXf &w, DMatrixHandle &x_handle,
  int ntree_limit, BoosterHandle &model, double &sens, double &spec) {
  // initialization
  Eigen::VectorXf yhat(y.size());

  // generate model predictions
  predict_xgboost_model(ntree_limit, model, x_handle, yhat);

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
    // if there are actual positives then calculate correctly
    // sens = true_positive / (true_positive + false_negative);
    sens = true_positive / total_positive;
  } else {
    // otherwise, if there are no positives then report if the predictions
    // got this right
    sens = static_cast<double>((yhat.array() >= 0.5).count() == 0);
  }

  // calculate model specificity
  if (std::abs(total_negative) > 1.0e-10) {
    // if there are actual negatives then calculate correctly
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
  std::vector<std::string> xgb_parameter_names,
  Rcpp::CharacterMatrix xgb_parameter_values,
  std::vector<std::size_t> n_xgb_rounds,
  std::vector<std::size_t> n_xgb_early_stopping_rounds,
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

  // prepare pu prediction data
  std::vector<MatrixXfRM> pu_predict_env_data(n_f_survey, pu_env_data);

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
  extract_k_fold_vector_data_from_indices(
    rij, xgb_train_folds2, survey_features_idx, train_y);
  extract_k_fold_vector_data_from_indices(
    wij, xgb_train_folds2, survey_features_idx, train_w);
  extract_k_fold_matrix_data_from_indices(
    pu_env_data, xgb_train_folds2, survey_features_idx, train_x);

  // prepare xgboost data structures for model evaluation
  std::vector<std::vector<Eigen::VectorXf>> test_y;
  std::vector<std::vector<Eigen::VectorXf>> test_w;
  std::vector<std::vector<MatrixXfRM>> test_x;
  extract_k_fold_vector_data_from_indices(
    rij, xgb_test_folds2, survey_features_idx, test_y);
  extract_k_fold_vector_data_from_indices(
    wij, xgb_test_folds2, survey_features_idx, test_w);
  extract_k_fold_matrix_data_from_indices(
    pu_env_data, xgb_test_folds2, survey_features_idx, test_x);

  // fit models
  fit_xgboost_models_and_assess_performance(
    rij, wij, survey_features_idx, feature_outcome_idx,
    xgb_parameter_names, xgb_parameter_values2,
    n_xgb_rounds, n_xgb_early_stopping_rounds,
    train_x, train_y, train_w, test_x, test_y, test_w, pu_predict_env_data,
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
