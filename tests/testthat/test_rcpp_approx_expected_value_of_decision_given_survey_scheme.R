context("rcpp_approx_expected_value_of_decision_given_survey_scheme")

test_that("correct result", {
  # data
  RandomFields::RFoptions(seed = 700)
  set.seed(500)
  n_f <- 3
  site_data <- simulate_site_data(n_sites = 7, n_features = n_f, 0.5)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  feature_data$preweight = round(feature_data$preweight)
  feature_data$postweight = round(feature_data$postweight)
  total_budget <- sum(site_data$management_cost * 0.8)
  site_data$survey <- FALSE
  site_data$survey[which(is.na(site_data$f1))[1:2]] <- TRUE
  n_reps <- 1
  n_outcomes_per_rep <- 2
  # prepare data
  site_occ_columns <- paste0("f", seq_len(n_f))
  site_prb_columns <- paste0("p", seq_len(n_f))
  site_env_columns <- c("e1", "e2", "e3")
  rij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_occ_columns])))
  pij <- prior_probability_matrix(
    site_data, feature_data, site_occ_columns, site_prb_columns,
    "survey_sensitivity", "survey_specificity",
    "model_sensitivity", "model_specificity")
  wij <- matrix(1, ncol = ncol(pij), nrow = nrow(pij))
  ejx <- as.matrix(sf::st_drop_geometry(site_data[, site_env_columns]))
  pu_model_prediction <- lapply(seq_len(nrow(feature_data)), function(i) {
    which(!site_data$survey & is.na(rij[i, ]))
  })
  # prepare xgboost inputs
  xgb_n_folds <- rep(5, n_f)
  xgb_parameters <-
    list(list(seed = "0", scale_pos_weight = "2",
              objective = "binary:logistic"))[rep(1, n_f)]
  ## folds for training and testing models
  xgb_folds <- lapply(seq_len(nrow(feature_data)), function(i) {
    pu_train_idx <- which(site_data$survey | !is.na(rij[i, ]))
    create_folds(unname(rij[i, pu_train_idx]), xgb_n_folds[i],
                 index = pu_train_idx,
                 na.fail = FALSE)
  })
  ## set NA rij values to -1
  rij[is.na(rij)] <- -1
  # calculations
  set.seed(1)
  r1 <- rcpp_approx_expected_value_of_decision_given_survey_scheme(
    rij = rij, pij = pij, wij = wij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction,
    pu_survey_costs = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = ejx,
    xgb_parameters = lapply(xgb_parameters, append, list(seed = "1")),
    n_xgb_nrounds = rep(8, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_preweight = feature_data$preweight,
    obj_fun_postweight = feature_data$postweight,
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    optim_gap = 0,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    method_approx_outcomes = "weighted_without_replacement",
    n_approx_states = 5)
  set.seed(1)
  r2 <- r_approx_expected_value_of_decision_given_survey_scheme(
    rij = rij, pij = pij, wij = wij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction,
    pu_survey_costs = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = ejx,
    xgb_parameters = lapply(xgb_parameters, append, list(seed = "1")),
    n_xgb_nrounds = rep(8, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_preweight = feature_data$preweight,
    obj_fun_postweight = feature_data$postweight,
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    optim_gap = 0,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    n_approx_states = 5)
  # tests
  expect_equal(r1, r2)
})
