context("evdsi")

test_that("equal weights", {
  # data
  RandomFields::RFoptions(seed = 800)
  set.seed(500)
  n_f <- 2
  site_data <- simulate_site_data(n_sites = 5, n_features = n_f, 0.5)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  total_budget <- sum(site_data$management_cost) * 0.4
  site_data$survey <- FALSE
  site_data$survey[which(is.na(site_data$f1))[1:2]] <- TRUE
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
  # prepare xgboost inputs
  xgb_n_folds <- rep(5, n_f)
  xgb_parameters <-
    list(list(objective = "binary:logistic"))[rep(1, n_f)]
  pu_predict_idx <- which(site_data$survey | !is.na(site_data$f1))
  xgb_folds <- lapply(seq_len(n_f), function(i) {
      withr::with_seed(1, {
        create_folds(unname(rij[i, pu_predict_idx]), xgb_n_folds[i],
                     index = pu_predict_idx,
                     na.fail = FALSE, seed = 1)
      })
  })
  ## set NA rij values to -1
  rij[is.na(rij)] <- -1
  # calculations
  r1 <- evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_env_vars_columns = site_env_columns,
    site_survey_scheme_column = "survey",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_parameters = lapply(xgb_parameters, append, list(nrounds = 8)),
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1)
  r2 <- r_expected_value_of_decision_given_survey_scheme(
    rij = rij, pij = pij, wij = wij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_survey_status = !is.na(site_data$f1),
    pu_survey_costs = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_env_data = ejx,
    xgb_parameters = lapply(xgb_parameters, append, list(seed = "1")),
    xgb_nrounds = rep(8, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_preweight = feature_data$preweight,
    obj_fun_postweight = feature_data$postweight,
    obj_fun_target = feature_data$target,
    n_approx_obj_fun_points = 1000,
    total_budget = total_budget,
    optim_gap = 0)
  # tests
  expect_equal(r1, r2)
})

test_that("variable weights", {
  # data
  RandomFields::RFoptions(seed = 501)
  set.seed(500)
  n_f <- 2
  site_data <- simulate_site_data(n_sites = 12, n_features = n_f, 0.5)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  total_budget <- sum(site_data$management_cost * 0.8)
  site_data$survey <- FALSE
  site_data$survey[which(is.na(site_data$f1))[1:2]] <- TRUE
  site_data$w1 <- runif(nrow(site_data)) + 1
  site_data$w2 <- runif(nrow(site_data)) + 1
  # prepare data
  site_occ_columns <- c("f1", "f2")
  site_prb_columns <- c("p1", "p2")
  site_env_columns <- c("e1", "e2", "e3")
  site_wgt_columns <- c("w1", "w2")
  rij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_occ_columns])))
  pij <- prior_probability_matrix(
    site_data, feature_data, site_occ_columns, site_prb_columns,
    "survey_sensitivity", "survey_specificity",
    "model_sensitivity", "model_specificity")
  wij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_wgt_columns])))
  ejx <- as.matrix(sf::st_drop_geometry(site_data[, site_env_columns]))
  # prepare xgboost inputs
  xgb_n_folds <- rep(5, n_f)
  xgb_parameters <-
    list(list(objective = "binary:logistic"))[rep(1, n_f)]
  pu_predict_idx <- which(site_data$survey | !is.na(site_data$f1))
  xgb_folds <- lapply(seq_len(n_f), function(i) {
      withr::with_seed(1, {
        create_folds(unname(rij[i, pu_predict_idx]), xgb_n_folds[i],
                     index = pu_predict_idx,
                     na.fail = FALSE, seed = 1)
      })
  })
  ## set NA rij values to -1
  rij[is.na(rij)] <- -1
  # calculations
  r1 <- evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_env_vars_columns = site_env_columns,
    site_weight_columns = site_wgt_columns,
    site_survey_scheme_column = "survey",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_parameters = lapply(xgb_parameters, append, list(nrounds = 8)),
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1)
  r2 <- r_expected_value_of_decision_given_survey_scheme(
    rij = rij, pij = pij, wij = wij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_survey_status = !is.na(site_data$f1),
    pu_survey_costs = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_env_data = ejx,
    xgb_parameters = lapply(xgb_parameters, append, list(seed = "1")),
    xgb_nrounds = rep(8, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_preweight = feature_data$preweight,
    obj_fun_postweight = feature_data$postweight,
    obj_fun_target = feature_data$target,
    n_approx_obj_fun_points = 1000,
    total_budget = total_budget,
    optim_gap = 0)
  # tests
  expect_equal(r1, r2)
})
