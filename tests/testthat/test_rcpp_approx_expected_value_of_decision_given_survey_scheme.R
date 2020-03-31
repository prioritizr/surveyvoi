context("rcpp_approx_expected_value_of_decision_given_survey_scheme_n_states")

test_that("correct result", {
  # data
  RandomFields::RFoptions(seed = 505)
  set.seed(500)
  n_f <- 2
  site_data <- simulate_site_data(n_sites = 8, n_features = n_f, 0.5)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  total_budget <- sum(site_data$management_cost * 0.8)
  site_data$survey <- FALSE
  site_data$survey[which(is.na(site_data$f1))[1:2]] <- TRUE
  n_reps <- 4
  n_states_per_rep <- 5
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
  set.seed(1)
  r1 <- r_approx_expected_value_of_decision_given_survey_scheme_n_states(
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
    obj_fun_alpha = feature_data$alpha,
    obj_fun_gamma = feature_data$gamma,
    n_approx_obj_fun_points = 1000,
    total_budget = total_budget,
    optim_gap = 0,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep)
  set.seed(1)
  r2 <- rcpp_approx_expected_value_of_decision_given_survey_scheme_n_states(
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
    n_xgb_nrounds = rep(8, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_alpha = feature_data$alpha,
    obj_fun_gamma = feature_data$gamma,
    n_approx_obj_fun_points = 1000,
    total_budget = total_budget,
    optim_gap = 0,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_states = "weighted_without_replacement")
  # tests
  expect_equal(r1, r2)
})
