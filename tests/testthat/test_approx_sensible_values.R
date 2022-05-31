context("approx sensible values")

test_that("lower voi when most of budget spent on surveys", {
  skip_if_not_installed("RandomFields")
  skip_on_os("Windows")
  # initialize rng
  set.seed(501)
  RandomFields::RFoptions(seed = 501)
  # data
  n_f <- 2
  site_data <- simulate_site_data(n_sites = 15, n_features = n_f, 0.1)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  feature_data$target <- c(2, 2)
  total_budget <- 5
  # manually define costs
  site_data$management_cost <- 0.1
  site_data$survey_cost <- c(4.8, rep(0.01, nrow(site_data) - 1))
  # create surveys
  site_data$survey1 <- FALSE
  site_data$survey1[1] <- TRUE
  site_data$survey2 <- FALSE
  site_data$survey2[c(3, 6)] <- TRUE
  # prepare data
  site_det_columns <- c("f1", "f2")
  site_n_columns <- c("n1", "n2")
  site_prb_columns <- c("p1", "p2")
  # calculations
  r1 <- approx_evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = site_det_columns,
    site_n_surveys_columns = site_n_columns,
    site_probability_columns = site_prb_columns,
    site_survey_scheme_column = "survey1",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    n_approx_replicates = 1,
    seed = 1)
  r2 <- approx_evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = site_det_columns,
    site_n_surveys_columns = site_n_columns,
    site_probability_columns = site_prb_columns,
    site_survey_scheme_column = "survey2",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    n_approx_replicates = 1,
    seed = 1)
  # tests
  expect_true(is.finite(r1))
  expect_true(is.finite(r2))
  expect_gt(r1, 0)
  expect_gt(r2, 0)
  expect_gt(r2, r1)
})

test_that("current == optimal info, when all pu selected", {
  skip_if_not_installed("gurobi")
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(7),
      y = x,
      f1 = c(0, 1, 1, 0, 0, 0, 0),
      f2 = c(0, 1, 0, 1, 0, 0, 0),
      n1 = c(1, 1, 1, 1, 0, 0, 0),
      n2 = c(1, 1, 1, 1, 0, 0, 0),
      p1 = c(0.99, 0.99, 0.99, 0.05, 0.99, 0.99, 0.6),
      p2 = c(0.05, 0.99, 0.05, 0.99, 0.05, 0.99, 0.4),
      survey_cost = c(1, 1, 1, 1, 5, 5, 1),
      management_cost = c(10, 10, 10, 10, 10, 10, 2),
      locked_in = FALSE),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:2],
    survey = rep(TRUE, 2),
    survey_sensitivity = rep(0.95, 2),
    survey_specificity = rep(0.9, 2),
    model_sensitivity = rep(0.8, 2),
    model_specificity = rep(0.85, 2),
    target = c(4, 4))
  # prepare data
  site_det_columns <- c("f1", "f2")
  site_n_columns <- c("n1", "n2")
  site_prb_columns <- c("p1", "p2")
  pm <- t(as.matrix(sf::st_drop_geometry(site_data)[, site_prb_columns]))
  # calculate expected values
  evd_current <- evdci(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = site_det_columns,
    site_n_surveys_columns = site_n_columns,
    site_probability_columns = site_prb_columns,
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = 100,
    site_management_locked_in_column = "locked_in",
    prior_matrix = pm)
  evd_ss <- approx_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = site_det_columns,
    site_n_surveys_columns = site_n_columns,
    site_probability_columns = site_prb_columns,
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = 100,
    survey_budget = 10,
    site_management_locked_in_column = "locked_in",
    n_approx_replicates = 1,
    prior_matrix = pm)
  # tests
  expect_equal(evd_current, max(attr(evd_ss, "ev")))
})

test_that("current < optimal info, some pu selected", {
  skip_if_not_installed("gurobi")
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(7),
      y = x,
      f1 = c(0, 1, 1, 0, 1, 1, 0),
      f2 = c(0, 1, 0, 1, 1, 0, 0),
      f3 = c(1, 0, 1, 0, 1, 1, 0),
      n1 = c(1, 1, 1, 1, 0, 1, 0),
      n2 = c(1, 1, 1, 1, 0, 1, 0),
      n3 = c(1, 1, 1, 1, 0, 1, 0),
      p1 = c(0.51, 0.99, 0.99, 0.05, 0.5, 0.99, 0.5),
      p2 = c(0.51, 0.99, 0.05, 0.99, 0.5, 0.05, 0.5),
      p3 = c(0.51, 0.05, 0.99, 0.99, 0.5, 0.99, 0.5),
      survey_cost = c(1, 1, 1, 1, 5, 100000, 8),
      management_cost = c(10, 10, 10, 10, 10, 10, 10),
      locked_in = FALSE),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    survey = rep(TRUE, 3),
    survey_sensitivity = rep(0.95, 3),
    survey_specificity = rep(0.9, 3),
    model_sensitivity = rep(0.8, 3),
    model_specificity = rep(0.85, 3),
    target = c(3, 3, 3))
  total_budget <- 57
  survey_budget <- 5
  # prepare data
  site_det_columns <- c("f1", "f2", "f3")
  site_n_columns <- c("n1", "n2", "n3")
  site_prb_columns <- c("p1", "p2", "p3")
  pm <- t(as.matrix(sf::st_drop_geometry(site_data)[, site_prb_columns]))
  # calculate expected values
  evd_current <- evdci(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = site_det_columns,
    site_n_surveys_columns = site_n_columns,
    site_probability_columns = site_prb_columns,
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    site_management_locked_in_column = "locked_in",
    prior_matrix = pm)
  evd_ss <- approx_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = site_det_columns,
    site_n_surveys_columns = site_n_columns,
    site_probability_columns = site_prb_columns,
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    survey_budget = survey_budget,
    site_management_locked_in_column = "locked_in",
    n_approx_replicates = 1,
    prior_matrix = pm)
  # tests
  expect_gt(max(attr(evd_ss, "ev")), evd_current)
})

test_that("locking out planning units lowers voi", {
  skip_if_not_installed("gurobi")
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(7),
      y = x,
      f1 = c(0, 1, 1, 0, 1, 1, 0),
      f2 = c(0, 1, 0, 1, 1, 0, 0),
      f3 = c(1, 0, 1, 0, 1, 1, 0),
      n1 = c(1, 1, 1, 1, 0, 1, 0),
      n2 = c(1, 1, 1, 1, 0, 1, 0),
      n3 = c(1, 1, 1, 1, 0, 1, 0),
      p1 = c(0.51, 0.99, 0.99, 0.05, 0.5, 0.99, 0.5),
      p2 = c(0.51, 0.99, 0.05, 0.99, 0.5, 0.05, 0.5),
      p3 = c(0.51, 0.05, 0.99, 0.99, 0.5, 0.99, 0.5),
      survey_cost = c(1, 1, 1, 1, 5, 100000, 8),
      management_cost = c(10, 10, 10, 10, 10, 10, 10),
      locked_in = FALSE,
      locked_out = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    survey = rep(TRUE, 3),
    survey_sensitivity = rep(0.95, 3),
    survey_specificity = rep(0.9, 3),
    model_sensitivity = rep(0.8, 3),
    model_specificity = rep(0.85, 3),
    target = c(3, 3, 3))
  xgb_parameters <-
    list(eta = c(0.1, 0.3, 0.5),
         lambda = c(0.01, 0.1, 0.5),
         objective = "binary:logistic")
  total_budget <- 57
  survey_budget <- 5
  # prepare data
  site_det_columns <- c("f1", "f2", "f3")
  site_n_columns <- c("n1", "n2", "n3")
  site_prb_columns <- c("p1", "p2", "p3")
  pm <- t(as.matrix(sf::st_drop_geometry(site_data)[, site_prb_columns]))
  # calculate expected values
  ## approx_evdsi
  evd_opt1 <- approx_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = site_det_columns,
    site_n_surveys_columns = site_n_columns,
    site_probability_columns = site_prb_columns,
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    survey_budget = survey_budget,
    site_management_locked_in_column = "locked_in",
    n_approx_replicates = 1,
    prior_matrix = pm)
  evd_opt2 <- approx_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = site_det_columns,
    site_n_surveys_columns = site_n_columns,
    site_probability_columns = site_prb_columns,
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    survey_budget = survey_budget,
    site_management_locked_in_column = "locked_in",
    site_management_locked_out_column = "locked_out",
    n_approx_replicates = 1,
    prior_matrix = pm)
  # tests
  expect_lt(max(attr(evd_opt2, "ev")), attr(evd_opt1, "ev"))
})

test_that("approx_evdsi >= evdci when solution is fixed", {
  skip_if_not_installed("gurobi")
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(7),
      y = x,
      f1 = c(0, 1, 1, 0, 0, 0, 0),
      f2 = c(0, 1, 0, 1, 0, 0, 0),
      n1 = c(1, 1, 1, 1, 0, 0, 0),
      n2 = c(1, 1, 1, 1, 0, 0, 0),
      p1 = c(0.99, 0.99, 0.99, 0.05, 0.99, 0.99, 0.6),
      p2 = c(0.05, 0.99, 0.05, 0.99, 0.05, 0.99, 0.4),
      survey_cost = c(1, 1, 1, 1, 5, 100000, 8),
      management_cost = c(10, 10, 10, 10, 10, 10, 2),
      locked_in = c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
      locked_out = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:2],
    survey = rep(TRUE, 2),
    survey_sensitivity = rep(0.95, 2),
    survey_specificity = rep(0.9, 2),
    model_sensitivity = rep(0.8, 2),
    model_specificity = rep(0.85, 2),
    target = c(4, 4))
  # prepare data
  site_det_columns <- c("f1", "f2")
  site_n_columns <- c("n1", "n2")
  site_prb_columns <- c("p1", "p2")
  pm <- t(as.matrix(sf::st_drop_geometry(site_data)[, site_prb_columns]))
  # calculate expected values
  evd_current <- evdci(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = site_det_columns,
    site_n_surveys_columns = site_n_columns,
    site_probability_columns = site_prb_columns,
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = 100,
    site_management_locked_in_column = "locked_in",
    site_management_locked_out_column = "locked_out",
    prior_matrix = pm)
  evd_ss <- approx_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = site_det_columns,
    site_n_surveys_columns = site_n_columns,
    site_probability_columns = site_prb_columns,
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = 100,
    survey_budget = 10,
    site_management_locked_in_column = "locked_in",
    site_management_locked_out_column = "locked_out",
    n_approx_replicates = 1,
    prior_matrix = pm)
  # tests
  expect_equal(evd_current, max(attr(evd_ss, "ev")))
})
