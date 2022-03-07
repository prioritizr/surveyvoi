context("fit_xgb_occupancy_models")

test_that("single species", {
  skip_on_cran()
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  n_pu <- 1000
  n_f <- 1
  n_vars <- 2
  x <- simulate_site_data(n_pu, n_f, 0.5, n_env_vars = n_vars)
  f <- simulate_feature_data(n_f = 1)
  f$survey_sensitivity <- 0.9
  f$survey_specificity <- 0.999
  tuning_parameters <-
    list(eta = c(0.1, 0.5), lambda = c(0.01, 0.05),
         objective = "binary:logistic", tree_method = "auto")
  # fit models
  suppressWarnings({
    r <- fit_xgb_occupancy_models(
      x, f, paste0("f", seq_len(n_f)), paste0("n", seq_len(n_f)),
      paste0("e", seq_len(n_vars)),
      "survey_sensitivity", "survey_specificity",
      xgb_tuning_parameters = tuning_parameters)
  })
  # tests
  y <- sf::st_drop_geometry(x)
  expect_is(r, "list")
  expect_is(r$parameters, "list")
  expect_is(r$predictions, "tbl_df")
  expect_equal(nrow(r$predictions), n_pu)
  expect_is(r$predictions$f1, "numeric")
  expect_lte(max(r$predictions$f1), 1)
  expect_lte(min(r$predictions$f1), 1)
  expect_gte(mean(round(as.matrix(r$predictions)) ==
                  round(as.matrix(y[, paste0("p", seq_len(n_f))]))), 0.8)
  expect_is(r$performance, "tbl_df")
  expect_equal(nrow(r$performance), n_f)
  expect_is(r$performance$feature, "character")
  expect_is(r$performance$train_tss_mean, "numeric")
  expect_is(r$performance$train_tss_std, "numeric")
  expect_is(r$performance$train_sensitivity_mean, "numeric")
  expect_is(r$performance$train_sensitivity_std, "numeric")
  expect_is(r$performance$train_specificity_mean, "numeric")
  expect_is(r$performance$train_specificity_std, "numeric")
  expect_is(r$performance$test_tss_mean, "numeric")
  expect_is(r$performance$test_tss_std, "numeric")
  expect_is(r$performance$test_sensitivity_mean, "numeric")
  expect_is(r$performance$test_sensitivity_std, "numeric")
  expect_is(r$performance$test_specificity_mean, "numeric")
  expect_is(r$performance$test_specificity_std, "numeric")
  expect_lte(max(r$performance$train_tss_mean), 1)
  expect_lte(max(r$performance$train_tss_std), 1)
  expect_lte(max(r$performance$train_sensitivity_mean), 1)
  expect_lte(max(r$performance$train_sensitivity_std), 1)
  expect_lte(max(r$performance$train_specificity_mean), 1)
  expect_lte(max(r$performance$train_specificity_std), 1)
  expect_lte(max(r$performance$test_tss_mean), 1)
  expect_lte(max(r$performance$test_tss_std), 1)
  expect_lte(max(r$performance$test_sensitivity_mean), 1)
  expect_lte(max(r$performance$test_sensitivity_std), 1)
  expect_lte(max(r$performance$test_specificity_mean), 1)
  expect_lte(max(r$performance$test_specificity_std), 1)
  expect_gte(min(r$performance$train_tss_mean), 0)
  expect_gte(min(r$performance$train_tss_std), 0)
  expect_gte(min(r$performance$train_sensitivity_mean), 0)
  expect_gte(min(r$performance$train_sensitivity_std), 0)
  expect_gte(min(r$performance$train_specificity_mean), 0)
  expect_gte(min(r$performance$train_specificity_std), 0)
  expect_gte(min(r$performance$test_tss_mean), 0)
  expect_gte(min(r$performance$test_tss_std), 0)
  expect_gte(min(r$performance$test_sensitivity_mean), 0)
  expect_gte(min(r$performance$test_sensitivity_std), 0)
  expect_gte(min(r$performance$test_specificity_mean), 0)
  expect_gte(min(r$performance$test_specificity_std), 0)
})

test_that("multiple species", {
  skip_on_cran()
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  n_pu <- 300
  n_f <- 3
  n_vars <- 2
  x <- simulate_site_data(n_pu, n_f, 0.5, n_env_vars = n_vars)
  f <- simulate_feature_data(n_f)
  tuning_parameters <-
    list(eta = c(0.1, 0.5), lambda = c(0.01, 0.05),
         objective = "binary:logistic")
  # fit models
  suppressWarnings({
    r <- fit_xgb_occupancy_models(
      x, f, paste0("f", seq_len(n_f)), paste0("n", seq_len(n_f)),
      paste0("e", seq_len(n_vars)),
      "survey_sensitivity", "survey_specificity",
      xgb_tuning_parameters = tuning_parameters,
      xgb_early_stopping_rounds = rep(5, n_f))
  })
  # tests
  y <- sf::st_drop_geometry(x)
  expect_is(r, "list")
  expect_is(r$parameters, "list")
  expect_is(r$predictions, "tbl_df")
  expect_equal(nrow(r$predictions), n_pu)
  expect_is(r$predictions$f1, "numeric")
  expect_lte(max(r$predictions$f1), 1)
  expect_lte(min(r$predictions$f1), 1)
  expect_gte(mean(round(as.matrix(r$predictions)) ==
                  round(as.matrix(y[, paste0("p", seq_len(n_f))]))), 0.8)
  expect_is(r$performance, "tbl_df")
  expect_equal(nrow(r$performance), n_f)
  expect_is(r$performance$feature, "character")
  expect_is(r$performance$train_tss_mean, "numeric")
  expect_is(r$performance$train_tss_std, "numeric")
  expect_is(r$performance$train_sensitivity_mean, "numeric")
  expect_is(r$performance$train_sensitivity_std, "numeric")
  expect_is(r$performance$train_specificity_mean, "numeric")
  expect_is(r$performance$train_specificity_std, "numeric")
  expect_is(r$performance$test_tss_mean, "numeric")
  expect_is(r$performance$test_tss_std, "numeric")
  expect_is(r$performance$test_sensitivity_mean, "numeric")
  expect_is(r$performance$test_sensitivity_std, "numeric")
  expect_is(r$performance$test_specificity_mean, "numeric")
  expect_is(r$performance$test_specificity_std, "numeric")
  expect_lte(max(r$performance$train_tss_mean), 1)
  expect_lte(max(r$performance$train_tss_std), 1)
  expect_lte(max(r$performance$train_sensitivity_mean), 1)
  expect_lte(max(r$performance$train_sensitivity_std), 1)
  expect_lte(max(r$performance$train_specificity_mean), 1)
  expect_lte(max(r$performance$train_specificity_std), 1)
  expect_lte(max(r$performance$test_tss_mean), 1)
  expect_lte(max(r$performance$test_tss_std), 1)
  expect_lte(max(r$performance$test_sensitivity_mean), 1)
  expect_lte(max(r$performance$test_sensitivity_std), 1)
  expect_lte(max(r$performance$test_specificity_mean), 1)
  expect_lte(max(r$performance$test_specificity_std), 1)
  expect_gte(min(r$performance$train_tss_mean), 0)
  expect_gte(min(r$performance$train_tss_std), 0)
  expect_gte(min(r$performance$train_sensitivity_mean), 0)
  expect_gte(min(r$performance$train_sensitivity_std), 0)
  expect_gte(min(r$performance$train_specificity_mean), 0)
  expect_gte(min(r$performance$train_specificity_std), 0)
  expect_gte(min(r$performance$test_tss_mean), 0)
  expect_gte(min(r$performance$test_tss_std), 0)
  expect_gte(min(r$performance$test_sensitivity_mean), 0)
  expect_gte(min(r$performance$test_sensitivity_std), 0)
  expect_gte(min(r$performance$test_specificity_mean), 0)
  expect_gte(min(r$performance$test_specificity_std), 0)
})

test_that("multiple species (sparse, multiple threads)", {
  skip_on_cran()
  # skip if not on Linux or package not installed
  skip_if(!requireNamespace("surveyvoi") ||
          !identical(.Platform$OS.type, "unix"))
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  n_pu <- 2000
  n_f <- 4
  n_vars <- 2
  x <- simulate_site_data(n_pu, n_f, 0.5, n_env_vars = n_vars)
  f <- simulate_feature_data(n_f)
  tuning_parameters <-
    list(eta = c(0.1, 0.5), lambda = c(0.01, 0.05),
         objective = "binary:logistic")
  # randomly set sites to have 0 surveys for certain species
  x2 <- x
  for (i in seq_len(n_f)) {
    fn <- paste0("f", i)
    nn <- paste0("n", i)
    non_zero <- which(x[[fn]] > 0)
    idx <- sample(non_zero, ceiling(length(non_zero) * 0.25))
    x2[[fn]][idx] <- 0
    x2[[nn]][idx] <- 0
    expect_gt(sum(x2[[fn]]), 0)
  }
  # fit models
  suppressWarnings({
    r <- fit_xgb_occupancy_models(
      x2, f, paste0("f", seq_len(n_f)), paste0("n", seq_len(n_f)),
      paste0("e", seq_len(n_vars)),
      "survey_sensitivity", "survey_specificity",
      xgb_tuning_parameters = tuning_parameters,
      xgb_early_stopping_rounds = rep(5, n_f),
      n_threads = 2)
  })
  # tests
  expect_is(r, "list")
  expect_is(r$parameters, "list")
  expect_is(r$predictions, "tbl_df")
  expect_equal(nrow(r$predictions), n_pu)
  expect_is(r$predictions$f1, "numeric")
  expect_lte(max(r$predictions$f1), 1)
  expect_lte(min(r$predictions$f1), 1)
  expect_gte(
    mean(round(as.matrix(r$predictions)) ==
    round(as.matrix(sf::st_drop_geometry(x)[,
      paste0("p", seq_len(n_f))]))), 0.8)
  expect_is(r$performance, "tbl_df")
  expect_equal(nrow(r$performance), n_f)
  expect_is(r$performance$feature, "character")
  expect_is(r$performance$train_tss_mean, "numeric")
  expect_is(r$performance$train_tss_std, "numeric")
  expect_is(r$performance$train_sensitivity_mean, "numeric")
  expect_is(r$performance$train_sensitivity_std, "numeric")
  expect_is(r$performance$train_specificity_mean, "numeric")
  expect_is(r$performance$train_specificity_std, "numeric")
  expect_is(r$performance$test_tss_mean, "numeric")
  expect_is(r$performance$test_tss_std, "numeric")
  expect_is(r$performance$test_sensitivity_mean, "numeric")
  expect_is(r$performance$test_sensitivity_std, "numeric")
  expect_is(r$performance$test_specificity_mean, "numeric")
  expect_is(r$performance$test_specificity_std, "numeric")
  expect_lte(max(r$performance$train_tss_mean), 1)
  expect_lte(max(r$performance$train_tss_std), 1)
  expect_lte(max(r$performance$train_sensitivity_mean), 1)
  expect_lte(max(r$performance$train_sensitivity_std), 1)
  expect_lte(max(r$performance$train_specificity_mean), 1)
  expect_lte(max(r$performance$train_specificity_std), 1)
  expect_lte(max(r$performance$test_tss_mean), 1)
  expect_lte(max(r$performance$test_tss_std), 1)
  expect_lte(max(r$performance$test_sensitivity_mean), 1)
  expect_lte(max(r$performance$test_sensitivity_std), 1)
  expect_lte(max(r$performance$test_specificity_mean), 1)
  expect_lte(max(r$performance$test_specificity_std), 1)
  expect_gte(min(r$performance$train_tss_mean), 0)
  expect_gte(min(r$performance$train_tss_std), 0)
  expect_gte(min(r$performance$train_sensitivity_mean), 0)
  expect_gte(min(r$performance$train_sensitivity_std), 0)
  expect_gte(min(r$performance$train_specificity_mean), 0)
  expect_gte(min(r$performance$train_specificity_std), 0)
  expect_gte(min(r$performance$test_tss_mean), 0)
  expect_gte(min(r$performance$test_tss_std), 0)
  expect_gte(min(r$performance$test_sensitivity_mean), 0)
  expect_gte(min(r$performance$test_sensitivity_std), 0)
  expect_gte(min(r$performance$test_specificity_mean), 0)
  expect_gte(min(r$performance$test_specificity_std), 0)
})
