context("fit_hglm_occupancy_models")

test_that("single species", {
  skip_on_cran()
  skip_if_not(suppressWarnings(is_jags_installed()))
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  n_pu <- 1000
  n_f <- 1
  n_vars <- 2
  x <- simulate_site_data(n_pu, n_f, 0.5, n_env_vars = n_vars)
  f <- simulate_feature_data(n_f = 1)
  # fit models
  suppressWarnings({
    r <- fit_hglm_occupancy_models(
      x, f, paste0("f", seq_len(n_f)), paste0("n", seq_len(n_f)),
      paste0("e", seq_len(n_vars)),
      "survey_sensitivity", "survey_specificity", n_folds = rep(2, n_f),
      jags_n_samples = rep(400, n_f), jags_n_burnin = rep(100, n_f),
      jags_n_thin = rep(1, n_f), jags_n_adapt = rep(50, n_f))
  })
  # tests
  y <- sf::st_drop_geometry(x)
  expect_is(r, "list")
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
  expect_is(r$performance$max_mpsrf, "numeric")
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
  expect_true(assertthat::noNA(r$performance$max_mpsrf))
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
  skip_if_not(suppressWarnings(is_jags_installed()))
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  n_pu <- 300
  n_f <- 3
  n_vars <- 2
  x <- simulate_site_data(n_pu, n_f, 0.5, n_env_vars = n_vars)
  f <- simulate_feature_data(n_f)
  # fit models
  suppressWarnings({
    r <- fit_hglm_occupancy_models(
      x, f, paste0("f", seq_len(n_f)), paste0("n", seq_len(n_f)),
      paste0("e", seq_len(n_vars)),
      "survey_sensitivity", "survey_specificity", n_folds = rep(2, n_f),
      jags_n_samples = rep(400, n_f), jags_n_burnin = rep(100, n_f),
      jags_n_thin = rep(1, n_f), jags_n_adapt = rep(50, n_f))
  })
  # tests
  y <- sf::st_drop_geometry(x)
  expect_is(r, "list")
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
  expect_is(r$performance$max_mpsrf, "numeric")
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
  expect_true(assertthat::noNA(r$performance$max_mpsrf))
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
  skip_if_not(is_jags_installed())
  # skip if using PSOCK cluster and package not installed
  skip_if_not(requireNamespace("surveyvoi") &&
              identical(.Platform$OS.type, "unix"))
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  n_pu <- 2000
  n_f <- 4
  n_vars <- 2
  x <- simulate_site_data(n_pu, n_f, 0.5, n_env_vars = n_vars)
  f <- simulate_feature_data(n_f)
  # randomly set sites to hvae 0 surveys for certain species
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
    r <- fit_hglm_occupancy_models(
      x2, f, paste0("f", seq_len(n_f)), paste0("n", seq_len(n_f)),
      paste0("e", seq_len(n_vars)),
      "survey_sensitivity", "survey_specificity", n_folds = rep(2, n_f),
      jags_n_samples = rep(400, n_f), jags_n_burnin = rep(100, n_f),
      jags_n_thin = rep(1, n_f), jags_n_adapt = rep(50, n_f),
      n_threads = 2)
  })
  # tests
  expect_is(r, "list")
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
  expect_is(r$performance$max_mpsrf, "numeric")
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
  expect_true(assertthat::noNA(r$performance$max_mpsrf))
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
