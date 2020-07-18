context("fit_occupancy_models")

test_that("single species", {
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  n_pu <- 10000
  n_f <- 1
  n_vars <- 2
  x <- simulate_site_data(n_pu, n_f, 0.5, n_env_vars = n_vars)
  y <- sf::st_drop_geometry(x)
  tuning_parameters <-
    list(eta = c(0.1, 0.5), lambda = c(0.01, 0.05),
         max_depth = c(1, 2, 3),
         colsample_bytree = c(0.1, 0.8), subsample = c(0.5, 0.8, 0.9),
         objective = "binary:logistic")
  # fit models
  suppressWarnings({
    r <- fit_occupancy_models(
      x, paste0("f", seq_len(n_f)), paste0("e", seq_len(n_vars)),
      tree_method = "hist",
      parameters = tuning_parameters, n_random_search_iterations = 20)
  })
  # tests
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
  expect_is(r$performance$train_auc_mean, "numeric")
  expect_is(r$performance$train_auc_std, "numeric")
  expect_is(r$performance$train_sensitivity_mean, "numeric")
  expect_is(r$performance$train_sensitivity_std, "numeric")
  expect_is(r$performance$train_specificity_mean, "numeric")
  expect_is(r$performance$train_specificity_std, "numeric")
  expect_is(r$performance$test_auc_mean, "numeric")
  expect_is(r$performance$test_auc_std, "numeric")
  expect_is(r$performance$test_sensitivity_mean, "numeric")
  expect_is(r$performance$test_sensitivity_std, "numeric")
  expect_is(r$performance$test_specificity_mean, "numeric")
  expect_is(r$performance$test_specificity_std, "numeric")
  expect_lte(max(r$performance$train_auc_mean), 1)
  expect_lte(max(r$performance$train_auc_std), 1)
  expect_lte(max(r$performance$train_sensitivity_mean), 1)
  expect_lte(max(r$performance$train_sensitivity_std), 1)
  expect_lte(max(r$performance$train_specificity_mean), 1)
  expect_lte(max(r$performance$train_specificity_std), 1)
  expect_lte(max(r$performance$test_auc_mean), 1)
  expect_lte(max(r$performance$test_auc_std), 1)
  expect_lte(max(r$performance$test_sensitivity_mean), 1)
  expect_lte(max(r$performance$test_sensitivity_std), 1)
  expect_lte(max(r$performance$test_specificity_mean), 1)
  expect_lte(max(r$performance$test_specificity_std), 1)
  expect_gte(min(r$performance$train_auc_mean), 0)
  expect_gte(min(r$performance$train_auc_std), 0)
  expect_gte(min(r$performance$train_sensitivity_mean), 0)
  expect_gte(min(r$performance$train_sensitivity_std), 0)
  expect_gte(min(r$performance$train_specificity_mean), 0)
  expect_gte(min(r$performance$train_specificity_std), 0)
  expect_gte(min(r$performance$test_auc_mean), 0)
  expect_gte(min(r$performance$test_auc_std), 0)
  expect_gte(min(r$performance$test_sensitivity_mean), 0)
  expect_gte(min(r$performance$test_sensitivity_std), 0)
  expect_gte(min(r$performance$test_specificity_mean), 0)
  expect_gte(min(r$performance$test_specificity_std), 0)
})

test_that("multiple species", {
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  n_pu <- 2000
  n_f <- 4
  n_vars <- 2
  x <- simulate_site_data(n_pu, n_f, 0.5, n_env_vars = n_vars)
  y <- sf::st_drop_geometry(x)
  tuning_parameters <-
    list(eta = c(0.1, 0.5), lambda = c(0.01, 0.05),
         max_depth = c(1, 2, 3),
         colsample_bytree = c(0.1, 0.8), subsample = c(0.5, 0.8, 0.9),
         objective = "binary:logistic")
  # fit models
  suppressWarnings({
    r <- fit_occupancy_models(
      x, paste0("f", seq_len(n_f)), paste0("e", seq_len(n_vars)),
      parameters = tuning_parameters, n_random_search_iterations = 5,
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
  expect_gte(mean(round(as.matrix(r$predictions)) ==
                  round(as.matrix(y[, paste0("p", seq_len(n_f))]))), 0.8)
  expect_is(r$performance, "tbl_df")
  expect_equal(nrow(r$performance), n_f)
  expect_is(r$performance$feature, "character")
  expect_is(r$performance$train_auc_mean, "numeric")
  expect_is(r$performance$train_auc_std, "numeric")
  expect_is(r$performance$train_sensitivity_mean, "numeric")
  expect_is(r$performance$train_sensitivity_std, "numeric")
  expect_is(r$performance$train_specificity_mean, "numeric")
  expect_is(r$performance$train_specificity_std, "numeric")
  expect_is(r$performance$test_auc_mean, "numeric")
  expect_is(r$performance$test_auc_std, "numeric")
  expect_is(r$performance$test_sensitivity_mean, "numeric")
  expect_is(r$performance$test_sensitivity_std, "numeric")
  expect_is(r$performance$test_specificity_mean, "numeric")
  expect_is(r$performance$test_specificity_std, "numeric")
  expect_lte(max(r$performance$train_auc_mean), 1)
  expect_lte(max(r$performance$train_auc_std), 1)
  expect_lte(max(r$performance$train_sensitivity_mean), 1)
  expect_lte(max(r$performance$train_sensitivity_std), 1)
  expect_lte(max(r$performance$train_specificity_mean), 1)
  expect_lte(max(r$performance$train_specificity_std), 1)
  expect_lte(max(r$performance$test_auc_mean), 1)
  expect_lte(max(r$performance$test_auc_std), 1)
  expect_lte(max(r$performance$test_sensitivity_mean), 1)
  expect_lte(max(r$performance$test_sensitivity_std), 1)
  expect_lte(max(r$performance$test_specificity_mean), 1)
  expect_lte(max(r$performance$test_specificity_std), 1)
  expect_gte(min(r$performance$train_auc_mean), 0)
  expect_gte(min(r$performance$train_auc_std), 0)
  expect_gte(min(r$performance$train_sensitivity_mean), 0)
  expect_gte(min(r$performance$train_sensitivity_std), 0)
  expect_gte(min(r$performance$train_specificity_mean), 0)
  expect_gte(min(r$performance$train_specificity_std), 0)
  expect_gte(min(r$performance$test_auc_mean), 0)
  expect_gte(min(r$performance$test_auc_std), 0)
  expect_gte(min(r$performance$test_sensitivity_mean), 0)
  expect_gte(min(r$performance$test_sensitivity_std), 0)
  expect_gte(min(r$performance$test_specificity_mean), 0)
  expect_gte(min(r$performance$test_specificity_std), 0)
})

test_that("multiple species (sparse)", {
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  n_pu <- 2000
  n_f <- 4
  n_vars <- 2
  x <- simulate_site_data(n_pu, n_f, 0.5, n_env_vars = n_vars)
  x2 <- x
  tuning_parameters <-
    list(eta = c(0.1, 0.5), lambda = c(0.01, 0.05),
         max_depth = c(1, 2, 3),
         colsample_bytree = c(0.1, 0.8), subsample = c(0.5, 0.8, 0.9),
         objective = "binary:logistic")
  # randomly add a missing values to each species
  for (i in paste0("f", n_f)) {
    non_na <- which(!is.na(x[[i]]))
    idx <- sample(non_na, ceiling(length(non_na) * 0.25))
    x2[[i]][idx] <- NA_real_
    expect_gt(sum(is.na(x2[[i]])), sum(is.na(x[[i]])))
  }
  # fit models
  suppressWarnings({
    r <- fit_occupancy_models(
      x2, paste0("f", seq_len(n_f)), paste0("e", seq_len(n_vars)),
      parameters = tuning_parameters, n_random_search_iterations = 5,
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
  expect_is(r$performance$train_auc_mean, "numeric")
  expect_is(r$performance$train_auc_std, "numeric")
  expect_is(r$performance$train_sensitivity_mean, "numeric")
  expect_is(r$performance$train_sensitivity_std, "numeric")
  expect_is(r$performance$train_specificity_mean, "numeric")
  expect_is(r$performance$train_specificity_std, "numeric")
  expect_is(r$performance$test_auc_mean, "numeric")
  expect_is(r$performance$test_auc_std, "numeric")
  expect_is(r$performance$test_sensitivity_mean, "numeric")
  expect_is(r$performance$test_sensitivity_std, "numeric")
  expect_is(r$performance$test_specificity_mean, "numeric")
  expect_is(r$performance$test_specificity_std, "numeric")
  expect_lte(max(r$performance$train_auc_mean), 1)
  expect_lte(max(r$performance$train_auc_std), 1)
  expect_lte(max(r$performance$train_sensitivity_mean), 1)
  expect_lte(max(r$performance$train_sensitivity_std), 1)
  expect_lte(max(r$performance$train_specificity_mean), 1)
  expect_lte(max(r$performance$train_specificity_std), 1)
  expect_lte(max(r$performance$test_auc_mean), 1)
  expect_lte(max(r$performance$test_auc_std), 1)
  expect_lte(max(r$performance$test_sensitivity_mean), 1)
  expect_lte(max(r$performance$test_sensitivity_std), 1)
  expect_lte(max(r$performance$test_specificity_mean), 1)
  expect_lte(max(r$performance$test_specificity_std), 1)
  expect_gte(min(r$performance$train_auc_mean), 0)
  expect_gte(min(r$performance$train_auc_std), 0)
  expect_gte(min(r$performance$train_sensitivity_mean), 0)
  expect_gte(min(r$performance$train_sensitivity_std), 0)
  expect_gte(min(r$performance$train_specificity_mean), 0)
  expect_gte(min(r$performance$train_specificity_std), 0)
  expect_gte(min(r$performance$test_auc_mean), 0)
  expect_gte(min(r$performance$test_auc_std), 0)
  expect_gte(min(r$performance$test_sensitivity_mean), 0)
  expect_gte(min(r$performance$test_sensitivity_std), 0)
  expect_gte(min(r$performance$test_specificity_mean), 0)
  expect_gte(min(r$performance$test_specificity_std), 0)
})
