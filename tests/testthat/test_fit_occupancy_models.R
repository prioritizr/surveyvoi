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
  tuning_parameters <- list(list(nrounds = 10, objective = "binary:logistic"))
  pu_not_surveyed <-
    rowSums(is.na(as.matrix(y[, paste0("f", seq_len(n_f))]))) > 0
  # fit models
  r <- fit_occupancy_models(x, paste0("f", seq_len(n_f)),
                            paste0("e", seq_len(n_f)),
                            parameters = tuning_parameters)
  # tests
  expect_is(r, "list")
  expect_is(r$predictions, "tbl_df")
  expect_is(r$performance, "tbl_df")
  expect_equal(nrow(r$performance), n_f)
  expect_is(r$performance$name, "character")
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
  expect_equal(r$performance$name, paste0("f", seq_len(n_f)))
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
  expect_equal(nrow(r$predictions), n_pu)
  expect_equal(ncol(r$predictions), n_f)
  expect_lte(max(as.matrix(r$predictions)), 1)
  expect_gte(min(as.matrix(r$predictions)), 0)
  expect_gte(mean(round(as.matrix(r$predictions)) ==
                  round(as.matrix(y[, paste0("p", seq_len(n_f))]))), 0.8)
})

test_that("multiple species", {
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  n_pu <- 10000
  n_f <- 3
  n_vars <- 5
  x <- simulate_site_data(n_pu, n_f, 0.5, n_env_vars = n_vars)
  y <- sf::st_drop_geometry(x)
  tuning_parameters <-
    list(list(eta = 0.3, nrounds = 10,
              objective = "binary:logistic"))[rep(1, n_f)]
  pu_not_surveyed <-
    rowSums(is.na(as.matrix(y[, paste0("f", seq_len(n_f))]))) > 0
  # fit models
  r <- fit_occupancy_models(x, paste0("f", seq_len(n_f)),
                            paste0("e", seq_len(n_f)),
                            parameters = tuning_parameters)
  # tests
  expect_is(r, "list")
  expect_is(r$predictions, "tbl_df")
  expect_is(r$performance, "tbl_df")
  expect_equal(nrow(r$performance), n_f)
  expect_is(r$performance$name, "character")
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
  expect_equal(r$performance$name, paste0("f", seq_len(n_f)))
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
  expect_equal(nrow(r$predictions), n_pu)
  expect_equal(ncol(r$predictions), n_f)
  expect_lte(max(as.matrix(r$predictions)), 1)
  expect_gte(min(as.matrix(r$predictions)), 0)
  expect_gte(mean(round(as.matrix(r$predictions)) ==
                  round(as.matrix(y[, paste0("p", seq_len(n_f))]))), 0.8)
})

test_that("species with few presences", {
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  n_pu <- 10000
  n_f <- 1
  n_vars <- 5
  x <- simulate_site_data(n_pu, n_f, 0.5, n_env_vars = n_vars)
  y <- sf::st_drop_geometry(x)
  tuning_parameters <-
    list(list(eta = 0.3, nrounds = 10,
              objective = "binary:logistic"))[rep(1, n_f)]
  pu_not_surveyed <-
    rowSums(is.na(as.matrix(y[, paste0("f", seq_len(n_f))]))) > 0
  # manually encoded one presence
  x$f1[!is.na(x$f1)] <- 0
  x$f1[sample(which(!is.na(x$f1)), 1)] <- 1
  # fit models
  r <- fit_occupancy_models(x, paste0("f", seq_len(n_f)),
                            paste0("e", seq_len(n_f)),
                            parameters = tuning_parameters)
  # tests
  expect_is(r, "list")
  expect_is(r$predictions, "tbl_df")
  expect_is(r$performance, "tbl_df")
  expect_equal(nrow(r$performance), n_f)
  expect_is(r$performance$name, "character")
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
  expect_equal(r$performance$name, paste0("f", seq_len(n_f)))
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

test_that("species with few absences", {
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  n_pu <- 10000
  n_f <- 1
  n_vars <- 5
  x <- simulate_site_data(n_pu, n_f, 0.5, n_env_vars = n_vars)
  y <- sf::st_drop_geometry(x)
  tuning_parameters <-
    list(list(eta = 0.3, nrounds = 10,
              objective = "binary:logistic"))[rep(1, n_f)]
  pu_not_surveyed <-
    rowSums(is.na(as.matrix(y[, paste0("f", seq_len(n_f))]))) > 0
  # manually encoded one presence
  x$f1[!is.na(x$f1)] <- 0
  x$f1[sample(which(!is.na(x$f1)), 1)] <- 1
  # fit models
  r <- fit_occupancy_models(x, paste0("f", seq_len(n_f)),
                            paste0("e", seq_len(n_f)),
                            parameters = tuning_parameters)
  # tests
  expect_is(r, "list")
  expect_is(r$predictions, "tbl_df")
  expect_is(r$performance, "tbl_df")
  expect_equal(nrow(r$performance), n_f)
  expect_is(r$performance$name, "character")
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
  expect_equal(r$performance$name, paste0("f", seq_len(n_f)))
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
  expect_equal(nrow(r$predictions), n_pu)
  expect_equal(ncol(r$predictions), n_f)
})
