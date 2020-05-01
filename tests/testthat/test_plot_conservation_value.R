context("plot_conservation_value")

test_that("plot = TRUE", {
  # initialize rng
  RandomFields::RFoptions(seed = 700)
  set.seed(500)
  # data
  n_f <- 4
  site_data <- simulate_site_data(n_sites = 100, n_features = n_f, 0.3)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  # make plot
  p <- plot_conservation_value(
  site_data = site_data,
  feature_data = feature_data,
  site_occupancy_columns = paste0("f", seq_len(n_f)),
  feature_preweight_column = "preweight",
  feature_postweight_column  = "postweight",
  feature_target_column = "target")
  # tests
  expect_is(p, "gg")
})

test_that("plot = FALSE", {
  # initialize rng
  RandomFields::RFoptions(seed = 700)
  set.seed(500)
  # data
  n_f <- 4
  site_data <- simulate_site_data(n_sites = 30, n_features = n_f, 0.8)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  # make plot
  p <- plot_conservation_value(
  site_data = site_data,
  feature_data = feature_data,
  site_occupancy_columns = paste0("f", seq_len(n_f)),
  feature_preweight_column = "preweight",
  feature_postweight_column  = "postweight",
  feature_target_column = "target",
  plot = FALSE)
  # tests
  expect_is(p, "tbl_df")
})
