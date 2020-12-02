context("cluster")

test_that("FORK", {
  skip_on_cran()
  skip_on_os("windows")
  # data
  x <- seq_len(3)
  f <- function(z) {
    # verify can access internal functions
    g <- rcpp_approx_expected_value_of_decision_given_survey_scheme
    # return output
    z + 0.01
  }
  # main processing
  cl <- start_cluster(2, "x", "FORK")
  r <- suppressWarnings(plyr::laply(x, f, .parallel = TRUE))
  cl <- stop_cluster(cl)
  # tests
  expect_equal(r, x + 0.01)
})

test_that("PSOCK", {
  skip_on_cran()
  # data
  x <- seq_len(3)
  f <- function(z) {
    # verify can access internal functions
    g <- rcpp_n_states
    # return output
    z + 0.01
  }
  # main processing
  cl <- start_cluster(2, c("x", "rcpp_n_states"), "PSOCK")
  r <- suppressWarnings(plyr::laply(x, f, .parallel = TRUE))
  cl <- stop_cluster(cl)
  # tests
  expect_equal(r, x + 0.01)
})
