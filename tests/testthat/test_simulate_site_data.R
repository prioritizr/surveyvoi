context("simulate_site_data")

test_that("simulate_site_data", {
  # data
  d <- simulate_site_data(n_sites = 1000, n_features = 2, prop = 0.5,
                          max_number_surveys_per_site = 10)
  # tests
  ## data
  expect_is(d, "sf")
  expect_equal(nrow(d), 1000)
  ## d1 column
  expect_is(d$f1, "numeric")
  expect_true(max(d$f1) <= 1)
  expect_true(min(d$f1) >= 0)
  ## d2 column
  expect_is(d$f2, "numeric")
  expect_true(max(d$f2) <= 1)
  expect_true(min(d$f2) >= 0)
  ## n1 column
  expect_is(d$n1, "numeric")
  expect_true(max(d$n1) <= 10)
  expect_true(min(d$n1) == 0)
  ## n2 column
  expect_is(d$n2, "numeric")
  expect_true(max(d$n2) <= 10)
  expect_true(min(d$n2) == 0)
  ## p1 column
  expect_is(d$p1, "numeric")
  expect_true(all(is.finite(d$p1)))
  expect_true(all(d$p1 >= 0))
  expect_true(all(d$p1 <= 1))
  ## p2 column
  expect_is(d$p2, "numeric")
  expect_true(all(is.finite(d$p2)))
  expect_true(all(d$p2 >= 0))
  expect_true(all(d$p2 <= 1))
  ## e1 column
  expect_is(d$e1, "numeric")
  expect_true(all(is.finite(d$e1)))
  ## e2 column
  expect_is(d$e2, "numeric")
  expect_true(all(is.finite(d$e2)))
  # survey_cost column
  expect_is(d$survey_cost, "numeric")
  expect_true(all(is.finite(d$survey_cost)))
  expect_true(all(d$survey_cost >= 0))
  # management_cost column
  expect_is(d$management_cost, "numeric")
  expect_true(all(is.finite(d$management_cost)))
  expect_true(all(d$management_cost >= 0))
})
