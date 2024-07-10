context("greedy_heuristic_prioritization")

test_that("expected result", {
  # data
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(3),
      y = x,
      management_cost = c(100, 100, 100),
      locked_in = c(FALSE, FALSE, FALSE),
      locked_out = c(FALSE, FALSE, FALSE),
      p1 = c(0.8, 0.8, 0.8),
      p2 = c(0.05, 0.8, 0.8)
    ),
    coords = c("x", "y")
  )
  feature_data <- tibble::tibble(
    name = letters[1:2],
    target = c(1, 1)
  )
  budget = 300
  # results
  x <- greedy_heuristic_prioritization(
    site_data,
    feature_data,
    site_probability_columns = c("p1", "p2"),
    site_management_cost_column = "management_cost",
    feature_target_column = "target",
    total_budget = budget,
    site_management_locked_in_column = "locked_in",
    site_management_locked_out_column = "locked_out"
  )
  # tests
  expect_equal(x$x, c(FALSE, TRUE, TRUE))
  expect_equal(
    x$objval,
    # sum of probability values that each species meets its target
    (1 - ((1 - 0.8) * (1 - 0.8))) * 2
  )
})

test_that("expected result (default locked in and locked out)", {
  # data
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(3),
      y = x,
      management_cost = c(100, 100, 100),
      locked_in = c(FALSE, FALSE, FALSE),
      locked_out = c(FALSE, FALSE, FALSE),
      p1 = c(0.8, 0.8, 0.8),
      p2 = c(0.05, 0.8, 0.8)
    ),
    coords = c("x", "y")
  )
  feature_data <- tibble::tibble(
    name = letters[1:2],
    target = c(1, 1)
  )
  budget = 300
  # results
  x <- greedy_heuristic_prioritization(
    site_data,
    feature_data,
    site_probability_columns = c("p1", "p2"),
    site_management_cost_column = "management_cost",
    feature_target_column = "target",
    total_budget = budget
  )
  # tests
  expect_equal(x$x, c(FALSE, TRUE, TRUE))
  expect_equal(
    x$objval,
    # sum of probability values that each species meets its target
    (1 - ((1 - 0.8) * (1 - 0.8))) * 2
  )
})
