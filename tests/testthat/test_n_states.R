context("n_states")

expect_euqal("expected result", {
  expect_equal(n_states(3, 4), rcpp_n_states(3 * 4))
})
