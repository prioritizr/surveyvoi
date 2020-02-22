context("rcpp_pmedian_constraint_matrix")

r_pmedian_constraint_matrix <- function(x, costs) {
  # initialization
  n <- nrow(x)
  col_names <- c(paste0("X_", seq_len(n)),
                 paste0("Y_", apply(expand.grid(n1 = seq_len(n),
                                                n2 = seq_len(n)), 1,
                                    function(x) paste(rev(x),
                                                      collapse = "_"))))
  # create constraints
  A <- Matrix::sparseMatrix(i = 1, j = 1, x = 0,
                            dims = c(1 + n + length(x), length(col_names)),
                            dimnames = list(NULL, col_names))
  A <- Matrix::drop0(A)
  # ensure that the selected sites is less than the budget
  A[1, seq_len(n)] <- costs
  # ensure each point is assigned to a single selected point
  for (i in seq_len(n)) {
    A[i + 1, n + ((i - 1) * n) + seq_len(n)] <- 1
  }
  # ensure that points can only be assigned to selected points
  counter <- n + 1
  for (i in seq_len(n)) {
    A[counter + seq_len(n), i] <- -1
    pos <- matrix(c(counter + seq_len(n), n + ((seq_len(n) - 1) * n) + i),
                  ncol = 2)
    A[pos] <- 1
    counter <- counter + n
  }
  # return result
  A
}

test_that("expected result", {
  # data
  set.seed(500)
  n <- 3
  pts <- matrix(runif(n * 2), ncol = 2)
  x <- as.matrix(dist(pts))
  costs <- runif(n)
  # computation
  r1 <- rcpp_pmedian_constraint_matrix(x, costs)
  r1 <- Matrix::sparseMatrix(i = r1$i, j = r1$j, x = r1$x, index1 = FALSE)
  r2 <- r_pmedian_constraint_matrix(x, costs)
  # tests
  expect_true(all(r1 == r2))
})
