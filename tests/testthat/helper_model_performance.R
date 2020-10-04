simulate_confusion_matrix <- function(
  n, p, reference_sensitivity, reference_specificity, test_sensitivity,
  test_specificity) {
  # assert that arguments are valid
  assert_that(is.count(n), is.number(p), p >= 0, p <= 1,
    is.number(reference_sensitivity), is.number(reference_specificity),
    reference_sensitivity >= 0, reference_sensitivity <= 1,
    reference_specificity >= 0, reference_specificity <= 1,
    is.number(test_sensitivity), is.number(test_specificity),
    test_sensitivity >= 0, test_specificity <= 1,
    test_sensitivity >= 0, test_specificity <= 1)
  # rename variables to follow Staquet (eqns 1--4)
  x <- n * p
  y <- n * (1 - p)
  sr <- reference_sensitivity
  sn <- test_sensitivity
  spr <- reference_specificity
  spn <- test_specificity
  ma <- (x * sr * sn) + (y * (1 - spn) * (1 - spr))
  mb <- (x * (1 - sn) * sr) + (y * spn * (1 - spr))
  mc <- (x * sn * (1 - sr)) + (y * spr * (1 - spn))
  md <- (x * (1 - sr) * (1 - sn)) + (y * spr * spn)
  # create confusion matrix following Staquet
  cm <- matrix(c(ma, mb, mc, md), ncol = 2,
         dimnames = list(
           c("test_presence", "test_absence"),
           c("reference_presence", "reference_absence")))
  # round values in matrix
  cm[] <- round(cm)
  # check for valid values
  assert_that(all(c(cm) >= 0), sum(cm) >= (0.9 * n), sum(cm) <= (1.1 * n))
  # return matrix
  cm
}

r_maxlik_sensitivity_and_specificity <- function(
  x, survey_sensitivity, survey_specificity) {
  # define negative log-likelihood function
  prec <- 1e+4
  nll <- function(y) {
    # calculate likelihood components
    t1 <-
      (survey_sensitivity * y[1] * y[3]) +
      ((1 - survey_specificity) * (1 - y[2]) * (1 - y[3]))
    t2 <-
      (survey_sensitivity * (1 - y[1]) * y[3]) +
      ((1 - survey_specificity) * y[2] * (1 - y[3]))
    t3 <-
      ((1 - survey_sensitivity) * y[1] * y[3]) +
      (survey_specificity * (1 - y[2]) * (1 - y[3]))
    t4 <-
      ((1 - survey_sensitivity) * (1 - y[1]) * y[3]) +
      (survey_specificity * y[2] * (1 - y[3]))
    # calculate negative log likelihood
   -as.numeric(sum(log(c(Rmpfr::mpfr(t1, prec) ^ Rmpfr::mpfr(x[1, 1], prec),
                         Rmpfr::mpfr(t2, prec) ^ Rmpfr::mpfr(x[2, 1], prec),
                         Rmpfr::mpfr(t3, prec) ^ Rmpfr::mpfr(x[1, 2], prec),
                         Rmpfr::mpfr(t4, prec) ^ Rmpfr::mpfr(x[2, 2], prec)))))
  }
  # estimate test sensitivity and specificity by minimizing negative loglik,
  # note parameters are bounded to avoid issues with log(0)
  res <- nloptr::bobyqa(c(0.9, 0.9, 0.5), nll,
                        lower = rep(1e-10, 3), upper = rep(1 - 1e-10, 3))
  # return result
  c(res$par[[1]], res$par[[2]])
}

r_formula_sensitivity_and_specificity <- function(
  x, survey_sensitivity, survey_specificity) {
  # calculate values using imperfect detection
  # formula from: https://doi.org/10.1111/j.1466-8238.2010.00605.x
  sens <-
    ((sum(x[1, ]) * survey_specificity) - x[1, 2]) /
    ((sum(x) * (survey_specificity - 1)) + sum(x[, 1]))
  spec <-
    (sum(x[2, ] * survey_sensitivity) - x[2, 1]) /
    ((sum(x) * survey_sensitivity) - sum(x[, 1]))
  c(sens, spec)
}

r_model_performance <- function(
  actual, predicted, weights, survey_sensitivity, survey_specificity,
  method = "auto") {
  # validate arguments
  assertthat::assert_that(
    is.numeric(actual), assertthat::noNA(actual),
    is.numeric(predicted), assertthat::noNA(predicted),
    identical(length(actual), length(predicted)),
    assertthat::is.number(survey_sensitivity),
    isTRUE(survey_sensitivity > 0),
    isTRUE(survey_sensitivity < 1),
    assertthat::is.number(survey_specificity),
    isTRUE(survey_specificity > 0),
    isTRUE(survey_specificity < 1),
    is.numeric(weights), assertthat::noNA(weights),
    identical(length(actual), length(weights)),
    min(actual) == 0, max(actual) == 1,
    assertthat::is.string(method))
  assertthat::assert_that(method %in% c("formula", "maxlik", "auto"))
  # normalize weights
  wts <- (weights / sum(weights)) * length(weights)
  # create confusion matrix
  cmx <- matrix(
    c(sum(wts[(actual >= 0.5) & (predicted >= 0.5)]),
      sum(wts[(actual >= 0.5) & (predicted < 0.5)]),
      sum(wts[(actual < 0.5) & (predicted >= 0.5)]),
      sum(wts[(actual < 0.5) & (predicted < 0.5)])),
    nrow = 2, ncol = 2)
  # if no predicted presences or absences then report random performance
  if ((cmx[1, 1] < 1e-3) || (cmx[2, 2] < 1e-3))
    return(c(tss = 0, sensitivity = 0.5, specificity = 0.5))
  # calculate parameter estimates
  if (method %in% c("formula", "auto")) {
    est <- r_formula_sensitivity_and_specificity(
      cmx, survey_sensitivity, survey_specificity)
  }
  if (method == "maxlik" ||
      ((method == "auto") && (any(est < 0) || any(est > 1)))) {
    est <- r_maxlik_sensitivity_and_specificity(
      cmx, survey_sensitivity, survey_specificity)
  }
  # return model performance parameters
  c(tss = sum(est) - 1, sensitivity = est[1], specificity = est[2])
}
