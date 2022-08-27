# load packages
library(testthat)
library(surveyvoi)

# enable parallel testing
Sys.unsetenv("R_TESTS")

# determine reporter
if (identical(Sys.getenv("CI"), "true")) {
  reporter = "progress"
} else {
  reporter = testthat::check_reporter()
}

# check if on Fedora
os_name <- utils::sessionInfo()$running
is_Fedora <- TRUE
if (
  is.character(os_name) &&
  identical(length(os_name), 1L) &&
  all(!is.na(os_name))
) {
  is_Fedora <- any(grepl("fedora", os_name, ignore.case = TRUE, fixed = TRUE))
}

# run tests (but not on Fedora systems)
if (isTRUE(is_Fedora)) {
  message("skipping tests on Fedora system")
} else {
  test_check("surveyvoi", reporter = reporter)
}
