# load packages
library(testthat)
library(surveyvoi)

# enable parallel testing
Sys.unsetenv("R_TESTS")

# run tests
test_check("surveyvoi")
