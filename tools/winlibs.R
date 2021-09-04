# Determine GMP version from command line arguments
VERSION <- commandArgs(TRUE)

# Verify R version validity
if (getRversion() < "3.3.0") {
  stop("R version too old. On Windows this package requires at least R-3.3")
}

# Download GMP from rwinlib (https://github.com/rwinlib/gmp)
test_file <- sprintf("../windows/gmp-%s/include/gmpxx.h", VERSION)
if (!file.exists(test_file)) {
  download.file(
    sprintf("https://github.com/rwinlib/gmp/archive/v%s.zip", VERSION),
    "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}
