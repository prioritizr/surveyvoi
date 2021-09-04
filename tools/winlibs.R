if (getRversion() < "3.3.0") {
  stop("R version too old. On Windows this package requires at least R-3.3")
}

# Download gmp-6.1.2 from rwinlib
if (!file.exists("../windows/gmp-6.1.2/include/gmpxx.hpp")) {
  download.file(
    "https://github.com/rwinlib/gmp/archive/v6.1.2.zip", "lib.zip",
    quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}
