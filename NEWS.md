# surveyvoi 1.0.3.5

- Autoconf is used for installation on Linux and macOS operating systems.

# surveyvoi 1.0.3.4

- Update `prior_probability_matrix()` to compute prior probabilities when
  no existing survey data are available.

# surveyvoi 1.0.3.3

- Fix issue with `fit_xgb_occupancy_models()` using more than specified number
  of threads for parallel processing.
- Ensure that PSOCK and FORK clusters used for parallel processing are
  terminated correctly, even when processing is interrupted.

# surveyvoi 1.0.3.2

- Fix compatibility issues with updates to the xgboost package (version 1.5.0).
- Fix parallel processing tests given updates to the testthat package
  (version 3.1.2).
- Fix tests for environmental and geographic survey schemes given updates to
  the gurobi package (version 9.5.0).

# surveyvoi 1.0.3.1

- GMP dependencies on Windows systems are now handled using RWinLib
  (see https://github.com/rwinlib/gmp).
- Package configuration now reports compilation variables
  (i.e. PKG_CPPFLAGS and PKG_LIBS variables).
- The Free Software Foundation is no longer listed as a contributor and
  copyright holder because GMP source files are no longer distributed
  with the package (because GMP dependencies are obtained via RWinLib).
- Configuration variables can now (optionally) be used to specify location of
  GMP and MPFR dependencies for package installation (i.e. GMP_INCLUDE_DIR,
  GMP_LIB_DIR, MPFR_INCLUDE_DIR and MPFR_LIB_DIR). Although the package
  configuration routine attempts to deduce these variables automatically,
  the variables can be used if installation with default settings fails.
  For example, the variables can be set using the following system command:
  ```
  R CMD INSTALL --configure-vars='GMP_INCLUDE_DIR=... GMP_LIB_DIR=... MPFR_INCLUDE_DIR=... MPFR_LIB_DIR=...'
  ```

# surveyvoi 1.0.3

- CRAN release.
- Fix issue with missing gmpxx file.
- Remove unused dependencies.

# surveyvoi 1.0.2

- CRAN release.
- Update CRAN comments.
- Update examples.
- Improve documentation for functions that depend on external software.

# surveyvoi 1.0.1

-  Fix typos in documentation (#38).
-  Add details for tuning xgboost models to documentation (#39).

# surveyvoi 1.0.0

- Refactor for official release.
- Add support for generating surveys with the Rsymphony package.

# surveyvoi 0.0.76

- All functions appear to work.

# surveyvoi 0.0.1

- Initial commit.
