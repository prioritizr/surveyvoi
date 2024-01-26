# surveyvoi 1.0.6

- Fix installation for Windows on arm64 (#50).
- Remove CXX specification in Makevars to avoid NOTEs in package checks.

# surveyvoi 1.0.5.1

- Update citation.

# surveyvoi 1.0.5

- CRAN release.
- Skip tests on Fedora systems because unable to reproduce error on CRAN.
- Fix compiler warning thrown during installation.
- Fix compatibility issues with upcoming _Matrix_ package (version >= 1.4-2) .
- Fix issues unloading package properly.
- Update citation information in README.

# surveyvoi 1.0.4

- CRAN release.

# surveyvoi 1.0.3.12

- Update CRAN comments document and remove pkg-config from system requirements.
- Fix compiler warning thrown during installation.

# surveyvoi 1.0.3.11

- Update `simulate_site_data()` function to remove _RandomFields_ package
  as a dependency.

# surveyvoi 1.0.3.10

- Fix CRAN note about _utils_ package dependency.
- Skip unit tests that require _RandomFields_ package on Windows to avoid
  spurious failures.

# surveyvoi 1.0.3.9

- Update documentation for new code repository location.
- Add remote for _RandomFields_ package to facilitate installation.

# surveyvoi 1.0.3.8

- Standardize spelling (en-US).

# surveyvoi 1.0.3.7

- Tweak package documentation.
- The _RandomFields_ package is now an optional dependency.

# surveyvoi 1.0.3.6

- Update README with system requirements for _PoissonBinomial_ package (#42).

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

- Fix compatibility issues with updated _xgboost_ package (version 1.5.0).
- Fix parallel processing tests given updates to the _testthat_ package
  (version 3.1.2).
- Fix tests for environmental and geographic survey schemes given updates to
  the _gurobi_ package (version 9.5.0).

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
- Add support for generating surveys with the _Rsymphony_ package.

# surveyvoi 0.0.76

- All functions appear to work.

# surveyvoi 0.0.1

- Initial commit.
