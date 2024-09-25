Dear CRAN volunteers,

Thank you for reviewing this submission. This submission contains minor updates for the _surveyvoi_ package to (i) improve algorithmic performance, (ii) fix a typo in the citation, and (iii) improve compatibility with updates to dependencies.

Cheers,

Jeff

# Test environments

* [Debian (testing), R-release](https://github.com/r-devel/rcheckserver) ([based on rcheckserver](https://statmath.wu.ac.at/AASC/debian/))
* [Fedora 33, clang, R-devel](https://github.com/prioritizr/surveyvoi/actions?query=workflow%3AFedora)
* [Ubuntu 22.04, R-release](https://github.com/prioritizr/surveyvoi/actions?query=workflow%3AUbuntu)
* [Ubuntu 22.04, R-devel](https://github.com/prioritizr/surveyvoi/actions?query=workflow%3AUbuntu)
* [macOS 10.15, R-release](https://github.com/prioritizr/surveyvoi/actions?query=workflow%3A%22Mac+OSX%22)
* [macOS 11.5.2 (arm64), R-release (macOS builder)](https://mac.r-project.org/macbuilder/submit.html)
* [Windows Server 2019, R-release](https://github.com/prioritizr/surveyvoi/actions?query=workflow%3AWindows)
* [Windows Server 2008 (x64), R-devel (Win-Builder)](https://win-builder.r-project.org/)

# R CMD check results

0 errors | 0 warnings | 3 notes

# Notes

* checking installed package size ... NOTE
  installed size is 14.4Mb
  sub-directories of 1Mb or more:
    libs  13.8Mb

    **The package makes extensive use of C++ code to reduce run time.**

* checking package dependencies ... NOTE
  Packages suggested but not available for checking: 'gurobi'

    **The package uses the gurobi R package that is distributed with Gurobi software suite (and not available on CRAN). The DESCRIPTION, README, and package documentation provide instructions for installing the gurobi R package.**

* found the following (possibly) invalid URLs:
  URL: https://support.gurobi.com/hc/en-us/articles/4534161999889-How-do-I-install-Gurobi-Optimizer
    From: man/approx_optimal_survey_scheme.Rd
          man/env_div_survey_scheme.Rd
          man/feasible_survey_schemes.Rd
          man/geo_cov_survey_scheme.Rd
          man/optimal_survey_scheme.Rd
          man/surveyvoi.Rd
          man/weighted_survey_scheme.Rd
          inst/doc/surveyvoi.html
    Status: 403
    Message: Forbidden

    **I have checked this URL and can confirm that it is valid.**

# Downstream dependencies

There are no existing packages that depend on this package.
