## Test environments

* [Ubuntu 20.04, R-release](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3AUbuntu)
* [Ubuntu 20.04, R-devel](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3AUbuntu)
* [Mac OSX 10.15, R-release](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3A%22Mac+OSX%22)
* [Windows Server 2019, R-release](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3AWindows)
* Windows Server 2008 (x64), R-devel (win-builder)

## R CMD check results

0 errors | 0 warnings | 2 notes

## Notes

* checking installed package size ... NOTE
  installed size is 14.4Mb
  sub-directories of 1Mb or more:
    libs  13.8Mb

    The package makes extensive use of C++ code to reduce run time.

* checking package dependencies ... NOTE
  Packages suggested but not available for checking: 'gurobi'

    The package uses the gurobi R package that is distributed with Gurobi
    software suite (and not available on CRAN). The DESCRIPTION, README, and
    package documentation contain information for installing the gurobi R
    package.

## Comments from CRAN maintainers

* Please unwrap the examples if they are executable in < 5 sec, or replace
  \dontrun{} with \donttest{}. Examples which depend on 'Gurobi' should be left
  in \dontrun{}.

    Thank you for this suggestion. I have removed \dontrun{} from several
    examples based on this feedback. For reference, the following functions
    depend on external software (i.e. Gurobi or JAGS) and so have examples
    containing \dontrun{}: feasible_survey_scheme, optimal_survey_scheme, and
    approx_optimal_survey_scheme. Additionally, the following
    funcitons have examples that are not executable within less than 5 seconds
    and so contain in \donttest{}: fit_xgb_occupancy_models.

* Please always add all authors, contributors and copyright holders in the
  Authors@R field with the appropriate roles.
  e.g.: Free Software Foundation, Inc
  Please explain in the submission comments what you did about this issue.

    I have added all authors, contributors, and copyright holders in the
    Authors@R field with the appropriate roles. With regards to the Free
    Software Foundation, the surveyvoi package is distributed with a file (i.e.
    the src/gmp/gmpxx.h file) for which the Free Software Foundation owns the
    copyright. As such, I have listed Free Software Foundation in the Authors@R
    field with a contributor and copyright role following R package conventions
    (e.g. see conventions in the leaflet R package). This information is
    detailed in the LICENSE.note file per CRAN policies.

## Downstream dependencies

There are no existing packages that depend on this package.
