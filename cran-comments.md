## R CMD check results

0 errors | 0 warnings | 3 notes

## Notes

* checking installed package size ... NOTE
  installed size is 14.4Mb
  sub-directories of 1Mb or more:
    libs  13.8Mb

    **The package makes extensive use of C++ code to reduce run time.**

* checking package dependencies ... NOTE
  Packages suggested but not available for checking: 'gurobi'

    **The package uses the gurobi R package that is distributed with Gurobi software suite (and not available on CRAN). The DESCRIPTION, README, and package documentation files contain information for installing the gurobi R package.**

* README.md:12:73 404: Not Found [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/surveyvoi)](https://CRAN.R-project.org/package=surveyvoi)

  **The README contains a link to the CRAN web page where the package will be available once the package is published on CRAN. Since the package is not yet available on CRAN, this link does not yet work. However, once the package is available on CRAN, the link will work successfully.**

## Test environments

* [Ubuntu 20.04, R-release](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3AUbuntu)
* [Ubuntu 20.04, R-devel](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3AUbuntu)
* [Mac OSX 10.15, R-release](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3A%22Mac+OSX%22)
* [Windows Server 2019, R-release](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3AWindows)
* [Windows Server 2019, R-devel](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3AWindows)
* Windows Server 2008 (x64), R-devel (win-builder)

## Downstream dependencies

There are no existing packages that depend on this package.
