Dear CRAN volunteers,

Thank you very much for reviewing this submission. I recognize that there were many issues with the previous submission of this package (back in May 2021), and I apologize for my carelessness. I am extremely grateful to CRAN volunteers for providing advice on addressing these issues. I have verified that the package passes CRAN package checks on both the Win-Builder and macOS platforms, along with several other environments using GitHub Actions and Docker (see Test environments below). I have also aimed to verify that the system requirements for this package are available on CRAN's check servers, and provide links to relevant files/documentation that indicate the presence of these requirements on CRAN's servers (see System requirements below).

Cheers,

Jeff

# Test environments

* [Debian (testing), R-release](https://github.com/r-devel/rcheckserver) ([based on rcheckserver](https://statmath.wu.ac.at/AASC/debian/))
* [Ubuntu 20.04, R-release](https://github.com/prioritizr/surveyvoi/actions?query=workflow%3AUbuntu)
* [Ubuntu 20.04, R-devel](https://github.com/prioritizr/surveyvoi/actions?query=workflow%3AUbuntu)
* [macOS 10.15, R-release](https://github.com/prioritizr/surveyvoi/actions?query=workflow%3A%22Mac+OSX%22)
* [macOS 11.5.2 (arm64), R-release (macOS builder)](https://mac.r-project.org/macbuilder/submit.html)
* [Windows Server 2019, R-release](https://github.com/prioritizr/surveyvoi/actions?query=workflow%3AWindows)
* [Windows Server 2008 (x64), R-devel (Win-Builder)](https://win-builder.r-project.org/)

# R CMD check results

0 errors | 0 warnings | 2 notes

# Notes

* checking installed package size ... NOTE
  installed size is 14.4Mb
  sub-directories of 1Mb or more:
    libs  13.8Mb

    **The package makes extensive use of C++ code to reduce run time.**

* checking package dependencies ... NOTE
  Packages suggested but not available for checking: 'gurobi'

    **The package uses the gurobi R package that is distributed with Gurobi software suite (and not available on CRAN). The DESCRIPTION, README, and package documentation provide instructions for installing the gurobi R package.**

# System requirements

The package has system requirements. Some of these requirements are mandatory -- and are required for successful installation -- and others are optional. To help ensure that these requirements are available on CRAN systems, I have checked that they are available under the Windows and macOS toolchains (i.e., [RTools](https://cran.r-project.org/bin/windows/Rtools/rtools40.html) and [macOS recipes](https://github.com/R-macos/recipes), respectively) and the rcheckserver Debian meta-package used by CRAN's Debian server. Below, I have provided information on whether each requirement is optional or mandatory, as well as a link indicating the requirement is available on CRAN's Debian, Windows, and macOS systems.

| Software | Debian | Windows | macOS |
|:--------|:---------:|:--------:|:------:|
| JAGS (>= 4.3.0) | Optional |Optional | Optional |
| gmp (>= 6.2.1) | [Mandatory (1)](https://statmath.wu.ac.at/AASC/debian/dists/stable/main/binary-amd64/Packages) | [Mandatory](https://github.com/r-windows/rtools-packages/blob/master/mingw-w64-gmp/PKGBUILD) | [Mandatory](https://github.com/r-windows/rtools-packages/blob/master/mingw-w64-gmp/PKGBUILD) | [Mandatory](https://github.com/R-macos/recipes/blob/master/recipes/gmp) |
| gmpxx (>= 6.2.1) | [Mandatory (1)](https://statmath.wu.ac.at/AASC/debian/dists/stable/main/binary-amd64/Packages) | [Mandatory](https://github.com/r-windows/rtools-packages/blob/master/mingw-w64-gmp/PKGBUILD) | [Mandatory](https://github.com/R-macos/recipes/blob/master/recipes/gmp) |
| mpfr (>= 4.1.0) | [Mandatory (2)](https://statmath.wu.ac.at/AASC/debian/dists/stable/main/binary-amd64/Packages) | [Mandatory](https://github.com/r-windows/rtools-packages/blob/master/mingw-w64-mpfr/PKGBUILD) | [Mandatory](https://github.com/R-macos/recipes/blob/master/recipes/mpfr) |
| pkgconfig (>= 0.29.2) | Optional | [Optional](https://github.com/r-windows/rtools-packages/blob/master/mingw-w64-mpfr/PKGBUILD) | [Optional](https://github.com/R-macos/recipes/blob/master/recipes/pkgconfig) |
| autoconf (>= 2.69) | Optional | Optional | [Mandatory](https://github.com/R-macos/recipes/blob/master/recipes/autoconf) |
| automake (>= 1.16.5) | Optinal | Optional | [Mandatory](https://github.com/R-macos/recipes/blob/master/recipes/automake) |

(1) This system requirement is available on CRAN's Debian server because it is available via the Debian libgmp3-dev package, and this package is a dependency of the rcheckserver Debian meta-package.

(2) This requirement is available on CRAN's Debian server because it is available via the Debian libmpfr-dev package, and this package is a dependency of the rcheckserver Debian meta-package.

# Downstream dependencies

There are no existing packages that depend on this package.

# Previous comments from CRAN volunteers

* Correct your SystemRequirements.

  **Thanks for catching this. I have added the C++ bindings for the GMP library (i.e., gmpxx) to the SystemRequirements field, and updated the minimum versions number for gmp and gmpxx.**

* Correct your configure script (see 'Writing R Extensions' for recommendations for using autoconf to write a *reliable* script, including logging what it does).

  **Thank you for the suggestion. I have updated the configure script to use autoconf and provide logging.**

* Inform the maintainers of binary packages, as this is a new
requirement on external software. It seems the macOS ones are not built
with C++ bindings (and you could and should have checked for yourself
from the published 'recipes' at https://github.com/R-macos/recipes).

  **Thank you very much for this advice. I am sorry that I did not check these requirements prior to my previous submission. After submitting a pull request to the maintainer of the macOS binary package (https://github.com/R-macos/recipes/pull/19), the macOS libraries now include C++ bindings for GMP. I have also confirmed with the maintainer of the Windows binary packages that the C++ bindings for GMP are available via Rtools (https://github.com/rwinlib/utils/issues/1#issuecomment-912371834). Additionally, to ensure that outdated versions of GMP do not cause issues on Windows systems, the package now uses the RWinLib infrastructure to obtain recent versions of gmp and gmpxx (https://github.com/rwinlib/gmp), see `tools/winlibs.R`). Since the package passes checks on the WinBuilder and macOS builder platforms, I am hopeful that these issues have been resolved.**

* Please always add all authors, contributors and copyright holders in the Authors@R field with the appropriate roles. e.g.: Free Software Foundation, Inc. Please explain in the submission comments what you did about this issue.

  **This comment was raised on a previous submission wherein the Free Software Foundation was listed as a copyright holder. I had previously listed the foundation as a copyright holder because the package contained a header file copied from the gmp library. Since the package now uses the RWinLib infrastructure to handle gmpxx dependencies (https://github.com/rwinlib/gmp), the header file has been removed from package and the Free Software Foundation removed as a copyright holder.**
