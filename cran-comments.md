Dear CRAN volunteers,

Thank you very much for reviewing this submission. This submission aims to fix the issues causing the package to fail CRAN checks. Specifically, it fixes the compiler warnings during package checks (e.g., on Debian-clang flavor), and address unit test errors (i.e., on Fedora flavors). I have also taken this opportunity to update the package to be compatible with the upcoming Matrix package (version >= 1.4-2), whilst maintaining backwards compatibility.

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

The package has system requirements. Some of these requirements are mandatory -- and are required for successful installation -- and others are optional. To ensure that all mandatory requirements are available on CRAN systems, I have checked the software installed on CRAN's various systems. Briefly, this information was obtained from the [Windows](https://github.com/r-windows/rtools-packages) and [macOS](https://github.com/R-macos/recipes) toolchains; the [rcheckserver Debian meta-package](https://statmath.wu.ac.at/AASC/debian/dists/stable/main/binary-amd64/Packages), and correspondence with CRAN volunteers.

Below, I have provided a markdown table detailing whether each system requirement is optional or mandatory under various operating systems. I have also included details and URLs verifying the availability of mandatory system requirements on CRAN's systems.

| Software | Debian | Fedora | Windows | macOS |
|:--------|:---------:|:--------:|:------:|:------:|
| JAGS (>= 4.3.0) | Optional | Optional |Optional | Optional |
| fftw3 (>= 3.3) | Mandatory (1) | Mandatory (1) | Mandatory (1) | Mandatory (1) |
| gmp (>= 6.2.1) | Mandatory (2) | Mandatory (3) | [Mandatory](https://statmath.wu.ac.at/AASC/debian/dists/stable/main/binary-amd64/Packages) | [Mandatory](https://github.com/r-windows/rtools-packages/blob/master/mingw-w64-gmp/PKGBUILD) | [Mandatory](https://github.com/r-windows/rtools-packages/blob/master/mingw-w64-gmp/PKGBUILD) | [Mandatory](https://github.com/R-macos/recipes/blob/master/recipes/gmp) |
| gmpxx (>= 6.2.1) | [Mandatory (2)](https://statmath.wu.ac.at/AASC/debian/dists/stable/main/binary-amd64/Packages) | Mandatory (3) | [Mandatory](https://github.com/r-windows/rtools-packages/blob/master/mingw-w64-gmp/PKGBUILD) | [Mandatory](https://github.com/R-macos/recipes/blob/master/recipes/gmp) |
| mpfr (>= 4.1.0) | [Mandatory (4)](https://statmath.wu.ac.at/AASC/debian/dists/stable/main/binary-amd64/Packages) | Mandatory (4) | [Mandatory](https://github.com/r-windows/rtools-packages/blob/master/mingw-w64-mpfr/PKGBUILD) | [Mandatory](https://github.com/R-macos/recipes/blob/master/recipes/mpfr) |
| autoconf (>= 2.69) | Optional | Optional | Optional | [Mandatory](https://github.com/R-macos/recipes/blob/master/recipes/autoconf) |
| automake (>= 1.16.5) | Optional | Optional | Optional | [Mandatory](https://github.com/R-macos/recipes/blob/master/recipes/automake) |

(1) This requirement is available on all of CRAN's servers because the PoissonBinomial R package relies on the fftw3 library as a system requirement and this package passes checks on all systems (see https://cran.r-project.org/package=PoissonBinomial).

(2) This requirement is available on CRAN's Debian server(s) because it is available via the Debian libgmp3-dev package, and this package is a dependency of the rcheckserver Debian meta-package.

(3) This requirement is available on CRAN's Fedora server(s) because Prof. Brian Ripely checked that that the gmp-devel RPM is available on each machine (i.e., https://fedora.pkgs.org/35/fedora-x86_64/gmp-devel-6.2.0-7.fc35.i686.rpm.html) (Prof Brian Ripley, personal communication, August 25, 2022).

(4) This requirement is available on all of CRAN's servers because the Rmpfr R package relies on the mpfr library as a system requirement and this package passes checks on all systems (see https://cran.r-project.org/package=Rmpfr).

# Downstream dependencies

There are no existing packages that depend on this package.

# Comments from CRAN volunteers on previous submissions

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

  **This comment was raised on a previous submission wherein the Free Software Foundation was listed as a copyright holder. Previously, the foundation was listed a copyright holder because the package contained a header file from the gmp library to help ensure that dependencies were available. Since the dependencies should now all be available on CRAN's systems, the header file has been removed and, as such, the foundation has been removed as a copyright holder.**
