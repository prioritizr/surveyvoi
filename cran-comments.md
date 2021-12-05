Dear CRAN volunteers,

Thank you very much for reviewing this submission. I recognize that there were many issues with the previous submission of this package (back in May 2021), and I apologize for my carelessness. I am extremely grateful to CRAN volunteers for providing advice on addressing these issues. I have verified that the package passes CRAN checks on both the Win-Builder and macOS platforms, along with several other environments using GitHub Actions.

Cheers,

Jeff

# Test environments

* [Ubuntu 20.04, R-release](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3AUbuntu)
* [Ubuntu 20.04, R-devel](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3AUbuntu)
* [Mac OSX 10.15, R-release](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3A%22Mac+OSX%22)
* [macOS 11.5.2 (arm64), R-release (macOS builder)](https://mac.r-project.org/macbuilder/submit.html)
* [Windows Server 2019, R-release](https://github.com/jeffreyhanson/surveyvoi/actions?query=workflow%3AWindows)
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

# Downstream dependencies

There are no existing packages that depend on this package.

# Previous comments from CRAN volunteers

* Correct your SystemRequirements.

  **Thanks for catching this. I have added the C++ bindings for the GMP library (GMPXX) to the SystemRequirements field, and updated the minimum version number for GMP and GMPXX.**

* Correct your configure script (see 'Writing R Extensions' for recommendations for using autoconf to write a *reliable* script, including logging what it does).

  **Thank you for the suggestion. I have updated the configure script to use autoconf and provide logging.**

* Inform the maintainers of binary packages, as this is a new
requirement on external software. It seems the macOS ones are not built
with C++ bindings (and you could and should have checked for yourself
from the published 'recipes' at https://github.com/R-macos/recipes).

  **Thank you very much for this advice. I am sorry that I did not check these requirements prior to my previous submission. After submitting a pull request to the maintainer of the macOS binary package (https://github.com/R-macos/recipes/pull/19), the macOS libraries now include C++ bindings for GMP. I have also confirmed with the maintainer of the Windows binary packages that the C++ bindings for GMP are available via Rtools (https://github.com/rwinlib/utils/issues/1#issuecomment-912371834). Additionally, to ensure that outdated versions of GMP do not cause issues on Windows systems, the package now uses the RWinLib infrastructure to obtain recent versions of GMP and GMPXX (https://github.com/rwinlib/gmp), see `tools/winlibs.R`). Since the package passes checks on the WinBuilder and macOS builder platforms, I am hopeful that these issues have been resolved.**

* Please always add all authors, contributors and copyright holders in the Authors@R field with the appropriate roles. e.g.: Free Software Foundation, Inc. Please explain in the submission comments what you did about this issue.

  **This comment was raised on a previous submission wherein the Free Software Foundation was listed as a copyright holder. I had previously listed the foundation as a copyright holder because the package contained a header file copied from the GMPXX library. Since the package now uses the RWinLib infrastructure to handle GMPXX dependencies (https://github.com/rwinlib/gmp), the header file has been removed from package and the Free Software Foundation removed as a copyright holder.**
