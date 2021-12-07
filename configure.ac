# Initialization
## setup configure
AC_PREREQ([2.69])
AC_INIT([surveyvoi], [see.DESCRIPTION], [https://github.com/jeffreyhanson/surveyvoi/issues])

## set GMP variables
GMP_NAME="gmp"
GMP_CONFIG_NAME="gmpxx"
GMP_DEB_NAME="libgmp3-dev"
GMP_RPM_NAME="gmp-devel"
GMP_CSW_NAME=""
GMP_BREW_NAME="gmp"
GMP_TEST_C_HEADER="<gmp.h>"
GMP_TEST_CXX_HEADER="<gmpxx.h>"

## set MFRP variables
MPFR_NAME="mpfr"
MPFR_CONFIG_NAME="mpfr"
MPFR_DEB_NAME="libmpfr-dev"
MPFR_RPM_NAME="mpfr-devel"
MPFR_CSW_NAME=""
MPFR_BREW_NAME="mpfr"
MPFR_TEST_HEADER="<mpfr.h>"

# Preliminary processing
## ensure C++ is set up as R expects
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi

## define compiler flags
if test [$RVER_MAJOR = "development"]; then
    CXX11=`"${RBIN}" CMD config CXX11`
    CXX11STD=`"${RBIN}" CMD config CXX11STD`
else
    if test [$RVER_MAJOR -lt 3] -o [$RVER_MAJOR -eq 3 -a $RVER_MINOR -lt 3]; then
        AC_MSG_ERROR([surveyvoi is not compatible with R versions before 3.3.0])
    else
        if test [$RVER_MAJOR -eq 3 -a $RVER_MINOR -eq 3]; then
            CXX11=`"${RBIN}" CMD config CXX1X`
            CXX11STD=`"${RBIN}" CMD config CXX1XSTD`
        else
            CXX11=`"${RBIN}" CMD config CXX11`
            CXX11STD=`"${RBIN}" CMD config CXX11STD`
        fi
    fi
fi

# pick all flags for testing from R
: ${CC=`"${RBIN}" CMD config CC`}
: ${CXX=${CXX11} ${CXX11STD}}
: ${CFLAGS=`"${RBIN}" CMD config CFLAGS`}
: ${CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`}
: ${CXXFLAGS=`"${RBIN}" CMD config CXXFLAGS`}
: ${LDFLAGS=`"${RBIN}" CMD config LDFLAGS`}

### provide compiler information
AC_MSG_NOTICE([CC: ${CC}])
AC_MSG_NOTICE([CXX: ${CXX}])

# Specify settings
## set up pkgconfig
AC_DEFUN([AC_PROG_PKGCONFIG], [AC_CHECK_PROG(PKGCONFIG,pkg-config,yes)])
AC_PROG_PKGCONFIG

## specify flags and library paths
if test x"${PKGCONFIG}" == x"yes"; then
  ### GMP settings
  if pkg-config --exists ${GMP_CONFIG_NAME}; then
    GMP_CFLAGS=`pkg-config --cflags --silence-errors ${GMP_CONFIG_NAME}`
    GMP_LIBS=`pkg-config --libs --silence-errors ${GMP_CONFIG_NAME}`
  else
    GMP_CFLAGS=""
    GMP_LIBS="-lgmpxx -lgmp"
  fi
  ### MPFR settings
  if pkg-config --exists ${MPFR_CONFIG_NAME}; then
    MPFR_CFLAGS=`pkg-config --cflags --silence-errors ${MPFR_CONFIG_NAME}`
    MPFR_LIBS=`pkg-config --libs --silence-errors ${MPFR_CONFIG_NAME}`
  else
    MPFR_CFLAGS=""
    MPFR_LIBS="-lmpfr -lgmp"
  fi
else
  ## manually specify defaults if pkg-config missing
  GMP_CFLAGS=""
  GMP_LIBS="-lgmpxx -lgmp"
  MPFR_CFLAGS=""
  MPFR_LIBS="-lmpfr -lgmp"
fi

## add flags for subsequent tests
PKG_CPPFLAGS="${CPPFLAGS} ${GMP_CFLAGS} ${MPFR_CFLAGS}"
PKG_LIBS="${LIBS} ${GMP_LIBS} ${MPFR_LIBS}"
CPPFLAGS="${PKG_CPPFLAGS}"
LIBS="${PKG_LIBS}"

# Verify GMP settings work
## display message
AC_MSG_CHECKING([if ${GMP_NAME} compiles])

## try compiling file with settings
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <stdlib.h>
#include ${GMP_TEST_C_HEADER}
#include ${GMP_TEST_CXX_HEADER}
#if defined(__GNU_MP_VERSION)
#define SUCCESS 1
#else
#define SUCCESS 0
#endif
int main() {
  if (SUCCESS > 0) {
    exit(0);
  } else {
    exit(1);
  }
}
]])],
[gmp_ok=yes],
[gmp_ok=no],
[gmp_ok=yes])
if test x"${gmp_ok}" == x"no"; then
    AC_MSG_RESULT([no])
else
    AC_MSG_RESULT([yes])
fi

## throw error if needed
if test x"${gmp_ok}" == x"no"; then
  AC_MSG_NOTICE([Configuration failed because ${GMP_NAME} was not found. Try installing:])
  AC_MSG_NOTICE([  * deb: ${GMP_DEB_NAME} (Debian, Ubuntu)])
  AC_MSG_NOTICE([  * rpm: ${GMP_RPM_NAME} (Fedora, EPEL)])
  AC_MSG_NOTICE([  * brew: ${GMP_BREW_NAME} (Mac OSX)])
  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([If ${GMP_NAME} is already installed, check that 'pkg-config' is in])
  AC_MSG_NOTICE([your PATH and PKG_CONFIG_PATH contains a ${GMP_CONFIG_NAME}.pc file.])
  AC_MSG_NOTICE([])
  AC_MSG_ERROR([ERROR: Installation failed])
fi

# Verify MFPR settings work
## display message
AC_MSG_CHECKING([if ${MPFR_NAME} compiles])

## try compiling file with settings
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <stdlib.h>
#include ${MPFR_TEST_HEADER}
#if defined(MPFR_VERSION_MAJOR)
#define SUCCESS 1
#else
#define SUCCESS 0
#endif
int main() {
  if (SUCCESS > 0) {
    exit(0);
  } else {
    exit(1);
  }
}
]])],
[mpfr_ok=yes],
[mpfr_ok=no],
[mpfr_ok=yes])
if test x"${mpfr_ok}" == x"no"; then
    AC_MSG_RESULT([no])
else
    AC_MSG_RESULT([yes])
fi

## throw error if needed
if test x"${mpfr_ok}" == x"no"; then
  AC_MSG_NOTICE([Configuration failed because ${MPFR_NAME} was not found. Try installing:])
  AC_MSG_NOTICE([  * deb: ${MPFR_DEB_NAME} (Debian, Ubuntu)])
  AC_MSG_NOTICE([  * rpm: ${MPFR_RPM_NAME} (Fedora, EPEL)])
  AC_MSG_NOTICE([  * brew: ${MPFR_BREW_NAME} (Mac OSX)])
  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([If ${MPFR_NAME} is already installed, check that 'pkg-config' is in])
  AC_MSG_NOTICE([your PATH and PKG_CONFIG_PATH contains a ${MPFR_CONFIG_NAME}.pc file.])
  AC_MSG_NOTICE([])
  AC_MSG_ERROR([ERROR: Installation failed])
fi

# Export settings
## save settings to Makevars
AC_SUBST([PKG_CPPFLAGS],["${PKG_CPPFLAGS}"])
AC_SUBST([PKG_LIBS],["${PKG_LIBS}"])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

## display settings
AC_MSG_NOTICE([Package CPPFLAGS: ${PKG_CPPFLAGS}])
AC_MSG_NOTICE([Package LIBS:     ${PKG_LIBS}])