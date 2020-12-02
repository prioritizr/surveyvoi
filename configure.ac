# Initialization
AC_PREREQ(2.61)
AC_INIT([surveyvoi], 1.0.0)

# Find the compiler and compiler flags used by R.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

# Set compiler flags
## Get R compilers and flags.
## Add user-supplied flags to flags specified by R if they have been set.
CC=$("${R_HOME}/bin/R" CMD config CC)
CFLAGS="$("${R_HOME}/bin/R" CMD config CFLAGS) ${CFLAGS}"
CPP="$("${R_HOME}/bin/R" CMD config CC) -E"
CXX=$("${R_HOME}/bin/R" CMD config CXX)
CXXFLAGS="$("${R_HOME}/bin/R" CMD config CXXFLAGS) ${CXXFLAGS}"
CXXCPP="$("${R_HOME}/bin/R" CMD config CXX) -E"

# Set C++ specific options
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

# Set defaults
## usr include

## gmp
have_gmp=no
GMP_CXXFLAGS=""
GMP_LIBS="-lgmp"

## gmpxx
have_gmp=no
GMPXX_CXXFLAGS=""
GMPXX_LIBS="-lgmpxx"

## mpfr
have_gmp=no
MPFR_CXXFLAGS=""
MPFR_LIBS="-lmpfr"

## nlopt
have_nlopt=no
NLOPT_CXXFLAGS=""
NLOPT_LIBS="-lnlopt -lm -here"

# Find configurations
PKG_PROG_PKG_CONFIG

## GMP
if test [ -n "$PKG_CONFIG" ] ; then
  AC_MSG_CHECKING([pkg-config for gmp])
  if $PKG_CONFIG --exists gmp; then
    have_gmp=yes
    GMP_CXXFLAGS=`"${PKG_CONFIG}" --cflags gmp`
    GMP_LIBS=`"${PKG_CONFIG}" --libs gmp`
  fi
  AC_MSG_RESULT([${have_gmp}])
fi

if test "x${have_gmp}" = xno; then
  AC_CHECK_HEADERS(gmp.h, [have_gmp=yes], [AC_ERROR(
    [Header file gmp.h not found]
  )])
fi

## GMPXX
if test [ -n "$PKG_CONFIG" ] ; then
  AC_MSG_CHECKING([pkg-config for gmpxx])
  if $PKG_CONFIG --exists gmpxx; then
    have_gmpxx=yes
    GMPXX_CXXFLAGS=`"${PKG_CONFIG}" --cflags gmpxx`
    GMPXX_LIBS=`"${PKG_CONFIG}" --libs gmpxx`
  fi
  AC_MSG_RESULT([${have_gmpxx}])
fi

if test "x${have_gmpxx}" = xno; then
  AC_CHECK_HEADERS(gmpxx.h, [have_gmpxx=yes], [AC_ERROR(
    [Header file gmpxx.h not found]
  )])
fi

## MPFR
if test [ -n "$PKG_CONFIG" ] ; then
  AC_MSG_CHECKING([pkg-config for mpfr])
  if $PKG_CONFIG --exists mpfr; then
    have_mpfr=yes
    MPFR_CXXFLAGS=`"${PKG_CONFIG}" --cflags mpfr`
    MPFR_LIBS=`"${PKG_CONFIG}" --libs mpfr`
  fi
  AC_MSG_RESULT([${have_mpfr}])
fi

if test "x${have_mpfr}" = xno; then
  AC_CHECK_HEADERS(mpfr.h, [have_mpfr=yes], [AC_ERROR(
    [Header file mpfr.h not found]
  )])
fi

## NLOPT
if test [ -n "$PKG_CONFIG" ] ; then
  AC_MSG_CHECKING([pkg-config for nlopt])
  if $PKG_CONFIG --exists nlopt; then
    have_nlopt=yes
    NLOPT_CXXFLAGS=`"${PKG_CONFIG}" --cflags nlopt`
    NLOPT_LIBS=`"${PKG_CONFIG}" --libs nlopt`
  fi
  AC_MSG_RESULT([${have_nlopt}])
fi

if test "x${have_nlopt}" = xno; then
  AC_CHECK_HEADERS(nlopt.hpp, [have_nlopt=yes], [AC_ERROR(
    [Header file nlopt.hpp not found]
  )])
fi

# Write the flags into the src/Makevars file.
AC_SUBST([PKG_CPPFLAGS], ["${PKG_CPPFLAGS} -I/usr/local/include/ ${GMP_CXXFLAGS} ${GMPXX_CXXFLAGS} ${MPFR_CXXFLAGS} ${NLOPT_CXXFLAGS}"])
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LDFLAGS)
AC_SUBST([PKG_LIBS], ["${LIBS} -L/usr/local/lib/ ${PKG_LIBS} ${GMP_LIBS} ${GMPXX_LIBS} ${MPFR_LIBS} ${NLOPT_LIBS}"])
AC_OUTPUT

echo "
  --------------------------------------------------
  Configuration for ${PACKAGE_NAME} ${PACKAGE_VERSION}

    cppflags: ${CPPFLAGS} ${PKG_CPPFLAGS}
    libs:     ${PKG_LIBS}

  --------------------------------------------------
"
