#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(wlc, 0.0.0, cortini@lptl.jussieu.fr)
AC_CONFIG_SRCDIR([src/wlc.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AM_INIT_AUTOMAKE
LT_INIT

# Checks for libraries.
AC_CHECK_LIB(c, main)
AC_CHECK_LIB(m, [sqrt])
CFLAGS="-Wall -g -O0 -Wextra -Wshadow"

# check for GSL
AC_ARG_WITH(gsl,
            [  --with-gsl=DIR        Directory where the GSL is installed (optional)],
            GSL_DIR="$withval")

GSL_CONFIG=
AC_PATH_PROG([GSL_CONFIG], [gsl-config],, [$GSL_DIR/bin:$PATH])

GSL_INCLUDES=
GSL_LIBS=

if test -x "$GSL_CONFIG"; then
  GSL_INCLUDES=`$GSL_CONFIG --cflags`
  GSL_LIBS=`$GSL_CONFIG --libs`
else
  AC_ERROR([GSL not found. You need gsl to compile this package.
      See http://www.gnu.org/software/gsl/])
fi

AC_SUBST(GSL_INCLUDES)
AC_SUBST(GSL_LIBS)

# Checks for header files.
AC_CHECK_HEADERS([stdio.h stdlib.h math.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC

# output
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
