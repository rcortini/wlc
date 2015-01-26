#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(wlc, 0.1.0, cortini@lptl.jussieu.fr)
AC_CONFIG_SRCDIR([src/wlc.c])
AC_CONFIG_HEADERS([config.h])

PACKAGE=wlc
VERSION=0.1.0

# Checks for programs.
AC_PROG_CC
AM_INIT_AUTOMAKE
LT_INIT
AM_PROG_CC_C_O

AC_SUBST(PACKAGE)
AC_SUBST(VERSION)

# Checks for libraries.
AC_CHECK_LIB(c, main)
AC_CHECK_LIB(m, [sqrt])

# Global CFLAGS: thanks to Svyatoslav Kondrat for this part of the code
WARN_FLAGS="-Wall -Wextra -Wshadow -Wno-variadic-macros --pedantic"
CFLAGS="$WARN_FLAGS -O2"
CFLAGS_save="$CFLAGS"

# enable "optimize option", which will compile using only the -O3 option
AC_ARG_ENABLE([optimize],
    [AC_HELP_STRING([--enable-optimize],
	    [disable debugging (removes -g flags), disable compiler warnings\
	     and set optimization level O3 (disabled)])],
    [CFLAGS="-O3"], [CFLAGS="$CFLAGS_save"])

# enable "more-debug", which will compile using extra debugging options
CFLAGS_save="$CFLAGS"
AC_ARG_ENABLE([more-debug],
    [AC_HELP_STRING([--enable-more-debug],
	    [enable more debugging information (disabled)])],
     [CFLAGS="$WARN_FLAGS -g3 -ggdb -O0"], [CFLAGS="$CFLAGS_save"])

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
WLC_CFLAGS="-I$includedir"
WLC_LIBS="-L$libdir -lwlc"

AC_SUBST(WLC_CFLAGS)
AC_SUBST(WLC_LIBS)


# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC

# output
AC_CONFIG_FILES([Makefile
		 src/Makefile
		 wlc/Makefile
		 wlc-config wlc.pc])
AC_OUTPUT
