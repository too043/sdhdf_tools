#                                           -*- Autoconf -*-
# This files generates the configure file
#
# If you change this file, you usually should run
# autoreconf -fvi
#
# Initialise
AC_PREREQ(2.61)
AC_INIT([SDHDFPROC], [1.0], [george.hobbs@csiro.au])
AC_CONFIG_MACRO_DIR(config)
AC_CONFIG_SRCDIR([src])
AC_DISABLE_SHARED
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_CONFIG_HEADER([config.h])

# Check for programs
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_CC
AC_PROG_F77
AC_PROG_LIBTOOL

# Check for libraries
LIB_PGPLOT

# Check for header files
AC_HEADER_STDC
#AC_CHECK_HEADER([sofa.h],
#                [AC_DEFINE([HAVE_SOFA_H], [1],
#                [Define to 1 if you have <sofa.h>.])],
#                [AC_MSG_ERROR([header sofa.h not found])])

AC_CHECK_HEADER([erfa.h],
                [AC_DEFINE([HAVE_ERFA_H], [1],
                [Define to 1 if you have <erfa.h>.])],
                [AC_MSG_ERROR([header erfa.h (SOFA replacement not found])])

AC_CHECK_HEADER([calceph.h],
                [AC_DEFINE([HAVE_CALEPH_H], [1],
                [Define to 1 if you have <calceph.h>.])],
                [AC_MSG_ERROR([header calceph.h not found])])

AC_CONFIG_FILES([
                 Makefile
                 src/Makefile
                ])
AC_OUTPUT

echo
echo "INSPECTA is ready to be compiled. Please run 'make'"
echo
