#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Prerequisites.
AC_PREREQ([2.59])
AC_INIT([NL2SOL],[2.0],[dakota@sandia.gov])
AC_CONFIG_SRCDIR([d1mach.c])
AC_CONFIG_HEADER([nl2sol_config.h])
AM_INIT_AUTOMAKE([1.9 foreign])

# Since we're using libtool, don't build shared libraries!
AC_DISABLE_SHARED

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL

# Tells autoconf that libtool will be used.
AC_LANG_PUSH([C])
AC_PROG_LIBTOOL
AC_LANG_POP([C])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

# Write the Makefile and exit.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
