#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Prerequisites.
AC_PREREQ([2.59])
AC_INIT([FSUDace],[1.0],[dakota@sandia.gov])
AC_CONFIG_SRCDIR([fsu.H])
AC_CONFIG_HEADER([fsudace_config.h])
AM_INIT_AUTOMAKE([1.9 foreign])
AC_LANG([C++])
# Future feature (for now, we rely on ACLOCAL_AMFLAGS in Makefile.am):
#AC_CONFIG_MACRO_DIR([../../m4])

# Since we're using libtool, don't build shared libraries!
AC_DISABLE_SHARED

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL

# Tells autoconf that libtool will be used.
AC_LANG_PUSH([C])
AC_PROG_LIBTOOL
AC_LANG_POP([C])

# Checks for libraries.
AC_CHECK_LIB([m],[pow])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# WJB - ToDo: clean this mess up prior to CMakeification
dnl AC_CXX_HAVE_STD
AC_DEFINE(HAVE_STD,,[Unconditional dependency on an ISO C++ standard library/compiler])

# Checks for library functions.
AC_CHECK_FUNCS([rand random])

# Write the Makefile and exit.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
