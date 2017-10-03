#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Prerequisites.
AC_PREREQ([2.60])
AC_INIT([CONMIN],[1.0],[dakota@sandia.gov])
AM_INIT_AUTOMAKE([1.9.6 foreign])
AC_LANG([Fortran 77])

# Since we're using libtool, don't build shared libraries!
AC_DISABLE_SHARED

# Checks for programs.
AC_LANG_PUSH([C])
AC_PROG_F77
AC_PROG_INSTALL

# Tells autoconf that libtool will be used.
AC_PROG_LIBTOOL
AC_LANG_POP([C])

# Write the Makefiles and exit.
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
