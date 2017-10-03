#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Prerequisites.
AC_PREREQ([2.60])
AC_INIT([dfftpack],[4.0],[dakota@sandia.gov])
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

# Check for BLAS using a test from the GNU autoconf macro archive.
dnl ACX_BLAS
AC_LANG_POP([C])

# Write the Makefile and exit.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
