#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Prerequisites.
AC_PREREQ([2.60])
AC_INIT([LHS],[1.0],[dakota@sandia.gov])
AM_INIT_AUTOMAKE([1.9.6 foreign])
AC_LANG([Fortran])

# Since we're using libtool, don't build shared libraries!
AC_DISABLE_SHARED

# Checks for programs.
AC_LANG_PUSH([C])
AC_PROG_FC([],[90])
AC_PROG_INSTALL

# Tells autoconf that libtool will be used.
AC_PROG_LIBTOOL
AC_LANG_POP([C])

# Checks for compiler characteristics.
AC_FC_SRCEXT([f90])
AX_FC_MODULE_FLAG
AC_FC_FIXEDFORM

# Write the Makefile and exit.
AC_CONFIG_FILES([Makefile
		 mods/Makefile
		 mods/mod/Makefile])
AC_OUTPUT
