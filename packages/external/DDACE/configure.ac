# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# Preliminaries.
AC_PREREQ([2.59])
AC_INIT([DDACE],[1.0],[pdhough@sandia.gov])
AC_CONFIG_SRCDIR([src/Analyzer/ColumnHeader.cpp])
AC_CONFIG_HEADERS([ddace_config.h])
AM_INIT_AUTOMAKE([1.9 foreign])

# Since we're using libtool, don't build shared libraries!
AC_DISABLE_SHARED

# Checks for programs.
AC_PROG_CC([xlc cc gcc])
AC_PROG_CXX([xlC cxx CC c++ g++])
AC_PROG_F77
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Tells autoconf that libtool will be used.
AC_PROG_LIBTOOL

# Enable DDACE standalone mode (adding FATools for MARS response
# surface fitting capability, and corresponding tests
# TestMarsAnalyzer.cpp, TestOneWayANOVA.cpp)
AC_ARG_ENABLE([standalone],
	      AC_HELP_STRING([--enable-standalone],
    		             [enable DDACE standalone library capability]),
	      [enable_standalone=$enableval],[enable_standalone=no])
AC_MSG_CHECKING([whether DDACE standalone is enabled])
if test x"$enable_standalone" = xyes; then
  AC_MSG_RESULT([yes])
  AC_DEFINE([HAVE_STANDALONE],[1],
            [Macro to handle code for the standalone library.])
else
  AC_MSG_RESULT([no])
fi
AM_CONDITIONAL(ENABLE_STANDALONE, test x"$enable_standalone" = xyes)

# GSL is only presently used for src/Analyzer/OneWayANOVA.cpp and related tests
DAC_GSL

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h stddef.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_CXX_CPPFLAGS_STD_LANG([CXXFLAGS])
AC_F77_WRAPPERS
AC_TYPE_SIZE_T
AC_HEADER_TIME

# Platform specific local template instantiation.
case "$ac_cv_cxx_compiler_vendor" in
    sgi)
        AX_CXXFLAGS_IRIX_OPTION([-ptused])
    ;;
    compaq)
        AX_CXXFLAGS_OSF_OPTION([-tlocal])
    ;;
    *)
    ;;
esac

# Checks for library functions.
AC_CHECK_FUNCS([getcwd gettimeofday memset strstr srand48 drand48])

# Output files
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/Analyzer/Makefile
                 src/Bose/Makefile
                 src/FATools/Makefile
                 src/Random/Makefile
                 src/Sampling/Makefile
                 tests/Makefile])
AC_OUTPUT
