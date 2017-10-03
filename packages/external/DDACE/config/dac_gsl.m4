dnl @synopsis DAC_GSL
dnl
dnl Process DDACE options for GSL
dnl
AC_DEFUN([DAC_GSL],[
  dnl GSL package checks.
  AC_ARG_WITH([gsl],
              AC_HELP_STRING([--with-gsl<=DIR>],
                             [use GPL package GSL (default no); optionally 
                              specify the root DIR for GSL include/lib]),
              [],[with_gsl="no"])

  acx_gsl_ok=disable
  case $with_gsl in
    no) ;;
    yes | "") acx_gsl_ok=want ;;
    *)
      if test -d "$withval" -a -d "$withval/include" -a -d "$withval/lib"; then
	acx_gsl_ok=want
        GSL_CPPFLAGS="-I$withval/include"
        GSL_LDFLAGS="-L$withval/lib"
      else
        AC_MSG_ERROR([specified GSL directory $withval must exist and contain 
                      include/ and lib/])
      fi
      ;;
  esac

  dnl prepend user-specified location then search for GSL
  dnl no granularity to notify user which used
  if test "x$acx_gsl_ok" = xwant; then

    if test -n $GSL_CPPFLAGS; then
      CPPFLAGS="$CPPFLAGS $GSL_CPPFLAGS"
      LDFLAGS="$LDFLAGS $GSL_LDFLAGS"
    fi

    AC_MSG_NOTICE([checking for GSL...])
    acx_gsl_ok=yes;
    AC_CHECK_LIB([m],[cos],,acx_gsl_ok=no)
    AC_CHECK_LIB([gslcblas],[cblas_dgemm],,acx_gsl_ok=no)
    AC_CHECK_LIB([gsl],[gsl_ran_fdist_pdf],,acx_gsl_ok=no)
    AC_CHECK_HEADERS([gsl/gsl_randist.h],,acx_gsl_ok=no)

    if test "x$acx_gsl_ok" = xyes; then
      AC_MSG_NOTICE([GNU GPL package GSL found])
      AC_MSG_NOTICE([NOTE: your build includes GNU GPL (binary) library GSL!])
      AC_DEFINE([HAVE_GSL],[1],[Macro to handle code which depends on GSL.])
      AC_SUBST(GSL_CPPFLAGS)
      AC_SUBST(GSL_LDFLAGS)
    else
      AC_MSG_ERROR([GSL requested but not found])
    fi

  fi
  AM_CONDITIONAL([WITH_GSL], [test "x$acx_gsl_ok" = xyes])

]) 
