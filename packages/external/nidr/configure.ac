# See notes in Makefile.am

# TODO
# * Probe for necessary headers/libs?
# * Manage optional libraries for dlopen
#   and also use a more generic name than dl_solver

AC_PREREQ([2.60])
AC_INIT([NIDR],[1.201006],[dakota-developers@development.sandia.gov])

AM_INIT_AUTOMAKE([1.9 foreign])
AM_MAINTAINER_MODE

AC_LANG([C])

dnl DL_SOLVER package check (toggles NIDR DL features and DL_SOLVER).
AC_ARG_WITH([dl_solver],
            AS_HELP_STRING([--with-dl_solver], [turn DL_SOLVER support on]),
	    [with_dl_solver=$withval], [with_dl_solver=no])
dnl if no dl_solver, -DNO_NIDR_DYNLIB set in src/Makefile.am
AM_CONDITIONAL([WITH_DL_SOLVER], [test "x$with_dl_solver" = xyes ])


dnl Specification maintenance mode; must be enabled to update NIDR
dnl itself or the specification of programs depending on it
AC_ARG_ENABLE([spec-maint],
              AS_HELP_STRING([--enable-spec-maint], 
                             [enable specification maintenance]),
              [enable_spec_maint=$enableval], [enable_spec_maint=no])
AM_CONDITIONAL([ENABLE_SPEC_MAINT], [test "x$enable_spec_maint" = xyes ])

AC_PROG_CC
AC_PROG_RANLIB

dnl Probe for programs needed to regenerate spec or NIDR itself
#if test "x$enable_spec_maint" = xyes; then
#  AC_PROG_LEX
#  AC_PROG_YACC
#  AC_PROG_SED
#fi
LEX=lex
AC_SUBST(LEX) 
YACC="bison -y"
AC_SUBST(YACC) 
SED=sed
AC_SUBST(SED) 


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
