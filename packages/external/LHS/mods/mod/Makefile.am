## Process this file with automake to produce 'Makefile.in'.

## Need the m4 directory for non-native macro tests.
ACLOCAL_AMFLAGS = -I ../../m4

## HACK NOTICE: Libtool 1.5.X doesn't provide tagging support for F90 sources.
## Until it does (Libtool 2.0.X), we need to use the F77 tag rules.
LTFCCOMPILE = $(LIBTOOL) --tag=F77 --mode=compile $(FC) $(AM_FCFLAGS) \
	      $(FCFLAGS) $(FCFLAGS_f90)
FCLINK = $(LIBTOOL) --tag=F77 --mode=link $(FCLD) $(AM_FCFLAGS) $(FCFLAGS) \
	 $(FCFLAGS_f90) $(AM_LDFLAGS) $(LDFLAGS) -o $@

## Set the convenience library name for the primary LHS module.
noinst_LTLIBRARIES = libmod.la

## First LHS module source (built before anything else).
libmod_la_SOURCES = Parms.f90

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure depcomp install-sh \
                       missing

## Remove the module files on 'make clean'.
CLEANFILES = *.mod

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations
