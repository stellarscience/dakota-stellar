## Process this file with automake to produce 'Makefile.in'.

## Need the m4 directory for non-native macro tests.
ACLOCAL_AMFLAGS = -I ../m4

## We must build the primary module before the rest of the modules.
SUBDIRS = mod

## HACK NOTICE: Libtool 1.5.X doesn't provide tagging support for F90 sources.
## Until it does (Libtool 2.0.X), we need to use the F77 tag rules.
LTFCCOMPILE = $(LIBTOOL) --tag=F77 --mode=compile $(FC) $(AM_FCFLAGS) \
	      $(FCFLAGS) $(FCFLAGS_f90) $(AX_FC_MODFLAG)mod
FCLINK = $(LIBTOOL) --tag=F77 --mode=link $(FCLD) $(AM_FCFLAGS) $(FCFLAGS) \
	 $(FCFLAGS_f90) $(AM_LDFLAGS) $(LDFLAGS) -o $@

## Set the convenience library name for the primary LHS modules.
noinst_LTLIBRARIES = libmods.la

## Define the sources for the primary LHS module.
LHS_MODS_SRC = Ccmatr.f90 Cparam.f90 Crank.f90   Csamp.f90 Cworkc.f90 \
               Cworkx.f90 Distnm.f90 Killfile.f90 Miscmod.f90 ModInByC.f90
libmods_la_SOURCES = $(LHS_MODS_SRC)

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure depcomp install-sh \
                       missing

## Remove the module files on 'make clean'.
CLEANFILES = *.mod

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations
