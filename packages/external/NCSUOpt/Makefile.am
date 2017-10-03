## Process this file with automake to produce 'Makefile.in'.

## Set the NCSUOpt library name.
lib_LIBRARIES = libncsuopt.a

## NCSUOpt core library sources.
libncsuopt_a_SOURCES = DIRect.f DIRbatch.f DIRsubrout.f

## Extra files to remove when doing 'make clean'.
CLEANFILES = *.OUT *.Out *.out fort.*

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure depcomp install-sh \
                       missing

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Other files which must be included in the distribution.
EXTRA_DIST = SIPRA.INI KEYWORD.DAT InpData.txt
