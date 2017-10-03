## Process this file with automake to produce 'Makefile.in'.

## Need the m4 directory for non-native macro tests.
ACLOCAL_AMFLAGS = -I m4

## We must build the module directory first.
SUBDIRS = mods

## Include the module directories from the build tree.
MODULES = $(AX_FC_MODFLAG)$(top_builddir)/mods \
	  $(AX_FC_MODFLAG)$(top_builddir)/mods/mod

## HACK NOTICE: Libtool 1.5.X doesn't provide tagging support for F90 sources.
## Until it does (Libtool 2.0.X), we need to use the F77 tag rules.
LTFCCOMPILE = $(LIBTOOL) --tag=F77 --mode=compile $(FC) $(AM_FCFLAGS) \
	      $(FCFLAGS) $(FCFLAGS_f90) $(MODULES)
FCLINK = $(LIBTOOL) --tag=F77 --mode=link $(FCLD) $(AM_FCFLAGS) $(FCFLAGS) \
	 $(FCFLAGS_f90) $(AM_LDFLAGS) $(LDFLAGS) -o $@
FCCOMPILE = $(FC) $(AM_FCFLAGS) $(FCFLAGS) $(FCFLAGS_f90) $(MODULES)

## Set the LHS library name.
lib_LTLIBRARIES = liblhs.la

## LHS library sources.
## -> LHS utility sources
LHS_UTIL_SRC = Banner.f90 Beta.f90 Betafn.f90 Betaic.f90 Betaln.f90 Binom.f90 \
               Chkdat.f90 Chkdim.f90 Chkemp.f90 Chkstr.f90 Chkzro.f90 \
               Chlsky.f90 Cmcrd.f90 Corcal.f90 Corout.f90 Cumulc.f90 \
               Cumuld.f90 C_wrap.f90 Datout.f90 Datsqz.f90 Dmfsd.f90 Dminmax.f90 \
               Dsinv.f90 Entrpf.f90 Entrpy.f90 Errchk.f90 Errget.f90 \
               Erstgt.f90 Erxset.f90 Expon.f90 Factor.f90 Factr2.f90 \
               Fileoc.f90 Findit.f90 Finvnor.f90 Gammab.f90 Gamma.f90 \
               Gammam.f90 Geom.f90 Hgeom.f90 Histo.f90 Hpsrt.f90 Hstout.f90 \
               Hypgeo.f90 Igaus1.f90 Igaus.f90 Igausf.f90 Imtql2.f90 \
               Interp.f90 Intrpd.f90 Ljust.f90 Matinv.f90 Mix.f90 Nbinom.f90 \
               Normal.f90 Outcrd.f90 Outdat.f90 Pareto.f90 Pmtrx.f90 \
               Poison.f90 Posdef.f90 Ranker.f90 Rierfc1.f90 Rmcnp2.f90 \
               Rmcnp.f90 Rmcnpi2.f90 Rmcnpi.f90 RnumLHS1.f90 RnumLHS2.f90 \
               Samout.f90 SamStor.f90 Setdef.f90 Sift.f90 Sspev.f90 Table.f90 \
               Tqlrat.f90 Trbak3.f90 Tred3.f90 Triang.f90 Unifrm.f90 Vif.f90 \
               Weibul.f90 Wrtcrd.f90 Wrtpar.f90 Gumbel.f90 Frechet.f90 \
		Lhssetseed.f90 RnumLHS10.f90 RnumLHS20.f90

## -> LHS core sources.
LHS_CORE_SRC = Lhs_cls.f90 Lhs_cnst.f90 Lhs_corr.f90 Lhs_cout.f90 Lhs_dist.f90 \
               Lhs_file.f90 Lhs_imem.f90 Lhs_init.f90 Lhs_optn.f90 Lhs_prep.f90 \
               Lhs_rsed.f90 Lhs_run.f90 Lhs_same.f90 Lhs_sdst.f90 Lhs_udst.f90

## Define the master LHS library sources.
liblhs_la_SOURCES = $(LHS_UTIL_SRC) $(LHS_CORE_SRC)

## Add the module convenience libraries to the master LHS library.
liblhs_la_LIBADD = $(top_builddir)/mods/mod/libmod.la \
		   $(top_builddir)/mods/libmods.la

## LHS standalone binary.
bin_PROGRAMS = lhsdrv
lhsdrv_SOURCES = fileleng.f90 flname.f90 lhs.f90 lhsdrv.f90 lread.f90 \
                 newcrd.f90 rdpar2.f90 rdpar.f90 read.f90
lhsdrv_LDADD = $(top_builddir)/liblhs.la

## Need to copy data files to the build directory for LHS to work.
all-local: $(srcdir)/InpData.txt $(srcdir)/KEYWORD.DAT $(srcdir)/SIPRA.INI
	for f in InpData.txt KEYWORD.DAT SIPRA.INI; do \
	  if test ! -e $${f}; then cp $(srcdir)/$${f} $${f}; fi; \
	done

## LHS test program.
check_PROGRAMS = TestRun
TestRun_SOURCES = TestRun.f90
TestRun_LDADD = $(top_builddir)/liblhs.la
TESTS = TestRun

## Extra files to remove when doing 'make clean'.
CLEANFILES = *.OUT *.Out *.out *.mod fort.*

## Extra files to remove when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 config.guess config.sub \
		       configure depcomp install-sh ltmain.sh missing

## Remove EXTRA_DIST files from the build directory if we're running with
## a separated build/source tree.
distclean-local: $(EXTRA_DIST)
	if test ! $(top_builddir) = $(top_srcdir); then rm -f $(EXTRA_DIST); fi

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Other files which must be included in the distribution.
EXTRA_DIST = SIPRA.INI KEYWORD.DAT InpData.txt
