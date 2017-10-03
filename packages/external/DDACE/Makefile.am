## Process this file with automake to produce 'Makefile.in'.

## Need the config directory for its m4 files.
ACLOCAL_AMFLAGS = -I config

## Source subdirectories to be recursed into during the build.
SUBDIRS = src tests

## Header files which need to be distributed.
include_HEADERS = include/ColumnHeader.h \
		  include/DataValue.h \
		  include/DDaceArraySampler.h \
		  include/DDaceBoxBehnkenSampler.h \
		  include/DDaceCentralCompositeSampler.h \
		  include/DDaceFactorialSampler.h \
		  include/DDaceLHSampler.h \
		  include/DDaceOALHSampler.h \
		  include/DDaceOASampler.h \
		  include/DDaceRandomSampler.h \
		  include/DDaceSamplePoint.h \
		  include/DDaceSampler.h \
		  include/DDaceUserInputSampler.h \
		  include/Distribution.h \
		  include/Factor.h \
		  include/FuncApproxBase.h \
		  include/FuncApprox.h \
		  include/MainEffectAnalyzer.h \
		  include/MainEffectsAnalyzer3.h \
		  include/MainEffectsConverter.h \
		  include/MainEffectsExcelOutput.h \
		  include/Mars.h \
		  include/Mean.h \
		  include/NormalDistribution.h \
		  include/OneWayANOVA.h \
		  include/PseudoRandomTestsOnly.h \
		  include/Response.h \
		  include/SmartPtr.h \
		  include/Statistics.h \
		  include/StdDeviation.h \
		  include/UniformDistribution.h \
		  include/ValueAndRowIndexAndColumnIndex.h \
		  include/VectorCountingNumbersAndCount.h \
		  include/cdflib.h \
		  include/galois.h \
		  include/oa.h \
		  include/xtndeclare.h \
		  include/xtndispatch.h \
		  include/xtnset.h \
		  tests/TestDDaceArraySampler.h \
		  tests/TestDDaceBoxBehnkenSampler.h \
		  tests/TestDDaceCentralCompositeSampler.h \
		  tests/TestDDaceFactorialSampler.h \
		  tests/TestDDaceLHSampler.h \
		  tests/TestDDaceOALHSampler.h \
		  tests/TestDDaceOASampler.h \
		  tests/TestDDaceRandomSampler.h \
		  tests/TestDDaceSamplePoint.h \
		  tests/TestDDaceSampler.h \
		  tests/TestDDaceUserInputSampler.h \
		  tests/TestDistribution.h \
		  tests/TestFactor.h \
		  tests/TestMainEffectsAnalyzer.h \
		  tests/TestMainEffectsConverter.h \
		  tests/TestMainEffectsExcelOutput.h \
		  tests/TestMarsAnalyzer.h \
		  tests/TestMean.h \
		  tests/TestNormalDistribution.h \
		  tests/TestOneWayANOVA.h \
		  tests/TestPseudoRandom.h \
		  tests/TestResponse.h \
		  tests/TestStdDeviation.h \
		  tests/TestUniformDistribution.h \
		  tests/arrcmp.h \
		  tests/Stack.h \
		  tests/suite.h \
		  tests/test.h

## Remove the autotools-generated files when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure ddace_config.h.in \
		       compile config.guess config.sub depcomp install-sh \
		       ltmain.sh missing

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Other files which must be included in the distribution.
EXTRA_DIST = COPYRIGHT GNU_LGPL INSTALL README
