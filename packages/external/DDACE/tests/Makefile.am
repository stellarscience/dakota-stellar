## Process this file with automake to produce Makefile.in

## Since these programs are defined in the TESTS variable, they will be run
## when 'make check' is invoked.
TESTS = AllTests
check_PROGRAMS = $(TESTS)

## AllTests is dependent on these files being available (in the event of a
## VPATH build, they should be copied to the build directory).
check_DATA = TestDDaceUserInputSamplerData

## AllTests is dependent upon these sources being compiled.
## Add TestMarsAnalyzer to the link object list iff ENABLE_STANDALONE has
## been set!
AllTests_SOURCES = AllTests.cpp \
                   arrcmp.cpp \
                   suite.cpp \
                   TestDDaceArraySampler.cpp \
                   TestDDaceBoxBehnkenSampler.cpp \
                   TestDDaceCentralCompositeSampler.cpp \
                   TestDDaceFactorialSampler.cpp \
                   TestDDaceLHSampler.cpp \
                   TestDDaceOALHSampler.cpp \
                   TestDDaceOASampler.cpp \
                   TestDDaceRandomSampler.cpp \
                   TestDDaceSamplePoint.cpp  \
                   TestDDaceSampler.cpp \
                   TestDDaceUserInputSampler.cpp \
                   TestDistribution.cpp \
                   TestFactor.cpp \
                   TestMainEffectsAnalyzer.cpp \
                   TestMainEffectsConverter.cpp \
                   TestMainEffectsExcelOutput.cpp
if ENABLE_STANDALONE
AllTests_SOURCES += TestMarsAnalyzer.cpp
endif
AllTests_SOURCES += TestMean.cpp \
		    TestNormalDistribution.cpp
if ENABLE_STANDALONE
AllTests_SOURCES += TestOneWayANOVA.cpp
endif
AllTests_SOURCES += TestPseudoRandom.cpp \
		    test.cpp \
		    TestResponse.cpp \
		    TestStdDeviation.cpp \
		    TestUniformDistribution.cpp

## We need to link against the DDACE and F77 libraries.
AllTests_LDADD = $(top_builddir)/src/libddace.la $(FLIBS)
if WITH_GSL
AllTests_LDADD += $(GSL_LDFLAGS)
endif

## Include directory for DDACE headers.
INCLUDES = -I$(top_srcdir)/include

## Remove any other files created during a build in this directory.
CLEANFILES = funcApprox.dat marsOutput.txt runTestLog.txt *_Log *.xml

## Remove the Makefile.in when doing 'make maintainer-clean'.
MAINTAINERCLEANFILES = Makefile.in

## make clean extension to handle compiler cache directories and the like.
clean-local:
	-rm -rf ii_files/ SunWS_cache/ so_locations

## Other files which must be included in the distribution.
EXTRA_DIST = TestDDaceUserInputSamplerData
