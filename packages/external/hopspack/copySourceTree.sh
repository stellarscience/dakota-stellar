#!/bin/sh

#-- THIS UTILITY COPIES SOURCE, EXAMPLES, TEST, AND DOCS INTO A TAR FILE.
#-- THE TAR CAN BE SENT AS AN ALTERNATIVE TO FULL SVN ACCESS.
#-- DIRECTORIES AND FILES ARE NAMED EXPLICITLY TO AVOID COPYING .svn.
#--
#-- BUILD A PDF USER MANUAL ON WINDOWS, COPY THE FILE TO THIS DIRECTORY AS
#--   HopspackUserManual_2_0_x.pdf
#-- AND SET _DOCNAME BELOW.
#--
#-- ARGUMENT -nodoc WILL SUPPRESS CREATION OF DOXYGEN FILES.

_RELNAME=hopspack-2.0.2-src
_DOCNAME=HopspackUserManual_2_0_2.pdf
_TOPDIR=..

_NEWROOT=$_TOPDIR/$_RELNAME

if [ -d $_NEWROOT ]; then
  echo "*** Directory already exists: " $_NEWROOT
  exit 1
fi
if [ ! -f $_DOCNAME ]; then
  echo "*** Cannot find user manual $_DOCNAME, continuing"
fi

mkdir $_NEWROOT
cp CMakeLists.txt    $_NEWROOT/.
cp Configure*.cmake  $_NEWROOT/.


#-- COPY ALL THE src SUBDIRECTORIES.

_NEXT=src
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt    $_NEWROOT/$_NEXT/.
cp $_NEXT/Doxyfile.txt      $_NEWROOT/$_NEXT/.
cp $_NEXT/LICENSE_HOPSPACK  $_NEWROOT/$_NEXT/.

_NEXT=src/src-citizens
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/*.hpp           $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=src/src-citizens/citizen-gss
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/*.hpp           $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.
cp $_NEXT/*.h             $_NEWROOT/$_NEXT/.
cp $_NEXT/*.c             $_NEWROOT/$_NEXT/.

_NEXT=src/src-citizens/citizen-gss/cddlib
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/*.h             $_NEWROOT/$_NEXT/.
cp $_NEXT/*.c             $_NEWROOT/$_NEXT/.
cp $_NEXT/LICENSE_CDDLIB  $_NEWROOT/$_NEXT/.

_NEXT=src/src-citizens/citizen-gss-common
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/*.hpp           $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=src/src-citizens/citizen-gss-ms
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/*.hpp           $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=src/src-citizens/citizen-gss-nlc
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/*.hpp           $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=src/src-evaluator
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/*.hpp           $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=src/src-executor
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/*.hpp           $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=src/src-framework
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/*.hpp           $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=src/src-main
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/*.hpp           $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=src/src-shared
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt                 $_NEWROOT/$_NEXT/.
cp $_NEXT/*.hpp                          $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp                          $_NEWROOT/$_NEXT/.
cp $_NEXT/HOPSPACK_common.hpp.cmake      $_NEWROOT/$_NEXT/.


#-- COPY ALL THE examples SUBDIRECTORIES.

_NEXT=examples
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/README.txt      $_NEWROOT/$_NEXT/.

_NEXT=examples/1-var-bnds-only
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/example*.txt    $_NEWROOT/$_NEXT/.
cp $_NEXT/*.c             $_NEWROOT/$_NEXT/.

_NEXT=examples/2-linear-constraints
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/example*.txt    $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=examples/3-degen-linear-constraints
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/example*.txt    $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=examples/4-nonlinear-constraints
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/example*.txt    $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=examples/5-multi-start
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/example*.txt    $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=examples/linked-evaluator-example
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/README*.txt     $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.
cp $_NEXT/*.hpp           $_NEWROOT/$_NEXT/.


#-- COPY ALL THE test SUBDIRECTORIES.

_NEXT=test
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=test/lapack-tests
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/README.txt      $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=test/penalty-tests
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/README.txt      $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.

_NEXT=test/startpoint-tests
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/CMakeLists.txt  $_NEWROOT/$_NEXT/.
cp $_NEXT/README.txt      $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp           $_NEWROOT/$_NEXT/.


#-- CREATE DOXYGEN DOCUMENTATION.
if [ "$1" != "-nodoc" ]; then
  mkdir $_NEWROOT/src/doc_doxygen
  _PWD=`pwd`
  cd $_NEWROOT/src
  doxygen Doxyfile.txt
  cd $_PWD
fi


#-- COPY THE USER MANUAL.
if [ ! -f $_DOCNAME ]; then
  echo "*** Need to add the User Manual under $_NEWROOT/doc/"
else
  mkdir $_NEWROOT/doc
  cp $_DOCNAME       $_NEWROOT/doc/.
chmod ugo-wx $_NEWROOT/doc/$_DOCNAME
fi


#-- ZIP IT UP.
`(cd $_TOPDIR; tar cf $_RELNAME.tar $_RELNAME; gzip $_RELNAME.tar)`

echo "  Created $_NEWROOT and $_TOPDIR/$_RELNAME.tar.gz"
