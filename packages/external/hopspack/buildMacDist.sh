#!/bin/sh

#-- THIS UTILITY ASSUMES SERIAL AND MULTITHREADED EXECUTABLES HAVE BEEN
#-- BUILT IN DIRECTORIES AT THE SAME LEVEL AS trunk.  COMMANDS ARE:
#--
#--   mkdir build_serial
#--   cd build_serial
#--   cmake ../hopspack-2.0 -DLAPACK_ADD_LIBS=/usr/lib/libblas.dylib -Ddebug=false
#--   gnumake
#--
#--   mkdir build_mt
#--   cd build_mt
#--   cmake ../hopspack-2.0 -DLAPACK_ADD_LIBS=/usr/lib/libblas.dylib -Ddebug=false -Dmt=yes
#--   gnumake
#--
#-- BUILD A PDF USER MANUAL ON WINDOWS (BECAUSE pdflatex ON LINUX HAS PROBLEMS
#-- WITH HYPERLINKS).  COPY THE FILE TO THIS DIRECTORY AS
#--   HopspackUserManual_2_0_x.pdf
#-- AND SET _DOCNAME BELOW.


_RELNAME=hopspack-2.0.2-mac32
_DOCNAME=HopspackUserManual_2_0_2.pdf
_TOPDIR=..

_NEWROOT=$_TOPDIR/$_RELNAME

if [ -d $_NEWROOT ]; then
  echo "*** Directory already exists: " $_NEWROOT
  exit 1
fi
if [ ! -d ../build_serial ]; then
  echo "*** Cannot find ../build_serial"
  exit 1
fi
if [ ! -d ../build_mt ]; then
  echo "*** Cannot find ../build_mt"
  exit 1
fi
if [ ! -f $_DOCNAME ]; then
  echo "*** Cannot find user manual $_DOCNAME"
  exit 1
fi


rm -f ../$_RELNAME.tar.gz

echo "Distribution will include executables from ../build_serial and ../build_mt"


mkdir $_NEWROOT

cp distREADME_mac.txt  $_NEWROOT/README_mac.txt
chmod ugo-w $_NEWROOT/README_mac.txt

cp src/LICENSE         $_NEWROOT/.
chmod ugo-w $_NEWROOT/LICENSE

cp ../build_serial/HOPSPACK_main_serial  $_NEWROOT/.
cp ../build_mt/HOPSPACK_main_threaded    $_NEWROOT/.
chmod ugo-w $_NEWROOT/HOPSPACK_*


#-- COPY THE examples SUBDIRECTORIES, EXCEPT linked-evaluator-example
#-- BECAUSE IT REQUIRES SOURCE CODE.

mkdir $_NEWROOT/examples

cp examples/README.txt    $_NEWROOT/examples/.
chmod ugo-w $_NEWROOT/examples/README.txt

_NEXT=examples/1-var-bnds-only
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/example*.txt          $_NEWROOT/$_NEXT/.
cp $_NEXT/*.c                   $_NEWROOT/$_NEXT/.
cp ../build_serial/$_NEXT/var*  $_NEWROOT/$_NEXT/.

_NEXT=examples/2-linear-constraints
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/example*.txt          $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp                 $_NEWROOT/$_NEXT/.
cp ../build_serial/$_NEXT/lin*  $_NEWROOT/$_NEXT/.

_NEXT=examples/3-degen-linear-constraints
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/example*.txt          $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp                 $_NEWROOT/$_NEXT/.
cp ../build_serial/$_NEXT/deg*  $_NEWROOT/$_NEXT/.

_NEXT=examples/4-nonlinear-constraints
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/example*.txt          $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp                 $_NEWROOT/$_NEXT/.
cp ../build_serial/$_NEXT/non*  $_NEWROOT/$_NEXT/.

_NEXT=examples/5-multi-start
mkdir $_NEWROOT/$_NEXT
cp $_NEXT/example*.txt          $_NEWROOT/$_NEXT/.
cp $_NEXT/*.cpp                 $_NEWROOT/$_NEXT/.
cp ../build_serial/$_NEXT/mul*  $_NEWROOT/$_NEXT/.


#-- COPY THE USER MANUAL.
mkdir $_NEWROOT/doc
cp $_DOCNAME                    $_NEWROOT/doc/.
chmod ugo-wx $_NEWROOT/doc/$_DOCNAME


#-- ZIP IT UP.
`(cd $_TOPDIR; tar cf $_RELNAME.tar $_RELNAME; gzip $_RELNAME.tar)`

echo "  Created $_NEWROOT and $_TOPDIR/$_RELNAME.tar.gz"
