#!/bin/csh
#
# A csh script to checkout OPT++ and create a source tar file for distribution.
# Arguments - tar file tag, optpp-version

if ($#argv < 2) then
  echo "Usage:   extract_src.csh file_tag [dakota-version]" 
  echo "Note:    use the string 'head' to extract the head version"
  echo "Example: extract_src.csh 3_1_2 head "
  exit
endif

# optpp Version
echo "Begin script to create optpp_$1.src.tar"
if ($2 == 'head') then
  cvs co -P OPT++
else
  cvs co -r -P $2 OPT++
endif

# ---------------------
# Miscellaneous Cleanup
# 08/10/2004 PJW This section needs to be revisited
# Want to be in top level directory?
# ---------------------
#cd ..
# remove all CVS directories
#find . -name CVS | xargs rm -rf
# remove other unneeded directories
#\rm -rf tests/sqa/v*
#cd ..

# -------------------
# Create src tar file
# -------------------
pwd
# use explict tar commands to tar only the files wanted
tar rvf optpp_$1.src.tar OPT++/config*
tar rvf optpp_$1.src.tar OPT++/docs
tar rvf optpp_$1.src.tar OPT++/include
tar rvf optpp_$1.src.tar OPT++/newmat11
tar rvf optpp_$1.src.tar OPT++/src
tar rvf optpp_$1.src.tar OPT++/tests
tar rvf optpp_$1.src.tar OPT++/INSTALL
tar rvf optpp_$1.src.tar OPT++/GNU_LGPL
tar rvf optpp_$1.src.tar OPT++/COPYRIGHT
tar rvf optpp_$1.src.tar OPT++/README
tar rvf optpp_$1.src.tar OPT++/configure.in
tar rvf optpp_$1.src.tar OPT++/Makefile.in

# Clean up this OPT++ area, since the tar file is now the master
# (and the OPT++ working directory and tar file will differ)
\rm -rf OPT++

echo "OPT++ tar file optpp_$1.src.tar complete"
