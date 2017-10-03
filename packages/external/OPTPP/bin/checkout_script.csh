#!/bin/csh
#
# csh script to checkout OPT++
#
# Usage:   checkout_script.csh root_directory base_directory optpp-version
# Example: checkout_script.csh /home/pwillia/regress_tests 
#                   /home/pwillia/regress_tests/test_linux head 

#read in arguments from the command line to use for checkouts
set ROOT = $1 # root directory, e.g. /home/pwillia/regress_tests
set BASE = $2 # base directory for builds, e.g. $ROOT/test_linux/
set DV   = $3 # OPT++ Version
set SQA  = $ROOT/OPT++/bin # sqa directory for script access
set CVSROOT = /var/cvs

# if base directory does not exist, create
if (-e $BASE) then
  echo "Base directory is $BASE"
else
  mkdir -p $BASE
  echo "Creating base directory $BASE"
endif

# Clean up previous checkouts
cd $BASE
if (-e OPT++) then # if OPT++ exist, remove it
  \rm -rf OPT++
endif

# Get the correct version of OPT++
if ($DV == 'head' || $DV == '') then
  echo "Checkout head version of OPT++"
  cvs co -P OPT++
else if ($DV == 'none') then
  echo "Bypass OPT++ checkout"
else
  echo "Checkout OPT++ -r $DV"
  cvs co -P -r $DV OPT++
endif

echo "Checkout Process Complete"
