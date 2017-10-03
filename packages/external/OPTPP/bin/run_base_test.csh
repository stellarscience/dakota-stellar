#!/bin/csh
#
# This is the base test script for OPT++.  The script is written to run as
# a csh script because that preserves csh environment variables needed for 
# making OPT++
#
# Usage: run_base_test.csh root_directory optpp_version
#
# The package versions are passed on to the checkout_script.
# The script can be started at the command line or by using crontab.

# set up the environment variables
set ROOT = $1
set SQA  = $ROOT/OPT++/bin
set UNAME = `uname`
if ( $UNAME == 'SunOS' ) then
#  if ($4 == 'tflop')
#    set PLATFORM = tflop
#  else
    set PLATFORM = solaris
#  endif
else if ( $UNAME == 'AIX' ) then
  set PLATFORM = aix
else if ( $UNAME == 'IRIX64') then
  set PLATFORM = irix
else if ( $UNAME == 'Linux' ) then
  set PLATFORM = linux
else if ( $UNAME == 'OSF1' ) then
  set PLATFORM = osf
endif
set BASE = $ROOT/test_$PLATFORM

# if base directory does not exist, create
#if (-e $BASE) then
#  echo "Base directory is $BASE"
#else
#  mkdir -p $BASE
#  echo "Creating base directory $BASE"
#endif

# Checkout the correct versions of DAKOTA and its packages.
$SQA/checkout_script.csh $ROOT $BASE 

# Run make script which configures (with hardwired arguments) and builds dakota.
setenv OPTPP  $BASE/OPT++
cd $OPTPP
$SQA/make_script.csh 
#$SQA/make_script.csh --with-mpi 

# create VOTD distributions
$SQA/create_votd.csh $ROOT $PLATFORM

# Check to see if there were any errors.
# If not, run test suites (native builds only).
cd $OPTPP
if (! -e make_optpp.err) then
  make tests
  echo "Tests completed"
endif 

# Finish
cd $ROOT
