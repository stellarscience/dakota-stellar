#!/bin/csh
#
# This script creates the version of the day (votd) tar files
#
################################################################
#
# Usage: create_votd_tar.csh root_directory platform
################################################################

# set base directory
set ROOT     = $1
set PLATFORM = $2
set BASE     = $ROOT/test_$PLATFORM
set SQA      = $ROOT/OPT++/bin

# Make sure config/make ran successfully
if (-e $BASE/OPT++/make_optpp.err) then
  echo "OPTPP build failed" 
else
  echo "OPTPP built successfully"

  # if votd directory does not exist, create
  if (-e $BASE/votd) then
    echo "VOTD directory is $BASE/votd"
  else
    mkdir -p $BASE/votd
    echo "Creating base directory $BASE/votd"
  endif

  # change to votd directory
  cd $BASE/votd
  \rm -f votd.log
  \rm -f *.gz
  if (-e optpp) then # if OPT++ exist, remove it
    \rm -rf optpp
  endif

  # EXTERNAL VOTD (Linux only)
  # call extract_src script with head OPT++, 
  if ( $PLATFORM == 'linux' ) then
    echo "External VOTD: extract src tar files" >> votd.log
    # Create optpp_votd.src.tar
    $SQA/extract_src.csh votd head >>& votd.log
    gzip optpp*.tar
    # Secure copy VOTD tar files to csmr and extract VOTD HTML docs
    scp optpp_votd.src.tar.gz csmr.ca.sandia.gov:/home/pwillia/votd/.
  endif
