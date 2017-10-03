#!/bin/csh
#
#  A csh script to make OPT++.  This is a simple shell script
#  which sets up the necessary symbolic links and runs the basic 
#  make and configure commands.
#
################################################################
#
#  Usage: make_script.csh configure_options
################################################################

setenv OPTPP $PWD
echo $OPTPP
set SQA = ../OPT++/bin
echo $SQA
unsetenv DAKOTA

# configue with arguments passed in and make clean
cd $OPTPP
./configure $*
make clean


# Do a full make
cd $OPTPP
if (-e make_optpp.out) then
  \rm -f make_optpp.out
endif
make >& $OPTPP/make_optpp.out
$SQA/grep_error.perl make_optpp.out
