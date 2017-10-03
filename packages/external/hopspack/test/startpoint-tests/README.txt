# $Id: README.txt 154 2010-01-05 18:19:29Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/test/startpoint-tests/README.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#               Copyright 2009-2010 Sandia Corporation
# ************************************************************************


ABOUT THE STARTPOINT TESTS

This directory provides an executable that tests the projection of start
points onto linear constraints in HOPSPACK.  The tests verify correct
projection for a number of situations.

The test provided runs on a single processor using a single thread.


HOW TO RUN THE STARTPOINT TESTS

o Build information is provided in the User Manual.  Briefly,
  - Install the CMake build tool
  - Build the source and test code with CMake commands

o To run the test execute

    HOPSPACK_main_startpoint_tests


DESCRIPTION OF THE OUTPUT

o The executable runs several separate tests.  If successful, then a
  single line is printed out.  All tests should pass.

o Example output:

    Passed - Bnds test 1
    Passed - Bnds test 2
    Passed - Eqs test 1
    ...
