# $Id: README.txt 149 2009-11-12 02:40:41Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/test/penalty-tests/README.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


ABOUT THE PENALTY TESTS

This directory provides an executable that tests penalty functions used
for nonlinear constraints in HOPSPACK.  The tests verify that certain
penalty functions compute correctly for simple test cases.

The penalty test provided runs on a single processor using a single thread.


HOW TO RUN THE PENALTY TESTS

o Build information is provided in the User Manual.  Briefly,
  - Install the CMake build tool
  - Build the source and test code with CMake commands

o To run the test execute

    HOPSPACK_main_penalty_tests


DESCRIPTION OF THE OUTPUT

o The executable runs several separate tests.  If successful, then a
  single line is printed out.  All tests should pass.

o Example output:

    Passed - L1 eq test
    Passed - L1 ineq test
    Passed - L1 smooth eq test
    Passed - L1 smooth ineq test
    Passed - Linf eq test
    Passed - Linf ineq test
    Passed - Linf smooth eq test
    Passed - Linf smooth ineq test
