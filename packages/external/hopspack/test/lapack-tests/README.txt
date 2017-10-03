# $Id: README.txt 149 2009-11-12 02:40:41Z tplante $
# $URL: https://software.sandia.gov/svn/hopspack/trunk/test/lapack-tests/README.txt $
#
# ************************************************************************
#         HOPSPACK: Hybrid Optimization Parallel Search Package
#                 Copyright 2009 Sandia Corporation
#
#   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#   the U.S. Government retains certain rights in this software.
# ************************************************************************


ABOUT THE LAPACK TESTS

This directory provides an executable that tests LAPACK subroutines needed
by HOPSPACK.  LAPACK subroutines must be supplied in an external library
(see the User Manual for information about finding an LAPACK library).
The tests verify that the library executes correctly for simple test
problems.  In many cases LAPACK libraries come with their own test set
that is more comprehensive than the tests here.

HOPSPACK only uses LAPACK on the master processor and "citizen slaves".
If MPI is used to distribute HOPSPACK functionality, then LAPACK does not
need to be installed on processors that perform only function evaluations
("evalution slaves").

The LAPACK test provided runs on a single processor using a single thread.
It does not use MPI.


HOW TO RUN THE LAPACK TESTS

o Build information is provided in the User Manual.  Briefly,
  - Install the CMake build tool
  - Build the source and test code with CMake commands

o To run the LAPACK test execute

    HOPSPACK_main_lapack_tests


DESCRIPTION OF THE OUTPUT

o Each LAPACK subroutine runs a separate test.  If successful, then a
  single line is printed out.  All tests should pass.

o Example output:

    Passed - ddot test
    Passed - dgemv test 'N'
    Passed - dgemv test 'T'
    Passed - dgemm test 'N'
    Passed - dgemm test 'T'
    Passed - dgesvd test
    Passed - dgglse test
    Passed - dgelqf test
