// $Id: cmake_test_link_lapack.cpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/test/cmake_test_link_lapack.cpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This file is part of HOPSPACK.
//
// HOPSPACK is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library.  If not, see http://www.gnu.org/licenses/.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov)
//                 or Todd Plantenga (tplante@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/*!
  @file cmake_test_link_lapack.cpp
  @brief Main program used by ConfigureLapack.cmake to test LAPACK linking.
*/


//--------------------------------------------------------------------
//! Main program used by ConfigureLapack.cmake to test LAPACK linking.
/*!
 *  This may or may not be compiled during CMake configuration, but
 *  is never executed.  The purpose is to test whether an LAPACK function
 *  can be linked correctly.
 */
//--------------------------------------------------------------------


extern "C"
{
    //---- EXTERNAL DEFINITIONS ARE INFERRED FROM FORTRAN API.
    double ddot_ (int    *  n,
                  double *  x,        //-- VECTOR
                  int    *  incx,
                  double *  y,        //-- VECTOR
                  int    *  incy);
    void  dgemv_ (char   *  trans,
                  int    *  m,
                  int    *  n,
                  double *  alpha,    //-- SCALAR
                  double *  A,        //-- MATRIX
                  int    *  ldA,
                  double *  x,        //-- VECTOR
                  int    *  incx,
                  double *  beta,     //-- SCALAR
                  double *  y,        //-- VECTOR
                  int    *  incy);
    void  dgemm_ (char   *  transA,
                  char   *  transB,
                  int    *  m,
                  int    *  n,
                  int    *  k,
                  double *  alpha,    //-- SCALAR
                  double *  A,        //-- MATRIX
                  int    *  ldA,
                  double *  B,        //-- MATRIX
                  int    *  ldB,
                  double *  beta,     //-- SCALAR
                  double *  C,        //-- MATRIX
                  int    *  ldC);
    void  dgesvd_ (char   *  jobU,
                   char   *  jobVT,
                   int    *  m,
                   int    *  n,
                   double *  A,       //-- MATRIX
                   int    *  ldA,
                   double *  Sigma,   //-- VECTOR
                   double *  U,       //-- MATRIX
                   int    *  ldU,
                   double *  VT,      //-- MATRIX
                   int    *  ldVT,
                   double *  work,    //-- VECTOR
                   int    *  workSize,
                   int    *  info);
    void  dgglse_ (int    *  m,
                   int    *  n,
                   int    *  p,
                   double *  A,       //-- MATRIX
                   int    *  ldA,
                   double *  B,       //-- MATRIX
                   int    *  ldB,
                   double *  c,       //-- VECTOR
                   double *  d,       //-- VECTOR
                   double *  x,       //-- VECTOR
                   double *  work,    //-- VECTOR
                   int    *  workSize,
                   int    *  info);
    void  dgelqf_ (int    *  m,
                   int    *  n,
                   double *  A,       //-- MATRIX
                   int    *  ldA,
                   double *  Tau,     //-- MATRIX
                   double *  work,    //-- VECTOR
                   int    *  workSize,
                   int    *  info);
}


/*
 *  LAPACK calls merely match the API.  They should compile but are not
 *  intended for execution.
 */
int  main (void)
{
    char    c = 'U';
    int     n = 1;
    double  d = 1.0;
    double  daA[1] = { 0.0 };

  #if defined(TEST_ddot)
    ddot_ (&n, daA, &n, daA, &n);
  #endif
  #if defined(TEST_dgemv)
    dgemv_ (&c, &n, &n, &d, daA, &n, daA, &n, &d, daA, &n);
  #endif
  #if defined(TEST_dgemm)
    dgemm_ (&c, &c, &n, &n, &n, &d, daA, &n, daA, &n, &d, daA, &n);
  #endif
  #if defined(TEST_dgesvd)
    dgesvd_ (&c, &c, &n, &n, daA, &n, daA, daA, &n, daA, &n, daA, &n, &n);
  #endif
  #if defined(TEST_dgglse)
    dgglse_ (&n, &n, &n, daA, &n, daA, &n, &d, &d, &d, daA, &n, &n);
  #endif
  #if defined(TEST_dgelqf)
    dgelqf_ (&n, &n, daA, &n, daA, &d, &n, &n);
  #endif

    return( 0 );
}
