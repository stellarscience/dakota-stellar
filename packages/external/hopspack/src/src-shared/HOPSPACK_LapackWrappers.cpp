// $Id: HOPSPACK_LapackWrappers.cpp 208 2012-08-01 22:59:33Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_LapackWrappers.cpp $

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
  @file HOPSPACK_LapackWrappers.cpp
  @brief Implement functions declared in HOPSPACK_LapackWrappers.hpp.
*/

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_LapackWrappers.hpp"
#ifdef HAVE_FORTRAN_COMPILER
  #include "HOPSPACK_f90_config.h"
#endif

extern "C"
{
    //---- EXTERNAL DEFINITIONS ARE INFERRED FROM FORTRAN API.

    #ifdef HAVE_FORTRAN_COMPILER
      #define ddot_   HOPSPACK_F90_GLOBAL(ddot, DDOT)
      #define dgemv_  HOPSPACK_F90_GLOBAL(dgemv, DGEMV)
      #define dgemm_  HOPSPACK_F90_GLOBAL(dgemm, DGEMM)
      #define dgesvd_ HOPSPACK_F90_GLOBAL(dgesvd, DGESVD)
      #define dgglse_ HOPSPACK_F90_GLOBAL(dgglse, DGGLSE)
      #define dgelss_ HOPSPACK_F90_GLOBAL(dgelss, DGELSS)
      #define dgelqf_ HOPSPACK_F90_GLOBAL(dgelqf, DGELQF)
    #endif

    #if defined(HAVE_BLAS_F2C_WRAPPERS)
      #define ddot_ f2c_ddot
      #define dgemv_ f2c_dgemv
      #define dgemm_ f2c_dgemm
    #endif

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
    void  dgelss_ (int    *  m,
                   int    *  n,
                   int    *  nrhs,
                   double *  A,       //-- MATRIX
                   int    *  ldA,
                   double *  B,       //-- MATRIX
                   int    *  ldB,
                   double *  s,       //-- VECTOR
                   double *  rcond,
                   int    *  rank,
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


namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Static method:  getTheInstance
//----------------------------------------------------------------------
LapackWrappers &  LapackWrappers::getTheInstance (void)
{
    static LapackWrappers  theInstance;
    return( theInstance );
}


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
LapackWrappers::LapackWrappers (void)
{
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
LapackWrappers::~LapackWrappers (void)
{
    return;
}


//----------------------------------------------------------------------
//  Method ddot
//----------------------------------------------------------------------
double LapackWrappers::ddot (const int             n,
                             const double * const  x,
                             const double * const  y) const
{
  #if defined(HAVE_LAPACK)

    //---- SHOULD COPY THE ARRAY ARGUMENTS TO GUARANTEE CONSTNESS,
    //---- BUT WE WILL TRUST THAT LAPACK LEAVES THEM UNTOUCHED.
    double *  daXrecast = const_cast< double * >(x);
    double *  daYrecast = const_cast< double * >(y);

    int  nCopy = n;
    int  inc1 = 1;
    double  dResult = ::ddot_ (&nCopy, daXrecast, &inc1, daYrecast, &inc1);

    return( dResult );

  #else
    //---- NO LAPACK, BUT SIMPLE ENOUGH TO COMPUTE EXPLICITLY.
    double  dResult = 0.0;
    for (int i = 0; i < n; i++)
        dResult += x[i] * y[i];
    return( dResult );

  #endif
}


//----------------------------------------------------------------------
//  Method dgemv
//----------------------------------------------------------------------
void  LapackWrappers::dgemv (const char            trans,
                             const int             m,
                             const int             n,
                             const double          alpha,
                             const double * const  A,
                             const double * const  x,
                             const double          beta,
                                   double * const  y) const
{
  #if defined(HAVE_LAPACK)

    //---- SHOULD COPY THE ARRAY ARGUMENTS TO GUARANTEE CONSTNESS,
    //---- BUT WE WILL TRUST THAT LAPACK LEAVES THEM UNTOUCHED.
    double *  daArecast = const_cast< double * >(A);
    double *  daXrecast = const_cast< double * >(x);

    char    cTransCopy = trans;
    int     mCopy = m;
    int     nCopy = n;
    double  dAlphaCopy = alpha;
    double  dBetaCopy = beta;
    int     inc1 = 1;
    ::dgemv_ (&cTransCopy, &mCopy, &nCopy, &dAlphaCopy, daArecast,
	      &mCopy, daXrecast, &inc1, &dBetaCopy, y, &inc1);

    return;

  #else
    cerr << "ERROR: Cannot call dgemv, need to build with LAPACK" << endl;
    throw LAPACK_BUILD_ERROR;

  #endif
}


//----------------------------------------------------------------------
//  Method dgemm
//----------------------------------------------------------------------
void  LapackWrappers::dgemm (const char            transA,
                             const char            transB,
                             const int             m,
                             const int             n,
                             const int             k,
                             const double          alpha,
                             const double * const  A,
                             const double * const  B,
                             const double          beta,
                                   double * const  C) const
{
  #if defined(HAVE_LAPACK)

    if (transA != 'T')
    {
        //---- THIS CASE IS NOT IMPLEMENTED, BUT IS PROBABLY EASY.
        cerr << "ERROR: Cannot call dgemm with A untransposed" << endl;
        throw LAPACK_BUILD_ERROR;
    }

    //---- SHOULD COPY THE ARRAY ARGUMENTS TO GUARANTEE CONSTNESS,
    //---- BUT WE WILL TRUST THAT LAPACK LEAVES THEM UNTOUCHED.
    double *  daArecast = const_cast< double * >(A);
    double *  daBrecast = const_cast< double * >(B);

    char    cTransACopy = transA;
    char    cTransBCopy = transB;
    int     mCopy = m;
    int     nCopy = n;
    int     kCopy = k;
    double  dAlphaCopy = alpha;
    double  dBetaCopy = beta;

    int  ldB;
    if (transB == 'T')
        ldB = n;
    else
        ldB = k;
    ::dgemm_ (&cTransACopy, &cTransBCopy,
	      &mCopy, &nCopy, &kCopy,
	      &dAlphaCopy,
	      daArecast, &kCopy,
	      daBrecast, &ldB,
	      &dBetaCopy,
	      C, &mCopy);

    return;

  #else
    cerr << "ERROR: Cannot call dgemm, need to build with LAPACK" << endl;
    throw LAPACK_BUILD_ERROR;

  #endif
}


//----------------------------------------------------------------------
//  Method dgesvd
//----------------------------------------------------------------------
bool  LapackWrappers::dgesvd (const char            jobU,
                              const char            jobVT,
                              const int             m,
                              const int             n,
                                    double * const  A,
                                    double * const  Sigma,
                                    double * const  U,
                              const int             ldU,
                                    double * const  VT,
                              const int             ldVT) const
{
  #if defined(HAVE_LAPACK)

    //---- GENERAL SITUATIONS HAVE NOT BEEN TESTED, SO ASSUME THEY
    //---- MIGHT NOT WORK.
    if ((jobU != 'A') && (jobVT != 'A'))
    {
        cerr << "ERROR: Cannot call dgesvd for untested job types" << endl;
        throw LAPACK_BUILD_ERROR;
    }

    char  cJobUCopy = jobU;
    char  cJobVTCopy = jobVT;
    int   mCopy = m;
    int   nCopy = n;

    //---- A DOUBLE PRECISION WORK ARRAY IS REQUIRED WITH
    //---- SIZE >= MAX (3 * MIN(m,n) + MAX(m,n), 5 * MIN(m,n)).
    //---- WORKS BETTER IF SLIGHTLY LARGER; HENCE, DOUBLE THE SIZE.
    int  nWorkSize = 2 * max (3 * min(m,n) + max(m,n), 5 * min(m,n));
    double *  daWork = new double[nWorkSize];

    int  info = -1;
    ::dgesvd_ (&cJobUCopy, &cJobVTCopy,
	       &mCopy, &nCopy,
	       A, &mCopy,
	       Sigma, U, &mCopy,
	       VT, &nCopy,
	       daWork, &nWorkSize,
	       &info);

    delete[]  daWork;

    if (info != 0)
    {
        cerr << "WARNING: Call to LAPACK dgesvd failed" << endl;
        return( false );
    }
    return( true );

  #else
    cerr << "ERROR: Cannot call dgesvd, need to build with LAPACK" << endl;
    throw LAPACK_BUILD_ERROR;

  #endif
}


//----------------------------------------------------------------------
//  Method dgglse
//----------------------------------------------------------------------
bool  LapackWrappers::dgglse (const int             m,
                              const int             n,
                              const int             p,
                                    double * const  A,
                                    double * const  B,
                                    double * const  c,
                                    double * const  d,
                                    double * const  x) const
{
  #if defined(HAVE_LAPACK)

    int   mCopy = m;
    int   nCopy = n;
    int   pCopy = p;

    //---- A DOUBLE PRECISION WORK ARRAY IS REQUIRED WITH SIZE >= m+n+p.
    //---- FOR OPTIMIUM PERFORMANCE CHOOSE SIZE >= p+min(m,n)+max(m,n)*NB,
    //---- WHERE NB IS AN UPPER BOUND FOR THE OPTIMAL BLOCKSIZE FOR
    //---- DGEQRF, SGERQF, DORMQR AND SORMRQ.
    //---- IF nWorkSize = -1, THEN A WORKSPACE QUERY IS ASSUMED.
    //----
    //---- NOT SURE WHERE THIS HEURISTIC VALUE CAME FROM.
    int  nWorkSize = max (n*(n + 2), m + n + p);
    double *  daWork = new double[nWorkSize];

    int  info = -1;
    ::dgglse_ (&mCopy, &nCopy, &pCopy,
	       A, &mCopy,
	       B, &pCopy,
	       c, d, x,
	       daWork, &nWorkSize,
	       &info);

    delete[]  daWork;

    if (info != 0)
    {
        cerr << "WARNING: Call to LAPACK dgglse failed" << endl;
        return( false );
    }
    for (int  i = 0; i < m; i++)
    {
        //---- THIS WAS OBSERVED WITH INCONSISTENT (AND RANK DEFICIENT)
        //---- CONSTRAINTS:  x1 + x2 = 1, x1 + x2 = 2.
        if (isDoubleValid (x[i]) == false)
        {
            cerr << "WARNING: Call to LAPACK dgglse returned NaN result" << endl;
            return( false );
        }
    }
    return( true );

  #else
    cerr << "ERROR: Cannot call dgglse, need to build with LAPACK" << endl;
    throw LAPACK_BUILD_ERROR;

  #endif
}


//----------------------------------------------------------------------
//  Method dgelss
//----------------------------------------------------------------------
bool  LapackWrappers::dgelss (const int             m,
                              const int             n,
                                    double * const  A,
                              const double * const  c,
                                    double * const  x) const
{
  #if defined(HAVE_LAPACK)

    if (m < n)
    {
        cerr << "ERROR: Cannot call dgelss for underdetermined systems" << endl;
        throw LAPACK_BUILD_ERROR;
    }

    int   mCopy = m;
    int   nCopy = n;
    int   nOne = 1;

    //---- A DOUBLE PRECISION WORK ARRAY IS REQUIRED WITH SIZE OF AT LEAST
    //---- 3*min(m,n) + max( 2*min(m,n), max(m,n) ).  MORE WORK SPACE
    //---- ALLOWS BETTER PERFORMANCE, SO TWICE AS MUCH IS USED (NO TESTING
    //---- WAS PERFORMED TO OBTAIN THIS).
    //---- IF nWorkSize = -1, THEN A WORKSPACE QUERY IS ASSUMED.
    int  nWorkSize = 2 * (3 * min (m, n) + max (2 * min (m, n), max (m, n)));
    double *  daWork = new double[nWorkSize];

    int       nTmp = min (m, n);
    double *  daS = new double[nTmp];
    double    dRcond = 1.0e-12;
    int       nRank;

    double *  daCopyC = new double[m];
    for (int  i = 0; i < m; i++)
        daCopyC[i] = c[i];

    int  info = -1;
    ::dgelss_ (&mCopy, &nCopy, &nOne,
	       A, &mCopy,
	       daCopyC, &mCopy,
	       daS, &dRcond, &nRank,
	       daWork, &nWorkSize,
	       &info);
    for (int  i = 0; i < n; i++)
        x[i] = daCopyC[i];

    delete[]  daCopyC;
    delete[]  daWork;
    delete[]  daS;

    if (info != 0)
    {
        cerr << "WARNING: Call to LAPACK dgelss failed" << endl;
        return( false );
    }
    return( true );

  #else
    cerr << "ERROR: Cannot call dgelss, need to build with LAPACK" << endl;
    throw LAPACK_BUILD_ERROR;

  #endif
}


//----------------------------------------------------------------------
//  Method dgelqf
//----------------------------------------------------------------------
bool  LapackWrappers::dgelqf (const int             m,
                              const int             n,
                                    double * const  A,
                                    double * const  tau) const
{
  #if defined(HAVE_LAPACK)

    int   mCopy = m;
    int   nCopy = n;

    //---- A DOUBLE PRECISION WORK ARRAY IS REQUIRED WITH SIZE >= m.
    //---- FOR OPTIMIUM PERFORMANCE CHOOSE SIZE >= m*NB
    //---- WHERE NB IS AN UPPER BOUND FOR THE OPTIMAL BLOCKSIZE.
    //---- IF nWorkSize = -1, THEN A WORKSPACE QUERY IS ASSUMED.
    //----
    //---- NOT SURE WHERE THIS HEURISTIC VALUE CAME FROM.
    int  nWorkSize = max (n*(n + 2), m);
    double *  daWork = new double[nWorkSize];

    int  info = -1;
    ::dgelqf_ (&mCopy, &nCopy,
	       A, &mCopy,
	       tau,
	       daWork, &nWorkSize,
	       &info);
    
    delete[]  daWork;

    if (info != 0)
    {
        cerr << "WARNING: Call to LAPACK dgelqf failed" << endl;
        return( false );
    }
    return( true );

  #else
    cerr << "ERROR: Cannot call dgelqf, need to build with LAPACK" << endl;
    throw LAPACK_BUILD_ERROR;

  #endif
}


}          //-- namespace HOPSPACK
