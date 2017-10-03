// $Id: HOPSPACK_float.cpp 177 2010-11-24 19:56:01Z briadam $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_float.cpp $ 

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
  \file HOPSPACK_float.cpp
  \brief Implement functions declared in HOPSPACK_float.hpp.
*/

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_float.hpp"
#if defined(WIN32)
#elif defined(__APPLE__)
  #include "math.h"          //-- FOR isinf AND isnan
#elif defined (__SUNPRO_CC)
  #include "ieeefp.h"        //-- FOR finite
#else
  #include "math.h"          //-- FOR isinf AND isnan
#endif


//---- MACHINE_PRECISION defines "machine epsilon" for double precision numbers.
#if defined(DBL_EPSILON)
  #define MACHINE_PRECISION  DBL_EPSILON
#else
  #define MACHINE_PRECISION  2.220446049250313e-16
#endif


static int  g_nRandomNumberBase = 1;


double  HOPSPACK::dne (void)
{
    return( HOPSPACK_DBL_MAX );
}


bool  HOPSPACK::exists (const double  value)
{
    return( value != HOPSPACK_DBL_MAX );
}


double  HOPSPACK::getMachineEpsilon (void)
{
    return( MACHINE_PRECISION );
}


//----------------------------------------------------------------------
//  Function genRanNum
//----------------------------------------------------------------------
/*  Random number generator providing a uniform distribution over [0..1).
 *
 *  The rand() function in older MSVC runtimes is extremely poor.
 *  Linux is probably OK, but for consistency use our own generator.
 *
 *  The following linear congruential method is from "Algorithms in C",
 *  Robert Sedgewick, Addison-Wesley, 1990, who in turn credits D.E. Knuth.
 *  The "seed" is the initial value of g_nRandomNumberBase.
 */
double  HOPSPACK::genRandomNumber (void)
{
    const int  nMx = 100000000;     /*-- 100,000,000 */
    const int  nM1 = 10000;         /*-- SQRT OF nMx */
    const int  nBs = 31415821;

    int  p0, p1, q0, q1;
    int  k;

    //---- MULTIPLY TWO LARGE INTEGERS WITHOUT OVERFLOW.
    p0 = g_nRandomNumberBase % nM1;
    p1 = g_nRandomNumberBase / nM1;
    q0 = nBs % nM1;
    q1 = nBs / nM1;
    k = (((p0*q1 + p1*q0) % nM1) * nM1 + (p0*q0)) % nMx;

    k = (k + 1) % nMx;
    g_nRandomNumberBase = k;

    return( ((double) k) / ((double) nMx) );
}


//----------------------------------------------------------------------
//  Function isDoubleValid
//----------------------------------------------------------------------
bool  HOPSPACK::isDoubleValid (const double  d)
{
    //---- CMAKE WAS NOT HELPFUL IN FINDING THESE FUNCTIONS, SO THEY
    //---- ARE SPECIFICALLY CODED FOR EACH PLATFORM.
    //---- IF SUITABLE PRIMITIVES CANNOT BE FOUND FOR A COMPILER,
    //---- THEN IT IS ACCEPTABLE TO ALWAYS RETURN true.

    #if defined(WIN32)
        if (_finite (d) == 0)
            return( false );
        if (d != d)
            return( false );

    #elif defined(__APPLE__)
        if (isinf (d) != 0)
            return( false );
        if (isnan (d) != 0)
            return( false );

    #elif defined (__SUNPRO_CC)
        if (finite (d) == 0)
	    return( false );
        if (d != d)
            return( false );

    #else
        if (isinf (d) != 0)
            return( false );
        if (isnan (d) != 0)
            return( false );

    #endif

    return( true );
}
