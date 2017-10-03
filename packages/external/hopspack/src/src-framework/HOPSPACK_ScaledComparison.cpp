// $Id: HOPSPACK_ScaledComparison.cpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-framework/HOPSPACK_ScaledComparison.cpp $ 

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
  @file HOPSPACK_ScaledComparison.cpp
  @brief Implement HOPSPACK::ScaledComparison.
*/

#include <math.h>     //-- FOR fabs

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_ScaledComparison.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Define static variables
//----------------------------------------------------------------------

bool    ScaledComparison::_bIsScalingDefined = false;
Vector  ScaledComparison::_cScalingFactors;

//---- TOLERANCE OF ZERO WILL CAUSE THE CACHE TO DECIDE TWO POINTS ARE
//---- DIFFERENT EVEN THOUGH THEIR COORDINATES AGREE TO MACHINE EPSILON.
double  ScaledComparison::_dToleranceTau = 2.0 * getMachineEpsilon();


//----------------------------------------------------------------------
//  Method setScaling
//----------------------------------------------------------------------
void  ScaledComparison::setScaling (const Vector &  cScaling)
{
    _cScalingFactors.resize (cScaling.size());
    for (int  i = 0; i < cScaling.size(); i++)
    {
        if (cScaling[i] <= 0.0)
        {
            cerr << "ERROR: Scaling vector elements must be positive"
                 << "  <ScaledComparison>"  << endl;
            throw INTERNAL_ERROR;
        }
        _cScalingFactors[i] = cScaling[i];
    }

    _bIsScalingDefined = true;
    return;
}


//----------------------------------------------------------------------
//  Method setTolerance
//----------------------------------------------------------------------
void  ScaledComparison::setTolerance (const double  dTolerance)
{
    if (dTolerance < 0.0)
    {
        cerr << "ERROR: Cache comparison tolerance cannot be negative"
             << "  <ScaledComparison>"  << endl;
        throw INTERNAL_ERROR;
    }
    _dToleranceTau = dTolerance;
    return;
}


//----------------------------------------------------------------------
//  Method isEqual
//----------------------------------------------------------------------
bool  ScaledComparison::isEqual (const Vector &  cPointA,
                                 const Vector &  cPointB)
{
    return( !isNotEqual (cPointA, cPointB) );
}


//----------------------------------------------------------------------
//  Method isNotEqual
//----------------------------------------------------------------------
bool  ScaledComparison::isNotEqual (const Vector &  cPointA,
                                    const Vector &  cPointB)
{
    checkSizes_ (cPointA, cPointB);

    for (int  i = 0; i < cPointA.size(); i++)
    {
        double  dTol = _dToleranceTau;
        if (_bIsScalingDefined == true)
            dTol = dTol * _cScalingFactors[i];

        if (fabs (cPointA[i] - cPointB[i]) > dTol)
        {
            //---- POINTS ARE NOT EQUAL.
            return( true );
        }
    }

    //---- POINTS ARE EQUAL.
    return( false );
}


//----------------------------------------------------------------------
//  Method isGreaterThan
//----------------------------------------------------------------------
bool  ScaledComparison::isGreaterThan (const Vector &  cPointA,
                                       const Vector &  cPointB)
{
    checkSizes_ (cPointA, cPointB);

    for (int  i = 0; i < cPointA.size(); i++)
    {
        double  dTol = _dToleranceTau;
        if (_bIsScalingDefined == true)
            dTol = dTol * _cScalingFactors[i];

        if (fabs (cPointA[i] - cPointB[i]) > dTol)
        {
            //---- POINTS ARE NOT EQUAL.
            if ((cPointA[i] - cPointB[i]) > dTol)
                //---- A > B.
                return( true );
            else
                //---- A < B.
                return( false );
        }
    }

    //---- A EQUALS B, IS NOT GREATER THAN.
    return( false );
}


//----------------------------------------------------------------------
//  Method isLessThan
//----------------------------------------------------------------------
bool  ScaledComparison::isLessThan (const Vector &  cPointA,
                                    const Vector &  cPointB)
{
    checkSizes_ (cPointA, cPointB);

    for (int  i = 0; i < cPointA.size(); i++)
    {
        double  dTol = _dToleranceTau;
        if (_bIsScalingDefined == true)
            dTol = dTol * _cScalingFactors[i];

        if (fabs (cPointA[i] - cPointB[i]) > dTol)
        {
            //---- POINTS ARE NOT EQUAL.
            if ((cPointB[i] - cPointA[i]) > dTol)
                //---- A < B.
                return( true );
            else
                //---- A > B.
                return( false );
        }
    }

    //---- A EQUALS B, IS NOT LESS THAN.
    return( false );
}


//----------------------------------------------------------------------
//  Method printDebugInfo
//----------------------------------------------------------------------
void  ScaledComparison::printDebugInfo (void)
{
    cout << "  HOPSPACK_ScaledComparison" << endl;
    cout << "    Tolerance (tau) = " << _dToleranceTau
         << " (Cache Comparison Tolerance)" << endl;
    if (_bIsScalingDefined == true)
    {
        for (int  i = 0; i < _cScalingFactors.size(); i++)
            cout << "    Scaling[" << i << "] = " << _cScalingFactors[i] << endl;
    }
    else
    {
        cout << "    Scaling factors all equal 1 (default)" << endl;
    }

    return;
}


//----------------------------------------------------------------------
//  Private Method checkSizes_
//----------------------------------------------------------------------
void  ScaledComparison::checkSizes_ (const Vector &  cPointA,
                                     const Vector &  cPointB)
{
    if (cPointA.size() != cPointB.size())
    {
        cerr << "ERROR: Cannot compare vectors of different sizes"
             << "  <ScaledComparison>"  << endl;
        throw INTERNAL_ERROR;
    }
    if (_bIsScalingDefined == true)
    {
        if (cPointA.size() != _cScalingFactors.size())
        {
            cerr << "ERROR: Cannot compare scale vector of different size"
                 << "  <ScaledComparison>"  << endl;
            throw INTERNAL_ERROR;
        }
    }

    return;
}


}     //-- namespace HOPSPACK
