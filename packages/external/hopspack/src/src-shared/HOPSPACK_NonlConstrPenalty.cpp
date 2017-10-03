// $Id: HOPSPACK_NonlConstrPenalty.cpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_NonlConstrPenalty.cpp $

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
  @file HOPSPACK_NonlConstrPenalty.cpp
  @brief Implement HOPSPACK::NonlConstrPenalty.
*/

#include <math.h>     //-- FOR fabs, sqrt

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_NonlConstrPenalty.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Declare internal enumerations and constants
//----------------------------------------------------------------------

enum {
    NONE = 0,
    L2_SQUARED,
    L1,
    L1_SMOOTHED,
    L2,
    L2_SMOOTHED,
    LINF,
    LINF_SMOOTHED
};

static const string  sL2_SQUARED_NAME    = "L2 Squared";
static const string  sL1_NAME            = "L1";
static const string  sL1_SMOOTHED_NAME   = "L1 (smoothed)";
static const string  sL2_NAME            = "L2";
static const string  sL2_SMOOTHED_NAME   = "L2 (smoothed)";
static const string  sLINF_NAME          = "L_inf";
static const string  sLINF_SMOOTHED_NAME = "L_inf (smoothed)";

//---- VALIDATED WITH test/penalty-tests.
static const double  dMIN_L1_ALPHA   = 1.0e-20;
static const double  dMIN_LINF_ALPHA = 1.0e-20;


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
NonlConstrPenalty::NonlConstrPenalty (void)
{
    _nPenType = NONE;
    _dPenCoef = 0.0;
    _dSmoothingFactor = 0.0;
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
NonlConstrPenalty::~NonlConstrPenalty (void)
{
    return;
}


//----------------------------------------------------------------------
//  Method defineFunction
//----------------------------------------------------------------------
bool  NonlConstrPenalty::defineFunction (const string &  sPenaltyName,
                                         const double    dPenaltyCoefficient,
                                         const double    dSmoothingFactor)
{
    if (sPenaltyName == sL2_SQUARED_NAME)
        _nPenType = L2_SQUARED;
    else if (sPenaltyName == sL1_NAME)
        _nPenType = L1;
    else if (sPenaltyName == sL1_SMOOTHED_NAME)
        _nPenType = L1_SMOOTHED;
    else if (sPenaltyName == sL2_NAME)
        _nPenType = L2;
    else if (sPenaltyName == sL2_SMOOTHED_NAME)
        _nPenType = L2_SMOOTHED;
    else if (sPenaltyName == sLINF_NAME)
        _nPenType = LINF;
    else if (sPenaltyName == sLINF_SMOOTHED_NAME)
        _nPenType = LINF_SMOOTHED;
    else
    {
        cerr << "ERROR: Unknown penalty function '" << sPenaltyName
             << "'" << endl;
        return( false );
    }

    _dPenCoef = dPenaltyCoefficient;
    if (exists (_dPenCoef) == false)
    {
        cerr << "ERROR: Penalty function coefficient does not exist" << endl;
        return( false );
    }

    updateSmoothing (dSmoothingFactor);

    return( true );
}


//----------------------------------------------------------------------
//  Method isDefined
//----------------------------------------------------------------------
bool  NonlConstrPenalty::isDefined (void) const
{
    return( _nPenType != NONE );
}


//----------------------------------------------------------------------
//  Method updateCoefficient
//----------------------------------------------------------------------
void  NonlConstrPenalty::updateCoefficient (const double  dPenaltyCoefficient)
{
    _dPenCoef = dPenaltyCoefficient;
    return;
}


//----------------------------------------------------------------------
//  Method updateSmoothing
//----------------------------------------------------------------------
void  NonlConstrPenalty::updateSmoothing (const double  dSmoothingFactor)
{
    _dSmoothingFactor = dSmoothingFactor;

    if (   (_nPenType != L1_SMOOTHED)
        && (_nPenType != L2_SMOOTHED)
        && (_nPenType != LINF_SMOOTHED) )
    {
        _dSmoothingFactor = 0.0;
    }

    if ((_nPenType == L1_SMOOTHED) && (_dSmoothingFactor < dMIN_L1_ALPHA))
    {
        _dSmoothingFactor = dMIN_L1_ALPHA;
        cerr << "WARNING: Smoothing parameter for '" << sL1_SMOOTHED_NAME
             << "' cannot be too close to zero" << endl;
        cerr << "         Changing smoothing parameter to "
             << _dSmoothingFactor << endl;
    }

    if ((_nPenType == LINF_SMOOTHED) && (_dSmoothingFactor < dMIN_LINF_ALPHA))
    {
        _dSmoothingFactor = dMIN_LINF_ALPHA;
        cerr << "WARNING: Smoothing parameter for '" << sLINF_SMOOTHED_NAME
             << "' cannot be too close to zero" << endl;
        cerr << "         Changing smoothing parameter to "
             << _dSmoothingFactor << endl;
    }

    return;
}


//----------------------------------------------------------------------
//  Method getPenaltyName
//----------------------------------------------------------------------
const string &  NonlConstrPenalty::getPenaltyName (void) const
{
    switch (_nPenType)
    {
    case L2_SQUARED:
        return( sL2_SQUARED_NAME );
    case L1:
        return( sL1_NAME );
    case L1_SMOOTHED:
        return( sL1_SMOOTHED_NAME );
    case L2:
        return( sL2_NAME );
    case L2_SMOOTHED:
        return( sL2_SMOOTHED_NAME );
    case LINF:
        return( sLINF_NAME );
    case LINF_SMOOTHED:
        return( sLINF_SMOOTHED_NAME );
    }

    cerr << "ERROR: Undefined penalty function type"
         << "  <HOPSPACK::NonlConstrPenalty>" << endl;
    throw INTERNAL_ERROR;
}


//----------------------------------------------------------------------
//  Method getCoefficient
//----------------------------------------------------------------------
double  NonlConstrPenalty::getCoefficient (void) const
{
    return( _dPenCoef );
}


//----------------------------------------------------------------------
//  Method getSmoothing
//----------------------------------------------------------------------
double  NonlConstrPenalty::getSmoothing (void) const
{
    return( _dSmoothingFactor );
}


//----------------------------------------------------------------------
//  Method computePenalty
//----------------------------------------------------------------------
double  NonlConstrPenalty::computePenalty (const Vector &  cEqs,
                                           const Vector &  cIneqs) const
{
    switch (_nPenType)
    {
    case L2_SQUARED:
        return( computeL2Sqrd_ (cEqs, cIneqs) );
    case L1:
        return( computeL1_ (cEqs, cIneqs) );
    case L1_SMOOTHED:
        return( computeL1Smoothed_ (cEqs, cIneqs) );
    case L2:
        return( computeL2_ (cEqs, cIneqs) );
    case L2_SMOOTHED:
        return( computeL2Smoothed_ (cEqs, cIneqs) );
    case LINF:
        return( computeLinf_ (cEqs, cIneqs) );
    case LINF_SMOOTHED:
        return( computeLinfSmoothed_ (cEqs, cIneqs) );
    }

    //---- THE INSTANCE TYPE MAY BE UNDEFINED IF A SUBPROBLEM DOES NOT
    //---- HAVE NONLINEAR CONSTRAINTS.
    return( 0.0 );
}


//----------------------------------------------------------------------
//  Method printDefinition
//----------------------------------------------------------------------
void  NonlConstrPenalty::printDefinition (void) const
{
    cout << "Nonlinear Constraint Penalty Function" << endl;

    if (isDefined() == false)
    {
        cout << "  Not defined" << endl;
        return;
    }

    cout << "  Type: " << getPenaltyName() << endl;

    cout << "  Penalty function weight = " << _dPenCoef << endl;
    cout << "  Smoothing factor        = " << _dSmoothingFactor << endl;

    return;
}


//----------------------------------------------------------------------
//  Private Method computeL2Sqrd_
//----------------------------------------------------------------------
double  NonlConstrPenalty::computeL2Sqrd_ (const Vector &  cEqs,
                                           const Vector &  cIneqs) const
{
    //---- L2^2 = rho * {  \sum_{i in eqs} c_i(x)^2
    //----               + \sum_{i in ineqs} max(-c_i(x), 0)^2 }

    return( _dPenCoef * computeSumSqs_ (cEqs, cIneqs) );
}


//----------------------------------------------------------------------
//  Private Method computeL1_
//----------------------------------------------------------------------
double  NonlConstrPenalty::computeL1_ (const Vector &  cEqs,
                                       const Vector &  cIneqs) const
{
    //---- L1 = rho * {  \sum_{i in eqs} |c_i(x)|
    //----             + \sum_{i in ineqs} max(-c_i(x), 0) }

    double  dNorm = 0.0;

    for (int  i = 0; i < cEqs.size(); i++)
        dNorm += fabs (cEqs[i]);

    for (int  i = 0; i < cIneqs.size(); i++)
        if (cIneqs[i] < 0.0)
            dNorm += fabs (cIneqs[i]);

    return( _dPenCoef * dNorm );
}


//----------------------------------------------------------------------
//  Private Method computeL1Smoothed_
//----------------------------------------------------------------------
double  NonlConstrPenalty::computeL1Smoothed_ (const Vector &  cEqs,
                                               const Vector &  cIneqs) const
{
    //---- MATHEMATICAL DEFINITION IS
    //----   L1sm =   \sum{i in eqs} theta (rho * c_i(x), alpha)
    //----          + \sum{i in ineqs} psi (-rho * c_i(x), alpha)
    //---- WHERE
    //----   theta (t,a) = a ln[2 + 2cosh(t/a)]
    //----   psi   (t,a) = a ln[1 + exp(t/a)]
    //----
    //---- THE SOFTWARE IMPLEMENTATION PROVIDES SAFEGUARDS WHEN alpha -> 0
    //---- (PARAMETER SATISFIES  0 < alpha <= 1).
    //----   CASE 1 (x >= 0):
    //----     ln(2 + 2cosh(x)) = ln[1 + exp(x) + 1 + exp(-x)]
    //----                      = ln[exp(x) (2*exp(-x) + 1 + exp(-2x))]
    //----                      = x + ln[(1 + exp(-x))^2]
    //----                      = x + 2 ln[1 + exp(-x)]
    //----     WE SEE THAT 1 < 1 + exp(-x) < 2 FOR x > 0.
    //----
    //----   CASE 2:
    //----     ln(1 + exp(x)) = ln[exp(|x|)exp(-|x|)(1 + exp(x))]
    //----                    = |x| + ln[exp(-|x|) + exp(x - |x|)]
    //----     WE SEE THAT exp(-|x|) + exp(x - |x|) < 2 FOR x > 0.

    double  dNorm = 0.0;

    for (int  i = 0; i < cEqs.size(); i++)
    {
        //---- x IS ALWAYS NONNEGATIVE.
        double  x = _dPenCoef * fabs (cEqs[i]) / _dSmoothingFactor;
        dNorm += _dSmoothingFactor * (x + (2.0 * log (1.0 + exp (-x))));
    }

    for (int  i = 0; i < cIneqs.size(); i++)
    {
        //---- x CAN BE POSITIVE (FEASIBLE) OR NEGATIVE (INFEASIBLE).
        double  x = - _dPenCoef * cIneqs[i] / _dSmoothingFactor;
        if (x <= 0.0)
        {
            dNorm += _dSmoothingFactor * log (1 + exp (x));
        }
        else
        {
            double  xAbs = fabs (x);
            dNorm += _dSmoothingFactor
                     * (xAbs + log (exp (-xAbs) + exp (x - xAbs)));
        }
    }

    return( dNorm );
}


//----------------------------------------------------------------------
//  Private Method computeL2_
//----------------------------------------------------------------------
double  NonlConstrPenalty::computeL2_ (const Vector &  cEqs,
                                       const Vector &  cIneqs) const
{
    //---- L2 = rho * sqrt{  \sum_{i in eqs} c_i(x)^2
    //----                 + \sum_{i in ineqs} max(-c_i(x), 0)^2 }

    double  dSumSqs = computeSumSqs_ (cEqs, cIneqs);

    if (dSumSqs <= 0.0)
        return( 0.0 );
    double  dNorm = sqrt (dSumSqs);

    return( _dPenCoef * dNorm );
}


//----------------------------------------------------------------------
//  Private Method computeL2Smoothed_
//----------------------------------------------------------------------
double  NonlConstrPenalty::computeL2Smoothed_ (const Vector &  cEqs,
                                               const Vector &  cIneqs) const
{
    //---- L2sm = rho * sqrt{  \sum_{i in eqs} c_i(x)^2
    //----                   + \sum_{i in ineqs} max{-c_i(x), 0}^2
    //----                   + (alpha / rho)^2 }

    if (_dPenCoef == 0.0)
        return( 0.0 );

    double  dSumSqs = computeSumSqs_ (cEqs, cIneqs);

    double  dTmp = _dSmoothingFactor / _dPenCoef;
    dSumSqs += dTmp * dTmp;

    if (dSumSqs <= 0.0)
        return( 0.0 );

    return( _dPenCoef * sqrt (dSumSqs) );
}


//----------------------------------------------------------------------
//  Private Method computeLinf_
//----------------------------------------------------------------------
double  NonlConstrPenalty::computeLinf_ (const Vector &  cEqs,
                                         const Vector &  cIneqs) const
{
    //---- Linf = rho * {  max( max_{i in eqs} |c_i(x)| ,
    //----                      max_{i in ineqs} max(-c_i(x), 0) ) }

    double  dNorm = 0.0;

    for (int  i = 0; i < cEqs.size(); i++)
        dNorm = max (dNorm, fabs (cEqs[i]));

    for (int  i = 0; i < cIneqs.size(); i++)
        if (cIneqs[i] < 0.0)
            dNorm = max (dNorm, fabs (cIneqs[i]));

    return( _dPenCoef * dNorm );
}


//----------------------------------------------------------------------
//  Private Method computeLinfSmoothed_
//----------------------------------------------------------------------
double  NonlConstrPenalty::computeLinfSmoothed_ (const Vector &  cEqs,
                                                 const Vector &  cIneqs) const
{
    //---- MATHEMATICAL DEFINITION IS
    //----   Linf sm = alpha ln[1 + \sum{i in eqs} 2cosh(rho * c_i(x) / alpha)
    //----                        + \sum{i in ineqs} exp(-rho * c_i(x) / alpha)]
    //----
    //---- THE SOFTWARE IMPLEMENTATION PROVIDES SAFEGUARDS WHEN alpha -> 0
    //---- (PARAMETER SATISFIES  0 < alpha <= 1).
    //----   LET  veq_i   = ceq_i(x)*rho/alpha
    //----        vineq_i = -cineq_i(x)*rho/alpha
    //----   SO THAT
    //----     Linf sm = alpha ln[1 + \sum{i in eqs} exp(veq_i) + exp(-veq_i)
    //----                          + \sum{i in ineqs} exp(vineq_i)]
    //----   DEFINE THE VECTOR
    //----        v = (veq_1, ..., vineq_1, ...)
    //----   LET  |v| = ||v||_inf
    //----   THEN
    //----     Linf sm
    //----     = alpha ln[exp(|v|) * exp(-|v|) *
    //----                (1 + \sum{i in eqs} exp(veq_i) + exp(-veq_i)
    //----                   + \sum{i in ineqs} exp(vineq_i) )]
    //----     = alpha * |v|
    //----       + alpha ln[exp(-|v|)
    //----                  + \sum{i in eqs} exp(veq_i - |v|) + exp(-veq_i - |v|)
    //----                  + \sum{i in ineqs} exp(vineq_i - |v|)]
    //----
    //----     NOW ALL EXPONENTS ARE LESS THAN OR EQUAL TO ZERO.

    double  cNorm = computeLinf_ (cEqs, cIneqs) / _dPenCoef;
    double  dFactor = _dPenCoef / _dSmoothingFactor;

    double  dInnerSum = exp (-cNorm * dFactor);

    for (int  i = 0; i < cEqs.size(); i++)
    {
        dInnerSum += exp (( cEqs[i] - cNorm) * dFactor);
        dInnerSum += exp ((-cEqs[i] - cNorm) * dFactor);
    }

    for (int  i = 0; i < cIneqs.size(); i++)
    {
        //---- cIneqs CAN BE POSITIVE (FEASIBLE) OR NEGATIVE (INFEASIBLE).
        dInnerSum += exp ((-cIneqs[i] - cNorm) * dFactor);
    }

    double  dNorm = _dPenCoef * cNorm + _dSmoothingFactor * (log (dInnerSum));
    return( dNorm );
}


//----------------------------------------------------------------------
//  Private Method computeSumSqs_
//----------------------------------------------------------------------
double  NonlConstrPenalty::computeSumSqs_ (const Vector &  cEqs,
                                           const Vector &  cIneqs) const
{
    double  dEqSumSqs = 0.0;
    for (int  i = 0; i < cEqs.size(); i++)
        dEqSumSqs += cEqs[i] * cEqs[i];

    double  dIneqSumSqs = 0.0;
    for (int  i = 0; i < cIneqs.size(); i++)
        if (cIneqs[i] < 0.0)
            dIneqSumSqs += cIneqs[i] * cIneqs[i];

    return( dEqSumSqs + dIneqSumSqs );
}


}     //-- namespace HOPSPACK
