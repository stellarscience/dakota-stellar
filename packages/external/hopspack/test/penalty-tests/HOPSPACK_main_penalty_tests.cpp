// $Id: HOPSPACK_main_penalty_tests.cpp 217 2013-11-25 21:59:49Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/test/penalty-tests/HOPSPACK_main_penalty_tests.cpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2013 Sandia Corporation
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
  @file HOPSPACK_main_penalty_tests.cpp
  @brief Main program that tests nonlinear constraint penalty functions.
*/

#include <iomanip>
#include <math.h>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_NonlConstrPenalty.hpp"
#include "HOPSPACK_Vector.hpp"


//---- SET THESE TO GET DETAILED OUTPUT FOR DEBUGGING.
static const bool  bDISPLAY_L1SM_EQ   = false;
static const bool  bDISPLAY_L1SM_INEQ = false;
static const bool  bDISPLAY_LINFSM_EQ   = false;
static const bool  bDISPLAY_LINFSM_INEQ = false;


//--------------------------------------------------------------------
//  Internal Function testL1_eq
//--------------------------------------------------------------------
// Test L1 penalty for a simple equality constraint.
static void  testL1_eq (void)
{
    using HOPSPACK::NonlConstrPenalty;
    using HOPSPACK::Vector;

    double  dCoef = 1.0;
    NonlConstrPenalty *  pPenalty = new NonlConstrPenalty();
    pPenalty->defineFunction ("L1", dCoef, 0.0);

    Vector  cEqs;
    cEqs.resize (1);
    Vector  cIneqs;
    for (int  nCoefCase = 0; nCoefCase < 3; nCoefCase++)
    {
        //---- CHECK VALUES AT A FEW POINTS OF c(x):  x-1 = 0.
        cEqs[0] = -1.0;     //-- x=0
        if (pPenalty->computePenalty (cEqs, cIneqs) != dCoef)
        {
            cout << "FAILED x=0 rho=" << dCoef << " - L1 eq test" << endl;
            throw "Exception from testL1_eq";
        }
        cEqs[0] =  0.0;     //-- x=1
        if (pPenalty->computePenalty (cEqs, cIneqs) != 0.0)
        {
            cout << "FAILED x=1 rho=" << dCoef << " - L1 eq test" << endl;
            throw "Exception from testL1_eq";
        }
        cEqs[0] =  1.0;     //-- x=2
        if (pPenalty->computePenalty (cEqs, cIneqs) != dCoef)
        {
            cout << "FAILED x=2 rho=" << dCoef << " - L1 eq test" << endl;
            throw "Exception from testL1_eq";
        }

        //---- UPDATE TO THE NEXT PENALTY COEFFICIENT.
        dCoef = dCoef * 100.0;
        pPenalty->updateCoefficient (dCoef);
    }
    cout << "Passed - L1 eq test" << endl;

    delete pPenalty;
    return;
}


//--------------------------------------------------------------------
//  Internal Function testL1_ineq
//--------------------------------------------------------------------
// Test L1 penalty for a simple inequality constraint.
static void  testL1_ineq (void)
{
    using HOPSPACK::NonlConstrPenalty;
    using HOPSPACK::Vector;

    double  dCoef = 1.0;
    NonlConstrPenalty *  pPenalty = new NonlConstrPenalty();
    pPenalty->defineFunction ("L1", dCoef, 0.0);

    Vector  cEqs;
    Vector  cIneqs;
    cIneqs.resize (1);
    for (int  nCoefCase = 0; nCoefCase < 3; nCoefCase++)
    {
        //---- CHECK VALUES AT A FEW POINTS OF c(x):  x-1 >= 0.
        cIneqs[0] = -1.0;     //-- x=0
        if (pPenalty->computePenalty (cEqs, cIneqs) != dCoef)
        {
            cout << "FAILED x=0 rho=" << dCoef << " - L1 ineq test" << endl;
            throw "Exception from testL1_ineq";
        }
        cIneqs[0] =  0.0;     //-- x=1
        if (pPenalty->computePenalty (cEqs, cIneqs) != 0.0)
        {
            cout << "FAILED x=1 rho=" << dCoef << " - L1 ineq test" << endl;
            throw "Exception from testL1_ineq";
        }
        cIneqs[0] =  1.0;     //-- x=2
        if (pPenalty->computePenalty (cEqs, cIneqs) != 0.0)
        {
            cout << "FAILED x=2 rho=" << dCoef << " - L1 ineq test" << endl;
            throw "Exception from testL1_ineq";
        }

        //---- UPDATE TO THE NEXT PENALTY COEFFICIENT.
        dCoef = dCoef * 100.0;
        pPenalty->updateCoefficient (dCoef);
    }
    cout << "Passed - L1 ineq test" << endl;

    delete pPenalty;
    return;
}


//--------------------------------------------------------------------
//  Internal Function computeFitToL1
//--------------------------------------------------------------------
// Return maximum and average error at a selected set of points.
static void  computeFitToL1 (const HOPSPACK::NonlConstrPenalty &  smoothedPen,
                             const double * const                 daCValues,
                             const bool                           bIsEq,
                             const int                            nNumCValues,
                                   double * const                 pMaxErr,
                                   double * const                 pAvgErr)
{
    using HOPSPACK::NonlConstrPenalty;
    using HOPSPACK::Vector;
    using HOPSPACK::isDoubleValid;

    double  dCoef = smoothedPen.getCoefficient();
    NonlConstrPenalty *  pPenalty = new NonlConstrPenalty();
    pPenalty->defineFunction ("L1", dCoef, 0.0);

    Vector  cConstraints;
    Vector  cDummy;
    cConstraints.resize (1);

    *pMaxErr = 0.0;
    double  dTotalErr = 0.0;
    double  dErr;

    for (int  i = 0; i < nNumCValues; i++)
    {
        cConstraints[0] = daCValues[i];

        double  d1;
        double  d2;
        if (bIsEq)
        {
            d1 = smoothedPen.computePenalty (cConstraints, cDummy);
            d2 = pPenalty->computePenalty (cConstraints, cDummy);
        }
        else
        {
            d1 = smoothedPen.computePenalty (cDummy, cConstraints);
            d2 = pPenalty->computePenalty (cDummy, cConstraints);
        }
        if (isDoubleValid (d1) == false)
        {
            cout << "  ***** Smoothed penalty value is infinite!" << endl;
            *pMaxErr = d1;
            *pAvgErr = d1;
            return;
        }
        dErr = d1 - d2;

        //---- SMOOTHED SHOULD ALWAYS BE GREATER THAN ACTUAL FOR THIS EXAMPLE.
        //---- ROUND-OFF ERROR CAN BE SIGNIFICANT FOR LARGE dCoef.
        if (dErr < -1.0e-9)
            cout << "Smoothed is not greater (" << dErr << ")"
                 << "  <computeFitToL1>" << endl;

        dTotalErr += fabs (dErr);
        if (dErr > *pMaxErr)
            *pMaxErr = dErr;
    }

    *pAvgErr = dTotalErr / ((double) nNumCValues);

    //---- MAXIMUM ERROR SHOULD BE AT c(x) = 0.
    //---- THIS TEST ASSUMES ONE OF daCValues IS ZERO.
    if (*pMaxErr > 1.0e-9)
    {
        cConstraints[0] = 0.0;
        if (bIsEq)
            dErr = smoothedPen.computePenalty (cConstraints, cDummy)
                   - pPenalty->computePenalty (cConstraints, cDummy);
        else
            dErr = smoothedPen.computePenalty (cDummy, cConstraints)
                   - pPenalty->computePenalty (cDummy, cConstraints);
        if (dErr != *pMaxErr)
            cout << "Smoothed max error not at [0] ("
                 << dErr << " vs " << *pMaxErr << ")"
                 << "  <computeFitToL1>" << endl;
    }

    return;
}


//--------------------------------------------------------------------
//  Internal Function testL1smooth_eq
//--------------------------------------------------------------------
// Test L1 smoothed penalty for a simple equality constraint.
// Primary concern is robustness for small values of alpha.
static void  testL1smooth_eq (void)
{
    using HOPSPACK::NonlConstrPenalty;
    using HOPSPACK::Vector;
    using HOPSPACK::isDoubleValid;

    //---- VALUES OF c(x) = x - 1 AT VARIOUS POINTS.
    double  daCValues[17] = { -2.0, -1.0, -0.5, -0.4, -0.3,
                              -0.2, -0.1, -0.05, 0.0,  0.05,
                               0.1,  0.2,  0.3,  0.4,  0.5,
                               1.0,  2.0 };

    double  dCoef = 1.0;
    double  dAlpha = 1.0;
    NonlConstrPenalty *  pPenalty = new NonlConstrPenalty();
    pPenalty->defineFunction ("L1 (smoothed)", dCoef, dAlpha);

    if (bDISPLAY_L1SM_EQ)
        cout << " Detailed results, L1 smoothed, for an equality" << endl;

    for (int  nCoefCase = 0; nCoefCase < 3; nCoefCase++)
    {
        pPenalty->updateCoefficient (dCoef);
        double  dAlpha = 1.0;
        pPenalty->updateSmoothing (dAlpha);

        for (int  nAlphaCase = 0; nAlphaCase < 20; nAlphaCase++)
        {
            double  dMaxErr;
            double  dAvgErr;
            computeFitToL1 (*pPenalty, daCValues, true, 17, &dMaxErr, &dAvgErr);
            if (bDISPLAY_L1SM_EQ)
            {
                streamsize  nCurrPrec = cout.precision (1);
                cout.setf (ios::scientific | ios::right);
                cout << "  r=" << dCoef << ", a=" << dAlpha;
                cout.precision (3);
                cout << " :  MaxErr = " << dMaxErr
                     << "  AvgErr = " << dAvgErr << endl;
                cout.precision (nCurrPrec);
            }

            if (isDoubleValid (dMaxErr) == false)
            {
                cout << "FAILED 1: rho=" << dCoef << ", alpha=" << dAlpha
                     << " - L1 smooth eq test" << endl;
                throw "Exception from testL1smooth_eq";
            }
            //---- EMPIRICAL TESTS.
            if ((dCoef == 1.0) && (dAlpha < 1.0) && (dMaxErr > 0.14))
            {
                cout << "FAILED 2: rho=" << dCoef << ", alpha=" << dAlpha
                     << " - L1 smooth eq test" << endl;
                throw "Exception from testL1smooth_eq";
            }
            if ((dAlpha <= 0.001) && (dAvgErr > 1.0e-4))
            {
                cout << "FAILED 3: rho=" << dCoef << ", alpha=" << dAlpha
                     << " - L1 smooth eq test" << endl;
                throw "Exception from testL1smooth_eq";
            }

            //---- UPDATE TO THE NEXT ALPHA.
            dAlpha = dAlpha / 10.0;
            pPenalty->updateSmoothing (dAlpha);
        }
        if (bDISPLAY_L1SM_EQ)
            cout << endl;

        //---- UPDATE TO THE NEXT PENALTY COEFFICIENT.
        dCoef = dCoef * 1000.0;
        pPenalty->updateCoefficient (dCoef);
    }
    cout << "Passed - L1 smooth eq test" << endl;

    delete pPenalty;
    return;
}


//--------------------------------------------------------------------
//  Internal Function testL1smooth_ineq
//--------------------------------------------------------------------
// Test L1 smoothed penalty for a simple inequality constraint.
// Primary concern is robustness for small values of alpha.
static void  testL1smooth_ineq (void)
{
    using HOPSPACK::NonlConstrPenalty;
    using HOPSPACK::Vector;
    using HOPSPACK::isDoubleValid;

    //---- VALUES OF c(x) = x - 1 AT VARIOUS POINTS.
    double  daCValues[17] = { -2.0, -1.0, -0.5, -0.4, -0.3,
                              -0.2, -0.1, -0.05, 0.0,  0.05,
                               0.1,  0.2,  0.3,  0.4,  0.5,
                               1.0,  2.0 };

    double  dCoef = 1.0;
    double  dAlpha = 1.0;
    NonlConstrPenalty *  pPenalty = new NonlConstrPenalty();
    pPenalty->defineFunction ("L1 (smoothed)", dCoef, dAlpha);

    if (bDISPLAY_L1SM_INEQ)
        cout << " Detailed results, L1 smoothed, for an inequality" << endl;

    for (int  nCoefCase = 0; nCoefCase < 3; nCoefCase++)
    {
        pPenalty->updateCoefficient (dCoef);
        double  dAlpha = 1.0;
        pPenalty->updateSmoothing (dAlpha);

        for (int  nAlphaCase = 0; nAlphaCase < 20; nAlphaCase++)
        {
            double  dMaxErr;
            double  dAvgErr;
            computeFitToL1 (*pPenalty, daCValues, false, 17, &dMaxErr, &dAvgErr);
            if (bDISPLAY_L1SM_INEQ)
            {
                streamsize  nCurrPrec = cout.precision (1);
                cout.setf (ios::scientific | ios::right);
                cout << "  r=" << dCoef << ", a=" << dAlpha;
                cout.precision (3);
                cout << " :  MaxErr = " << dMaxErr
                     << "  AvgErr = " << dAvgErr << endl;
                cout.precision (nCurrPrec);
            }

            if (isDoubleValid (dMaxErr) == false)
            {
                cout << "FAILED 1: rho=" << dCoef << ", alpha=" << dAlpha
                     << " - L1 smooth eq test" << endl;
                throw "Exception from testL1smooth_ineq";
            }
            //---- EMPIRICAL TESTS.
            if (dMaxErr > 0.7)
            {
                cout << "FAILED 2: rho=" << dCoef << ", alpha=" << dAlpha
                     << " - L1 smooth ineq test" << endl;
                throw "Exception from testL1smooth_ineq";
            }
            if ((dAlpha <= 0.001) && (dAvgErr > 1.0e-4))
            {
                cout << "FAILED 3: rho=" << dCoef << ", alpha=" << dAlpha
                     << " - L1 smooth ineq test" << endl;
                throw "Exception from testL1smooth_ineq";
            }

            //---- UPDATE TO THE NEXT ALPHA.
            dAlpha = dAlpha / 10.0;
            pPenalty->updateSmoothing (dAlpha);
        }
        if (bDISPLAY_L1SM_INEQ)
            cout << endl;

        //---- UPDATE TO THE NEXT PENALTY COEFFICIENT.
        dCoef = dCoef * 1000.0;
        pPenalty->updateCoefficient (dCoef);
    }
    cout << "Passed - L1 smooth ineq test" << endl;

    delete pPenalty;
    return;
}


//--------------------------------------------------------------------
//  Internal Function testLinf_eq
//--------------------------------------------------------------------
// Test Linf penalty for the case of two simple equality constraints:
//   c1(x,y) = x - 1 = 0  and  c2(x,y) = y - 2 = 0
// Test at points (-1, 1), (0,1), (1, 1), (2, 1), (3, 0), (2, 3), (1, 2).
static void  testLinf_eq (void)
{
    using HOPSPACK::NonlConstrPenalty;
    using HOPSPACK::Vector;

    //---- VALUES OF CONSTRAINTS AND EXPECTED RESULT AT TEST POINTS.
    double  daC1Values[7] = { -2.0, -1.0,  0.0,  1.0,  2.0,  1.0,  0.0 };
    double  daC2Values[7] = { -1.0, -1.0, -1.0, -1.0, -2.0,  1.0,  0.0 };
    double  daExpNorm[7]  = {  2.0,  1.0,  1.0,  1.0,  2.0,  1.0,  0.0 };

    double  dCoef = 1.0;
    NonlConstrPenalty *  pPenalty = new NonlConstrPenalty();
    pPenalty->defineFunction ("L_inf", dCoef, 0.0);

    Vector  cEqs;
    Vector  cIneqs;
    cEqs.resize (2);
    for (int  nCoefCase = 0; nCoefCase < 3; nCoefCase++)
    {
        for (int  i = 0; i < 7; i++)
        {
            cEqs[0] = daC1Values[i];
            cEqs[1] = daC2Values[i];
            double  dTmp = pPenalty->computePenalty (cEqs, cIneqs);
            if (dTmp != (daExpNorm[i] * dCoef))
            {
                cout << "FAILED rho=" << dCoef
                     << " ([" << i << "] "
                     << dTmp << " vs " << daExpNorm[i] << ")"
                     << " - Linf eq test" << endl;
                throw "Exception from testLinf_eq";
            }
        }

        //---- UPDATE TO THE NEXT PENALTY COEFFICIENT.
        dCoef = dCoef * 100.0;
        pPenalty->updateCoefficient (dCoef);
    }
    cout << "Passed - Linf eq test" << endl;

    delete pPenalty;
    return;
}


//--------------------------------------------------------------------
//  Internal Function testLinf_ineq
//--------------------------------------------------------------------
// Test Linf penalty for the case of two simple inequality constraints:
//   c1(x,y) = x - 1 >= 0  and  c2(x,y) = y - 2 >= 0
// Test at points (-1, 1), (0,1), (1, 1), (2, 1), (3, 0), (2, 3), (1, 2).
static void  testLinf_ineq (void)
{
    using HOPSPACK::NonlConstrPenalty;
    using HOPSPACK::Vector;

    //---- VALUES OF CONSTRAINTS AND EXPECTED RESULT AT TEST POINTS.
    double  daC1Values[7] = { -2.0, -1.0,  0.0,  1.0,  2.0,  1.0,  0.0 };
    double  daC2Values[7] = { -1.0, -1.0, -1.0, -1.0, -2.0,  1.0,  0.0 };
    double  daExpNorm[7]  = {  2.0,  1.0,  1.0,  1.0,  2.0,  0.0,  0.0 };

    double  dCoef = 1.0;
    NonlConstrPenalty *  pPenalty = new NonlConstrPenalty();
    pPenalty->defineFunction ("L_inf", dCoef, 0.0);

    Vector  cEqs;
    Vector  cIneqs;
    cIneqs.resize (2);
    for (int  nCoefCase = 0; nCoefCase < 3; nCoefCase++)
    {
        for (int  i = 0; i < 7; i++)
        {
            cIneqs[0] = daC1Values[i];
            cIneqs[1] = daC2Values[i];
            double  dTmp = pPenalty->computePenalty (cEqs, cIneqs);
            if (dTmp != (daExpNorm[i] * dCoef))
            {
                cout << "FAILED rho=" << dCoef
                     << " ([" << i << "] "
                     << dTmp << " vs " << daExpNorm[i] << ")"
                     << " - Linf eq test" << endl;
                throw "Exception from testLinf_ineq";
            }
        }

        //---- UPDATE TO THE NEXT PENALTY COEFFICIENT.
        dCoef = dCoef * 100.0;
        pPenalty->updateCoefficient (dCoef);
    }
    cout << "Passed - Linf ineq test" << endl;

    delete pPenalty;
    return;
}




//--------------------------------------------------------------------
//  Internal Function computeFitToLinf
//--------------------------------------------------------------------
// Return maximum and average error at a selected set of points.
static void  computeFitToLinf (const HOPSPACK::NonlConstrPenalty &  smoothedPen,
                               const double * const                 daC1Values,
                               const double * const                 daC2Values,
                               const bool                           bIsEq,
                               const int                            nNumCValues,
                                     double * const                 pMaxErr,
                                     double * const                 pAvgErr)
{
    using HOPSPACK::NonlConstrPenalty;
    using HOPSPACK::Vector;
    using HOPSPACK::isDoubleValid;

    double  dCoef = smoothedPen.getCoefficient();
    NonlConstrPenalty *  pPenalty = new NonlConstrPenalty();
    pPenalty->defineFunction ("L_inf", dCoef, 0.0);

    Vector  cConstraints;
    Vector  cDummy;
    cConstraints.resize (2);

    *pMaxErr = 0.0;
    double  dTotalErr = 0.0;
    double  dErr;

    for (int  i = 0; i < nNumCValues; i++)
    {
        cConstraints[0] = daC1Values[i];
        cConstraints[1] = daC2Values[i];

        double  d1;
        double  d2;
        if (bIsEq)
        {
            d1 = smoothedPen.computePenalty (cConstraints, cDummy);
            d2 = pPenalty->computePenalty (cConstraints, cDummy);
        }
        else
        {
            d1 = smoothedPen.computePenalty (cDummy, cConstraints);
            d2 = pPenalty->computePenalty (cDummy, cConstraints);
        }
        if (isDoubleValid (d1) == false)
        {
            cout << "  ***** Smoothed penalty value is infinite!" << endl;
            *pMaxErr = d1;
            *pAvgErr = d1;
            return;
        }
        dErr = d1 - d2;

        //---- SMOOTHED SHOULD ALWAYS BE GREATER THAN ACTUAL FOR THIS EXAMPLE.
        //---- ROUND-OFF ERROR CAN BE SIGNIFICANT FOR LARGE dCoef.
        if (dErr < -1.0e-9)
            cout << "Smoothed is not greater (" << dErr << ")"
                 << "  <computeFitToLinf>" << endl;

        dTotalErr += fabs (dErr);
        if (dErr > *pMaxErr)
            *pMaxErr = dErr;
    }

    *pAvgErr = dTotalErr / ((double) nNumCValues);

    //---- MAXIMUM ERROR SHOULD BE AT c(x) = 0.
    //---- THIS TEST ASSUMES ONE OF daCValues IS ZERO.
    if (*pMaxErr > 1.0e-9)
    {
        cConstraints[0] = 0.0;
        if (bIsEq)
            dErr = smoothedPen.computePenalty (cConstraints, cDummy)
                   - pPenalty->computePenalty (cConstraints, cDummy);
        else
            dErr = smoothedPen.computePenalty (cDummy, cConstraints)
                   - pPenalty->computePenalty (cDummy, cConstraints);
        if (dErr != *pMaxErr)
            cout << "Smoothed max error not at [0] ("
                 << dErr << " vs " << *pMaxErr << ")"
                 << "  <computeFitToLinf>" << endl;
    }

    return;
}


//--------------------------------------------------------------------
//  Internal Function testLinfSmooth_eq
//--------------------------------------------------------------------
// Test Linf smoothed penalty for two simple equality constraints:
//   c1(x,y) = x - 1 = 0  and  c2(x,y) = y - 2 = 0
// Test at points (-1, 1), (0,1), (1, 1), (2, 1), (3, 0), (2, 3), (1, 2).
// Primary concern is robustness for small values of alpha.
static void  testLinfSmooth_eq (void)
{
    using HOPSPACK::NonlConstrPenalty;
    using HOPSPACK::Vector;
    using HOPSPACK::isDoubleValid;

    //---- VALUES OF CONSTRAINTS AND EXPECTED RESULT AT TEST POINTS.
    double  daC1Values[7] = { -2.0, -1.0,  0.0,  1.0,  2.0,  1.0,  0.0 };
    double  daC2Values[7] = { -1.0, -1.0, -1.0, -1.0, -2.0,  1.0,  0.0 };

    double  dCoef = 1.0;
    double  dAlpha = 1.0;
    NonlConstrPenalty *  pPenalty = new NonlConstrPenalty();
    pPenalty->defineFunction ("L_inf (smoothed)", dCoef, dAlpha);

    if (bDISPLAY_LINFSM_EQ)
        cout << " Detailed results, Linf smoothed, for equalities" << endl;

    for (int  nCoefCase = 0; nCoefCase < 3; nCoefCase++)
    {
        pPenalty->updateCoefficient (dCoef);
        double  dAlpha = 1.0;
        pPenalty->updateSmoothing (dAlpha);

        for (int  nAlphaCase = 0; nAlphaCase < 20; nAlphaCase++)
        {
            double  dMaxErr;
            double  dAvgErr;
            computeFitToLinf (*pPenalty, daC1Values, daC2Values, true, 7,
                              &dMaxErr, &dAvgErr);
            if (bDISPLAY_LINFSM_EQ)
            {
                streamsize  nCurrPrec = cout.precision (1);
                cout.setf (ios::scientific | ios::right);
                cout << "  r=" << dCoef << ", a=" << dAlpha;
                cout.precision (3);
                cout << " :  MaxErr = " << dMaxErr
                     << "  AvgErr = " << dAvgErr << endl;
                cout.precision (nCurrPrec);
            }

            if (isDoubleValid (dMaxErr) == false)
            {
                cout << "FAILED 1: rho=" << dCoef << ", alpha=" << dAlpha
                     << " - Linf smooth eq test" << endl;
                throw "Exception from testLinfSmooth_eq";
            }
            //---- EMPIRICAL TESTS.
            if (dMaxErr > (2.0 * dAlpha))
            {
                cout << "FAILED 2: rho=" << dCoef << ", alpha=" << dAlpha
                     << " - Linf smooth eq test" << endl;
                throw "Exception from testLinfSmooth_ineq";
            }
            if ((dAlpha < 1.0) && (dAvgErr > (0.8 * dAlpha)))
            {
                cout << "FAILED 3: rho=" << dCoef << ", alpha=" << dAlpha
                     << " - Linf smooth eq test" << endl;
                throw "Exception from testLinfSmooth_ineq";
            }

            //---- UPDATE TO THE NEXT ALPHA.
            dAlpha = dAlpha / 10.0;
            pPenalty->updateSmoothing (dAlpha);
        }
        if (bDISPLAY_LINFSM_EQ)
            cout << endl;

        //---- UPDATE TO THE NEXT PENALTY COEFFICIENT.
        dCoef = dCoef * 1000.0;
        pPenalty->updateCoefficient (dCoef);
    }
    cout << "Passed - Linf smooth eq test" << endl;

    delete pPenalty;
    return;
}


//--------------------------------------------------------------------
//  Internal Function testLinfSmooth_ineq
//--------------------------------------------------------------------
// Test Linf smoothed penalty for two simple inequality constraints:
//   c1(x,y) = x - 1 >= 0  and  c2(x,y) = y - 2 >= 0
// Test at points (-1, 1), (0,1), (1, 1), (2, 1), (3, 0), (2, 3), (1, 2).
// Primary concern is robustness for small values of alpha.
static void  testLinfSmooth_ineq (void)
{
    using HOPSPACK::NonlConstrPenalty;
    using HOPSPACK::Vector;
    using HOPSPACK::isDoubleValid;

    //---- VALUES OF CONSTRAINTS AND EXPECTED RESULT AT TEST POINTS.
    double  daC1Values[7] = { -2.0, -1.0,  0.0,  1.0,  2.0,  1.0,  0.0 };
    double  daC2Values[7] = { -1.0, -1.0, -1.0, -1.0, -2.0,  1.0,  0.0 };

    double  dCoef = 1.0;
    double  dAlpha = 1.0;
    NonlConstrPenalty *  pPenalty = new NonlConstrPenalty();
    pPenalty->defineFunction ("L_inf (smoothed)", dCoef, dAlpha);

    if (bDISPLAY_LINFSM_INEQ)
        cout << " Detailed results, Linf smoothed, for inequalities" << endl;

    for (int  nCoefCase = 0; nCoefCase < 3; nCoefCase++)
    {
        pPenalty->updateCoefficient (dCoef);
        double  dAlpha = 1.0;
        pPenalty->updateSmoothing (dAlpha);

        for (int  nAlphaCase = 0; nAlphaCase < 20; nAlphaCase++)
        {
            double  dMaxErr;
            double  dAvgErr;
            computeFitToLinf (*pPenalty, daC1Values, daC2Values, false, 7,
                              &dMaxErr, &dAvgErr);
            if (bDISPLAY_LINFSM_INEQ)
            {
                streamsize  nCurrPrec = cout.precision (1);
                cout.setf (ios::scientific | ios::right);
                cout << "  r=" << dCoef << ", a=" << dAlpha;
                cout.precision (3);
                cout << " :  MaxErr = " << dMaxErr
                     << "  AvgErr = " << dAvgErr << endl;
                cout.precision (nCurrPrec);
            }

            if (isDoubleValid (dMaxErr) == false)
            {
                cout << "FAILED 1: rho=" << dCoef << ", alpha=" << dAlpha
                     << " - Linf smooth ineq test" << endl;
                throw "Exception from testLinfSmooth_ineq";
            }
            //---- EMPIRICAL TESTS.
            if (dMaxErr > (2.0 * dAlpha))
            {
                cout << "FAILED 2: rho=" << dCoef << ", alpha=" << dAlpha
                     << " - Linf smooth ineq test" << endl;
                throw "Exception from testLinfSmooth_ineq";
            }
            if ((dAlpha < 1.0) && (dAvgErr > (0.5 * dAlpha)))
            {
                cout << "FAILED 3: rho=" << dCoef << ", alpha=" << dAlpha
                     << " - Linf smooth ineq test" << endl;
                throw "Exception from testLinfSmooth_ineq";
            }

            //---- UPDATE TO THE NEXT ALPHA.
            dAlpha = dAlpha / 10.0;
            pPenalty->updateSmoothing (dAlpha);
        }
        if (bDISPLAY_LINFSM_INEQ)
            cout << endl;

        //---- UPDATE TO THE NEXT PENALTY COEFFICIENT.
        dCoef = dCoef * 1000.0;
        pPenalty->updateCoefficient (dCoef);
    }
    cout << "Passed - Linf smooth eq test" << endl;

    delete pPenalty;
    return;
}


//--------------------------------------------------------------------
//! Main routine for test executable.
/*!
 *  Return 0 if successful, nonzero value if not (for CTest automatic testing).
 *  Brief results for each test are printed to the screen.
 */
//--------------------------------------------------------------------
int  main (void)
{
    try
    {
        testL1_eq();
        testL1_ineq();
        testL1smooth_eq();
        testL1smooth_ineq();
        testLinf_eq();
        testLinf_ineq();
        testLinfSmooth_eq();
        testLinfSmooth_ineq();
    }
    catch (...)
    {
        cout << "EXITING: Exception caught from a test." << endl;
        return( 1 );
    }

    return( 0 );
}
