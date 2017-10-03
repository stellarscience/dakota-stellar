// $Id: HOPSPACK_main_startpoint_tests.cpp 167 2010-03-24 00:07:00Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/test/startpoint-tests/HOPSPACK_main_startpoint_tests.cpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2010 Sandia Corporation
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
  @file HOPSPACK_main_startpoint_tests.cpp
  @brief Main program that tests start points with linear constraints.
*/

#include <math.h>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_Matrix.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_ProblemDef.hpp"
#include "HOPSPACK_Vector.hpp"


//---- FOR DETAILED DEBUGGING OUTPUT, SET bDEBUG = true IN
//---- src/src_shared/HOPSPACK_SolveLinConstrProj.cpp AND RECOMPILE.


//--------------------------------------------------------------------
//  Internal Function check_2var_result
//--------------------------------------------------------------------
static void  check_2var_result (const HOPSPACK::LinConstr &  linConstr,
                                      HOPSPACK::Vector    &  x,
                                const double                 dSol_0,
                                const double                 dSol_1,
                                const string              &  sTestName)
{
    if (linConstr.projectToFeasibility (x) == false)
    {
        cout << "FAILED 1 - " << sTestName << endl;
        throw ("Exception from " + sTestName);
    }

    if (   (fabs (x[0] - dSol_0) > 4.0e-15)
        || (fabs (x[1] - dSol_1) > 4.0e-15) )
    {
        double  dMaxErr = max (fabs (x[0] - dSol_0), fabs (x[1] - dSol_1));
        cout << "FAILED 2: x[0]=" << x[0] << ", x[1]=" << x[1]
             << ", maxerr=" << dMaxErr
             << " - " << sTestName << endl;
        throw ("Exception from " + sTestName);
    }

    cout << "Passed - " << sTestName << endl;
    return;
}


//--------------------------------------------------------------------
//  Internal Function test_bnds_1
//--------------------------------------------------------------------
// Test projection of a problem with variable bounds only.
static void  test_bnds_1 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vUpBnds (2, 10.0);
    probDefList.setParameter ("Upper Bounds", vUpBnds);
    Vector  vLoBnds (2, -10.0);
    probDefList.setParameter ("Lower Bounds", vLoBnds);
    Vector  vScaling(2);
    vScaling[0] = 1.0;
    vScaling[1] = 2.0;
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] =  35.0;
    x[1] = -11.7;
    check_2var_result (linConstr, x, 10.0, -10.0, "Bnds test 1");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_bnds_2
//--------------------------------------------------------------------
// Test projection of a problem with variable bounds only, and one variable
// fixed by the bounds.
static void  test_bnds_2 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vUpBnds (2, 10.0);
    probDefList.setParameter ("Upper Bounds", vUpBnds);
    Vector  vLoBnds (2);
    vLoBnds[0] = -10.0;
    vLoBnds[1] =  10.0;
    probDefList.setParameter ("Lower Bounds", vLoBnds);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] =  3.5;
    x[1] = -11.7;
    check_2var_result (linConstr, x, 3.5, 10.0, "Bnds test 2");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_eqs_1
//--------------------------------------------------------------------
// Test projection of a problem with one equality constraint only.
static void  test_eqs_1 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons(2);
    vCons[0] = 1.0;
    vCons[1] = 2.0;
    Matrix  mEqs;
    mEqs.addRow (vCons);
    linConstrList.setParameter ("Equality Matrix", mEqs);
    Vector  vRhs (1, 1.0);
    linConstrList.setParameter ("Equality Bounds", vRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 0.5;
    x[1] = 1.5;
    check_2var_result (linConstr, x, 0.0, 0.5, "Eqs test 1");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_eqs_2
//--------------------------------------------------------------------
// Test projection of a problem with two equality constraints only.
static void  test_eqs_2 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons1(2);
    vCons1[0] = 1.0;
    vCons1[1] = 2.0;
    Vector  vCons2(2);
    vCons2[0] =  1.0;
    vCons2[1] = -1.0;
    Matrix  mEqs;
    mEqs.addRow (vCons1);
    mEqs.addRow (vCons2);
    linConstrList.setParameter ("Equality Matrix", mEqs);
    Vector  vRhs(2);
    vRhs[0] = 1.0;
    vRhs[1] = 0.0;
    linConstrList.setParameter ("Equality Bounds", vRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 0.5;
    x[1] = 1.5;
    check_2var_result (linConstr, x, (1.0 / 3.0), (1.0 / 3.0), "Eqs test 2");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_eqs_3
//--------------------------------------------------------------------
// Test projection of a problem with inconsistent equality constraints.
static void  test_eqs_3 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    cout << "----- vvvvv --------------------" << endl;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 3);
    Vector  vScaling (3, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons1(3);
    vCons1[0] = 1.0;
    vCons1[1] = 2.0;
    vCons1[2] = 0.0;
    Vector  vCons2(3);
    vCons2[0] = 1.0;
    vCons2[1] = 0.0;
    vCons2[2] = 0.0;
    Vector  vCons3(3);
    vCons3[0] = 0.0;
    vCons3[1] = 1.0;
    vCons3[2] = 0.0;
    Matrix  mEqs;
    mEqs.addRow (vCons1);
    mEqs.addRow (vCons2);
    mEqs.addRow (vCons3);
    linConstrList.setParameter ("Equality Matrix", mEqs);
    Vector  vRhs(3);
    vRhs[0] = 1.0;
    vRhs[1] = 1.0;
    vRhs[2] = 1.0;
    linConstrList.setParameter ("Equality Bounds", vRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(3);
    x[0] = 2.0;
    x[1] = 3.0;
    x[2] = 4.0;
    //---- PROJECTION SHOULD FAIL BECAUSE CONSTRAINTS ARE INCONSISTENT.
    if (linConstr.projectToFeasibility (x) == true)
    {
        cout << "FAILED 1 - Eqs test 3" << endl;
        throw "Exception from Eqs test 3";
    }
    cout << "Passed - Eqs test 3 (errors from dgglse are expected)" << endl;
    cout << "----- ^^^^^ --------------------" << endl;

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_eqs_4
//--------------------------------------------------------------------
// Test projection of a problem with redundant equality constraints.
static void  test_eqs_4 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons1(2);
    vCons1[0] = 1.0;
    vCons1[1] = 2.0;
    Vector  vCons2(2);
    vCons2[0] = 1.0;
    vCons2[1] = 2.0;
    Matrix  mEqs;
    mEqs.addRow (vCons1);
    mEqs.addRow (vCons2);
    linConstrList.setParameter ("Equality Matrix", mEqs);
    Vector  vRhs(2);
    vRhs[0] = 1.0;
    vRhs[1] = 1.0;
    linConstrList.setParameter ("Equality Bounds", vRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 0.5;
    x[1] = 1.5;

    //---- SOME LAPACK ROUTINES RETURN NaN, OTHERS RETURN BAD NUMBERS.
    if (linConstr.projectToFeasibility (x) == false)
    {
        cout << "Passed - Eqs test 4 (errors from dgglse are expected)" << endl;
        return;
    }
    //---- PROJECTION SHOULD BE INCORRECT BECAUSE CONSTRAINTS ARE REDUNDANT.
    //---- UNFORTUNATELY, NO GOOD WAY TO DETECT THE ERROR.
    if (   (fabs (x[0] - 0.0) > 1.0e-15)
        || (fabs (x[1] - 0.5) > 1.0e-15) )
    {
        cout << "Passed - Eqs test 4 (incorrect projection is expected)" << endl;
        cout << "         x=(" << x[0] << ", " << x[1] << ")"
             << "  vs  (0.0, 0.5)" << endl;
    }
    else
    {
        cout << "Passed - Eqs test 4 (correct projection not expected!)" << endl;
    }

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_bndseqs_1
//--------------------------------------------------------------------
// Test projection of a problem with an equality constraint to the
// interior of the variable bounds.
static void  test_bndseqs_1 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vUpBnds (2, 5.0);
    probDefList.setParameter ("Upper Bounds", vUpBnds);
    Vector  vLoBnds (2, -5.0);
    probDefList.setParameter ("Lower Bounds", vLoBnds);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons(2);
    vCons[0] = 1.0;
    vCons[1] = 2.0;
    Matrix  mEqs;
    mEqs.addRow (vCons);
    linConstrList.setParameter ("Equality Matrix", mEqs);
    Vector  vRhs (1, 1.0);
    linConstrList.setParameter ("Equality Bounds", vRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 0.0;
    x[1] = 0.0;
    check_2var_result (linConstr, x, 0.2, 0.4, "BndsEqs test 1");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_bndseqs_2
//--------------------------------------------------------------------
// Test projection of a problem with one variable fixed by an equality to
// a variable bound.
static void  test_bndseqs_2 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vUpBnds (2, 5.0);
    probDefList.setParameter ("Upper Bounds", vUpBnds);
    Vector  vLoBnds (2, 3.0);
    probDefList.setParameter ("Lower Bounds", vLoBnds);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons (2);
    vCons[0] = 0.0;
    vCons[1] = 1.0;
    Matrix  mEqs;
    mEqs.addRow (vCons);
    linConstrList.setParameter ("Equality Matrix", mEqs);
    Vector  vRhs (1, 3.0);
    linConstrList.setParameter ("Equality Bounds", vRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 0.0;
    x[1] = 0.0;
    check_2var_result (linConstr, x, 3.0, 3.0, "BndsEqs test 2");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_bndseqs_3
//--------------------------------------------------------------------
// Test projection of a problem with one variable fixed by an equality,
// and both bounds active.  There is a redundant constraint, but the
// active set method never adds the extra constraint.
static void  test_bndseqs_3 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vUpBnds (2);
    vUpBnds[0] = -8.0;
    vUpBnds[1] = -9.0;
    probDefList.setParameter ("Upper Bounds", vUpBnds);
    Vector  vLoBnds (2);
    vLoBnds[0] = -9.0;
    vLoBnds[1] = -9.0;
    probDefList.setParameter ("Lower Bounds", vLoBnds);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons (2);
    vCons[0] = 0.0;
    vCons[1] = 1.0;
    Matrix  mEqs;
    mEqs.addRow (vCons);
    linConstrList.setParameter ("Equality Matrix", mEqs);
    Vector  vRhs (1, -9.0);
    linConstrList.setParameter ("Equality Bounds", vRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 3.0;
    x[1] = 0.0;
    check_2var_result (linConstr, x, -8.0, -9.0, "BndsEqs test 3");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_bndseqs_4
//--------------------------------------------------------------------
// Test projection of a problem with an equality constraint outside
// the variable bounds:
//     x + 0.5y = 3, 0 <= x <= 1, 0 <= y <= 1
//
// In this example the active set method moves to the [1,1] corner of
// the bounds box and fails.  If the equality is changed to  y = 3,
// then the method makes y = 1 an active inequality and fails because
// subproblem equalities are dependent.
static void  test_bndseqs_4 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    cout << "----- vvvvv --------------------" << endl;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vUpBnds (2, 1.0);
    probDefList.setParameter ("Upper Bounds", vUpBnds);
    Vector  vLoBnds (2, 0.0);
    probDefList.setParameter ("Lower Bounds", vLoBnds);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons(2);
    vCons[0] = 0.5;
    vCons[1] = 1.0;
    Matrix  mEqs;
    mEqs.addRow (vCons);
    linConstrList.setParameter ("Equality Matrix", mEqs);
    Vector  vRhs (1, 3.0);
    linConstrList.setParameter ("Equality Bounds", vRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 0.0;
    x[1] = 0.0;
    //---- PROJECTION SHOULD FAIL BECAUSE CONSTRAINTS ARE INCONSISTENT.
    if (linConstr.projectToFeasibility (x) == true)
    {
        cout << "FAILED 1 - BndsEqs test 4" << endl;
        throw "Exception from BndsEqs test 4";
    }
    cout << "Passed - BndsEqs test 4 (error 'overdetermined system' is expected)"
         << endl;
    cout << "----- ^^^^^ --------------------" << endl;

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_ineqs_1
//--------------------------------------------------------------------
// Test projection of a problem to single-sided and double-sided
// inequality constraints that are all inactive.  The initial point
// is feasible and should not move.
static void  test_ineqs_1 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons1(2);
    vCons1[0] = 2.0;
    vCons1[1] = 0.0;
    Vector  vCons2(2);
    vCons2[0] = 2.0;
    vCons2[1] = 2.0;
    Matrix  mIneqs;
    mIneqs.addRow (vCons1);
    mIneqs.addRow (vCons2);
    linConstrList.setParameter ("Inequality Matrix", mIneqs);
    Vector  vLoRhs (2);
    vLoRhs[0] = 0.0;
    vLoRhs[1] = 0.0;
    linConstrList.setParameter ("Inequality Lower", vLoRhs);
    Vector  vUpRhs (2);
    vUpRhs[0] = HOPSPACK::dne();
    vUpRhs[1] = 6.0;
    linConstrList.setParameter ("Inequality Upper", vUpRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 1.0;
    x[1] = 2.0;
    check_2var_result (linConstr, x, 1.0, 2.0, "Ineqs test 1");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_ineqs_2
//--------------------------------------------------------------------
// Test projection of a problem to single-sided inequality constraints
// with an active lower bound but inactive upper bound.
//   x+y >= 1  x+2y <= 2  starting from (0,0)
// answer is (0.5,0.5)
static void  test_ineqs_2 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons1(2);
    vCons1[0] = 1.0;
    vCons1[1] = 1.0;
    Vector  vCons2(2);
    vCons2[0] = 1.0;
    vCons2[1] = 2.0;
    Matrix  mIneqs;
    mIneqs.addRow (vCons1);
    mIneqs.addRow (vCons2);
    linConstrList.setParameter ("Inequality Matrix", mIneqs);
    Vector  vLoRhs (2);
    vLoRhs[0] = 1.0;
    vLoRhs[1] = HOPSPACK::dne();
    linConstrList.setParameter ("Inequality Lower", vLoRhs);
    Vector  vUpRhs (2);
    vUpRhs[0] = HOPSPACK::dne();
    vUpRhs[1] = 2.0;
    linConstrList.setParameter ("Inequality Upper", vUpRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 0.0;
    x[1] = 0.0;
    check_2var_result (linConstr, x, 0.5, 0.5, "Ineqs test 2");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_ineqs_3
//--------------------------------------------------------------------
// Test projection of a problem to single-sided inequality constraints
// with an active upper bound but inactive lower bound.
//   x+y >= 1  x+2y <= 2  starting from (2,2.5)
// answer is (1,0.5)
static void  test_ineqs_3 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons1(2);
    vCons1[0] = 1.0;
    vCons1[1] = 1.0;
    Vector  vCons2(2);
    vCons2[0] = 1.0;
    vCons2[1] = 2.0;
    Matrix  mIneqs;
    mIneqs.addRow (vCons1);
    mIneqs.addRow (vCons2);
    linConstrList.setParameter ("Inequality Matrix", mIneqs);
    Vector  vLoRhs (2);
    vLoRhs[0] = 1.0;
    vLoRhs[1] = HOPSPACK::dne();
    linConstrList.setParameter ("Inequality Lower", vLoRhs);
    Vector  vUpRhs (2);
    vUpRhs[0] = HOPSPACK::dne();
    vUpRhs[1] = 2.0;
    linConstrList.setParameter ("Inequality Upper", vUpRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 2.0;
    x[1] = 2.5;
    check_2var_result (linConstr, x, 1.0, 0.5, "Ineqs test 3");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_ineqs_4
//--------------------------------------------------------------------
// Test projection of a problem to double-sided inequality constraints.
//   2 >= x+y >= 1  starting from (3,3)
// answer is (1,1)
static void  test_ineqs_4 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons(2);
    vCons[0] = 1.0;
    vCons[1] = 1.0;
    Matrix  mIneqs;
    mIneqs.addRow (vCons);
    linConstrList.setParameter ("Inequality Matrix", mIneqs);
    Vector  vLoRhs (1);
    vLoRhs[0] = 1.0;
    linConstrList.setParameter ("Inequality Lower", vLoRhs);
    Vector  vUpRhs (1);
    vUpRhs[0] = 2.0;
    linConstrList.setParameter ("Inequality Upper", vUpRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 3.0;
    x[1] = 3.0;
    check_2var_result (linConstr, x, 1.0, 1.0, "Ineqs test 4");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_ineqs_5
//--------------------------------------------------------------------
// Test projection of an infeasible problem with parallel constraints.
//   1 >= x+y  x+y >= 2
// An active set subproblem receives both constraints and is therefore
// hopelessly ill-conditioned.
static void  test_ineqs_5 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    cout << "----- vvvvv --------------------" << endl;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons(2);
    vCons[0] = 1.0;
    vCons[1] = 1.0;
    Matrix  mIneqs;
    mIneqs.addRow (vCons);
    mIneqs.addRow (vCons);
    linConstrList.setParameter ("Inequality Matrix", mIneqs);
    Vector  vLoRhs (2);
    vLoRhs[0] = HOPSPACK::dne();
    vLoRhs[1] = 2.0;
    linConstrList.setParameter ("Inequality Lower", vLoRhs);
    Vector  vUpRhs (2);
    vUpRhs[0] = 1.0;
    vUpRhs[1] = HOPSPACK::dne();
    linConstrList.setParameter ("Inequality Upper", vUpRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 0.0;
    x[1] = 0.0;
    //---- PROJECTION SHOULD FAIL BECAUSE CONSTRAINTS ARE INCONSISTENT.
    if (linConstr.projectToFeasibility (x) == true)
    {
        cout << "FAILED 1 - Ineqs test 5" << endl;
        throw "Exception from Ineqs test 5";
    }
    cout << "Passed - Ineqs test 5 (error 'no feasible start point' is expected)"
         << endl;
    cout << "----- ^^^^^ --------------------" << endl;

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_ineqs_6
//--------------------------------------------------------------------
// Test projection of an infeasible problem with nonparallel constraints.
//   x <= 0  y <= 0  x+y >= 1
static void  test_ineqs_6 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    cout << "----- vvvvv --------------------" << endl;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons1(2);
    vCons1[0] = 1.0;
    vCons1[1] = 0.0;
    Vector  vCons2(2);
    vCons2[0] = 0.0;
    vCons2[1] = 1.0;
    Vector  vCons3(2);
    vCons3[0] = 1.0;
    vCons3[1] = 1.0;
    Matrix  mIneqs;
    mIneqs.addRow (vCons1);
    mIneqs.addRow (vCons2);
    mIneqs.addRow (vCons3);
    linConstrList.setParameter ("Inequality Matrix", mIneqs);
    Vector  vLoRhs (3);
    vLoRhs[0] = HOPSPACK::dne();
    vLoRhs[1] = HOPSPACK::dne();
    vLoRhs[2] = 1.0;
    linConstrList.setParameter ("Inequality Lower", vLoRhs);
    Vector  vUpRhs (3);
    vUpRhs[0] = 0.0;
    vUpRhs[1] = 0.0;
    vUpRhs[2] = HOPSPACK::dne();
    linConstrList.setParameter ("Inequality Upper", vUpRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 0.0;
    x[1] = 0.0;
    //---- PROJECTION SHOULD FAIL BECAUSE CONSTRAINTS ARE INCONSISTENT.
    if (linConstr.projectToFeasibility (x) == true)
    {
        cout << "FAILED 1 - Ineqs test 6" << endl;
        throw "Exception from Ineqs test 6";
    }
    cout << "Passed - Ineqs test 6 (error 'constraints may be inconsistent' is expected)"
         << endl;
    cout << "----- ^^^^^ --------------------" << endl;

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_ineqs_7
//--------------------------------------------------------------------
// Test projection that goes to one inequality but has to drop it to
// reach another more binding one.
//   x + y >= 1  x + 0.5y >= 1
static void  test_ineqs_7 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons1(2);
    vCons1[0] = 1.0;
    vCons1[1] = 1.0;
    Vector  vCons2(2);
    vCons2[0] = 1.0;
    vCons2[1] = 0.5;
    Matrix  mIneqs;
    mIneqs.addRow (vCons1);
    mIneqs.addRow (vCons2);
    linConstrList.setParameter ("Inequality Matrix", mIneqs);
    Vector  vLoRhs (2);
    vLoRhs[0] = 1.0;
    vLoRhs[1] = 1.0;
    linConstrList.setParameter ("Inequality Lower", vLoRhs);
    Vector  vUpRhs (2);
    vUpRhs[0] = HOPSPACK::dne();
    vUpRhs[1] = HOPSPACK::dne();
    linConstrList.setParameter ("Inequality Upper", vUpRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = -0.1;
    x[1] =  1.0;
    check_2var_result (linConstr, x, 0.38, 1.24, "Ineqs test 7");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_ineqs_8
//--------------------------------------------------------------------
// Test projection that goes to one inequality but has to drop it to
// reach another more binding one, and the two are parallel.
//   x + y >= 1  x + y >= 2
static void  test_ineqs_8 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons1(2);
    vCons1[0] = 1.0;
    vCons1[1] = 1.0;
    Vector  vCons2(2);
    vCons2[0] = 1.0;
    vCons2[1] = 1.0;
    Matrix  mIneqs;
    mIneqs.addRow (vCons1);
    mIneqs.addRow (vCons2);
    linConstrList.setParameter ("Inequality Matrix", mIneqs);
    Vector  vLoRhs (2);
    vLoRhs[0] = 1.0;
    vLoRhs[1] = 2.0;
    linConstrList.setParameter ("Inequality Lower", vLoRhs);
    Vector  vUpRhs (2);
    vUpRhs[0] = HOPSPACK::dne();
    vUpRhs[1] = HOPSPACK::dne();
    linConstrList.setParameter ("Inequality Upper", vUpRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 0.0;
    x[1] = 0.0;
    check_2var_result (linConstr, x, 1.0, 1.0, "Ineqs test 8");

    return;
}


//--------------------------------------------------------------------
//  Internal Function test_ineqs_9
//--------------------------------------------------------------------
// Test projection to a pair of inequalities that act like an equality.
//   x + y >= 1  x + y <= 1
static void  test_ineqs_9 (void)
{
    using HOPSPACK::LinConstr;
    using HOPSPACK::Matrix;
    using HOPSPACK::ParameterList;
    using HOPSPACK::ProblemDef;
    using HOPSPACK::Vector;

    ParameterList  probDefList;
    probDefList.setParameter ("Number Unknowns", 2);
    Vector  vScaling (2, 1.0);
    probDefList.setParameter ("Scaling", vScaling);
    ProblemDef  probDef;
    probDef.initialize (probDefList);

    ParameterList  linConstrList;
    Vector  vCons1(2);
    vCons1[0] = 1.0;
    vCons1[1] = 1.0;
    Vector  vCons2(2);
    vCons2[0] = 1.0;
    vCons2[1] = 1.0;
    Matrix  mIneqs;
    mIneqs.addRow (vCons1);
    mIneqs.addRow (vCons2);
    linConstrList.setParameter ("Inequality Matrix", mIneqs);
    Vector  vLoRhs (2);
    vLoRhs[0] = 1.0;
    vLoRhs[1] = HOPSPACK::dne();
    linConstrList.setParameter ("Inequality Lower", vLoRhs);
    Vector  vUpRhs (2);
    vUpRhs[0] = HOPSPACK::dne();
    vUpRhs[1] = 1.0;
    linConstrList.setParameter ("Inequality Upper", vUpRhs);
    LinConstr  linConstr (probDef);
    linConstr.initialize (linConstrList);

    Vector  x(2);
    x[0] = 0.8;
    x[1] = 0.0;
    check_2var_result (linConstr, x, 0.9, 0.1, "Ineqs test 9");

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
        test_bnds_1();
        test_bnds_2();
        test_eqs_1();
        test_eqs_2();
        test_eqs_3();
        test_eqs_4();
        test_bndseqs_1();
        test_bndseqs_2();
        test_bndseqs_3();
        test_bndseqs_4();
        test_ineqs_1();
        test_ineqs_2();
        test_ineqs_3();
        test_ineqs_4();
        test_ineqs_5();
        test_ineqs_6();
        test_ineqs_7();
        test_ineqs_8();
        test_ineqs_9();
    }
    catch (...)
    {
        cout << "EXITING: Exception caught from a test." << endl;
        return( 1 );
    }

    return( 0 );
}
