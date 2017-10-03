// $Id: HOPSPACK_main_lapack_tests.cpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/test/lapack-tests/HOPSPACK_main_lapack_tests.cpp $

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
  @file HOPSPACK_main_lapack_tests.cpp
  @brief Main program that tests LAPACK functionality for HOPSPACK.
*/

#include <math.h>     //-- FOR fabs

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_LapackWrappers.hpp"
#include "HOPSPACK_Matrix.hpp"
#include "HOPSPACK_Vector.hpp"


//--------------------------------------------------------------------
//  Internal Function testDdot
//--------------------------------------------------------------------
// Directly test ddot in LAPACK.
static void  testDdot (HOPSPACK::LapackWrappers &  lapack)
{
    double  daX[3] = {0.0, 1.0, 2.0};

    double  dResult = lapack.ddot (3, daX, daX);
    if (dResult != 5.0)
    {
        cout << "FAILED - " << dResult << " instead of 5.0"
             << " - ddot test" << endl;
        throw "Exception from testDgemm";
    }
    cout << "Passed - ddot test" << endl;

    return;
}


//--------------------------------------------------------------------
//  Internal Function testDgemv
//--------------------------------------------------------------------
// Indirectly test dgemv with a matrix-vector multiply.
static void  testDgemv (HOPSPACK::LapackWrappers &  lapack)
{
    using HOPSPACK::Matrix;
    using HOPSPACK::Vector;

    double  daX[3] = {1.0, 2.0, 3.0};
    double  daY[2] = {1.0, 1.0};
    double  daA_row1[3] = { 6.0, 5.0, 4.0 };
    double  daA_row2[3] = { 3.0, 2.0, 1.0 };
    Matrix  A;
    Vector  aRow1 (3, daA_row1);
    Vector  aRow2 (3, daA_row2);
    A.addRow (aRow1);
    A.addRow (aRow2);
    Vector  x (3, daX);

    Vector  y (2, daY);
    A.multVec (x, y, Matrix::NO_TRANSPOSE);
    if ((y[0] == 28.0) && (y[1] == 10.0))
        cout << "Passed - dgemv test 'N'" << endl;
    else
    {
        cout << "FAILED - " << y[0] << "," << y[1]
             << " instead of 28,10"
             << " - dgemv test 'N'" << endl;
        throw "Exception from testDgemm";
    }

    y[0] = 1.0;
    y[1] = 1.0;
    A.multVec (y, x, Matrix::TRANSPOSE);
    if ((x[0] == 9.0) && (x[1] == 7.0) && (x[2] == 5.0))
        cout << "Passed - dgemv test 'T'" << endl;
    else
    {
        cout << "FAILED - " << x[0] << "," << x[1] << "," << x[2]
             << " instead of 9,7,5"
             << " - dgemv test 'T'" << endl;
        throw "Exception from testDgemm";
    }

    return;
}


//--------------------------------------------------------------------
//  Internal Function testDgemm
//--------------------------------------------------------------------
// Indirectly test dgemm with a matrix-matrix multiply.
static void  testDgemm (HOPSPACK::LapackWrappers &  lapack)
{
    using HOPSPACK::Matrix;
    using HOPSPACK::Vector;

    double  daA_row1[3] = { 6.0, 5.0, 4.0 };
    double  daA_row2[3] = { 3.0, 2.0, 1.0 };
    Matrix  A;
    Vector  aRow1 (3, daA_row1);
    Vector  aRow2 (3, daA_row2);
    A.addRow (aRow1);
    A.addRow (aRow2);

    double  daB_row1[2] = {  1.0,  0.0 };
    double  daB_row2[2] = {  2.0,  1.0 };
    double  daB_row3[2] = { -1.0, -2.0 };
    Matrix  B;
    Vector  bRow1 (2, daB_row1);
    Vector  bRow2 (2, daB_row2);
    Vector  bRow3 (2, daB_row3);
    B.addRow (bRow1);
    B.addRow (bRow2);
    B.addRow (bRow3);

    Matrix  C;
    A.multMat (B, C, Matrix::NO_TRANSPOSE);
    if ((C.getNrows() != 2) || (C.getNcols() != 2))
    {
        cout << "FAILED - wrong matrix size returned"
             << " - dgemm test 'N'" << endl;
        throw "Exception from testDgemm";
    }
    if (   (C.getRow(0)[0] == 12.0) && (C.getRow(0)[1] == -3.0)
        && (C.getRow(1)[0] ==  6.0) && (C.getRow(1)[1] ==  0.0))
    {
        cout << "Passed - dgemm test 'N'" << endl;
    }
    else
    {
        cout << "FAILED - " << C.getRow(0)[0] << "," << C.getRow(0)[1]
             << C.getRow(1)[0] << "," << C.getRow(1)[1]
             << " instead of 12,-3,6,0"
             << " - dgemm test 'N'" << endl;
        throw "Exception from testDgemm";
    }

    //---- NOW TEST THE TRANSPOSE OPTION.

    double  daBtrans_row1[3] = { -1.0, 2.0, 2.0 };
    Matrix  Btrans;
    Vector  btransRow1 (3, daBtrans_row1);
    Btrans.addRow (btransRow1);
    Matrix  Ctest2;
    A.multMat (Btrans, C, Matrix::TRANSPOSE);
    if ((C.getNrows() != 2) || (C.getNcols() != 1))
    {
        cout << "FAILED - wrong matrix size returned"
             << " - dgemm test 'T'" << endl;
        throw "Exception from testDgemm";
    }
    if ((C.getRow(0)[0] == 12.0) && (C.getRow(1)[0] == 3.0))
    {
        cout << "Passed - dgemm test 'T'" << endl;
    }
    else
    {
        cout << "FAILED - " << C.getRow(0)[0] << "," << C.getRow(1)[0]
             << " instead of 12,3"
             << " - dgemm test 'T'" << endl;
        throw "Exception from testDgemm";
    }

    return;
}


//--------------------------------------------------------------------
//  Internal Function testDgesvd
//--------------------------------------------------------------------
// Indirectly test dgesvd with a call to Matrix::svd.
static void  testDgesvd (HOPSPACK::LapackWrappers &  lapack)
{
    using HOPSPACK::Matrix;
    using HOPSPACK::Vector;

    double  daA_row1[3] = { 3.0, 1.0, 0.0 };
    double  daA_row2[3] = { 0.0, 2.0, 0.0 };
    Matrix  A;
    Vector  aRow1 (3, daA_row1);
    Vector  aRow2 (3, daA_row2);
    A.addRow (aRow1);
    A.addRow (aRow2);

    Matrix  U;
    Matrix  VT;
    Vector  Sigma;

    A.svd (U, Sigma, VT);
    if ((U.getNrows() != A.getNrows()) || (U.getNcols() != A.getNrows()))
    {
        cout << "FAILED - wrong matrix size returned"
             << " - dgesvd matrix U" << endl;
        throw "Exception from testDgesvd";
    }
    if ((VT.getNrows() != A.getNcols()) || (VT.getNcols() != A.getNcols()))
    {
        cout << "FAILED - wrong matrix size returned"
             << " - dgesvd matrix VT" << endl;
        throw "Exception from testDgesvd";
    }

    //---- CHECK THE ANSWER BY MULTIPLYING IT OUT.
    Matrix  SigmaTimesVT;
    for (int  i = 0; i < Sigma.size(); i++)
    {
        Vector  tmpRow;
        for (int  j = 0; j < Sigma.size(); j++)
            tmpRow.push_back (Sigma[i] * (VT.getRow(i))[j]);
        for (int  j = Sigma.size(); j < VT.getNrows(); j++)
            tmpRow.push_back (0.0);
        SigmaTimesVT.addRow (tmpRow);
    }
    Matrix  USVT;
    U.multMat (SigmaTimesVT, USVT);
    for (int  i = 0; i < A.getNrows(); i++)
    {
        const Vector &  nextArow = A.getRow (i);
        const Vector &  nextResRow = USVT.getRow (i);
        for (int  j = 0; j < nextArow.size(); j++)
        {
            if (fabs (nextArow[j] - nextResRow[j]) > 1.0e-12)
            {
                cout << "FAILED - reconstructed element [" << i << "][" << j
                     << "] differs by too much: " << nextArow[j]
                     << " vs " << nextResRow[j]
                     << " - dgesvd test" << endl;
                throw "Exception from testDgesvd";
            }
        }
    }
    cout << "Passed - dgesvd test" << endl;

    return;
}


//--------------------------------------------------------------------
//  Internal Function testDgglse
//--------------------------------------------------------------------
// Directly test dgglse in LAPACK.
/*
 * The problem is
 *    min || Ax - c ||_2  s.t.  Bx = d
 */
static void  testDgglse (HOPSPACK::LapackWrappers &  lapack)
{
    int  m = 2;
    int  n = 2;
    int  p = 1;

    double  daA[4] = {1.0, 0.0, -1.0, 1.0};     //-- m x n
    double  daB[2] = {1.0, 0.0};                //-- p x n
    double  daC[2] = {2.0, 0.0};                //-- n x 1
    double  daD[1] = {4.0};                     //-- p x 1
    double  daX[2] = {0.0, 0.0};                //-- n x 1 SOLUTION

    bool  bSuccess = lapack.dgglse (m, n, p, daA, daB, daC, daD, daX);
    if (bSuccess == false)
    {
        cout << "FAILED - call to dgglse returned false - dgglse test" << endl;
        throw "Exception from testDgglse";
    }
    if ((daX[0] != 4.0) && (daX[0] != 1.0))
    {
        cout << "FAILED - " << daX[0] << "," << daX[1] << " instead of 4,1"
             << " - dgglse test" << endl;
        throw "Exception from testDgglse";
    }
    cout << "Passed - dgglse test" << endl;

    return;
}


//--------------------------------------------------------------------
//  Internal Function testDepRows
//--------------------------------------------------------------------
// Indirectly test dgelqf in LAPACK.
static void  testDepRows (HOPSPACK::LapackWrappers &  lapack)
{
    using HOPSPACK::Matrix;
    using HOPSPACK::Vector;

    //---- THIS SIMPLE CASE DEMONSTRATES
    //----   |  1 0        |
    //----   |  1 10^(-k)  |
    //---- IS DEPENDENT IF THE THRESHOLD IS BIGGER THAN 10^(-k).
    double  daA_row1[2] = { 1.0, 0.0    };
    double  daA_row2[2] = { 1.0, 1.0e-8 };
    Matrix  A;
    Vector  aRow1 (2, daA_row1);
    Vector  aRow2 (2, daA_row2);
    A.addRow (aRow1);
    A.addRow (aRow2);
    Vector  b (2, 1.0);
    A.pruneDependentRows (b, 1.0e-9);
    if (A.getNrows() != 2)
    {
        cout << "FAILED - should not have pruned a row - dgelqf test" << endl;
        throw "Exception from testDepRows";
    }
    A.pruneDependentRows (b, 1.0e-7);
    if (A.getNrows() != 1)
    {
        cout << "FAILED - should have pruned a row - dgelqf test" << endl;
        throw "Exception from testDepRows";
    }
    cout << "Passed - dgelqf test" << endl;

    return;
}


//--------------------------------------------------------------------
//! Main routine for test executable.
/*!
 *  Return 0 if successful, nonzero value if not (for CTest automatic testing).
 *  Brief results for each test are printed to the screen.
 *  Failure usually means something is wrong with the LAPACK library.
 */
//--------------------------------------------------------------------
int  main (void)
{
    using HOPSPACK::LapackWrappers;

    try
    {
        LapackWrappers &  lapack = LapackWrappers::getTheInstance();

        testDdot (lapack);
        testDgemv (lapack);
        testDgemm (lapack);
        testDgesvd (lapack);
        testDgglse (lapack);
        testDepRows (lapack);
    }
    catch (...)
    {
        cout << "EXITING: Exception caught from a test." << endl;
        return( 1 );
    }

    return( 0 );
}
