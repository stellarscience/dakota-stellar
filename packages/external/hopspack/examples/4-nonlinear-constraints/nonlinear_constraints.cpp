// $Id: nonlinear_constraints.cpp 166 2010-03-22 19:58:07Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/examples/4-nonlinear-constraints/nonlinear_constraints.cpp $

/* @HEADER
 * ************************************************************************
 *
 *         HOPSPACK: Hybrid Optimization Parallel Search Package
 *                 Copyright 2009-2010 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * This file is part of HOPSPACK.
 *
 * HOPSPACK is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library.  If not, see http://www.gnu.org/licenses/.
 *
 * Questions? Contact Tammy Kolda (tgkolda@sandia.gov)
 *                 or Todd Plantenga (tplante@sandia.gov)
 *
 * ************************************************************************
 * @HEADER
 */

/*!
  @file nonlinear_constraints.cpp
  @brief Evaluates function and nonlinear constraint values for example 4.

  This example provides a C++ language evaluation program for an optimization
  problem with a maximization goal and one nonlinear inequality constraint.
  The problem is:

    maximize    f(x,y) = x + 2y

    subject to  x + y >= 1
                1 - (x - 1)^2 - (y - 1)^2 >= 0

  The problem constrains points to lie within a circle centered at (1,1).
  HOPSPACK expects nonlinear inequalities to return a negative value if in
  violation of the constraint.  A nonnegative value means there is no violation.

  The exact solution point is x = 1.44721 (1 +  sqrt(5)/5),
                              y = 1.89443 (1 + 2sqrt(5)/5),
                        where f = 5.23607
  Only the nonlinear inequality is active at the solution, with sensitivity
  value  1.11803 (sqrt(5)/2).

  The input file has format:
    FC
    2
    1.000e-01
    2.000e-01
  Linear constraints and bounds are defined in the parameters file.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


// --------------------------------------------------------------------
//  Internal Function evaluateF
// --------------------------------------------------------------------
static double  evaluateF (const vector< double > &  x)
{
    double  f = x[0] + (2 * x[1]);
    return( f );
}


// --------------------------------------------------------------------
//  Internal Function evaluateCineqs
// --------------------------------------------------------------------
static void  evaluateCineqs (const vector< double > &  x,
                                   vector< double > &  cInequalities)
{
    cInequalities.resize (1);
    cInequalities[0] = 1 - (x[0] - 1)*(x[0] - 1) - (x[1] - 1)*(x[1] - 1);
    return;
}


// --------------------------------------------------------------------
//  Internal Function readInputFile
// --------------------------------------------------------------------
static bool  readInputFile (const string &            szExeName,
                            const string &            szFileName,
                                  vector< double > &  x)
{
    const int  REQBUFLEN = 10;
    char       reqBuf[REQBUFLEN + 1];

    ifstream  fp;

    //---- OPEN THE INPUT FILE.
    fp.open (szFileName.c_str());
    if (!fp)
    {
        cerr << szExeName << " - ERROR opening input file '"
             << szFileName << "'." << endl;
        return( false );
    }

    //---- READ AND VERIFY THE INPUT REQUEST TYPE.
    fp >> reqBuf;
    for (int  i = 0; i < REQBUFLEN + 1; i++)
    {
        if ((reqBuf[i] == '\n') || (reqBuf[i] == '\r'))
        {
            reqBuf[i] = 0;
            break;
        }
    }
    if ((reqBuf[0] != 'F') || (reqBuf[1] != 'C') || (reqBuf[2] != 0))
    {
        cerr << szExeName << " - ERROR reading request type from '"
             << szFileName << "'." << endl;
        cerr << "  Read '" << reqBuf << "'" << endl;
        fp.close();
        return( false );
    }

    //---- READ THE LENGTH OF x.
    int  n;
    fp >> n;
    if (n != 2)
    {
        cerr << szExeName << " - ERROR reading n from '"
             << szFileName << "'." << endl;
        cerr << "  Read " << n << ", but problem size should be 2" << endl;
        fp.close();
        return( false );
    }

    x.resize (n);

    //---- READ x.
    for (int  i = 0; i < n; i++)
        fp >> x[i];

    fp.close();

    //---- RETURN SUCCESS.
    return( true );
}


// --------------------------------------------------------------------
//  Internal Function writeOutputFile
// --------------------------------------------------------------------
static bool  writeOutputFile (const string           &  szExeName,
                              const string           &  szFileName,
                              const double              f,
                              const vector< double > &  cInequalities)
{
    ofstream  fp;

    fp.open (szFileName.c_str(), ios::trunc);
    if (!fp)
    {
        cerr << szExeName << " - ERROR opening output file '"
             << szFileName << "'." << endl;
        return( false );
    }

    //-- WRITE ALL DIGITS USING SCIENTIFIC NOTATION.
    fp.setf (ios::scientific);
    fp.precision (15);

    //---- WRITE THE NUMBER OF OBJECTIVES AND THEIR VALUES TO THE OUTPUT.
    fp << "1" << endl;
    fp << f << endl;

    //---- WRITE THE NUMBER OF NONLINEAR EQUALITIES.
    fp << "0" << endl;

    //---- WRITE THE NUMBER OF NONLINEAR INEQUALITIES AND THEIR VALUES.
    fp << cInequalities.size() << endl;
    for (int  i = 0; i < (int) cInequalities.size(); i++)
        fp << cInequalities[i] << endl;

    fp.close();

    return( true );
}


/* --------------------------------------------------------------------
 *  Main routine for evaluation executable.
 *
 *  Each execution reads two variables from a file, evaluates the objective
 *  function and nonlinear constraints, and writes the function values
 *  to an output file.
 *  HOPSPACK generates the file names and decides when to invoke this
 *  executable.  If there is an error, either create no output file,
 *  or create an output file with an error message.
 *
 *  @param argc  Number of command line arguments.
 *  @param argv  Command line arguments (input and output file names).
 *  @return      0 if successful.
 * --------------------------------------------------------------------
 */
int  main (const int           argc,
           const char * const  argv[])
{
    vector< double >  x;
    double            f;
    vector< double >  cIneqs;


    //---- CHECK THE COMMAND LINE ARGUMENTS.
    //---- IN THIS EXAMPLE ONLY THE INPUT AND OUTPUT FILE NAMES ARE USED.
    if (argc != 5)
    {
        cerr << "usage: " << argv[0]
             << " <input file> <output file> <tag> <type>" << endl;
        return( -1 );
    }

    if (readInputFile (argv[0], argv[1], x) == false)
        return( -2 );

    f = evaluateF (x);
    evaluateCineqs (x, cIneqs);

    if (writeOutputFile (argv[0], argv[2], f, cIneqs) == false)
        return( -3 );

    return( 0 );
}
