// $Id: degen_linear_constraints.cpp 166 2010-03-22 19:58:07Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/examples/3-degen-linear-constraints/degen_linear_constraints.cpp $

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
  @file degen_linear_constraints.cpp
  @brief Evaluates function values for degen_linear_constraints example 3.

  This example is the same as example 2, except for an additional linear
  constraint that is degenerate (redundant) at the solution point.  Since the
  objectives are identical, this source code file is identical to the
  one for example 2.

  The file provides a C++ language evaluation program for an optimization
  problem with one objective, bound constraints, and linear constraints.
  The problem is:

    minimize  f(x[1],x[2],x[3],x[4]) = \sum (x[i] - 10)^2

    subject to  - x[1] - x[2] - x[3] - x[4] >= -10
                  x[1] - x[2] + x[3] - x[4] <=  -1
                 2x[1]        +2x[3]        <=   9
                 2x[1]        +2x[3] -7x[4] =    3
                -10 <= x[i] <= 10

  The exact solution point is (2.25, 4.6429, 2.25, 0.8571), where f = 232.416.
  All linear inequalities are active at the solution, but the 3rd inequality
  is weakly active.

  The input file has format:
    F
    4
    1.000e-01
    2.000e-01
    3.000e-01
    4.000e-01
  Linear constraints and bounds are defined in the parameters file.
  This program does not evaluate or check constraints, since the objective
  is well-defined for any point.
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
    double  f = 0.0;
    for (int  i = 0; i < (int) x.size(); i++)
        f += (x[i] - 10.0)*(x[i] - 10.0);

    return( f );
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
    if ((reqBuf[0] != 'F') || (reqBuf[1] != 0))
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
    if (n != 4)
    {
        cerr << szExeName << " - ERROR reading n from '"
             << szFileName << "'." << endl;
        cerr << "  Read " << n << ", but problem size should be 4" << endl;
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
static bool  writeOutputFile (const string &  szExeName,
                              const string &  szFileName,
                              const double    f)
{
    ofstream  fp;

    fp.open (szFileName.c_str(), ios::trunc);
    if (!fp)
    {
        cerr << szExeName << " - ERROR opening output file '"
             << szFileName << "'." << endl;
        return( false );
    }

    //---- WRITE THE NUMBER OF OBJECTIVES AND THEIR VALUES TO THE OUTPUT.
    fp << "1" << endl;
    fp.precision (15);     //-- WRITE ALL DIGITS
    fp << f << endl;

    fp.close();

    return( true );
}


/* --------------------------------------------------------------------
 *  Main routine for evaluation executable.
 *
 *  Each execution reads variables x[0] thru x[3] from a file, evaluates
 *  the objective function, and writes the function value to an output file.
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

    if (writeOutputFile (argv[0], argv[2], f) == false)
        return( -3 );

    return( 0 );
}
