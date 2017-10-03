// $Id: multi_start.cpp 166 2010-03-22 19:58:07Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/examples/5-multi-start/multi_start.cpp $

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
  @file multi_start.cpp
  @brief Evaluates function values for multi-start example 5.

  This example provides a C++ language evaluation program for a standard
  global optimization problem known as the "six-humped camel back problem".   
  The problem is:

    minimize    f(x,y) = x^2(4 - 2.1x^2 + x^4/3) + xy + 4y^2(y^2 - 1)

    subject to  -3 <= x <= 3
                -2 <= y <= 3

  The problem has two global minima and four local minima:
    x =  0.089842     x =  1.70361      x =  1.60711
    y = -0.712656     y = -0.796084     y =  0.568651
    f = -1.03163      f = -0.215464     f =  2.10425

    x = -0.089842     x = -1.70361      x = -1.60711
    y =  0.712656     y =  0.796084     y = -0.568651
    f = -1.03163      f = -0.215464     f =  2.10425

  The input file has format:
    FC
    2
    0.000e-01
    0.000e-01
  Variable bounds are defined in the parameters file.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


// --------------------------------------------------------------------
//  Internal Function evaluateF
// --------------------------------------------------------------------
static double  evaluateF (const vector< double > &  cVectorX)
{
    double  x = cVectorX[0];
    double  y = cVectorX[1];
    double  xSqd = x*x;
    double  ySqd = y*y;
    double  f =   xSqd * (4.0 - (2.1 * xSqd) + (xSqd * xSqd / 3.0))
                + x * y
                + 4.0 * ySqd * (ySqd - 1.0);
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

    //-- WRITE ALL DIGITS USING SCIENTIFIC NOTATION.
    fp.setf (ios::scientific);
    fp.precision (15);

    //---- WRITE THE NUMBER OF OBJECTIVES AND THEIR VALUES TO THE OUTPUT.
    fp << "1" << endl;
    fp << f << endl;

    fp.close();

    return( true );
}


/* --------------------------------------------------------------------
 *  Main routine for evaluation executable.
 *
 *  Each execution reads two variables from a file, evaluates the objective
 *  function, and writes the function values to an output file.
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
