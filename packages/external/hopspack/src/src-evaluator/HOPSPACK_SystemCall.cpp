// $Id: HOPSPACK_SystemCall.cpp 220 2014-01-02 21:24:59Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-evaluator/HOPSPACK_SystemCall.cpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2014 Sandia Corporation
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
  @file HOPSPACK_SystemCall.cpp
  @brief Implement HOPSPACK::SystemCall.
*/

#include <stdio.h>     //-- FOR remove()
#include <stdlib.h>    //-- FOR system()
#include <sstream>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Print.hpp"
#include "HOPSPACK_SystemCall.hpp"
#include "HOPSPACK_utils.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
SystemCall::SystemCall (const ParameterList &  cEvalParams)
{
    _sExecutableName  = cEvalParams.getParameter ("Executable Name", "a.out");
    _sInputPrefix     = cEvalParams.getParameter ("Input Prefix", "input");
    _sOutputPrefix    = cEvalParams.getParameter ("Output Prefix", "output");
    _bSaveIOFiles     = cEvalParams.getParameter ("Save IO Files", false);
    _bDebug           = cEvalParams.getParameter ("Debug Eval Worker", false);

    _nPrecisionDigits = cEvalParams.getParameter ("File Precision", 14);
    if (_nPrecisionDigits < 0)
    {
        cerr << "WARNING: Illegal 'File Precision' value in 'Evaluator' sublist"
             << endl;
        cerr << "         Changing 'File Precision' to zero" << endl;
        _nPrecisionDigits = 0;
    }

    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
SystemCall::~SystemCall (void)
{
    return;
}


//----------------------------------------------------------------------
//  Method evalF
//----------------------------------------------------------------------
void  SystemCall::evalF (const int       nTag,
                         const Vector &  cX,
                               Vector &  cFns,
                               string &  sMsg)
{
    string  sReqType = "F";

    //---- CREATE FILE NAMES FOR I/O.
    string  sInputFileName;
    string  sOutputFileName;
    string  sSysCall;
    generateStrings_ (nTag, sReqType,
                      sInputFileName, sOutputFileName, sSysCall);

    //---- WRITE THE INPUT FILE.
    if (writeInputFile_ (sInputFileName, sReqType, cX) == false)
    {
        //---- RETURN AN EMPTY cFns.  THE MESSAGE NEEDS TO BE GENERIC
        //---- AND SHOULD NOT INCLUDE THE UNIQUE FILE NAME.
        cFns.resize (0);
        sMsg = "Could not write input file.";
        return;
    }

    //---- RUN THE EXECUTABLE, BLOCKING UNTIL COMPLETE.
    //---- INTERPRET OUTPUT FILE CONTENTS TO DECIDE IF IT WORKED.
    if (_bDebug)
    {
        cout << "  SystemCall::evalF calling '"
             << sSysCall << "'" << endl;
    }
    int  nSysReturnCode = system (sSysCall.c_str());
    if (nSysReturnCode != 0)
    {
        //---- RETURN AN EMPTY cFns.  THE MESSAGE NEEDS TO BE GENERIC
        //---- AND SHOULD NOT INCLUDE THE UNIQUE FILE NAME.
        cerr << "ERROR: Call failed: '" << sSysCall << "'"
             << " <SystemCall>" << endl;
        if (_bDebug)
            cerr << "  Return code = " << nSysReturnCode << endl;
        cFns.resize (0);
        sMsg = "Eval call failed.";

        deleteIOFile_ (sInputFileName);
        return;
    }

    //---- READ THE OUTPUT FILE.  RETURN AN EMPTY cFns ON FAILURE.
    ifstream  fptr;
    fptr.open (sOutputFileName.c_str(), ios::in);
    if (!fptr)
    {
        //---- RETURN AN EMPTY cFns.  THE MESSAGE NEEDS TO BE GENERIC
        //---- AND SHOULD NOT INCLUDE THE UNIQUE FILE NAME.
        cerr << "ERROR: Could not find file '" << sOutputFileName << "'"
             << " <SystemCall>" << endl;
        cFns.resize (0);
        sMsg = "No output file found.";

        deleteIOFile_ (sInputFileName);
        return;
    }
    if (readVector_ (fptr, sOutputFileName, cFns, sMsg) == false)
    {
        if (Print::doPrint (Print::EVALUATED_POINTS))
            cerr << "WARNING: Function evaluation returned an error for tag "
                 << nTag << endl;
        cFns.resize (0);
    }
    else
    {
        sMsg = "Success";
    }
    fptr.close();

    //---- DELETE THE TWO FILES.
    deleteIOFile_ (sInputFileName);
    deleteIOFile_ (sOutputFileName);

    return;
}


//----------------------------------------------------------------------
//  Method evalFC
//----------------------------------------------------------------------
void  SystemCall::evalFC (const int       nTag,
                          const Vector &  cX,
                                Vector &  cFns,
                                Vector &  cEqs,
                                Vector &  cIneqs,
                                string &  sMsg)
{
    string  sReqType = "FC";

    //---- CREATE FILE NAMES FOR I/O.
    string  sInputFileName;
    string  sOutputFileName;
    string  sSysCall;
    generateStrings_ (nTag, sReqType,
                      sInputFileName, sOutputFileName, sSysCall);

    //---- WRITE THE INPUT FILE.
    if (writeInputFile_ (sInputFileName, sReqType, cX) == false)
    {
        //---- RETURN AN EMPTY cFns.  THE MESSAGE NEEDS TO BE GENERIC
        //---- AND SHOULD NOT INCLUDE THE UNIQUE FILE NAME.
        cFns.resize (0);
        sMsg = "Could not write input file.";
        return;
    }

    //---- RUN THE EXECUTABLE, BLOCKING UNTIL COMPLETE.
    //---- INTERPRET OUTPUT FILE CONTENTS TO DECIDE IF IT WORKED.
    if (_bDebug)
    {
        cout << "  SystemCall::evalFC calling '"
             << sSysCall << "'" << endl;
    }
    system (sSysCall.c_str());

    //---- READ THE OUTPUT FILE.  RETURN AN EMPTY cFns ON FAILURE.
    ifstream  fptr;
    fptr.open (sOutputFileName.c_str(), ios::in);
    if (!fptr)
    {
        //---- RETURN AN EMPTY cFns.  THE MESSAGE NEEDS TO BE GENERIC
        //---- AND SHOULD NOT INCLUDE THE UNIQUE FILE NAME.
        cerr << "ERROR: Could not find file '" << sOutputFileName << "'"
             << " <SystemCall>" << endl;
        cFns.resize (0);
        cEqs.resize (0);
        cIneqs.resize (0);
        sMsg = "No output file found.";

        deleteIOFile_ (sInputFileName);
        return;
    }
    string  sDummyMsg;
    if (readVector_ (fptr, sOutputFileName, cFns, sMsg) == false)
    {
        if (Print::doPrint (Print::EVALUATED_POINTS))
            cerr << "WARNING: Function evaluation returned an error for tag "
                 << nTag << endl;
        cFns.resize (0);
        cEqs.resize (0);
        cIneqs.resize (0);
    }
    else if (readVector_ (fptr, sOutputFileName, cEqs, sDummyMsg) == false)
    {
        if (Print::doPrint (Print::EVALUATED_POINTS))
            cerr << "WARNING: Nonlinear equalities evaluation returned"
                 << " an error for tag " << nTag << endl;
        cFns.resize (0);
        cEqs.resize (0);
        cIneqs.resize (0);
    }
    else if (readVector_ (fptr, sOutputFileName, cIneqs, sDummyMsg) == false)
    {
        if (Print::doPrint (Print::EVALUATED_POINTS))
            cerr << "WARNING: Nonlinear inequalities evaluation returned"
                 << " an error for tag " << nTag << endl;
        cFns.resize (0);
        cEqs.resize (0);
        cIneqs.resize (0);
    }
    {
        sMsg = "Success";
    }
    fptr.close();

    //---- DELETE THE TWO FILES.
    deleteIOFile_ (sInputFileName);
    deleteIOFile_ (sOutputFileName);

    return;
}


//----------------------------------------------------------------------
//  Method getEvaluatorType
//----------------------------------------------------------------------
string  SystemCall::getEvaluatorType (void) const
{
    return( "System Call" );
}


//----------------------------------------------------------------------
//  Method printDebugInfo
//----------------------------------------------------------------------
void  SystemCall::printDebugInfo (void) const
{
    cout << "  HOPSPACK_SystemCall --"
         << " make a system call for evaluations" << endl;
    cout << "    Executable name:    " << _sExecutableName << endl;
    cout << "    Input file prefix:  " << _sInputPrefix << endl;
    cout << "    Output file prefix: " << _sOutputPrefix << endl;
    cout << "    File Precision:     " << _nPrecisionDigits << endl;
    cout << "    Save IO Files:      ";
    if (_bSaveIOFiles)
        cout << "true" << endl;
    else
        cout << "false" << endl;

    return;
}


//----------------------------------------------------------------------
//  Private Method generateStrings_
//----------------------------------------------------------------------
void  SystemCall::generateStrings_
          (const int       nTag,
           const string &  sReqType,
                 string &  sInputFileName,
                 string &  sOutputFileName,
                 string &  sSysCall) const
{
    stringstream  ssTag;
    ssTag << nTag;
    sInputFileName = _sInputPrefix + "." + ssTag.str() + "_" + sReqType;
    sOutputFileName = _sOutputPrefix + "." + ssTag.str() + "_" + sReqType;

    sSysCall = _sExecutableName + " " + sInputFileName
                                + " " + sOutputFileName
                                + " " + ssTag.str()
                                + " " + sReqType;

    return;
}


//----------------------------------------------------------------------
//  Private Method writeInputFile_
//----------------------------------------------------------------------
bool  SystemCall::writeInputFile_
          (const string &  sInputFileName,
           const string &  sReqType,
           const Vector &  cX) const
{
    ofstream  fptr;
    fptr.open (sInputFileName.c_str(), ios::out | ios::trunc);
    if (!fptr)
    {
        cerr << "ERROR: Could not open file '" << sInputFileName << "'"
             << " <SystemCall>" << endl;
        return( false );
    }

    //---- FIRST LINE IS THE REQUEST TYPE.
    fptr << sReqType << endl;

    //---- SECOND LINE IS THE NUMBER OF VARIABLES TO FOLLOW.
    fptr << cX.size() << endl;

    //---- ONE VARIABLE PER LINE.
    fptr.setf (ios::scientific);
    fptr.precision (_nPrecisionDigits);
    for (int  i = 0; i < cX.size(); i++)
        fptr << cX[i] << endl;

    fptr.close();
    return( true );
}


//----------------------------------------------------------------------
//  Private Method readVector_
//----------------------------------------------------------------------
bool  SystemCall::readVector_ (      ifstream &  fptr,
                                        const string   &  sFileName,
                                              Vector   &  cV,
                                              string   &  sMsg) const
{
    //---- WIPE OUT ANY EXISTING VECTOR CONTENTS.
    cV.resize (0);

    if (fptr.eof())
    {
        cerr << "ERROR: Unexpected end of file '" << sFileName << "'"
             << " <SystemCall>" << endl;
        return( false );
    }

    string  sNextLine;
    string::size_type  nPos;

    //---- READ THE SIZE OF THE VECTOR, AN INTEGER ON ONE LINE,
    //---- OR CHECK FOR A POSSIBLE ERROR MESSAGE.
    getline (fptr, sNextLine);
    nPos = 0;
    int  nVectorSize;
    if (HOPSPACK::getNextInt (sNextLine, nPos, nVectorSize) == false)
    {
        //---- ASSUME THIS IS THE FIRST LINE OF THE FILE AND THE CONTENTS
        //---- ARE AN EVALUATION ERROR MESSAGE.
        sMsg = sNextLine;
        return( false );
    }

    //---- READ THE VECTOR.
    for (int  i = 0; i < nVectorSize; i++)
    {
        if (fptr.eof())
        {
            cerr << "ERROR: Not enough vector components in file '"
                 << sFileName << "' <SystemCall>" << endl;
            sMsg = "Not enough vector components found.";
            return( false );
        }
        getline (fptr, sNextLine);
        string  sTmp;
        HOPSPACK::getNextString (sNextLine, nPos, sTmp);
        if (sTmp == "DNE")
            cV.push_back (dne());
        else
        {
            double  dTmp;
          #if defined(HAVE_MSVC_SECURE_STRING_FNS)
            if (sscanf_s (sTmp.c_str(), "%le", &dTmp) == 1)
          #else
            if (sscanf (sTmp.c_str(), "%le", &dTmp) == 1)
          #endif
                cV.push_back (dTmp);
            else
            {
                cerr << "ERROR: Expected a number, found '" << sTmp
                     << "' in file '" << sFileName
                     << "' <SystemCall>" << endl;
                sMsg = "Bad number found while reading vector components.";
                return( false );
            }
        }
    }

    return( true );
}


//----------------------------------------------------------------------
//  Private Method deleteIOFile_
//----------------------------------------------------------------------
void  SystemCall::deleteIOFile_ (const string &  sFileName) const
{
    if (_bSaveIOFiles == false)
    {
        remove (sFileName.c_str());
    }
    return;
}

    
}     //-- namespace HOPSPACK
