// $Id: HOPSPACK_ExecutorMpi.cpp 220 2014-01-02 21:24:59Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-executor/HOPSPACK_ExecutorMpi.cpp $

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
  @file HOPSPACK_ExecutorMpi.cpp
  @brief Implement HOPSPACK::ExecutorMpi.
*/

#include <iomanip>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_ExecutorMpi.hpp"
#include "HOPSPACK_GenProcComm.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Print.hpp"
#include "HOPSPACK_SystemTimer.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Static Definitions
//----------------------------------------------------------------------

const static int  nWORKER_IDLE = -9999;


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
ExecutorMpi::ExecutorMpi (void)
    :
    _bIsInitialized (false),
    _pTimers (NULL),
    _naWorkerIDs(),
    _naWorkerTag()
{
    //---- DO THE WORK IN "initialize()" SO IT CAN HANDLE ERRORS.
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
ExecutorMpi::~ExecutorMpi (void)
{
    if (_pTimers != NULL)
        delete _pTimers;
    return;
}


//----------------------------------------------------------------------
//  Method initialize
//----------------------------------------------------------------------
bool  ExecutorMpi::initialize (const vector< int > &  naWorkerIDs)
{
    if (_bIsInitialized == true)
        //---- ONLY INITIALIZE ONCE.
        return( true );

    for (int  i = 0; i < (int) naWorkerIDs.size(); i++)
    {
        _naWorkerIDs.push_back (naWorkerIDs[i]);
        _naWorkerTag.push_back (nWORKER_IDLE);
    }

    _pTimers = new SystemTimer (1 + naWorkerIDs.size());
    _pTimers->start (0);

    _bIsInitialized = true;
    return( true );
}


//----------------------------------------------------------------------
//  Method submit
//----------------------------------------------------------------------
bool  ExecutorMpi::submit (const int              nTag,
                           const Vector&          cX,
                           const EvalRequestType  nRequest)
{
    if (_bIsInitialized == false)
    {
        cerr << "ERROR: Must call initialize() before submit()"
             << "  <ExecutorMpi::submit>" << endl;
        throw INTERNAL_ERROR;
    }

    for (int  i = 0; i < (int) _naWorkerTag.size(); i++)
    {
        if (_naWorkerTag[i] == nWORKER_IDLE)
        {
            GenProcComm &  cGPC = GenProcComm::getTheInstance();
            cGPC.initSend();

            if (   (nRequest == EVALREQTYPE_F)
                || (nRequest == EVALREQTYPE_FC) )
            {
                cGPC.pack (nTag);
                cGPC.pack (cX);
                cGPC.pack (nRequest);
            }
            else
            {
                cerr << "ERROR: Evaluator request type " << nRequest
                     << " not implemented <ExecutorMpi::submit>" << endl;
                throw INTERNAL_ERROR;
            }

            _naWorkerTag[i] = nTag;

            int  nWorkerID = _naWorkerIDs[i];
            if (Print::doPrint (Print::MOST_VERBOSE))
                cout << "ExecutorMpi calling evaluator "
                     << nWorkerID << " for tag " << nTag << endl;

            _pTimers->start (1 + i);
            cGPC.send (GenProcComm::EVALUATE, nWorkerID);
            return( true );
        }
    }

    return( false );
}


//----------------------------------------------------------------------
//  Method isReadyForWork
//----------------------------------------------------------------------
bool  ExecutorMpi::isReadyForWork (void) const
{
    for (int  i = 0; i < (int) _naWorkerTag.size(); i++)
        if (_naWorkerTag[i] == nWORKER_IDLE)
            return( true );
    return( false );
}


//----------------------------------------------------------------------
//  Method recv
//----------------------------------------------------------------------
int  ExecutorMpi::recv (int    &  nTag,
                        Vector &  cF,
                        Vector &  cEqs,
                        Vector &  cIneqs,
                        string &  sMsg)
{
    if (_bIsInitialized == false)
    {
        cerr << "ERROR: Must call initialize() before recv()"
             << "  <ExecutorMpi::recv>" << endl;
        throw INTERNAL_ERROR;
    }

    GenProcComm &  cGPC = GenProcComm::getTheInstance();

    //---- IF THERE ARE IDLE WORKERS BUT NO REPLIES ARE IN YET,
    //---- THEN RETURN IMMEDIATELY SO THE CONVEYOR CAN SCHEDULE MORE WORKERS.
    //---- OTHERWISE, CONTINUE TO THE BLOCKING recv CALL.
    if (isReadyForWork() == true)
    {
        if (cGPC.probe (GenProcComm::EVALUATE) == 0)
            return( 0 );
    }

    int  nTmp;
    int  nWorkerID;
    cGPC.recv (GenProcComm::EVALUATE);
    cGPC.bufinfo (nTmp, nWorkerID);
    _pTimers->stop (nWorkerID);

    cGPC.unpack (nTag);
    cGPC.unpack (cF);
    cGPC.unpack (cEqs);
    cGPC.unpack (cIneqs);
    cGPC.unpack (sMsg);

    int  nWorkerIndex = -1;
    for (int  i = 0; i < (int) _naWorkerIDs.size(); i++)
        if (_naWorkerIDs[i] == nWorkerID)
        {
            nWorkerIndex = i;
            break;
        }
    if (nWorkerIndex == -1)
    {
        cerr << "ERROR: Evaluated result from unknown worker ID " << nWorkerID
             << endl;
        throw INTERNAL_ERROR;
    }
    _naWorkerTag[nWorkerIndex] = nWORKER_IDLE;

    if (Print::doPrint (Print::MOST_VERBOSE))
        cout << "ExecutorMpi recv result from evaluator "
             << nWorkerID << " for tag " << nTag << endl;
    return( nWorkerID );
}


//----------------------------------------------------------------------
//  Method getEvaluatorType
//----------------------------------------------------------------------
string  ExecutorMpi::getEvaluatorType (void) const
{
    //TBD...MPI master figures this out from eval params, not by
    //      constructing an evaluator (only workers construct one)
    return( "TBD" );
}


//----------------------------------------------------------------------
//  Method printDebugInfo
//----------------------------------------------------------------------
void  ExecutorMpi::printDebugInfo (void) const
{
    cout << "  HOPSPACK_ExecutorMpi -- number evaluation workers = "
         << _naWorkerTag.size() << endl;
    cout << "    isReadyForWork() returns = " << isReadyForWork() << endl;

    int  nNumIdle = 0;
    for (int  i = 0; i < (int) _naWorkerTag.size(); i++)
        if (_naWorkerTag[i] == nWORKER_IDLE)
            nNumIdle++;
    cout << "    number of currently idle workers = " << nNumIdle << endl;

    for (int  i = 0; i < (int) _naWorkerTag.size(); i++)
        cout << "    worker[" << i << "]:  MPI rank = " << _naWorkerIDs[i]
             << ", current tag = " << _naWorkerTag[i] << endl;

    return;
}


//----------------------------------------------------------------------
//  Method printTimingInfo
//----------------------------------------------------------------------
void  ExecutorMpi::printTimingInfo (void) const
{
    _pTimers->stop (0);
    cout << "Total wall clock time in MPI Executor: "
         << _pTimers->getTotalTime (0) << " secs" << endl;

    int  nCurrPrec  = cout.precision (3);
    cout.setf (ios::fixed | ios::right);
    for (int  i = 0; i < (int) _naWorkerTag.size(); i++)
    {
        cout << "  Eval worker "
             << setw (3) << (1 + i) << "         "
             << setw (8) << _pTimers->getTotalTime (1 + i)
             << "   (" << setw (5) << _pTimers->getNumStarts (1 + i)
             << " calls)" << endl;
    }
    cout.precision (nCurrPrec);

    return;
}

    
}     //-- namespace HOPSPACK
