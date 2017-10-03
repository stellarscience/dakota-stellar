// $Id: HOPSPACK_ExecutorMultiThreaded.cpp 220 2014-01-02 21:24:59Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-executor/HOPSPACK_ExecutorMultiThreaded.cpp $

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
  @file HOPSPACK_ExecutorMultiThreaded.cpp
  @brief Implement HOPSPACK::ExecutorMultiThreaded.
*/

#include <iomanip>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Evaluator.hpp"
#include "HOPSPACK_EvaluatorFactory.hpp"
#include "HOPSPACK_ExecutorMultiThreaded.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Print.hpp"
#include "HOPSPACK_SystemTimer.hpp"
#include "HOPSPACK_ThreadedEvalWorker.hpp"
#include "HOPSPACK_ThreadSimpleLock.hpp"
#include "HOPSPACK_ThreadSynchObject.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
ExecutorMultiThreaded::ExecutorMultiThreaded (void)
    :
    _bIsInitialized (false)
{
    //---- DO MOST OF THE WORK IN "initialize()" SO IT CAN HANDLE ERRORS.

    _naWorkersDoneAndWaiting.clear();
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
ExecutorMultiThreaded::~ExecutorMultiThreaded (void)
{
    if (_bIsInitialized == false)
        return;

    delete _pTimers;

    delete[] _baIsWorkerIdle;
    delete[] _baIsWorkerDone;
    delete[] _baExecNeedsToSynch;
    delete[] _baIsWorkerShutdown;

    delete[] _naSubmitTags;
    delete[] _caSubmitPoints;
    delete[] _naSubmitRequestTypes;
    delete[] _caRecvF;
    delete[] _caRecvEqs;
    delete[] _caRecvIneqs;
    delete[] _saRecvMsgs;

    for (int  i = 0; i < _nNumWorkers; i++)
    {
        delete _paLockTsoNewPointLoaded[i];
        delete _paTsoNewPointLoaded[i];
    }
    delete[] _paLockTsoNewPointLoaded;
    delete[] _paTsoNewPointLoaded;
    delete[] _baWrkrNeedsToSynch;

    delete _pLockIsWorkerDoneFlags;

    for (int  i = 0; i < _nNumWorkers; i++)
    {
        delete _paEvaluators[i];
        delete _paWorkers[i];
    }
    delete[] _paEvaluators;
    delete[] _paWorkers;

    //---- WORKERS DELETE THE CONTENTS OF THIS ARRAY.
    delete[] _paTsoWorkerDone;

    return;
}


//----------------------------------------------------------------------
//  Method initialize
//----------------------------------------------------------------------
bool  ExecutorMultiThreaded::initialize (const int              nNumWorkers,
                                         const ParameterList &  cEvalParams)
{
    if (_bIsInitialized == true)
        //---- ONLY INITIALIZE ONCE.
        return( true );

    if (nNumWorkers <= 0)
    {
        cerr << "ERROR: ExecutorMultiThreaded called with num workers = "
             << nNumWorkers << endl;
        throw INTERNAL_ERROR;
    }
    _nNumWorkers = nNumWorkers;

    _pTimers = new SystemTimer (1 + _nNumWorkers);
    _pTimers->start (0);

    //---- CREATE FLAGS FOR WORKER AND EXECUTOR COORDINATION.
    _baIsWorkerIdle = new bool [_nNumWorkers];
    _baIsWorkerDone = new bool [_nNumWorkers];
    _baExecNeedsToSynch = new bool [_nNumWorkers];
    _baIsWorkerShutdown = new bool [_nNumWorkers];
    for (int  i = 0; i < _nNumWorkers; i++)
    {
        _baIsWorkerIdle[i] = true;
        _baIsWorkerDone[i] = true;         //-- ALSO SIGNALS WORKER STARTED
        _baExecNeedsToSynch[i] = false;
        _baIsWorkerShutdown[i] = false;
    }

    //---- CREATE STORAGE FOR INPUTS AND OUTPUTS TO EACH WORKER.
    _naSubmitTags = new int [_nNumWorkers];
    _caSubmitPoints = new Vector [_nNumWorkers];
    _naSubmitRequestTypes = new EvalRequestType [_nNumWorkers];
    _caRecvF = new Vector [_nNumWorkers];
    _caRecvEqs = new Vector [_nNumWorkers];
    _caRecvIneqs = new Vector [_nNumWorkers];
    _saRecvMsgs = new string [_nNumWorkers];
    for (int  i = 0; i < _nNumWorkers; i++)
    {
        _naSubmitTags[i] = -1;
    }

    //---- CREATE SYNCHRONIZATION OBJECTS FOR EACH WORKER.
    _paLockTsoNewPointLoaded = new ThreadSimpleLock * [_nNumWorkers];
    _paTsoNewPointLoaded = new ThreadSynchObject * [_nNumWorkers];
    _baWrkrNeedsToSynch = new bool [_nNumWorkers];
    for (int  i = 0; i < _nNumWorkers; i++)
    {
        //---- EXECUTOR OWNS THE OBJECTS INITIALLY.
        ThreadSimpleLock *  pNextLock = new ThreadSimpleLock();
        pNextLock->releaseLock();
        _paLockTsoNewPointLoaded[i] = pNextLock;
        _baWrkrNeedsToSynch[i] = false;
        ThreadSynchObject *  pNext
            = new ThreadSynchObject (*pNextLock, &(_baWrkrNeedsToSynch[i]));
        _paTsoNewPointLoaded[i] = pNext;
    }

    _pLockIsWorkerDoneFlags = new ThreadSimpleLock();
    _pLockIsWorkerDoneFlags->releaseLock();

    _bIsTimeToShutdown = false;

    //---- CREATE A POOL OF THREADED WORKERS AND START THEM.
    _paEvaluators = new Evaluator * [_nNumWorkers];
    _paWorkers = new ThreadedEvalWorker * [_nNumWorkers];
    _paTsoWorkerDone = new ThreadSynchObject * [_nNumWorkers];
    for (int  i = 0; i < _nNumWorkers; i++)
    {
        //---- CREATE AN EVALUATOR INSTANCE FOR EACH WORKER.
        _paEvaluators[i] = EvaluatorFactory::newInstance (cEvalParams);
        if (_paEvaluators[i] == NULL)
        {
            cerr << "ERROR: Could not construct evaluator." << endl;
            return( false );
        }

        //---- CREATE THE WORKER.
        _paTsoWorkerDone[i] = NULL;
        ThreadedEvalWorker *  pNextWorker
            = new ThreadedEvalWorker (i + 1,
                                      _paEvaluators[i],
                                      &(_naSubmitTags[i]),
                                      &(_caSubmitPoints[i]),
                                      &(_naSubmitRequestTypes[i]),
                                      &(_caRecvF[i]),
                                      &(_caRecvEqs[i]),
                                      &(_caRecvIneqs[i]),
                                      &(_saRecvMsgs[i]),
                                      _paTsoNewPointLoaded[i],
                                      &(_baIsWorkerIdle[i]),
                                      &(_baIsWorkerDone[i]),
                                      _pLockIsWorkerDoneFlags,
                                      &(_paTsoWorkerDone[i]),
                                      &(_baExecNeedsToSynch[i]),
                                      &(_baWrkrNeedsToSynch[i]),
                                      &_bIsTimeToShutdown,
                                      &(_baIsWorkerShutdown[i]) );
        _paWorkers[i] = pNextWorker;

        //---- START THE WORKER IN ITS OWN THREAD.
        _paWorkers[i]->start();

        //---- WAIT UNTIL THE WORKER HAS CONSTRUCTED _paTsoWorkerDone[i].
        //---- AFTER THAT, THE WORKER WAITS FOR NOTIFICATION VIA
        //----_paTsoNewPointLoaded.
        while (_baIsWorkerDone[i] == true)
        {
            SystemTimer::sleepMilliSecs (1);
        }
    }

    _bIsInitialized = true;
    return( true );
}


//----------------------------------------------------------------------
//  Method submit
//----------------------------------------------------------------------
bool  ExecutorMultiThreaded::submit (const int              nTag,
                                     const Vector &         cX,
                                     const EvalRequestType  nRequest)
{
    if (_bIsInitialized == false)
    {
        cerr << "ERROR: Must call initialize() before submit()"
             << "  <ExecutorMultiThreaded::submit>" << endl;
        throw INTERNAL_ERROR;
    }

    if ((nRequest != EVALREQTYPE_F) && (nRequest != EVALREQTYPE_FC))
    {
        cerr << "ERROR: Evaluator request type " << nRequest
             << " not implemented <ExecutorMultiThreaded::submit>" << endl;
        throw INTERNAL_ERROR;
    }

    //---- FIND A WORKER AND ASSIGN THE POINT.
    for (int  i = 0; i < _nNumWorkers; i++)
    {
        if (_baIsWorkerIdle[i] == true)
        {
            //---- COPY THE INPUTS TO MEMORY SHARED WITH THE CHOSEN WORKER.
            _naSubmitTags[i] = nTag;
            _caSubmitPoints[i] = cX;
            _naSubmitRequestTypes[i] = nRequest;

            if (Print::doPrint (Print::MOST_VERBOSE))
                cout << "ExecutorMultiThreaded calling worker "
                     << (i + 1) << " for tag " << nTag << endl;

            //---- NOTIFY THE WORKER.
            //---- TELL IT THE EXECUTOR NEEDS TO CROSS THE SYNCHRONIZATION
            //---- BARRIER BEFORE THE WORKER NOTIFIES BACK.
            _baIsWorkerIdle[i] = false;
            _baExecNeedsToSynch[i] = true;
            _paTsoNewPointLoaded[i]->notify();

            _pTimers->start (1 + i);

            //---- SYNCHRONIZE BEFORE PROCEEDING.  THIS ENSURES THE EXECUTOR DOES
            //---- NOT LOOP THRU ITS WAIT BEFORE THE WORKER CAN NOTIFY IT.
            _paTsoWorkerDone[i]->synchronizeLock();
            _baExecNeedsToSynch[i] = false;

            return( true );
        }
    }

    //---- ALL WORKERS ARE BUSY.
    return( false );
}


//----------------------------------------------------------------------
//  Method isReadyForWork
//----------------------------------------------------------------------
bool  ExecutorMultiThreaded::isReadyForWork (void) const
{
    for (int  i = 0; i < _nNumWorkers; i++)
        if (_baIsWorkerIdle[i])
            return( true );

    return( false );
}


//----------------------------------------------------------------------
//  Method recv
//----------------------------------------------------------------------
int  ExecutorMultiThreaded::recv (int    &  nTag,
                                  Vector &  cF,
                                  Vector &  cEqs,
                                  Vector &  cIneqs,
                                  string &  sMsg)
{
    if (_bIsInitialized == false)
    {
        cerr << "ERROR: Must call initialize() before recv()"
             << "  <ExecutorMultiThreaded::recv>" << endl;
        throw INTERNAL_ERROR;
    }

    //---- IF THERE ARE IDLE WORKERS BUT NO REPLIES ARE IN YET,
    //---- THEN RETURN IMMEDIATELY SO THE CONVEYOR CAN SCHEDULE MORE WORKERS.
    if ((isReadyForWork() == true) && (_naWorkersDoneAndWaiting.size() == 0))
    {
        //---- PEEK AT WORKER STATUS FLAGS TO SEE IF THERE ARE ANY REPLIES.
        //---- IT'S POSSIBLE FOR A WORKER TO FINISH WHILE THIS LOOP RUNS,
        //---- BUT IT JUST MEANS QUICK ITERATION IN THE CONVEYOR AND THEN
        //---- THE REPLY WILL BE DISCOVERED.
        bool  bIsAnyWorkerDone = false;
        for (int  i = 0; i < _nNumWorkers; i++)
            if (_baIsWorkerDone[i])
            {
                bIsAnyWorkerDone = true;
                break;
            }
        if (bIsAnyWorkerDone == false)
            return( 0 );
    }

    //---- IF THERE ARE NO QUEUED RESULTS, THEN BLOCK UNTIL ONE OR MORE
    //---- WORKERS NOTIFY.
    if (_naWorkersDoneAndWaiting.size() == 0)
    {
        //---- THIS IS A POLLING LOOP THAT WAITS FOR AT LEAST ONE WORKER
        //---- TO RESPOND.  THEORETICALLY, AN OPERATING SYSTEM CALL SUCH AS
        //---- WINDOWS waitForMultipleObjects WILL DO THIS WITHOUT SPINNING,
        //---- BUT THE POTENTIALS FOR DEADLOCK ARE HARD TO ERADICATE.
        while (_naWorkersDoneAndWaiting.size() == 0)
        {
            _pLockIsWorkerDoneFlags->waitForLock();
            for (int  i = 0; i < _nNumWorkers; i++)
            {
                if (_baIsWorkerDone[i])
                {
                    _naWorkersDoneAndWaiting.push_back (i);
                    _baIsWorkerDone[i] = false;
                    _pTimers->stop (1 + i);
                }
            }
            _pLockIsWorkerDoneFlags->releaseLock();

            //---- WASTE SOME TIME SO WORKERS HAVE A CHANCE TO GET THE LOCK.
            SystemTimer::sleepMilliSecs (1);
        }
    }

    //---- CHOOSE ONE OF THE FINISHED WORKERS.
    int  k = _naWorkersDoneAndWaiting.back();

    //---- COPY ITS EVALUATION RESULTS.
    _naWorkersDoneAndWaiting.pop_back();
    nTag   = _naSubmitTags[k];
    cF     = _caRecvF[k];
    cEqs   = _caRecvEqs[k];
    cIneqs = _caRecvIneqs[k];
    sMsg   = _saRecvMsgs[k];

    //---- THIS WILL WORK IMMEDIATELY BECAUSE THE WORKER ALREADY NOTIFIES
    //---- BEFORE SETTING THE FLAG THAT IT IS DONE.
    _paTsoWorkerDone[k]->waitForNotify();

    //---- RELEASE THE WORKER TO ACCEPT MORE WORK.
    _baIsWorkerDone[k] = false;
    _baIsWorkerIdle[k] = true;

    //---- THIS IS A GOOD PLACE TO TEST FOR DEADLOCK.
    // SystemTimer::sleepMilliSecs (800);

    _paTsoNewPointLoaded[k]->takeOwnership();

    return( k + 1 );
}


//----------------------------------------------------------------------
//  Method shutdown
//----------------------------------------------------------------------
bool  ExecutorMultiThreaded::shutdown (void)
{
    if (_bIsInitialized == false)
    {
        cerr << "ERROR: Must call initialize() before shutdown()"
             << "  <ExecutorMultiThreaded::shutdown>" << endl;
        throw INTERNAL_ERROR;
    }

    //---- VERIFY ALL WORKERS ARE IDLE.  IT IS OK IF A WORKER'S LAST
    //---- EVALUATION HAS NOT BEEN PICKED UP BY A CALL TO recv().
    for (int  i = 0; i < _nNumWorkers; i++)
        if (_baIsWorkerIdle[i] == false)
        {
            cerr << "ERROR: Cannot shut down Executor,"
                 << " workers are still busy." << endl;
            return( false );
        }

    //---- NOTIFY ALL THE WORKERS THAT THEY SHOULD SHUT DOWN.
    _bIsTimeToShutdown = true;
    for (int  i = 0; i < _nNumWorkers; i++)
        _paTsoNewPointLoaded[i]->notify();

    //---- WAIT AND TEST THAT WORKERS ARE SHUT DOWN.
    SystemTimer::sleepMilliSecs (20);
    int  nNumLoops = 0;
    int  nMaxNumLoops = 10;     //-- 500 MSEC PER LOOP.
    for (; nNumLoops < nMaxNumLoops; nNumLoops++)
    {
        bool  bAllWorkersShutdown = true;
        for (int  i = 0; i < _nNumWorkers; i++)
            if (_baIsWorkerShutdown[i] == false)
                bAllWorkersShutdown = false;
        if (bAllWorkersShutdown)
            break;
        SystemTimer::sleepMilliSecs (500);
    }
    if (nNumLoops >= nMaxNumLoops)
    {
        cerr << "ERROR: Workers notified to shutdown but not responding"
             << "  <ExecutorMultiThreaded::shutdown>" << endl;
        return( false );
    }

    return( true );
}


//----------------------------------------------------------------------
//  Method getEvaluatorType
//----------------------------------------------------------------------
string  ExecutorMultiThreaded::getEvaluatorType (void) const
{
    if ((_bIsInitialized == false) || (_nNumWorkers == 0))
        return( "undefined" );
    return( _paEvaluators[0]->getEvaluatorType() );
}


//----------------------------------------------------------------------
//  Method printDebugInfo
//----------------------------------------------------------------------
void  ExecutorMultiThreaded::printDebugInfo (void) const
{
    cout << "  HOPSPACK_ExecutorMultiThreaded -- number evaluation workers = "
         << _nNumWorkers << endl;
    cout << "    isReadyForWork() returns = " << isReadyForWork() << endl;

    int  nNumIdle = 0;
    for (int  i = 0; i < _nNumWorkers; i++)
        if (_baIsWorkerIdle[i])
            nNumIdle++;
    cout << "    number of currently idle workers = " << nNumIdle << endl;

    for (int  i = 0; i < _nNumWorkers; i++)
    {
        cout << "    worker[" << i << "]:  is idle = " << _baIsWorkerIdle[i]
             << ", current tag = " << _naSubmitTags[i] << endl;
        _paEvaluators[i]->printDebugInfo();
    }

    return;
}


//----------------------------------------------------------------------
//  Method printTimingInfo
//----------------------------------------------------------------------
void  ExecutorMultiThreaded::printTimingInfo (void) const
{
    _pTimers->stop (0);
    cout << "Total wall clock time in Multithreaded Executor: "
         << _pTimers->getTotalTime (0) << " secs" << endl;

    streamsize  nCurrPrec  = cout.precision (3);
    cout.setf (ios::fixed | ios::right);
    for (int  i = 0; i < _nNumWorkers; i++)
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
