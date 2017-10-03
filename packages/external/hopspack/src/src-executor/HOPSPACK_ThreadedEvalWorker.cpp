// $Id: HOPSPACK_ThreadedEvalWorker.cpp 169 2010-04-21 19:35:53Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-executor/HOPSPACK_ThreadedEvalWorker.cpp $

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
  @file HOPSPACK_ThreadedEvalWorker.cpp
  @brief Implement HOPSPACK::ThreadedEvalWorker.
*/

#include <sstream>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Evaluator.hpp"
#include "HOPSPACK_Print.hpp"
#include "HOPSPACK_SystemTimer.hpp"
#include "HOPSPACK_ThreadedEvalWorker.hpp"
#include "HOPSPACK_ThreadRunnable.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{

//---- ENABLE THIS TO HELP DEBUG THREADING PROBLEMS.
static const bool  bDEBUG_THRDNG = false;


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
ThreadedEvalWorker::ThreadedEvalWorker
    (const int                         nWorkerNumber,
           Evaluator         *         pEvaluator,
     const int               *  const  pnSubmitTag,
     const Vector            *  const  pSubmitPointX,
     const EvalRequestType   *  const  pnSubmitRequestType,
           Vector            *  const  pResultF,
           Vector            *  const  pResultEqs,
           Vector            *  const  pResultIneqs,
           string            *  const  pResultMsg,
           ThreadSynchObject *  const  pTsoNewPointLoaded,
     const bool              *  const  pbIsWorkerIdle,
           bool              *  const  pbIsWorkerDone,
           ThreadSimpleLock  *  const  pLockIsWorkerDoneFlags,
           ThreadSynchObject ** const  ppTsoWorkerDone,
           bool              *  const  pbExecNeedsToSynch,
           bool              *  const  pbWrkrNeedsToSynch,
     const bool              *  const  pbIsTimeToShutdown,
           bool              *  const  pbIsWorkerShutdown)
    :
    _nWorkerNumber (nWorkerNumber),
    _pEvaluator (pEvaluator),
    _pnSubmitTag (pnSubmitTag),
    _pSubmitPointX (pSubmitPointX),
    _pnSubmitRequestType (pnSubmitRequestType),
    _pResultF (pResultF),
    _pResultEqs (pResultEqs),
    _pResultIneqs (pResultIneqs),
    _pResultMsg (pResultMsg),
    _pTsoNewPointLoaded (pTsoNewPointLoaded),
    _pbIsWorkerIdle (pbIsWorkerIdle),
    _pbIsWorkerDone (pbIsWorkerDone),
    _pLockIsWorkerDoneFlags (pLockIsWorkerDoneFlags),
    _ppTsoWorkerDone (ppTsoWorkerDone),
    _pbExecNeedsToSynch (pbExecNeedsToSynch),
    _pbWrkrNeedsToSynch (pbWrkrNeedsToSynch),
    _pbIsTimeToShutdown (pbIsTimeToShutdown),
    _pbIsWorkerShutdown (pbIsWorkerShutdown)
{
    //---- NO OPERATIONS ALLOWED ON SYNCHRONIZATION OBJECTS IN THE CONSTRUCTOR
    //---- BECAUSE IT IS STILL RUNNING IN THE MAIN THREAD.
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
ThreadedEvalWorker::~ThreadedEvalWorker (void)
{
    delete _pLockTsoWorkerDone;
    delete *_ppTsoWorkerDone;

    return;
}


//----------------------------------------------------------------------
//  Method run
//----------------------------------------------------------------------
void  ThreadedEvalWorker::run (void)
{
    //---- CONSTRUCT THE SYNCHRONIZING OBJECT SHARED WITH THE EXECUTOR.
    _pLockTsoWorkerDone = new ThreadSimpleLock();
    _pLockTsoWorkerDone->releaseLock();
    *_ppTsoWorkerDone = new ThreadSynchObject (*_pLockTsoWorkerDone,
                                               _pbExecNeedsToSynch);
    (*_ppTsoWorkerDone)->notify();

    //---- SET THIS SO THE EXECUTOR KNOWS THE ThreadSynchObject IS CONSTRUCTED.
    *_pbIsWorkerDone = false;

    while (true)
    {
        //---- WAIT FOR THE EXECUTOR TO NOTIFY THERE IS WORK.
        if (bDEBUG_THRDNG)
            cout << "ThDbg worker ["<<_nWorkerNumber<<"] is waiting...\n";
        _pTsoNewPointLoaded->waitForNotify();
        if (*_pbIsTimeToShutdown == true)
            break;
        if ((*_pbIsWorkerIdle == true) || (*_pbIsWorkerDone == true))
        {
            cerr << "ERROR: Worker " << _nWorkerNumber
                 << " notified but flag still shows idle." << endl;
            break;
        }

        if (Print::doPrint (Print::MOST_VERBOSE))
            cout << "ThreadedEvalWorker calling Evaluator for tag "
                 << *_pnSubmitTag
                 << " [TID=" << ThreadRunnable::getThreadID() << "]" << endl;

        //---- ACQUIRE THE MEANS TO SIGNAL BACK WHEN DONE.
        (*_ppTsoWorkerDone)->takeOwnership();

        //---- ACCESS THE POINT FROM SHARED MEMORY AND EVALUATE IT.
        if (*_pnSubmitRequestType == EVALREQTYPE_F)
        {
            _pResultF->resize (0);
            _pEvaluator->evalF (*_pnSubmitTag, *_pSubmitPointX,
                                *_pResultF, *_pResultMsg);
        }
        else if (*_pnSubmitRequestType == EVALREQTYPE_FC)
        {
            _pResultF->resize (0);
            _pResultEqs->resize (0);
            _pResultIneqs->resize (0);
            _pEvaluator->evalFC (*_pnSubmitTag, *_pSubmitPointX,
                                 *_pResultF, *_pResultEqs, *_pResultIneqs,
                                 *_pResultMsg);
        }
        else
        {
            cerr << "ERROR: Evaluator request type " << *_pnSubmitRequestType
                 << " not implemented <ThreadedEvalWorker>" << endl;
            break;
        }
        if (bDEBUG_THRDNG)
        {
            cout << "ThDbg worker will sleep...\n";
            SystemTimer::sleepMilliSecs (10 + (500 * _nWorkerNumber));
        }

        //---- TELL THE EXECUTOR THIS EVALUATION IS COMPLETE.
        *_pbWrkrNeedsToSynch = true;
        (*_ppTsoWorkerDone)->notify();
        if (bDEBUG_THRDNG)
        {
            cout << "ThDbg worker just notified it is done, sleeping 500ms\n";
            SystemTimer::sleepMilliSecs (500);
        }
        _pLockIsWorkerDoneFlags->waitForLock();
        *_pbIsWorkerDone = true;
        _pLockIsWorkerDoneFlags->releaseLock();

        //---- SYNCHRONIZE BEFORE PROCEEDING.  THIS ENSURES THE WORKER DOES
        //---- NOT LOOP THRU ITS WAIT BEFORE THE EXECUTOR CAN NOTIFY IT.
        _pTsoNewPointLoaded->synchronizeLock();
        *_pbWrkrNeedsToSynch = false;
        if (bDEBUG_THRDNG)
        {
            //---- THIS CAUSES A CRASH AT SHUTDOWN IF DESTRUCTORS ARE WRONG;
            //---- USEFUL FOR TESTING.
            // SystemTimer::sleepMilliSecs (500);
        }

    }     //-- WHILE WORKER STILL ACTIVE

    //---- END EXECUTION OF THE THREAD
    if (*_pbIsTimeToShutdown == true)
    {
        *_pbIsWorkerShutdown = true;
    }

    return;
}


}     //-- namespace HOPSPACK
