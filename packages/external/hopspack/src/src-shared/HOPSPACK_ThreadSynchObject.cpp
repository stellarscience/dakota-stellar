// $Id: HOPSPACK_ThreadSynchObject.cpp 169 2010-04-21 19:35:53Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_ThreadSynchObject.cpp $

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
  @file HOPSPACK_ThreadSynchObject.cpp
  @brief Implement HOPSPACK::ThreadSynchObject
*/

#include "HOPSPACK_common.hpp"

//---- THIS CLASS IS DEFINED ONLY FOR MULTI-THREADED BUILDS.
#if defined(HAVE_MT)

#if defined(_WIN32)
  #include <windows.h>
#else
  #include <pthread.h>
#endif

#include "HOPSPACK_SystemTimer.hpp"
#include "HOPSPACK_ThreadRunnable.hpp"
#include "HOPSPACK_ThreadSimpleLock.hpp"
#include "HOPSPACK_ThreadSynchObject.hpp"

namespace HOPSPACK
{

//---- ENABLE THIS TO HELP DEBUG THREADING PROBLEMS.
static const bool  bDEBUG_THRDNG = false;


//----------------------------------------------------------------------
//  Private structures
//----------------------------------------------------------------------

//---- HIDE OPERATING SYSTEM DETAILS IN THESE STRUCTURES.
struct _systemMonitorType
{
    #if defined(_WIN32)
        HANDLE  cWinMonitor;
    #else
        pthread_cond_t   cCondVar;
        bool             bNeedToWait;
        pthread_mutex_t  cCondMutex;
    #endif
};


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
ThreadSynchObject::ThreadSynchObject (ThreadSimpleLock &  cLock,
                                      bool * const        pbSynchFlag)
    :
    _cAssociatedLock (cLock),
    _pbOtherNeedsToSynch (pbSynchFlag)
{
    _pMonitor = new _systemMonitorType;
    #if defined(_WIN32)
        //---- WINDOWS EVENTS ARE EITHER SET ("SIGNALED")
        //---- OR RESET ("UNSIGNALED").
        //---- A WAIT CALL RETURNS AFTER ITS EVENT IS SET.
        //---- 
        //---- CONSTRUCT THE EVENT TO AUTOMATICALLY RESET AFTER WAITING,
        //---- AND TO BE INITIALLY UNSIGNALED (UNOWNED).
        _pMonitor->cWinMonitor = CreateEvent (NULL, FALSE, FALSE, NULL);
        if (_pMonitor->cWinMonitor == NULL)
        {
            cerr << "ERROR: Could not create Windows event." << endl;
            throw INTERNAL_ERROR;
        }
        if (bDEBUG_THRDNG)
            printf ("ThDSyn constructed %p\n",&(_pMonitor->cWinMonitor));

    #else
        //---- UNIX CONDITION VARIABLES NEED AN ASSOCIATED BOOLEAN
        //---- TO TELL IF THEY WERE SIGNALED BEFORE A WAIT WAS CALLED.
        if (pthread_mutex_init (&(_pMonitor->cCondMutex), NULL) != 0)
        {
            cerr << "ERROR: Could not create Unix mutex." << endl;
            throw INTERNAL_ERROR;
        }
        if (pthread_cond_init (&(_pMonitor->cCondVar), NULL) != 0)
        {
            cerr << "ERROR: Could not create Unix condition variable." << endl;
            throw INTERNAL_ERROR;
        }
        _pMonitor->bNeedToWait = true;
        if (bDEBUG_THRDNG)
            printf ("ThDSyn constructed %p\n",&(_pMonitor->cCondVar));
    #endif

    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
ThreadSynchObject::~ThreadSynchObject (void)
{
    #if defined(_WIN32)
        //---- WINDOWS DOES NOT REQUIRE OWNERSHIP TO DELETE.  THERE WILL BE
        //---- PROBLEMS IF OTHER THREADS ARE STILL WAITING ON THIS INSTANCE.
        CloseHandle (_pMonitor->cWinMonitor);

    #else
        if (pthread_mutex_destroy (&(_pMonitor->cCondMutex)) != 0)
        {
            cerr << "ERROR: Could not destroy Unix mutex." << endl;
            throw INTERNAL_ERROR;
        }
        if (pthread_cond_destroy (&(_pMonitor->cCondVar)) != 0)
        {
            cerr << "ERROR: Could not destroy Unix condition variable." << endl;
            throw INTERNAL_ERROR;
        }
    #endif

    delete _pMonitor;

    return;
}


//----------------------------------------------------------------------
//  Method synchronizeLock
//----------------------------------------------------------------------
void  ThreadSynchObject::synchronizeLock (void)
{
    _cAssociatedLock.waitForLock();
    _cAssociatedLock.releaseLock();

    return;
}


//----------------------------------------------------------------------
//  Method waitForNotify
//----------------------------------------------------------------------
void  ThreadSynchObject::waitForNotify (void)
{
    #if defined(_WIN32)
        if (bDEBUG_THRDNG)
            printf ("ThDSyn waiting %p\n",&(_pMonitor->cWinMonitor));
        //---- WINDOWS WILL RETURN WHEN THE EVENT IS SET, BY ANY THREAD.
        DWORD  dwTmp = WaitForSingleObject (_pMonitor->cWinMonitor, INFINITE);
        if (dwTmp != WAIT_OBJECT_0)
        {
            cerr << "ERROR: Windows WaitForSingleObject returned " << dwTmp
                 << " <ThreadSynchObject::waitForNotify>" << endl;
            throw INTERNAL_ERROR;
        }
        if (bDEBUG_THRDNG)
            printf ("ThDSyn done waiting %p\n",&(_pMonitor->cWinMonitor));
        //---- IMMEDIATELY SET THE CAPTURED EVENT SO THAT IT IS READY
        //---- FOR SOME OTHER THREAD TO TAKE OWNERSHIP.
        if (SetEvent (_pMonitor->cWinMonitor) == 0)
        {
            cerr << "ERROR: Windows SetEvent failed"
                 << " <ThreadSynchObject::waitForNotify>" << endl;
            throw INTERNAL_ERROR;
        }

    #else
        pthread_mutex_lock (&(_pMonitor->cCondMutex));
        if (bDEBUG_THRDNG)
            printf ("ThDSyn waiting %p\n",&(_pMonitor->cCondVar));
        if (_pMonitor->bNeedToWait)
        {
            int  nTmp = pthread_cond_wait (&(_pMonitor->cCondVar),
                                           &(_pMonitor->cCondMutex));
            if (nTmp != 0)
            {
                cerr << "ERROR: Unix pthread_cond_wait returned " << nTmp
                     << " <ThreadSynchObject::waitForNotify>" << endl;
                throw INTERNAL_ERROR;
            }
        }
        _pMonitor->bNeedToWait = false;
        if (bDEBUG_THRDNG)
            printf ("ThDSyn done waiting %p\n",&(_pMonitor->cCondVar));
        pthread_mutex_unlock (&(_pMonitor->cCondMutex));

    #endif

    return;
}


//----------------------------------------------------------------------
//  Method takeOwnership
//----------------------------------------------------------------------
void  ThreadSynchObject::takeOwnership (void)
{
    #if defined(_WIN32)
        if (bDEBUG_THRDNG)
            printf ("ThDSyn taking ownship %p\n",&(_pMonitor->cWinMonitor));
        DWORD  dwTmp = WaitForSingleObject (_pMonitor->cWinMonitor, INFINITE);
        if (dwTmp != WAIT_OBJECT_0)
        {
            cerr << "ERROR: Windows WaitForSingleObject returned " << dwTmp
                 << " <ThreadSynchObject::takeOwnership>" << endl;
            throw INTERNAL_ERROR;
        }
        if (bDEBUG_THRDNG)
            printf ("ThDSyn done taking ownship %p\n",&(_pMonitor->cWinMonitor));

    #else
        if (bDEBUG_THRDNG)
            printf ("ThDSyn taking ownship %p\n",&(_pMonitor->cCondVar));
        pthread_mutex_lock (&(_pMonitor->cCondMutex));
        if (_pMonitor->bNeedToWait)
        {
            int  nTmp = pthread_cond_wait (&(_pMonitor->cCondVar),
                                           &(_pMonitor->cCondMutex));
            if (nTmp != 0)
            {
                cerr << "ERROR: Unix pthread_cond_wait returned " << nTmp
                     << " <ThreadSynchObject::takeOwnership>" << endl;
                throw INTERNAL_ERROR;
            }
        }
        _pMonitor->bNeedToWait = true;
        pthread_mutex_unlock (&(_pMonitor->cCondMutex));
 
        if (bDEBUG_THRDNG)
            printf ("ThDSyn done taking ownship %p\n",&(_pMonitor->cCondVar));

    #endif

    _cAssociatedLock.releaseLock();

    return;
}


//----------------------------------------------------------------------
//  Method notify
//----------------------------------------------------------------------
void  ThreadSynchObject::notify (void)
{
    //---- TO PRECISELY TIGHTEN ALL DEADLOCK WINDOWS, TWO SYNCHRONIZING
    //---- THREADS MUST AGREE TO USE THE LOCK AS A BARRIER.
    //---- IN PARTICULAR, DO NOT NOTIFY UNTIL THE OTHER THREAD IS BEYOND
    //---- THE BARRIER.
    //----
    //---- A PROPERLY CODED WORKER SHOULD CROSS ITS BARRIER IMMEDIATELY AFTER
    //---- NOTIFYING, SO THE WHILE LOOP AT WORST IS VERY BRIEF.
    //---- IT APPEARS THE MSVC 8.0 COMPILER OPTIMIZES OUT A SIMPLE SPIN
    //---- INSTRUCTION LIKE   while (*_pbOtherNeedsToSynch == true);
    //---- THEREFORE, SLEEP IS USED.
    //----
    //---- NOTE THAT SOME OPERATING SYSTEMS PROVIDE COUNTABLE SPIN LOCKS
    //---- FOR THIS PURPOSE.
    while (*_pbOtherNeedsToSynch == true)
    {
        SystemTimer::sleepMilliSecs (1);
    }
    *_pbOtherNeedsToSynch = true;              //-- RESET THE BARRIER

    _cAssociatedLock.waitForLock();

    #if defined(_WIN32)
        if (bDEBUG_THRDNG)
            printf ("ThDSyn notifying %p\n",&(_pMonitor->cWinMonitor));
        if (SetEvent (_pMonitor->cWinMonitor) == 0)
        {
            cerr << "ERROR: Windows SetEvent failed"
                 << " <ThreadSynchObject::notify>" << endl;
            throw INTERNAL_ERROR;
        }

    #else
        if (bDEBUG_THRDNG)
            printf ("ThDSyn notifying %p\n",&(_pMonitor->cCondVar));
        pthread_mutex_lock (&(_pMonitor->cCondMutex));
        _pMonitor->bNeedToWait = false;
        int  nTmp = pthread_cond_signal (&(_pMonitor->cCondVar));
        if (nTmp != 0)
        {
            cerr << "ERROR: Unix pthread_cond_signal returned " << nTmp
                 << " <ThreadSynchObject::notify>" << endl;
            throw INTERNAL_ERROR;
        }
        pthread_mutex_unlock (&(_pMonitor->cCondMutex));

    #endif

    return;
}


}     //-- namespace HOPSPACK

#endif     //-- HAVE_MT
