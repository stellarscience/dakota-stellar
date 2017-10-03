// $Id: HOPSPACK_ThreadSimpleLock.cpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_ThreadSimpleLock.cpp $

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
  @file HOPSPACK_ThreadSimpleLock.cpp
  @brief Implement HOPSPACK::ThreadSimpleLock
*/

#include "HOPSPACK_common.hpp"

//---- THIS CLASS IS DEFINED ONLY FOR MULTI-THREADED BUILDS.
#if defined(HAVE_MT)

#if defined(_WIN32)
  #include <windows.h>
#else
  #include <errno.h>
  #include <pthread.h>
#endif

#include "HOPSPACK_ThreadSimpleLock.hpp"

namespace HOPSPACK
{

//---- ENABLE THIS TO HELP DEBUG THREADING PROBLEMS.
static const bool  bDEBUG_THRDNG = false;


//----------------------------------------------------------------------
//  Private structures
//----------------------------------------------------------------------

//---- HIDE OPERATING SYSTEM DETAILS IN THESE STRUCTURES.
struct _systemLockType
{
    #if defined(_WIN32)
        CRITICAL_SECTION  cCritSctn;
    #else
        pthread_mutex_t  cMutex;
    #endif
};


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
ThreadSimpleLock::ThreadSimpleLock (void)
{
    _pLock = new _systemLockType;
    #if defined(_WIN32)
        InitializeCriticalSection (&(_pLock->cCritSctn));

        //---- THE CONSTRUCTING THREAD SHOULD INITIALLY OWN THE LOCK.
        waitForLock();

    #else
        if (pthread_mutex_init (&(_pLock->cMutex), NULL) != 0)
        {
            cerr << "ERROR: Could not create Unix mutex." << endl;
            throw INTERNAL_ERROR;
        }

        //---- THE CONSTRUCTING THREAD SHOULD INITIALLY OWN THE LOCK.
        waitForLock();
    #endif

    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
ThreadSimpleLock::~ThreadSimpleLock (void)
{
    #if defined(_WIN32)
        DeleteCriticalSection (&(_pLock->cCritSctn));

    #else
        //---- EBUSY MEANS THE MUTEX IS LOCKED OR PART OF A WAIT.
        int  nStatus = pthread_mutex_destroy (&(_pLock->cMutex));
        if ((nStatus != 0) && (nStatus != EBUSY))
        {
            cerr << "ERROR: Could not destroy Unix mutex, error = "
                 << nStatus << endl;
            throw INTERNAL_ERROR;
        }

    #endif

    delete _pLock;

    return;
}


//----------------------------------------------------------------------
//  Method waitForLock
//----------------------------------------------------------------------
void  ThreadSimpleLock::waitForLock (void)
{
    if (bDEBUG_THRDNG)
        printf ("ThDLck enter waitForLock %p\n", this);

    #if defined(_WIN32)
        EnterCriticalSection (&(_pLock->cCritSctn));

    #else
        pthread_mutex_lock (&(_pLock->cMutex));

    #endif

    return;
}


//----------------------------------------------------------------------
//  Method releaseLock
//----------------------------------------------------------------------
void  ThreadSimpleLock::releaseLock (void)
{
    if (bDEBUG_THRDNG)
        printf ("ThDLck enter releaseLock %p\n", this);

    #if defined(_WIN32)
        LeaveCriticalSection (&(_pLock->cCritSctn));

    #else
        pthread_mutex_unlock (&(_pLock->cMutex));

    #endif

    return;
}


}     //-- namespace HOPSPACK

#endif     //-- HAVE_MT
