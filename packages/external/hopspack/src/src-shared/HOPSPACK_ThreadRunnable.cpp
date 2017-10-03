// $Id: HOPSPACK_ThreadRunnable.cpp 169 2010-04-21 19:35:53Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_ThreadRunnable.cpp $

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
  @file HOPSPACK_ThreadRunnable.cpp
  @brief Implement HOPSPACK::ThreadRunnable
*/

#include "HOPSPACK_common.hpp"

//---- THIS CLASS IS DEFINED ONLY FOR MULTI-THREADED BUILDS.
#if defined(HAVE_MT)

#if defined(_WIN32)
  #include <windows.h>
#else
  #include <pthread.h>
  #include <time.h>
#endif

#include "HOPSPACK_ThreadRunnable.hpp"

namespace HOPSPACK
{

//---- ENABLE THIS TO HELP DEBUG THREADING PROBLEMS.
static const bool  bDEBUG_THRDNG = false;


//----------------------------------------------------------------------
//  Private structures
//----------------------------------------------------------------------

//---- HIDE OPERATING SYSTEM DETAILS IN THESE STRUCTURES.
struct _systemThreadType
{
    bool  bIsThreadDefined;
    #if defined(_WIN32)
        HANDLE  hThreadHandle;
    #else
        pthread_t  pThreadID;
    #endif
};


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
ThreadRunnable::ThreadRunnable (void)
{
    _bIsStarted = false;
    _nnThreadID = 0;

    _pSystemThread = new _systemThreadType;
    _pSystemThread->bIsThreadDefined = false;

    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
ThreadRunnable::~ThreadRunnable (void)
{
    if (_pSystemThread->bIsThreadDefined == true)
    {
    #if defined(_WIN32)
        if (_pSystemThread->hThreadHandle != NULL)
        {
            //---- THIS IS HOW WINDOWS "JOINS" THE THREAD.
            WaitForSingleObject (_pSystemThread->hThreadHandle, INFINITE);
            CloseHandle (_pSystemThread->hThreadHandle);
        }
    #else
        //---- POSIX IMPLICITLY CALLS pthread_exit() WHEN THE
        //---- EXECUTING THREAD COMPLETES, SO "JOIN" SHOULD BE IMMEDIATE.
        void *  pExecStatus;
        pthread_join (_pSystemThread->pThreadID, &pExecStatus);
    #endif
    }

    delete _pSystemThread;

    return;
}


//----------------------------------------------------------------------
//  Local Method for starting a Unix thread
//----------------------------------------------------------------------
#if !defined(_WIN32)
static void *  localUnixThreadFunction (void *  pParam)
{
    ThreadRunnable *  pThis = (ThreadRunnable *) pParam;
    pThis->run();
    return( NULL );
}
#endif


//----------------------------------------------------------------------
//  Local Method for starting a Windows thread
//----------------------------------------------------------------------
#if defined(_WIN32)
static DWORD WINAPI  localWinThreadFunction (LPVOID  lpParam)
{
    ThreadRunnable *  pThis = (ThreadRunnable *) lpParam;
    pThis->run();
    return( 0 );
}
#endif


//----------------------------------------------------------------------
//  Method start
//----------------------------------------------------------------------
void  ThreadRunnable::start (void)
{
    if (_bIsStarted)
        return;

    //---- ALL OPERATING SYSTEMS CREATE A NEW THREAD THAT STARTS RUNNING A
    //---- LOCALLY DEFINED C FUNCTION.  THE FUNCTION CALLS THE run() METHOD
    //---- OF THE SUBCLASS THAT INVOKED THIS METHOD.  WHEN run() RETURNS,
    //---- THREAD EXECUTION ENDS; HOWEVER, THE DESTRUCTOR SHOULD BE CALLED
    //---- FROM THE MAIN THREAD TO CLEAN UP.

    #if defined(_WIN32)
        int  nStackSize = 0;  //-- DEFAULT IS 1 MByte
        DWORD  dwThreadID;
        _pSystemThread->hThreadHandle = CreateThread (NULL,
                                                      nStackSize,
                                                      localWinThreadFunction,
                                                      this,
                                                      0,
                                                      &dwThreadID);
        if (_pSystemThread->hThreadHandle == NULL)
        {
            cerr << "ERROR: Could not start Windows thread." << endl;
            throw INTERNAL_ERROR;
        }
        _pSystemThread->bIsThreadDefined = true;

        _nnThreadID = (unsigned long int) dwThreadID;

    #else
        //---- POSIX DOES NOT SPECIFY A DEFAULT STACK SIZE, SO MAKE SURE
        //---- IT IS AT LEAST 1 MBYTE.
        const size_t    nTARGET_STACK_SIZE = 1024 * 1024;
        pthread_attr_t  cThreadAttr;
        size_t          nStackSize;
        pthread_attr_init (&cThreadAttr);
        pthread_attr_getstacksize (&cThreadAttr, &nStackSize);
        if (nStackSize < nTARGET_STACK_SIZE)
            nStackSize = nTARGET_STACK_SIZE;
        pthread_attr_setstacksize (&cThreadAttr, nStackSize);

        int  nStatus = pthread_create (&(_pSystemThread->pThreadID),
                                       &cThreadAttr,
                                       localUnixThreadFunction,
                                       (void *) this);
        if (nStatus != 0)
        {
            cerr << "ERROR: Could not start Unix thread." << endl;
            throw INTERNAL_ERROR;
        }
        _pSystemThread->bIsThreadDefined = true;
        pthread_attr_destroy (&cThreadAttr);

        _nnThreadID = (unsigned long int) _pSystemThread->pThreadID;
    #endif

    if (bDEBUG_THRDNG)
        cout << "ThDRun started thread (" << getThreadID() << ")\n";

    return;
}


//----------------------------------------------------------------------
//  Method getThreadID
//----------------------------------------------------------------------
unsigned long int  ThreadRunnable::getThreadID (void) const
{
    return( _nnThreadID );
}


}     //-- namespace HOPSPACK

#endif     //-- HAVE_MT
