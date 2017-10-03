// $Id: HOPSPACK_SystemTimer.cpp 169 2010-04-21 19:35:53Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_SystemTimer.cpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2010 Sandia Corporation
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
  @file HOPSPACK_SystemTimer.cpp
  @brief Implement HOPSPACK::SystemTimer.
*/

#if defined(_WIN32)
  #include <windows.h>
#else
  #include <time.h>
  #include <sys/time.h>
#endif

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_SystemTimer.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Private structures
//----------------------------------------------------------------------

#if defined(HAVE_REALTIME_CLOCK)
//---- HIDE OPERATING SYSTEM DETAILS IN THESE STRUCTURES.
struct _systemTimerRealType
{
    #if defined(_WIN32)
        DWORD  nnRealTime;      //-- PRIMITIVE RETURNED BY WINDOWS TIMER CALL
    #else
        timeval  cRealTime;     //-- PRIMITIVE RETURNED BY UNIX TIMER CALL
    #endif
};
#endif


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
SystemTimer::SystemTimer (const int  nNumTimers)
{
#if defined(HAVE_REALTIME_CLOCK)
    if (nNumTimers <= 0)
    {
        _nNumTimers = 0;
        return;
    }
    _nNumTimers = nNumTimers;

    _baIsStarted = new bool[_nNumTimers];
    _daCumTimes = new double[_nNumTimers];
    _naNumCalls = new int[_nNumTimers];
    _taStartTimes = new _systemTimerRealType[_nNumTimers];

    for (int  i = 0; i < _nNumTimers; i++)
    {
        _baIsStarted[i] = false;
        _daCumTimes[i] = 0.0;
        _naNumCalls[i] = 0;
    }

    //---- WINDOWS timeBeginPeriod IS CALLED IN HOPSPACK_Hopspack.cpp.

#else     //-- HAVE_REALTIME_CLOCK IS UNDEFINED
    _nNumTimers = 0;
#endif

    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
SystemTimer::~SystemTimer (void)
{
    if (_nNumTimers == 0)
        return;

    delete[] _baIsStarted;
    delete[] _daCumTimes;
    delete[] _naNumCalls;

#if defined(HAVE_REALTIME_CLOCK)
    delete[] _taStartTimes;
#endif

    return;
}


//----------------------------------------------------------------------
//  Method getDateTime
//----------------------------------------------------------------------
void  SystemTimer::getDateTime (string &  sCurrentDateTime)
{
    #if defined(_WIN32)
        SYSTEMTIME  tNow;
        GetLocalTime (&tNow);
        char  szTmp[25];
        #if defined(HAVE_MSVC_SECURE_STRING_FNS)
            sprintf_s (szTmp, 25, "%2d/%02d/%4d %02d:%02d:%02d",
                       (int) tNow.wMonth,
                       (int) tNow.wDay,
                       (int) tNow.wYear,
                       (int) tNow.wHour,
                       (int) tNow.wMinute,
                       (int) tNow.wSecond);
        #else
            sprintf (szTmp, "%2d/%02d/%4d %02d:%02d:%02d",
                     (int) tNow.wMonth,
                     (int) tNow.wDay,
                     (int) tNow.wYear,
                     (int) tNow.wHour,
                     (int) tNow.wMinute,
                     (int) tNow.wSecond);
        #endif
        sCurrentDateTime = szTmp;

    #else
        time_t  tNow = time (NULL);
        struct tm  tConvertedNow;
        if (localtime_r (&tNow, &tConvertedNow) == NULL)
        {
            sCurrentDateTime = "Error getting time";
        }
        else
        {
            char  szTmp[25];
            sprintf (szTmp, "%2d/%02d/%4d %02d:%02d:%02d",
                     tConvertedNow.tm_mon + 1,
                     tConvertedNow.tm_mday,
                     tConvertedNow.tm_year + 1900,
                     tConvertedNow.tm_hour,
                     tConvertedNow.tm_min,
                     tConvertedNow.tm_sec);
            sCurrentDateTime = szTmp;
        }
    #endif

    return;
}


//----------------------------------------------------------------------
//  Method sleepMilliSecs
//----------------------------------------------------------------------
void  SystemTimer::sleepMilliSecs (const int  nMilliSecs)
{
    #if defined(_WIN32)
        //---- NOTE, WINDOWS DOES SOMETHING EVEN IF ZERO.
        Sleep ((DWORD) nMilliSecs);
    #else
        struct timespec  cTS;
        cTS.tv_sec  = nMilliSecs / 1000;
        cTS.tv_nsec = (nMilliSecs % 1000) * 1000000;
        nanosleep (&cTS, NULL);
    #endif

    return;
}


//-----------------------------------------------------------------------------
//  Method start
//-----------------------------------------------------------------------------
bool  SystemTimer::start (const int  nTimerID)
{
    if ((nTimerID < 0) || (nTimerID >= _nNumTimers))
        return( false );

#if defined(HAVE_REALTIME_CLOCK)
    //---- READ AND STORE THE CURRENT TIME.
    #if defined(_WIN32)
        _taStartTimes[nTimerID].nnRealTime = timeGetTime();

    #else
        gettimeofday (&(_taStartTimes[nTimerID].cRealTime), NULL);

    #endif

    _baIsStarted[nTimerID] = true;
    return( true );
#else
    return( false );
#endif
}


//-----------------------------------------------------------------------------
//  Method stop
//-----------------------------------------------------------------------------
bool  SystemTimer::stop (const int  nTimerID)
{
    if ((nTimerID < 0) || (nTimerID >= _nNumTimers))
        return( false );

    if (_baIsStarted[nTimerID] == false)
        return( false );

    //---- ADD ELAPSED TIME SINCE THE LAST CALL TO start().
    _daCumTimes[nTimerID] += getTimeSinceLastStart_ (nTimerID);
    _baIsStarted[nTimerID] = false;
    _naNumCalls[nTimerID]++;

    return( true );
}


//-----------------------------------------------------------------------------
//  Method getTotalTime
//-----------------------------------------------------------------------------
double  SystemTimer::getTotalTime (const int  nTimerID) const
{
    if ((nTimerID < 0) || (nTimerID >= _nNumTimers))
        return( -1.0 );

    if ((getNumStarts (nTimerID) == 0) && (_baIsStarted[nTimerID] == false))
        return( 0.0 );

    double  dResult = _daCumTimes[nTimerID];
    if (_baIsStarted[nTimerID] == true)
        dResult += getTimeSinceLastStart_ (nTimerID);
    return( dResult );
}


//-----------------------------------------------------------------------------
//  Method getNumStarts
//-----------------------------------------------------------------------------
int  SystemTimer::getNumStarts (const int  nTimerID) const
{
    if ((nTimerID < 0) || (nTimerID >= _nNumTimers))
        return( -1 );

    return( _naNumCalls[nTimerID] );
}


//-----------------------------------------------------------------------------
//  Method getAvgTime
//-----------------------------------------------------------------------------
double  SystemTimer::getAvgTime (const int  nTimerID) const
{
    if ((nTimerID < 0) || (nTimerID >= _nNumTimers))
        return( -1.0 );

    if (getNumStarts (nTimerID) == 0)
        return( 0.0 );

    return( _daCumTimes[nTimerID] / ((double) getNumStarts (nTimerID)) );
}


//-----------------------------------------------------------------------------
//  Method reset
//-----------------------------------------------------------------------------
void  SystemTimer::reset (const int  nTimerID)
{
    if ((nTimerID < 0) || (nTimerID >= _nNumTimers))
        return;

    _daCumTimes[nTimerID] = 0.0;
    _baIsStarted[nTimerID] = false;
    _naNumCalls[nTimerID] = 0;
    return;
}


//-----------------------------------------------------------------------------
//  Private Method getTimeSinceLastStart_
//-----------------------------------------------------------------------------
double  SystemTimer::getTimeSinceLastStart_ (const int  nTimerID) const
{
#if defined(HAVE_REALTIME_CLOCK)
    //---- READ THE CURRENT TIME AND COMPUTE ELAPSED TIME.
    #if defined(_WIN32)
        //---- WINDOWS TIMER RETURNS MILLISECOND TICKS.
        DWORD  nnNow = timeGetTime();
        DWORD  nnElapsed = nnNow - _taStartTimes[nTimerID].nnRealTime;
        return( ((double) nnElapsed) / 1000.0 );

    #else
        timeval  cNow;
        gettimeofday (&cNow, NULL);
        time_t  nNowSecs = cNow.tv_sec;
        time_t  nStartSecs = _taStartTimes[nTimerID].cRealTime.tv_sec;
        long int  nnNowUsecs = cNow.tv_usec;
        long int  nnStartUsecs = _taStartTimes[nTimerID].cRealTime.tv_usec;

        time_t  nDiffSecs = nNowSecs - nStartSecs;
        long int  nnDiffUsecs;
        if (nnNowUsecs >= nnStartUsecs)
                nnDiffUsecs = nnNowUsecs - nnStartUsecs;
        else
        {
            nDiffSecs--;
            nnDiffUsecs = 1000000 - (nnStartUsecs - nnNowUsecs);
        }
        double  dDiff =   ((double) nDiffSecs)
                        + ((double) nnDiffUsecs) * 1.0e-6;
        return( dDiff );

    #endif

#else     //-- HAVE_REALTIME_CLOCK IS UNDEFINED
    return( 0.0 );
#endif
}


}     //-- namespace HOPSPACK
