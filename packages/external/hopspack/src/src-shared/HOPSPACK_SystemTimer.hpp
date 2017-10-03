// $Id: HOPSPACK_SystemTimer.hpp 169 2010-04-21 19:35:53Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_SystemTimer.hpp $

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
  \file HOPSPACK_SystemTimer.hpp
  \brief Class declaration of HOPSPACK::SystemTimer.
*/

#ifndef HOPSPACK_SYSTEMTIMER_HPP
#define HOPSPACK_SYSTEMTIMER_HPP

#include "HOPSPACK_common.hpp"

namespace HOPSPACK 
{


//----------------------------------------------------------------------
//! Defines system-dependent timer for wall clock time and date.
/*!
 *  On Windows the wall clock timer has resolution of 1 millisecond,
 *  on Unix the resolution is 1 microsecond.  The CPU time required to
 *  query a timer is about 100 nanoseconds on Windows, 1.5 microseconds
 *  on Linux (of course it depends on the processor).
 *
 *  Example use:
 *    SystemTimer  timers (5);     //-- DEFINES TIMERS 0..4
 *    timers.start (0);
 *    ...
 *    timers.stop (0);
 *    double  dTimeInSec = timers.getTotalTime (0);
 */
//----------------------------------------------------------------------
class SystemTimer
{
  public:

    //! Constructor.
    /*!
     *  @param[in] nNumTimers  Number of timers to allocate.
     */
    SystemTimer (const int  nNumTimers);

    //! Destructor.
    ~SystemTimer (void);

    //! Static method to return current date and time.
    /*!
     *  @param[out]  sCurrentDateTime  Format will be MM/DD/YYYY HH:MM:SS.
     */
    static void  getDateTime (string &  sCurrentDateTime);

    //! Static method to put the current thread to sleep.
    /*!
     *  @param[in] nMilliSecs  Number of milliseconds to sleep.  If the number
     *                         is zero or negative, then no sleep occurs.
     *                         Exact sleep time depends on the operating system,
     *                         but is typically very close to the target.
     */
    static void  sleepMilliSecs (const int  nMilliSecs);

    //! Start a particular timer.  Should not be called twice in a row.
    /*!
     *  @param[in] nTimerID  The timer of interest (>= 0).
     *  @return              True if successful.
     */
    bool  start (const int  nTimerID);

    //! Stop a particular timer.  Should not be called twice in a row.
    /*!
     *  @param[in] nTimerID  The timer of interest (>= 0).
     *  @return              True if successful.
     */
    bool  stop (const int  nTimerID);

    //! Return the total "start" to "stop" duration for the timer.
    /*!
     *  @param[in] nTimerID  The timer of interest (>= 0).
     *  @return              Total elapsed time between all calls to start()
     *                       and stop() since the class instance was constructed,
     *                       or since reset() was called.  If there is currently
     *                       a start() with no matching stop(), then the total
     *                       also includes all time since the last start().
     */
    double  getTotalTime (const int  nTimerID) const;

    //! Return the number of times start() was called for the timer.
    int  getNumStarts (const int  nTimerID) const;

    //! Return the average "start" to "stop" duration for the timer.
    /*!
     *  @param[in] nTimerID  The timer of interest (>= 0).
     *  @return              Average elapsed time between all calls to start()
     *                       and stop() since the class instance was constructed,
     *                       or since reset() was called.  Unlike getTotalTime(),
     *                       if there is currently a start() with no matching
     *                       stop(), then the total ignores time since this
     *                       last start().
     *
     *  Result is the same as computing getTotalTime (n) / getNumStarts (n).
     */
    double  getAvgTime (const int  nTimerID) const;

    //! Reset the timer so it can be started again from zero.
    void  reset (const int  nTimerID);


  private:

    //! By design, there is no copy constructor.
    SystemTimer (const SystemTimer &);
    //! By design, there is no assignment operator.
    SystemTimer & operator= (const SystemTimer &);

    //! Return the time in seconds since the last call to start().
    double  getTimeSinceLastStart_ (const int  nTimerID) const;


    //! Hide operating system details in this private type.
    typedef struct _systemTimerRealType _systemTimerRealType;
    _systemTimerRealType *  _taStartTimes;

    int       _nNumTimers;
    bool *    _baIsStarted;
    double *  _daCumTimes;
    int *     _naNumCalls;
};

}          //-- namespace HOPSPACK

#endif     //-- SYSTEMTIMER_HPP
