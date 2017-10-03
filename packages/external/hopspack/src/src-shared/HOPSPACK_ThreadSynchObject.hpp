// $Id: HOPSPACK_ThreadSynchObject.hpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_ThreadSynchObject.hpp $

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
  @file HOPSPACK_ThreadSynchObject.hpp
  @brief Class declaration for HOPSPACK::ThreadSynchObject.
*/
#ifndef HOPSPACK_THREADSYNCHOBJECT_HPP
#define HOPSPACK_THREADSYNCHOBJECT_HPP

#include "HOPSPACK_common.hpp"

//---- THIS CLASS IS DEFINED ONLY FOR MULTI-THREADED BUILDS.
#if defined(HAVE_MT)

#include "HOPSPACK_ThreadSimpleLock.hpp"


namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Defines an object for synchronizing execution between threads.
/*!
 *  Threads synchronize by acquiring and releasing ownership of an instance of
 *  this class.  A call to waitForOwnership() will acquire the underlying
 *  monitor (typically a mutex).  Calling notify() will release the monitor
 *  lock and allow a waiting thread to acquire ownership.  To synchronize safely,
 *  the calls use a simple lock with method synchronizeLock(), and a boolean.
 *
 *  Example usage, signaling from thread 1 to thread 2:
 *    ThreadSimpleLock   cLock
 *    cLock.releaseLock
 *    bool  b2needsToSynch = false
 *    ThreadSynchObject  cSynchObj (cLock, *b2needsToSynch)
 *    thread 1:
 *      loop
 *        ...set up for message to thread 2
 *        cSynchObj.notify
 *        ...possibly wait for a synchronizing reply
 *        cSynchObj.takeOwnership
 *    thread 2 (cSynchObj):
 *      loop
 *        cSynchObj.waitForNotify
 *        b2needsToSynch = true
 *        ...process message from thread 1
 *        cSynchObj.synchronizeLock
 *        b2needsToSynch = false
 *  In the example, the lock ensures thread 1 will have ownership before
 *  thread 2 calls waitForOwnership.  The boolean ensures thread 1 does not
 *  notify until thread 2 calls synchronizeLock.
 *
 *  The Windows implementation uses an interprocess WINAPI Event object.
 *  Since there is only one process, the lightweight WINAPI Critical Section
 *  was considered; however, Critical Sections do not allow a single atomic
 *  wait over a set of objects, and do not support waiting with a timeout.
 *
 *  The Unix implementation uses a condition variable with associated mutex
 *  and boolean.  The mutex serves only to safeguard updating of the boolean.
 *  Other external locks are used for thread synchronization.
 *  I must admit there might be a simpler implementation of thread
 *  synchronization (on both platforms) using external condition variables.
 */
//----------------------------------------------------------------------
class ThreadSynchObject
{
  public:

    //! Constructor.
    /*!
     *  @param[in] synchLock    Simple lock that guards access to the monitor
     *                          element of the instance.  There should be one
     *                          lock per instance.
     *  @param[in] pbSynchFlag  Read-only flag that will be set by another thread
     *                          to enable synchronizeLock().
     *
     *  Initial ownership belongs to the calling thread.  The next call should
     *  therefore be notify() to release ownership.
     *  Initial ownership of the lock is decided by the calling thread;
     *  one reason it is passed instead of hidden internally.
     */
    ThreadSynchObject (ThreadSimpleLock &  synchLock,
                       bool * const        pbSynchFlag);

    //! Destructor.
    /*!
     *  Be careful not to destroy an instance if any threads are still waiting.
     */
    ~ThreadSynchObject (void);

    //! Wait until the simple lock can be acquired and released.
    /*!
     *  This call effects a synchronizing barrier with the simple lock, and
     *  is usually called just before waitForOwnership().  It is made available
     *  as a separate method for flexibility in initializing wait loops.
     *  The calling code must ensure it does not already have ownership.
     */
    void  synchronizeLock (void);

    //! Block the calling thread until the monitor is notified.
    /*!
     *  The calling code must ensure it does not already have ownership,
     *  or no other thread will be able to notify it.  Upon return, the
     *  caller does not have ownership.
     */
    void  waitForNotify (void);

    //! Take ownership of the monitor.
    /*!
     *  The calling code must ensure it does not already have ownership.
     *  Underlying code blocks if necessary to take ownership and then unlocks
     *  the synchronizing simple lock.
     */
    void  takeOwnership (void);

    //! Release ownership of the monitor immediately and return.
    /*!
     *  The calling code must ensure it already has ownership.
     *  Underlying code first waits for the synchronizing boolean to be true;
     *  ie, for the other thread to signal that it has acquired and released
     *  the synchronizing simple lock.  Then it obtains control of the
     *  synchronizing simple lock.
     */
    void  notify (void);


  private:

    //! By design, there is no copy constructor.
    ThreadSynchObject (const ThreadSynchObject &);
    //! By design, there is no assignment operator.
    ThreadSynchObject & operator= (const ThreadSynchObject &);


    //! Hide operating system details in this private type.
    typedef struct _systemMonitorType _systemMonitorType;
    _systemMonitorType *  _pMonitor;

    ThreadSimpleLock &  _cAssociatedLock;
    bool * const        _pbOtherNeedsToSynch;
};

}          //-- namespace HOPSPACK

#endif     //-- HAVE_MT

#endif     //-- HOPSPACK_THREADSYNCHOBJECT_HPP
