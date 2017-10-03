// $Id: HOPSPACK_ThreadSimpleLock.hpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_ThreadSimpleLock.hpp $

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
  @file HOPSPACK_ThreadSimpleLock.hpp
  @brief Class declaration for HOPSPACK::ThreadSimpleLock.
*/
#ifndef HOPSPACK_THREADSIMPLELOCK_HPP
#define HOPSPACK_THREADSIMPLELOCK_HPP

#include "HOPSPACK_common.hpp"

//---- THIS CLASS IS DEFINED ONLY FOR MULTI-THREADED BUILDS.
#if defined(HAVE_MT)

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Defines a simple lock/barrier for synchronization between two threads.
/*!
 *  The constructing thread owns the original lock.  Two threads should
 *  take turns releasing and waiting for the lock.
 *
 *  The Windows implementation uses the lightweight WINAPI Critical Section.
 *  It is simple and fast, but does not provide a single atomic wait over
 *  a set of lock, or waiting with a timeout.
 *  Every entry into a Critical Section must be paired with a leave; ie,
 *  the operating system counts each entry and leave.
 */
//----------------------------------------------------------------------
class ThreadSimpleLock
{
  public:

    //! Constructor.  The lock is owned by the calling thread after construction.
    ThreadSimpleLock (void);

    //! Destructor.
    /*!
     *  Be careful not to destroy an instance if some thread is still waiting.
     */
    ~ThreadSimpleLock (void);

    //! Block the calling thread until the lock is acquired.
    /*!
     *  The calling code must ensure it does not already own the lock.
     */
    void  waitForLock (void);

    //! Release ownership of the lock.
    /*!
     *  The calling code must ensure it already owns the lock.
     */
    void  releaseLock (void);


  private:

    //! By design, there is no copy constructor.
    ThreadSimpleLock (const ThreadSimpleLock &);
    //! By design, there is no assignment operator.
    ThreadSimpleLock & operator= (const ThreadSimpleLock &);


    //! Hide operating system details in this private type.
    typedef struct _systemLockType _systemLockType;
    _systemLockType *  _pLock;
};

}          //-- namespace HOPSPACK

#endif     //-- HAVE_MT

#endif     //-- HOPSPACK_THREADSIMPLELOCK_HPP
