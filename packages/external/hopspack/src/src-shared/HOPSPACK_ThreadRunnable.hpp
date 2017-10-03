// $Id: HOPSPACK_ThreadRunnable.hpp 169 2010-04-21 19:35:53Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_ThreadRunnable.hpp $

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
  @file HOPSPACK_ThreadRunnable.hpp
  @brief Class declaration for HOPSPACK::ThreadRunnable.
*/
#ifndef HOPSPACK_THREADRUNNABLE_HPP
#define HOPSPACK_THREADRUNNABLE_HPP

#include "HOPSPACK_common.hpp"

//---- THIS CLASS IS DEFINED ONLY FOR MULTI-THREADED BUILDS.
#if defined(HAVE_MT)

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Defines an interface for classes that run in their own thread of execution.
/*! This simple interface defines an abstract run() method which should be
 *  overloaded by subclasses to execute their code.  Subclasses typically
 *  include a constructor with arguments for synchronizing the thread.
 *  Thread execution ends when the run() method returns.
 */
//----------------------------------------------------------------------
class ThreadRunnable
{
  public:

    //! Create the thread and call run() to execute code in it.
    /*!
     *  Subclasses call this method to start the thread.  Subclasses should
     *  not overload start(); instead, overload run(), which is called by start().
     *  The method returns as soon as the new thread begins executing.
     *  Do not call this method more than once.
     */
    void  start (void);

    //! Execute code that runs in a separate thread.
    /*!
     *  Subclasses must overload this method to execute code that runs
     *  in the separate thread.  Do not call this method directly; instead,
     *  invoke indirectly by calling start().
     */
    virtual void  run (void) = 0;

    //! Return the thread ID.
    unsigned long int  getThreadID (void) const;


  protected:

    //! Constructor.
    ThreadRunnable (void);

    //! Destructor.
    ~ThreadRunnable (void);


  private:

    //! By design, there is no copy constructor.
    ThreadRunnable (const ThreadRunnable &);
    //! By design, there is no assignment operator.
    ThreadRunnable & operator= (const ThreadRunnable &);


    //! Hide operating system details in this private type.
    typedef struct _systemThreadType _systemThreadType;
    _systemThreadType *  _pSystemThread;

    bool               _bIsStarted;
    unsigned long int  _nnThreadID;
};

}          //-- namespace HOPSPACK

#endif     //-- HAVE_MT

#endif     //-- HOPSPACK_THREADRUNNABLE_HPP
