// $Id: HOPSPACK_ThreadedEvalWorker.hpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-executor/HOPSPACK_ThreadedEvalWorker.hpp $

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
  @file HOPSPACK_ThreadedEvalWorker.hpp
  @brief Declaration for HOPSPACK::ThreadedEvalWorker.
*/

#ifndef HOPSPACK_THREADEDEVALWORKER_HPP
#define HOPSPACK_THREADEDEVALWORKER_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Evaluator.hpp"
#include "HOPSPACK_ThreadRunnable.hpp"
#include "HOPSPACK_ThreadSynchObject.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Implements an Evaluator worker that runs in a separate thread.
//----------------------------------------------------------------------
class ThreadedEvalWorker : public ThreadRunnable
{
  public:

    //! Constructor.
    /*!
     *  @param[in] nWorkerNumber  Unique number of the worker, starting at 1.
     *  @param[in] pEvaluator     Unique evaluator instance for the worker.
     *  @param[in] pnSubmitTag    Pointer to storage for tag
     *                            of a point submitted for evaluation
     *  @param[in] pSubmitPointX  Pointer to storage for X position
     *                            of a point submitted for evaluation
     *  @param[in] pnSubmitRequestType  Pointer to storage for request type
     *                            of a point submitted for evaluation
     *  @param[in] pResultF       Pointer to storage for F results
     *                            of a point that was evaluated.
     *  @param[in] pResultEqs     Pointer to storage for equality constraints
     *                            of a point that was evaluated.
     *  @param[in] pResultIneqs   Pointer to storage for inequality constraints
     *                            of a point that was evaluated.
     *  @param[in] pResultMsg     Pointer to storage for message
     *                            of a point that was evaluated.
     *  @param[in] pTsoNewPointLoaded      Synch object used by Executor to
     *                                     notify the worker of a new point.
     *  @param[in] pbIsWorkerIdle          Synch flag set by Executor and
     *                                     tested by the worker.
     *  @param[in] pbIsWorkerDone          Synch flag modified by both Executor
     *                                     and worker.
     *  @param[in] pLockIsWorkerDoneFlags  Synch lock that guards access to
     *                                     pbIsWorkerDone.
     *  @param[in] pTsoWorkerDone          Synch object used by worker to
     *                                     notify the Executor evaluation is done.
     *  @param[in] pbExecNeedsToSynch      Synch flag for barrier.
     *  @param[in] pbWrkrNeedsToSynch      Synch flag for barrier.
     *  @param[in] pbIsTimeToShutdown      Synch flag set by Executor and
     *                                     tested by the worker.
     *  @param[in] pbIsWorkerShutdown      Synch flag set by the worker.
     */
    ThreadedEvalWorker (const int                         nWorkerNumber,
                              Evaluator         *  const  pEvaluator,
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
                              ThreadSynchObject ** const  pTsoWorkerDone,
                              bool              *  const  pbExecNeedsToSynch,
                              bool              *  const  pbWrkrNeedsToSynch,
                        const bool              *  const  pbIsTimeToShutdown,
                              bool              *  const  pbIsWorkerShutdown);
    //! Destructor.
    ~ThreadedEvalWorker (void);

    //! Execute code that runs in a separate thread.
    void  run (void);


  private:

    //! By design, there is no copy constructor.
    ThreadedEvalWorker (const ThreadedEvalWorker &);
    //! By design, there is no assignment operator.
    ThreadedEvalWorker & operator= (const ThreadedEvalWorker &);


    const int    _nWorkerNumber;
    Evaluator *  _pEvaluator;

    //---- STORAGE SHARED WITH THE EXECUTOR.
    const int             * const  _pnSubmitTag;
    const Vector          * const  _pSubmitPointX;
    const EvalRequestType * const  _pnSubmitRequestType;
          Vector          * const  _pResultF;
          Vector          * const  _pResultEqs;
          Vector          * const  _pResultIneqs;
          string          * const  _pResultMsg;

    //---- SYNCHRONIZATION ITEMS SHARED WITH THE EXECUTOR.
    ThreadSynchObject *  const  _pTsoNewPointLoaded;
    const bool        *  const  _pbIsWorkerIdle;
    bool              *  const  _pbIsWorkerDone;
    ThreadSimpleLock  *         _pLockIsWorkerDoneFlags;
    ThreadSimpleLock  *         _pLockTsoWorkerDone;
    ThreadSynchObject **        _ppTsoWorkerDone;
    bool              *  const  _pbExecNeedsToSynch;
    bool              *  const  _pbWrkrNeedsToSynch;
    const bool        *  const  _pbIsTimeToShutdown;
    bool              *  const  _pbIsWorkerShutdown;
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_THREADEDEVALWORKER_HPP
