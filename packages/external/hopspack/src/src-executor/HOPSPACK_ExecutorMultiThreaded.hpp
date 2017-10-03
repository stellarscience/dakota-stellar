// $Id: HOPSPACK_ExecutorMultiThreaded.hpp 220 2014-01-02 21:24:59Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-executor/HOPSPACK_ExecutorMultiThreaded.hpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2014 Sandia Corporation
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
  @file HOPSPACK_ExecutorMultiThreaded.hpp
  @brief Declaration for HOPSPACK::ExecutorMultiThreaded, subclass of Executor.
*/

#ifndef HOPSPACK_EXECUTORMULTITHREADED_HPP
#define HOPSPACK_EXECUTORMULTITHREADED_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Evaluator.hpp"
#include "HOPSPACK_Executor.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_SystemTimer.hpp"
#include "HOPSPACK_ThreadedEvalWorker.hpp"
#include "HOPSPACK_ThreadSimpleLock.hpp"
#include "HOPSPACK_ThreadSynchObject.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Implements Executor as a thread pool of Evaluator workers.
/*! The multi-threaded implementation pushes work to Evaluators when submit()
 *  is called.
 */
//----------------------------------------------------------------------
class ExecutorMultiThreaded : public Executor
{
  public:

    //! Constructor.
    ExecutorMultiThreaded (void);

    //! Destructor.
    ~ExecutorMultiThreaded (void);

    //! Initialize the instance.  Print error messages if necessary.
    /*!
     *  @param[in] nNumWorkers  Number of worker threads to start.
     *  @param[in] cEvalParams  Parameters in the "Evaluator" sublist.
     *  @return true            If initialized successfully.
     */
    bool  initialize (const int              nNumWorkers,
                      const ParameterList &  cEvalParams);

    //! Submit ("spawn") a point for evaluation.
    /*!
     *  Clients should first call isReadyForWork() to verify that the executor
     *  is ready; otherwise, they may be blocked until resources for accepting
     *  the point are available.
     *
     *  The multi-threaded executor finds the first available worker and signals
     *  it to begin parallel evaluation.  If no workers are available, then data
     *  is ignored and the method immediately returns false.
     *
     *  @param[in] nTag      Contains a unique tag for the evaluation which
     *                       can be used to name files, etc.
     *  @param[in] cX        The point at which to evaluate the function(s).
     *  @param[in] nRequest  Type of evaluation information requested.
     *  @return true         If the point is accepted.
     *  @throws FATAL_ERROR  If nRequest is undefined.
     */
    bool  submit (const int                nTag,
                  const Vector          &  cX,
                  const EvalRequestType    nRequest);

    //! Return true if the executor is ready for more points to be submitted.
    /*!
     *  The multi-threaded executor returns false if all workers are busy.
     */
    bool  isReadyForWork (void) const;

    //! Return any completed objective and constraint evaluations.
    /*!
     *  If an evaluation has completed, fill in the references that are passed
     *  to the method with results.  Do not block the caller while waiting for
     *  results if the instance can accept more work.
     *
     *  The multi-threaded executor returns results from the first evaluation
     *  result it finds, which is determined by the operating system if multiple
     *  workers are finished.  If no messages are pending and isReadyForWork()
     *  is false, then the method waits until a worker signals completion.
     *
     *  @param[out] nTag    Tag corresponding to the point that was evaluated.
     *  @param[out] cF      Vector of returned objective function values.
     *                      The vector may be empty, or specific elements may
     *                      contain the value HOPSPACK::dne() to indicate an
     *                      objective does not exist at the point.
     *  @param[out] cEqs    Vector of returned nonlinear equality constraint
     *                      values.
     *                      The vector may be empty, or specific elements may
     *                      contain the value HOPSPACK::dne() to indicate an
     *                      objective does not exist at the point.
     *  @param[out] cIneqs  Vector of returned nonlinear inequality constraint
     *                      values.
     *                      The vector may be empty, or specific elements may
     *                      contain the value HOPSPACK::dne() to indicate an
     *                      objective does not exist at the point.
     *  @param[out] sMsg    Message returned by the evaluation routine.
     *  @return             0 (zero) if no evaluation result was received;
     *                      otherwise, the worker ID (not the thread ID).
     */
    int  recv (int    &  nTag,
               Vector &  cF,
               Vector &  cEqs,
               Vector &  cIneqs,
               string &  sMsg);

    //! Stop all worker threads so Executor can be deleted.
    /*!
     *  @return  True if successful.
     *
     *  In the current implementation, a worker cannot be stopped if it is busy
     *  with an evaluation.  The method returns false if some worker is busy.
     */
    bool  shutdown (void);

    //! Return the string from parameter 'Evaluator Type'.
    string  getEvaluatorType (void) const;

    //! Print debug information about the Executor instance.
    void  printDebugInfo (void) const;

    //! Print timing information, usually when the Executor is finished.
    void  printTimingInfo (void) const;


  private:

    //! By design, there is no copy constructor.
    ExecutorMultiThreaded (const ExecutorMultiThreaded &);
    //! By design, there is no assignment operator.
    ExecutorMultiThreaded & operator= (const ExecutorMultiThreaded &);


    bool           _bIsInitialized;
    int            _nNumWorkers;
    SystemTimer *  _pTimers;

    //! List of workers that have evaluation results waiting for recv().
    vector< int >  _naWorkersDoneAndWaiting;

    //! Thread pool of evaluation workers.
    ThreadedEvalWorker **  _paWorkers;

    //! Evaluators, one for each worker.
    Evaluator **  _paEvaluators;

    //! The executor thread waits for this.
    ThreadSynchObject **  _paTsoWorkerDone;

    //! Workers wait for this.
    ThreadSimpleLock  **  _paLockTsoNewPointLoaded;
    ThreadSynchObject **  _paTsoNewPointLoaded;

    //! Synchronization barrier flags.
    bool *  _baExecNeedsToSynch;
    bool *  _baWrkrNeedsToSynch;

    //! If true, then the worker is idle.
    bool *  _baIsWorkerIdle;

    //! If true, then the worker is done and has notified the executor thread.
    bool *  _baIsWorkerDone;

    //! Guards access to _baIsWorkerDone.
    ThreadSimpleLock *  _pLockIsWorkerDoneFlags;

    //! Synchronization shutdown flags.
    bool    _bIsTimeToShutdown;
    bool *  _baIsWorkerShutdown;

    //! Storage for information passed to an evaluator.
    int             *  _naSubmitTags;
    Vector          *  _caSubmitPoints;
    EvalRequestType *  _naSubmitRequestTypes;

    //! Storage for information returned by an evaluator.
    Vector          *  _caRecvF;
    Vector          *  _caRecvEqs;
    Vector          *  _caRecvIneqs;
    string          *  _saRecvMsgs;
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_EXECUTORMULTITHREADED_HPP
