// $Id: HOPSPACK_ExecutorMpi.hpp 220 2014-01-02 21:24:59Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-executor/HOPSPACK_ExecutorMpi.hpp $

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
  @file HOPSPACK_ExecutorMpi.hpp
  @brief Declaration for HOPSPACK::ExecutorMpi, subclass of Executor.
*/

#ifndef HOPSPACK_EXECUTORMPI_HPP
#define HOPSPACK_EXECUTORMPI_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Executor.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_SystemTimer.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Implements Executor as a single thread that tasks multiple processes.
/*! The MPI implementation sends messages to other processes to perform
 *  evaluations when submit() is called.
 */
//----------------------------------------------------------------------
class ExecutorMpi : public Executor
{
  public:

    //! Constructor.
    ExecutorMpi (void);

    //! Destructor.
    ~ExecutorMpi (void);

    //! Initialize the instance.  Print error messages if necessary.
    /*!
     *  @param[in] naWorkerIDs  List of MPI process rank numbers corresponding
     *                          to evaluation workers.  The executor will
     *                          communicate with these MPI IDs.
     *  @return true            If initialized successfully.
     */
    bool  initialize (const vector< int > &  naWorkerIDs);

    //! Submit ("spawn") a point for evaluation.
    /*!
     *  Clients should first call isReadyForWork() to verify that the executor
     *  is ready; otherwise, they may be blocked until resources for accepting
     *  the point are available.
     *
     *  The MPI executor finds the first available worker and sends it data
     *  for parallel evaluation.  If no workers are available, then data is
     *  not sent and the method immediately returns false.
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
     *  The MPI executor returns false if all workers are busy.
     */
    bool  isReadyForWork (void) const;

    //! Return any completed objective and constraint evaluations.
    /*!
     *  If an evaluation has completed, fill in the references that are passed
     *  to the method with results.  Do not block the caller while waiting for
     *  results if the instance can accept more work.
     *
     *  The MPI executor returns results from the first evaluation message
     *  it finds.  If no messages are pending and isReadyForWork() is false,
     *  then the method waits until a message is received.
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
     *                      otherwise, the worker ID (always 1).
     */
    int  recv (int    &  nTag,
               Vector &  cF,
               Vector &  cEqs,
               Vector &  cIneqs,
               string &  sMsg);

    //! Return the string from parameter 'Evaluator Type'.
    string  getEvaluatorType (void) const;

    //! Print debug information about the Executor instance.
    void  printDebugInfo (void) const;

    //! Print timing information, usually when the Executor is finished.
    void  printTimingInfo (void) const;


  private:

    //! By design, there is no copy constructor.
    ExecutorMpi (const ExecutorMpi &);
    //! By design, there is no assignment operator.
    ExecutorMpi & operator= (const ExecutorMpi &);


    bool           _bIsInitialized;
    SystemTimer *  _pTimers;

    vector< int >  _naWorkerIDs;      //-- WORKER MPI PROCESSOR RANK IDS
    vector< int >  _naWorkerTag;      //-- WORKER'S CURRENTLY ASSIGNED TAG

};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_EXECUTORMPI_HPP
