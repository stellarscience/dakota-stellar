// $Id: HOPSPACK_Executor.hpp 220 2014-01-02 21:24:59Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-executor/HOPSPACK_Executor.hpp $ 

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
  @file HOPSPACK_Executor.hpp
  @brief Interface declaration for HOPSPACK::Executor.
*/

#ifndef HOPSPACK_EXECUTOR_HPP
#define HOPSPACK_EXECUTOR_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Interface that coordinates execution of evaluations.
/*!
 *  Declare an abstract class interface for execution of evaluations.
 *  Implementations will make optimal use of computing resources to evaluate
 *  many points.  Clients can call submit() any number of times, and interlace
 *  calls to recv() at any time.  The submit() call may block, but not if a
 *  prior call to isReadyForWork() returned true.  The recv() call may block
 *  waiting for evaluation results, but not if isReadyForWork() returns true.
 *
 *  An interface design makes it easy to extend HOPSPACK with an Executor
 *  implementation specialized for a particular computing architecture.
 *  New implementations must subclass Executor, implement all methods, and
 *  be used in a new "main()" program that constructs and uses an instance
 *  of HOPSPACK::Hopspack.
 */
//----------------------------------------------------------------------
class Executor
{
  public:

    //! Submit ("spawn") a point for evaluation.
    /*!
     *  Clients should first call isReadyForWork() to verify that the executor
     *  is ready; otherwise, they may be blocked until resources for accepting
     *  the point are available.
     *
     *  @param[in] nTag      Contains a unique tag for the evaluation which
     *                       can be used to name files, etc.
     *  @param[in] cX        The point at which to evaluate the function(s).
     *  @param[in] nRequest  Type of evaluation information requested.
     *  @return true         If the point is accepted, false if unable to
     *                       accept more work.
     *  @throws FATAL_ERROR  If nRequest is undefined.
     */
    virtual bool  submit (const int                nTag,
                          const Vector          &  cX,
                          const EvalRequestType    nRequest) = 0;

    //! Return true if the executor is ready for more points to be submitted.
    virtual bool  isReadyForWork (void) const = 0;

    //! Return any completed objective and constraint evaluations.
    /*!
     *  If an evaluation has completed, fill in the references that are passed
     *  to the method with results.  Do not block the caller while waiting for
     *  results if the instance can accept more work.
     *
     *  @param[out] nTag    Tag corresponding to the point that was evaluated.
     *  @param[out] cF      Vector of returned objective function values.
     *                      The vector may be empty, or specific elements may
     *                      contain the value HOPSPACK::dne() to indicate an
     *                      objective does not exist at the point.
     *  @param[out] cEqs    Vector of returned nonlinear equality constraint
     *                      values.
     *                      The vector may be empty, or specific elements may
     *                      contain the value HOPSPACK::dne() to indicate a
     *                      constraint could not be evaluated at the point.
     *  @param[out] cIneqs  Vector of returned nonlinear inequality constraint
     *                      values.
     *                      The vector may be empty, or specific elements may
     *                      contain the value HOPSPACK::dne() to indicate a
     *                      constraint could not be evaluated at the point.
     *  @param[out] sMsg    Message returned by the evaluation routine.
     *  @return             0 (zero) if no evaluation result was received;
     *                      otherwise, the worker ID (a positive integer).
     */
    virtual int  recv (int    &  nTag,
                       Vector &  cF,
                       Vector &  cEqs,
                       Vector &  cIneqs,
                       string &  sMsg) = 0;

    //! Return the string from parameter 'Evaluator Type'.
    virtual string  getEvaluatorType (void) const = 0;

    //! Print debug information about the Executor instance.
    virtual void  printDebugInfo (void) const = 0;

    //! Print timing information, usually when the Executor is finished.
    virtual void  printTimingInfo (void) const = 0;


  protected:

    //! Abstract class requires a virtual destructor.
    virtual ~Executor (void) { };
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_EXECUTOR_HPP
