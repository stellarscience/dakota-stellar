// $Id: HOPSPACK_Evaluator.hpp 220 2014-01-02 21:24:59Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-evaluator/HOPSPACK_Evaluator.hpp $ 

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
  @file HOPSPACK_Evaluator.hpp
  @brief Interface declaration for HOPSPACK::Evaluator.
*/

#ifndef HOPSPACK_EVALUATOR_HPP
#define HOPSPACK_EVALUATOR_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Interface that evaluates functions and nonlinear constraints.
/*! Declare an abstract class interface that evaluates functions and constraints
    at a single point.  An Evaluator instance evaluates at just one point,
    while an Executor instance coordinates the parallel execution of evaluations
    at many points.

    An interface design makes it easy for users to extend HOPSPACK with their
    own Evaluator implementation.  New implementations will typically subclass
    this interface and modify source code in EvaluatorFactory.
 */
//----------------------------------------------------------------------
class Evaluator
{
  public:

    //! Evaluate the objective function(s) at a point x.
    /*!
     *  @param[in] nTag   Contains a unique tag for the evaluation which can be
     *                    used to name files, etc.
     *  @param[in] cX     The point at which to evaluate the function(s).
     *  @param[out] cFns  On output, contains a vector of objective function
     *                    values computed at X.  Multiple objectives are allowed.
     *                    If an evaluation failed, return an empty vector or set
     *                    individual elements of the vector to HOPSPACK::dne().
     *  @param[out] sMsg  On output, contains a message about the evaluation;
     *                    typically the word "Success" or an error message.
     */
    virtual void  evalF (const int       nTag,
                         const Vector &  cX,
                               Vector &  cFns,
                               string &  sMsg) = 0;

    //! Evaluate the objective functions and nonlinear constraints at a point x.
    /*!
     *  @param[in] nTag     Contains a unique tag for the evaluation which can be
     *                      used to name files, etc.
     *  @param[in] cX       The point at which to evaluate the function(s).
     *  @param[out] cFns    On output, contains a vector of objective function
     *                      values computed at X.  Multiple objectives are
     *                      allowed.  If an evaluation failed, return an empty
     *                      vector or set individual function elements of the
     *                      vector to HOPSPACK::dne().
     *  @param[out] cEqs    On output, contains a vector of nonlinear equality
     *                      constraint function values computed at X.  If an
     *                      evaluation failed, return an empty vector or set
     *                      individual elements of the vector to HOPSPACK::dne().
     *  @param[out] cIneqs  On output, contains a vector of nonlinear inequality
     *                      constraint function values computed at X.  If an
     *                      evaluation failed, return an empty vector or set
     *                      individual elements of the vector to HOPSPACK::dne().
     *  @param[out] sMsg    On output, contains a message about the evaluation;
     *                      typically the word "Success" or an error message.
     */
    virtual void  evalFC (const int       nTag,
                          const Vector &  cX,
                                Vector &  cFns,
                                Vector &  cEqs,
                                Vector &  cIneqs,
                                string &  sMsg) = 0;

    //! Evaluate the nonlinear constraints at a point x.
    /*!
     *  @param[in] nTag     Contains a unique tag for the evaluation which can be
     *                      used to name files, etc.
     *  @param[in] cX       The point at which to evaluate the function(s).
     *  @param[out] cEqs    On output, contains a vector of nonlinear equality
     *                      constraint function values computed at X.  If an
     *                      evaluation failed, return an empty vector or set
     *                      individual elements of the vector to HOPSPACK::dne().
     *  @param[out] cIneqs  On output, contains a vector of nonlinear inequality
     *                      constraint function values computed at X.  If an
     *                      evaluation failed, return an empty vector or set
     *                      individual elements of the vector to HOPSPACK::dne().
     *  @param[out] sMsg    On output, contains a message about the evaluation;
     *                      typically the word "Success" or an error message.
     */
    /*TBD...not used by GSS solver, but other solvers might want this
    virtual void  evalC (const int       nTag,
                         const Vector &  cX,
                               Vector &  cEqs,
                               Vector &  cIneqs,
                               string &  sMsg) = 0;
    */

    //! Return the string from parameter 'Evaluator Type'.
    virtual string  getEvaluatorType (void) const = 0;

    //! Print debug information about the Evaluator instance.
    virtual void  printDebugInfo (void) const = 0;

};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_EVALUATOR_HPP
