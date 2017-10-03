// $Id: TBD $
// $URL: svn+ssh://software.sandia.gov/svn/private/hopspack/trunk/src/src-evaluator/HOPSPACK_EvaluatorFactory.hpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2014 Sandia Corporation
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
  @file HOPSPACK_EvaluatorFactory.hpp
  @brief Declaration for HOPSPACK::EvaluatorFactory.
*/

#ifndef HOPSPACK_EVALUATORFACTORY_HPP
#define HOPSPACK_EVALUATORFACTORY_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Evaluator.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Factory class to construct instances of an evaluator.
/*!
 *  The factory allows an evaluator to be determined at run time based on
 *  the parameter "Evaluator Type".
 */
//----------------------------------------------------------------------
class EvaluatorFactory
{
  public:

    //! Parameterized factory method that returns a new thread-safe instance.
    /*!
     *  Any returned instance can be used with the serial and MPI versions
     *  of HOPSPACK.  If the returned instance evaluates independently of
     *  other instances and is thread-safe, then it can be used with the
     *  multithreaded version of HOPSPACK.
     *
     *  @param[in] cEvalParams  Parameters in the "Evaluator" sublist.
     *                          Parameter value "Evaluator Type" determines
     *                          the particular implementation.
     *  @return                 The Evaluator instance if successful, else NULL.
     *                          Client must delete the instance when finished.
     */
    static Evaluator *  newInstance (const ParameterList &  cEvalParams);

    //! Destructor.
    ~EvaluatorFactory (void);

  protected:

    //! Default constructor is hidden, as only static methods are available.
    EvaluatorFactory (void);

};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_EVALUATORFACTORY_HPP
