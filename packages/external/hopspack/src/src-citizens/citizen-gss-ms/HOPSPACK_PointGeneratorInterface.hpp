// $Id: HOPSPACK_PointGeneratorInterface.hpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss-ms/HOPSPACK_PointGeneratorInterface.hpp $

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
  @file HOPSPACK_PointGeneratorInterface.hpp
  @brief Interface declaration for HOPSPACK::PointGenerator.
*/

#ifndef HOPSPACK_POINTGENERATORINTERFACE_HPP
#define HOPSPACK_POINTGENERATORINTERFACE_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_DataPoint.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Interface for algorithms that generate a set of start points.
/*! Generator algorithms effectively control the multi-start process by
 *  deciding on a sequence of start points.
 *
 *  When the GSS-MS citizen calls a point generator method, all main thread
 *  processing by HOPSPACK waits for the method to return.  Therefore, a
 *  point generator requiring excessive CPU time should offload the work
 *  using an asynchronous citizen worker.
 *
 *  An interface design makes it easy for users to extend HOPSPACK with their
 *  own multi-start generator.  New implementations must subclass this class
 *  and become registered by modifying the "newInstance" method.
 */
//----------------------------------------------------------------------
class PointGenerator
{
  public:

    //! Possible state of a result after optimizing from a start point.
    enum ResultState
    {
        //! Result is feasible and declared locally optimal.
        OPTIMAL = 1,

        //! Result is infeasible with respect to nonlinear constraints.
        INFEASIBLE,

        //! Result is feasible but not declared locally optimal.
        FEASIBLE_NOT_OPTIMAL
    };


    //! Factory method for constructing a point generator instance.
    /*!
     *  @param[in] sName           Name of the generator type to construct.
     *  @param[in] nNumToGenerate  Total number of start points to generate.
     *                             Implementations are free to ignore this,
     *                             but some algorithms need to know in advance
     *                             how many times they will be called upon.
     *  @param[in] cProbDef        Problem definition.
     *  @param[in] cLinConstr      Linear constraints definition.
     *  @return                    New instance if successful, else NULL.
     *                             Caller must delete the instance when finished.
     */
    static PointGenerator *  newInstance (const string     &  sName,
                                          const int           nNumToGenerate,
                                          const ProblemDef &  cProbDef,
                                          const LinConstr  &  cLinConstr);

    //! Destructor.
    virtual ~PointGenerator (void) = 0;


    //! Generate the next start point, which will be linearly feasible.
    /*!
     *  @param[out] cStartLocation  Filled with the position of the next start
     *                              point.  If cEvalList is non-empty, then this
     *                              argument should be ignored.
     *  @param[out] cEvalList       Filled with a list of points to evaluate,
     *                              which should then be returned by calling
     *                              addResultPair().  Caller is responsible
     *                              for deleting the points.
     *                              If empty, then cStartLocation is defined.
     *  @return                     False if there are no more start points.
     *
     *  Exactly one of the two output arguments are filled.  The caller should
     *  then either solve a subproblem from the new start location, or
     *  evaluate the list of data points.  In either case, call addResultPair()
     *  with at least one result before calling getNextPoint() again.
     *
     *  All returned points are guaranteed to be feasible with respect to
     *  variable bounds and linear constraints.
     */
    virtual bool  getNextPoint (Vector &                       cStartLocation,
                                vector< const DataPoint * > &  cEvalList) = 0;

    //! Add the result for a previously generated start point.
    /*!
     *  @param[in] cStartPoint   Start point, whose getX() method must return
     *                           a position previously generated by
     *                           getNextPoint().
     *  @param[in] cResultPoint  Result after optimization from the start point.
     *  @param[in] nResultState  State of the result point.
     *
     *  Implementations are free to ignore result information, but some
     *  algorithms generate new start points based on the results from
     *  previous start points.
     */
    virtual void  addResultPair (const DataPoint &  cStartPoint,
                                 const DataPoint &  cResultPoint,
                                 const ResultState  nResultState) = 0;

    //! Print information about the instance to standard output.
    virtual void  printDefinition (void) const = 0;


  protected:

    //! Subclasses should call this base class constructor.
    PointGenerator (void);
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_POINTGENERATORINTERFACE_HPP
