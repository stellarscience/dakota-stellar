// $Id: HOPSPACK_MultiStartRepository.hpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss-ms/HOPSPACK_MultiStartRepository.hpp $

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
  @file HOPSPACK_MultiStartRepository.hpp
  @brief Class declaration for HOPSPACK::MultiStartRepository
*/

#ifndef HOPSPACK_MULTISTARTREPOSITORY_HPP
#define HOPSPACK_MULTISTARTREPOSITORY_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_DataPoint.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Stores results from multi-start searches.
/*! An instance stores points resulting from a search over multiple start
 *  points.  The results are made available in ranked order.
 */
//----------------------------------------------------------------------
class MultiStartRepository
{
  public:

    //! Constructor.
    /*!
     *  @param[in] cProbDef        Problem definition.
     *  @param[in] cLinConstr      Linear constraints definition.
     *  @param[in] dComparisonTol  Tolerance for declaring two points "equal".
     */
    MultiStartRepository (const ProblemDef &  cProbDef,
                          const LinConstr  &  cLinConstr,
                          const double        dComparisonTol);

    //! Destructor.
    ~MultiStartRepository (void);


    //! Add a new result point to the repository.
    /*!
     *  The point is added if distinct, or grouped if not distinct.
     *  Two points \f$ x \f$ and \f$ y \f$ are not distinct if their
     *  scaled component-wise distances are within the comparison tolerance
     *  \f$ \tau \f$; i.e.,
     *  \f[
     *    |x_i - y_i| \leq \tau s_i \; \mbox{  for all } i = 1, \dots, n
     *  \f]
     *  where \f$ s \f$ is the scaling vector for the problem.
     *  If two points are indistinct, then only the one with the best
     *  objective is kept.
     */
    void  addResult (const DataPoint &  cNewResult);

    //! Get the best feasible result, or return false if none.
    /*!
     *  @param[out] cBestResult  Filled in with the best feasible point.
     *                           The result will be feasible with respect to
     *                           linear and nonlinear constraints.
     *  @return                  false if no feasible result is known
     */
    bool  getBestResult (DataPoint &  cBestResult) const;

    //! Get the n best feasible results, ranked from best to worst.
    /*!
     *  @param[in] nMaxNumResults  Maximum number of results to return.
     *  @param[out] cResults       Filled in with the best results, in order
     *                             beginning with the best.  All results are
     *                             feasible with respect to linear and nonlinear
     *                             constraints.  The number of points may be
     *                             less than nMaxNumResults.
     *
     *  Caller is responsible for freeing the DataPoint instances returned
     *  in cResults.
     */
    void  getBestResultList (const int                    nMaxNumResults,
                                   vector< DataPoint * >  cResults) const;


  private:

    //! By design, there is no copy constructor.
    MultiStartRepository (const MultiStartRepository &);
    //! By design, there is no assignment operator.
    MultiStartRepository & operator= (const MultiStartRepository &);


    //! Reference to the problem definition.
    const ProblemDef &  _cProbDef;

    //! Reference to the linear constraints.
    const LinConstr &  _cLinConstr;

    //! Tolerance for determining if two points are distinct.
    double  _dComparisonTol;
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_MULTISTARTREPOSITORY_HPP
