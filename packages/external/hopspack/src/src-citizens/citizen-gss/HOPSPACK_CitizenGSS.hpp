// $Id: HOPSPACK_CitizenGSS.hpp 183 2010-12-15 18:22:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss/HOPSPACK_CitizenGSS.hpp $ 

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
  @file HOPSPACK_CitizenGSS.hpp
  @brief Class declaration for HOPSPACK::CitizenGSS.
*/

#ifndef HOPSPACK_CITIZENGSS_HPP
#define HOPSPACK_CITIZENGSS_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Citizen.hpp"
#include "HOPSPACK_DataPoint.hpp"
#include "HOPSPACK_GssIterator.hpp"
#include "HOPSPACK_GssList.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_NonlConstrPenalty.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Implements Citizen using Generating Set Search (GSS).
/*! This Citizen solves optimization problems with continuous variables,
 *  variable bounds, and linear constraints.  It does not require or use
 *  derivatives.  The GSS algorithm is guaranteed to find a local stationary
 *  point under mild problem assumptions.
 *
 *  Starting from an initial feasible point, GSS produces a set of feasible
 *  trial points within a "step size" of the current best point.  Evaluated
 *  points can be accepted in any order, and the algorithm is well-suited
 *  for asynchronous parallel architectures.  GSS can accept "better"
 *  points found by other algorithms and will adjust its search accordingly.
 */
//----------------------------------------------------------------------
class CitizenGSS : public Citizen
{
  public:

    //! Constructor.
    /*!
     *  @param[in] nIdNumber   Unique HOPSPACK ID number of the citizen.
     *  @param[in] sName       Name of the citizen.
     *  @param[in] cParams     Sublist of user input parameters for this citizen.
     *  @param[in] cProbDef    Problem definition.
     *  @param[in] cLinConstr  Linear constraints definition.
     *  @param[in] pParent     Pointer to parent citizen for calling back,
     *                         or NULL if there is no parent.
     *  @throws "GSS Error" if parameters are not accepted.
     */
    CitizenGSS (const int                    nIdNumber,
                const string        &        sName,
                const ParameterList &        cParams,
                const ProblemDef    &        cProbDef,
                const LinConstr     &        cLinConstr,
                      Citizen       * const  pParent);

    //! Destructor.
    ~CitizenGSS (void);


    //! Called once by the Mediator prior to processing points by any citizen.
    void  preProcess (void);

    //! Called once by the Mediator after processing points by all citizens.
    /*!
     *  Typically, a citizen will display its state and the best point found.
     */
    void  postProcess (void);

    //! Exchange a list of evaluated points for a list of new trial points.
    /*!
     *  The citizen sees all new result points since the last call to this
     *  method.  The "exchange" happens implicitly:  after all citizens have
     *  viewed the newly evaluated points, then the Mediator erases them.
     *  The cOwnerTags parameter allows the citizen to find the subset
     *  of results matching previous requests that it made.
     *
     *  GSS performs the following tasks:
     *  <ul>
     *  <li> The cReturnList is first copied to an internal exchangeList.
     *       If a point is from another citizen, then the copy appears as though
     *       it has no GSS parent.
     *  <li> Infeasible points are removed from the internal list
     *       (for efficiency, only infeasible points with a better objective
     *       value are removed).
     *  <li> The internal list is then sent to GssIterator, which replaces
     *       the newly evaluated points with new trial points.  These are
     *       copied to the cWaitList.
     *  <li> The internal exchangeList is pruned (all points deleted).
     *  </ul>
     *
     *  @param[in] cReturnList  Input list of newly evaluated points.
     *  @param[in] cOwnerTags   List containing the tags of points in cReturnList
     *                          that this citizen owns.
     *  @param[out] cWaitList   Output list of trial points that will be added
     *                          to the wait queue.  Points will be evaluated
     *                          in this order, starting from the back.
     */
    void  exchange (const list< DataPoint * > &  cReturnList,
                    const list< int >         &  cOwnerTags,
                          list< DataPoint * > &  cWaitList);

    //! Mark the child to exit as soon as possible.
    /*!
     *  Called by the Mediator when it plans to ignore trial points.
     *  The Mediator may make one more call to the exchange() method,
     *  but the citizen should not waste effort adding trial points,
     *  because they will not be evaluated.
     */
    void  setEarlyExit (void);

    //! Return the unique ID number of the citizen.
    int  getIdNumber (void) const;

    //! Return the unique name of the citizen.
    const string &  getName (void) const;

    //! Return the current state of the citizen.
    State  getState (void) const;


  private:

    //! By design, there is no copy constructor.
    CitizenGSS (const CitizenGSS &);
    //! By design, there is no assignment operator.
    CitizenGSS & operator= (const CitizenGSS &);

    //! Print debug information before generating new points.
    void  printPreDiagnostics_ (void) const;

    //! Print debug information after generating new points.
    void  printPostDiagnostics_ (bool  bFoundNewBest) const;

    //! Remove "best" infeasible points until finding a feasible one.
    void  popBestInfeasiblePoints_ (void);


    //! ID of the citizen.
    int  _nID;

    //! Name of the citizen.
    string  _sName;

    //! Reference to the problem definition.
    const ProblemDef &  _cProbDef;

    //! Reference to the linear constraints.
    const LinConstr &  _cLinConstr;

    //! Input parameters for the citizen, which may grow with appended defaults.
    ParameterList  _cGssParams;

    //! Penalty function for nonlinear constraints.
    NonlConstrPenalty *  _pPenalty;

    //! Pointer to parent for making a callback when finished, or NULL if none.
    Citizen * const  _pParent;

    //! Internal exchange list.
    GssList  _cExchangeList;

    //! Object that does the work of a GSS iteration.
    GssIterator *  _pGssIterator;

    //! True if Mediator intends to exit the citizen.
    bool  _bWillExitEarly;

    //! Maximum old point requests to keep after a new best point is found.
    int  _nMaxToKeepAfterNewBest;

    //! Display level, set by user parameter.
    /*!
     *  0 - display nothing from GSS citizen
     *  1 - display new best points, initial params, and a final summary
     *  2 - additionally, display exchanged points
     *  3 - additionally, display generated directions
     */
    int  _nDisplayLevel;
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_CITIZENGSS_HPP
