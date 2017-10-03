// $Id: HOPSPACK_CitizenGssNlc.hpp 183 2010-12-15 18:22:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss-nlc/HOPSPACK_CitizenGssNlc.hpp $

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
  @file HOPSPACK_CitizenGssNlc.hpp
  @brief Class declaration for HOPSPACK::CitizenGssNlc.
*/

#ifndef HOPSPACK_CITIZENGSSNLC_HPP
#define HOPSPACK_CITIZENGSSNLC_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_CallbackToMediator.hpp"
#include "HOPSPACK_Citizen.hpp"
#include "HOPSPACK_DataPoint.hpp"
#include "HOPSPACK_gssChildReturnCodes.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_NonlConstrPenalty.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Implements Citizen using GSS for nonlinear constraints.
/*! This Citizen solves optimization problems with continuous variables
 *  and general constraints, which includes variable bounds, linear constraints,
 *  and nonlinear constraints.  It does not require or use derivatives.
 *
 *  Starting from an initial linearly feasible point, GSS-NLC treats nonlinear
 *  constraints with penalty terms in the objective.  A series of subproblems,
 *  each with only linear constraints, are solved by GSS child citizens.
 *  The penalty parameter is modified for successive subproblems to drive
 *  convergence towards an optimal solution.
 *
 *  For more information, read:
 *  "Nonlinearly-Constrained Optimization Using Asynchronous Parallel
 *   Generating Set Search",
 *  JD Griffin and TG Kolda, Sandia National Laboratories, 2007,
 *  Sandia Report SAND2007-3257.
 */
//----------------------------------------------------------------------
class CitizenGssNlc : public Citizen
{
  public:

    //! Constructor.
    /*!
     *  @param[in] nIdNumber   Unique HOPSPACK ID number of the citizen.
     *  @param[in] sName       Name of the citizen.
     *  @param[in] cParams     Sublist of user input parameters for this citizen.
     *  @param[in] cProbDef    Problem definition.
     *  @param[in] cLinConstr  Linear constraints definition.
     *  @param[in] pCallback   Interface for dynamically adding child citizens.
     *  @param[in] pParent     Pointer to parent for calling back, or NULL
     *                         if this is not a child citizen.
     *  @throws "GSSNLC Error" if parameters are not accepted.
     */
    CitizenGssNlc (const int                         nIdNumber,
                   const string             &        sName,
                   const ParameterList      &        cParams,
                   const ProblemDef         &        cProbDef,
                   const LinConstr          &        cLinConstr,
                         CallbackToMediator * const  pCallback,
                         Citizen            * const  pParent);

    //! Destructor.
    ~CitizenGssNlc (void);


    //! Called once by the Mediator prior to processing points by any citizen.
    /*!
     *  GSS-NLC initializes the penalty function and launches the first GSS
     *  subproblem solver.
     */
    void  preProcess (void);

    //! Called once by the Mediator after processing points by all citizens.
    /*!
     *  Typically, a citizen will display its state and the best point found.
     */
    void  postProcess (void);

    //! Exchange a list of evaluated points for a list of new trial points.
    /*!
     *  The GSS-NLC citizen does not use this method.  All points are exchanged
     *  by subproblems.
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

    //! Child citizen calls back to its parent citizen using this method.
    /*!
     *  GSS subproblem solvers call this when finished.
     *  GSS-NLC examines the result and either halts or starts a new
     *  subproblem solver.
     *
     *  @param[in] nIdNumber    Unique HOPSPACK ID number of the child instance
     *                          that is calling back.
     *  @param[in] nReturnCode  The parent expects one of the GSS return codes
     *                          from GssChildReturnCodesType.
     *  @param[in] cFinalPoint  Solution point returned by the child.
     *  @param[in] nTotalEvals  Evaluations made by the child citizen.
     *                          The number does not include points that were
     *                          "evaluated" from the cache.
     */
    void  callbackFromChild (const int          nIdNumber,
                             const int          nReturnCode,
                             const DataPoint &  cFinalPoint,
                             const int          nTotalEvals);

  private:

    //! By design, there is no copy constructor.
    CitizenGssNlc (const CitizenGssNlc &);
    //! By design, there is no assignment operator.
    CitizenGssNlc & operator= (const CitizenGssNlc &);

    //! Parse parent parameters, returning the remainder.
    /*!
     *  @param[in,out] cParams  All GSS-NLC configuration parameters.
     *                          Additional parameters with default values are
     *                          added for easy debugging.
     *  @param[out] cRemainder  New copy of remaining parameters after GSS-NLC
     *                          has extracted those it knows about.
     *  @return                 False if there was a problem with an
     *                          extracted parameter.
     *
     *  As a side effect, most extracted parameters are saved in private members.
     */
    bool  extractParameters_ (ParameterList &  cParams,
                              ParameterList &  cRemainder);

    //! Create a new subproblem citizen.
    /*!
     *  @param[in] cChildParams  Full set of parameters for the GSS child
     *                           except penalty parameters, which are added
     *                           by the method.
     *  @param[in] cProbDef      Problem definition for the GSS child.
     *  @param[in] cPenalty      Penalty function for the GSS child.
     *  @return                  Child citizen ID number.
     */
    int  createNewChildCitizen_ (      ParameterList     &  cChildParams,
                                 const ProblemDef        &  cProbDef,
                                 const NonlConstrPenalty &  cPenalty);

    //! Return true if the nonlinear solve is finished.
    /*!
     *  @param[in] nReturnCode  Return code from the child, excluding the
     *                          case of REASON_ERROR.
     *  @param[in] cTestPoint   Solution point returned by the child.
     */
    bool  isTimeToStop_ (const GssChildReturnCodesType  nReturnCode,
                         const DataPoint &              cTestPoint);

    //! Update the penalty function for the next subproblem.
    /*!
     *  @param[in] cSubprobSolution  Solution point from the latest subproblem.
     */
    void  updatePenalty_ (const DataPoint &  cSubprobSolution);

    //! Print a point with its penalty term.
    void  printPointWithPen_ (const DataPoint &  cPoint) const;


    //! ID of the citizen.
    int  _nID;

    //! ID of the current child citizen.
    int  _nChildID;

    //! Name of the citizen.
    string  _sName;

    //! Current state of the citizen.
    State  _nState;

    //! Reference to the problem definition.
    const ProblemDef &  _cProbDef;

    //! Total number of nonlinear equality and inequality constraints.
    int  _nM;

    //! Reference to the linear constraints.
    const LinConstr &  _cLinConstr;

    //! Input parameters augmented with default values.
    ParameterList  _cParentParams;

    //! Unused input parameters sent to GSS child citizens.
    ParameterList  _cSubprobParams;

    //! Pointer to callbacks implemented by the Mediator.
    CallbackToMediator * const   _pCallbackToMed;

    //! Pointer to parent for making a callback when finished, or NULL if none.
    Citizen * const  _pParent;

    //! Parameters used by the currently running child.
    ParameterList *  _pChildParams;

    //! Problem definition used by the currently running child.
    ProblemDef *  _pChildProbDef;

    //! Most recent subproblem solution.
    DataPoint *  _pLatestSubprobSol;

    //! True if Mediator intends to exit the citizen.
    bool  _bWillExitEarly;

    //! Total point evaluations by this citizen and its children.
    /*!
     *  The number does not include points that were "evaluated" from the cache.
     */
    int  _nTotalEvals;

    //! Maximum number of evaluations for the citizen.
    int  _nMaxGssNlcEvals;

    //! Maximum number of evaluations per subproblem.
    int  _nMaxSubprobEvals;

    //! Display level for parent citizen, set by user parameter.
    /*!
     *  0 - display nothing from GSS-NLC citizen
     *  1 - display initial params and a final summary
     *  2 - additionally, display subproblem solutions
     */
    int  _nDisplayLevel;

    //! Display level for child citizens, set by user parameter.
    int  _nDisplaySubLevel;

    //! Reason for finishing (see enumeration in definition code).
    int  _nFinishedReason;

    //! True if GSS-NLC and its children ignore points from other citizens.
    bool  _bIgnoreOtherPoints;

    //! Penalty function for nonlinear constraints.
    NonlConstrPenalty  _cPenalty;

    //! Initial value for the penalty parameter.
    double  _dInitialPenalty;

    //! Maximum value for the penalty parameter.
    double  _dMaxPenalty;

    //! Factor to increase the penalty parameter by.
    double  _dPenaltyIncrease;

    //! Initial step length for first subproblem.
    double  _dInitialStepLength;

    //! Current stopping tolerance for step length.
    double  _dCurrentStepTol;

    //! Stopping tolerance for step length.
    double  _dFinalStepTol;

    //! Factor to decrease the subproblem step tolerance by.
    double  _dStepTolDecrease;

    //! Initial smoothing factor for nonsmooth penalty functions.
    double  _dInitialSmoothing;

    //! Minimum smoothing factor for nonsmooth penalty functions.
    double  _dMinSmoothing;

    //! Factor to decrease smoothing factor by.
    double  _dSmoothingDecrease;
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_CITIZENGSSNLC_HPP
