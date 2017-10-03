// $Id: HOPSPACK_CitizenGssMS.hpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss-ms/HOPSPACK_CitizenGssMS.hpp $

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
  @file HOPSPACK_CitizenGssMS.hpp
  @brief Class declaration for HOPSPACK::CitizenGssMS.
*/

#ifndef HOPSPACK_CITIZENGSSMS_HPP
#define HOPSPACK_CITIZENGSSMS_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_CallbackToMediator.hpp"
#include "HOPSPACK_Citizen.hpp"
#include "HOPSPACK_gssChildReturnCodes.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_PointGeneratorInterface.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Implements Citizen that drives GSS subproblems with multiple start points.
/*! This Citizen solves optimization problems with general variables
 *  (continuous and integer-valued) and general constraints (variable bounds,
 *  linear constraints, and nonlinear constraints).
 *  It does not require or use derivatives.
 *
 *  The citizen generates multiple starting points and invokes a GSS
 *  subproblem to find a local solution from each one.
 */
//----------------------------------------------------------------------
class CitizenGssMS : public Citizen
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
     *  @throws "GSSMS Error" if parameters are not accepted.
     */
    CitizenGssMS (const int                         nIdNumber,
                  const string             &        sName,
                  const ParameterList      &        cParams,
                  const ProblemDef         &        cProbDef,
                  const LinConstr          &        cLinConstr,
                        CallbackToMediator * const  pCallback);

    //! Destructor.
    ~CitizenGssMS (void);


    //! Called once by the Mediator prior to processing points by any citizen.
    /*!
     *  TBD.
     */
    void  preProcess (void);

    //! Called once by the Mediator after processing points by all citizens.
    /*!
     *  Typically, a citizen will display its state and the best point found.
     */
    void  postProcess (void);

    //! Exchange a list of evaluated points for a list of new trial points.
    /*!
     *  The GSS-MS citizen does not use this method.  All points are exchanged
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

    //! Return the unique ID number of the citizen.
    int  getIdNumber (void) const;

    //! Return the unique name of the citizen.
    const string &  getName (void) const;

    //! Return the current state of the citizen.
    State  getState (void) const;

    //! Child citizen calls back to its parent citizen using this method.
    /*!
     *  GSS subproblem solvers call this when finished.
     *  GSS-MS examines the result and either halts or starts a new
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
    CitizenGssMS (const CitizenGssMS &);
    //! By design, there is no assignment operator.
    CitizenGssMS & operator= (const CitizenGssMS &);

    //! Parse parent parameters, returning the remainder.
    /*!
     *  @param[in,out] cParams  All GSS-MS configuration parameters.
     *                          Additional parameters with default values are
     *                          added for easy debugging.
     *  @param[out] cRemainder  New copy of remaining parameters after GSS-MS
     *                          has extracted those it knows about.
     *  @return                 False if there was a problem with an
     *                          extracted parameter.
     *
     *  As a side effect, most extracted parameters are saved in private members.
     */
    bool  extractParameters_ (ParameterList &  cParams,
                              ParameterList &  cRemainder);

    //! Perform the next iteration of multi-start.
    /*!
     *  @return  False if there was an error creating a child citizen.
     *
     *  An iteration gets new start points from the point generator until one
     *  of the following is true:
     *  <ul>
     *  <li> Generator returns a list of points to be evaluated.
     *  <li> Full number of concurrent child citizens are started.
     *  <li> Generator has no more start points.
     *  <li> Mediator signals that it is time to stop.
     *  <li> A child citizen fails to start.
     *  </ul>
     */
    bool  nextIteration_ (void);

    //! Return true if the candidate start point is linearly feasible.
    bool  isStartPointOK_ (const Vector &  cStartPoint) const;

    //! Return true if the multi-start solve is finished.
    /*!
     *  @param[in] nReturnCode  Return code from the child, excluding the
     *                          case of REASON_ERROR.
     *  @param[in] cTestPoint   Solution point returned by the child.
     */
    bool  isTimeToStop_ (const GssChildReturnCodesType  nReturnCode,
                         const DataPoint &              cTestPoint);


    //! ID of the citizen.
    int  _nID;

    //! Name of the citizen.
    string  _sName;

    //! Current state of the citizen.
    State  _nState;

    //! Reference to the problem definition.
    const ProblemDef &  _cProbDef;

    //! Reference to the linear constraints.
    const LinConstr &  _cLinConstr;

    //! Input parameters augmented with default values.
    ParameterList  _cParentParams;

    //! Unused input parameters sent to child citizens.
    ParameterList  _cSubprobParams;

    //! Pointer to callbacks implemented by the Mediator.
    CallbackToMediator * const   _pCallbackToMed;

    //! Display level for parent citizen, set by user parameter.
    /*!
     *  0 - display nothing from GSS-MS citizen
     *  1 - display initial params and a final summary
     *  2 - additionally, display subproblem solutions
     */
    int  _nDisplayLevel;

    //! Display level for child citizens, set by user parameter.
    int  _nDisplaySubLevel;

    //! Reason for finishing (see enumeration in definition code).
    int  _nFinishedReason;

    //! Maximum number of evaluations per subproblem.
    int  _nMaxSubprobEvals;

    //! Total number of start points to generate.
    int  _nTotalStartPoints;

    //! Current tally of start points that have been submitted.
    int  _nCurrentNumStartPoints;

    //! Number of subproblems to run at the same time.
    int  _nNumConcurrentSubprobs;

    //! Point generator.
    PointGenerator *  _pGenerator;

    //! Total number of evaluations by this citizen and its children.
    int  _nTotalEvals;


    //! Structure to hold child citizen information.
    typedef struct
    {
        ProblemDef *  pProbDef;    //-- PROBLEM DEFINITION
        int           nChildID;    //-- CHILD CITIZEN ID
    }
    ChildCtznInfoBlockType;

    //! Container for currently active child problem definitions.
    vector< ChildCtznInfoBlockType * >  _cChildren;
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_CITIZENGSSMS_HPP
