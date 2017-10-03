// $Id: HOPSPACK_GssIterator.hpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss/HOPSPACK_GssIterator.hpp $

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
  \file HOPSPACK_GssIterator.hpp
  \brief Class declaration for HOPSPACK::GssIterator.
*/

#ifndef HOPSPACK_GSSITERATOR_HPP
#define HOPSPACK_GSSITERATOR_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_GssDirections.hpp"
#include "HOPSPACK_GssList.hpp"
#include "HOPSPACK_GssPoint.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_NonlConstrPenalty.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//! Performs a single GSS iteration using reverse communication.
/*!
  The constructor reads the following parameters from the user inputs.
  <ul>
  <li>"Initial Step"
  <li>"Snap To Boundary"
  <li>"Snap Distance"
  <li>"Sufficient Improvement Factor"
  <li>"Use Random Order"
  </ul>

  The class GssIterator is capable of performing a single GSS iteration.
  Communication with an external driving class is performed via the method
  pointExchange().  Essentially, pointExchange() inputs a list of evaluated
  points and outputs a list of points it wishes evaluated.  This process
  continues until the GssIterator "converges" or is stopped externally.
*/
class GssIterator
{
  public:

    //! Constructor.
    /*!
     *  @param[in] cProbDef     Problem definition.
     *  @param[in] cLinConstr   Linear constraints definition.
     *  @param[in] pPenalty     Penalty function for nonlinear constraints.
     *  @param[in,out] cParams  Sublist of user parameters for this citizen.
     */
    GssIterator (const ProblemDef        &        cProbDef,
                 const LinConstr         &        cLinConstr,
                 const NonlConstrPenalty * const  pPenalty,
                       ParameterList     &        cParams);

    //! Destructor
    ~GssIterator (void);

    // ----------------------------------------------------------------------

    //@{ \name Accessors

    //! Return true if GssIterator has finished running.
    bool  isFinished (void) const;

    //! Return a read-only reference to the best point.
    const GssPoint &  getBestPoint (void) const;

    //! Return number of points submitted by this citizen that were evaluated.
    /*!
     *  The number does not include points that were "evaluated" from the cache.
     */
    int  getNumGssEvals (void) const;

    //! Return user parameter for initial step.
    double  getInitialStep (void) const;

    //! Return true if the iterator finished because it converged.
    /*!
     *  Converged might mean the step length tolerance was met, or the objective
     *  target was met.
     */
    bool  hasStoppedAndConverged (void) const;

    //! Return true if the iterator finished because it ran out of evaluations.
    bool  hasStoppedOutOfEvals (void) const;

    //@}

    // ----------------------------------------------------------------------

    //@{ \name Manipulators

    //! Performs a single GSS iteration, exchanging evaluated for trial points.
    /*!
     *  @param[in,out] newList    Input a list of evaluated points, and
     *                            output a list of points to be evaluated.
     *  @param[in] bShouldIgnore  If true, then points in newList that are
     *                            from other citizens should be ignored.
     *  @param[in] bPrintDetails  If true, then details of point construction
     *                            are displayed.
     *  @return true              If a new best point was found.
    */
    bool pointExchange (      GssList &  newList,
                        const bool       bShouldIgnore,
                        const bool       bPrintDetails = false);

    //@}

    //@{ \name Print Functions

    //! Prints solver parameter list.
    void printInitializationInformation() const;

    //! Prints the reason for stopping.
    void printStopReason() const;

    //! Prints current set of search directions.
    void printDirections(const string label = "") const;

    //@}


  private:

    //! By design, there is no copy constructor.
    GssIterator (const GssIterator &);
    //! By design, there is no assignment operator.
    GssIterator & operator= (const GssIterator &);

    //! Create the initial best point from problem data.
    GssPoint *  initializeBestPointPtr (const ProblemDef &  cProbDef,
                                        const LinConstr  &  cLinConstr) const;

    //! Process a new best point, replacing any old best point.
    /*!
      If the function tolerance convergence test is being employed,
      check for convergence and return if convergence is detected.
      Otherwise, generate new search directions, re-initialize the step array
      and reset various bookkeeping variables.

      All step lengths start the same, and they are calculated as follows:
      \f[
        {\rm step}_{\rm new} = \max \{ {\rm step}_{\rm best},
                                       2 * {\rm stepTolerance} \}
      \f]
    */
    void  processNewBestPoint (GssPoint *  newBestPointPtr = NULL);

    /*! Trial points are generated corresponding to vectors in _cDirections 
    that are not converged and do not already have an associated trial point.

    For each direction \f$d_i\f$ in _cDirections, a new trial point \f$ x\f$ is
    computed satisfying
    \f[
    x = x_{\rm best} + \Delta_{\rm max} d_i,
    \f] 
    where \f$\Delta_{\rm max}\f$ denotes
    the maximum scalar in interval \f$[0,\Delta_i]\f$ such that
    \f$ x_{\rm best} + \Delta_{\rm max} d_i,\f$ is feasible.  
    - \f$ \Delta_i\f$ denotes the corresponding step value to \f$ d_i \f$
    whose value is determined via a call to GssDirections::getStep().
    - \f$ \Delta_{\rm max} \f$ is determined by a call to
    Constraints::Linear::maxStep().

    If the user enables useSnapTo, then a call is made to determine if \f$ x \f$
    is near a boundary point.  If a nearby boundary point \f$ x_b \f$ is found,
    then the trial point is modified:
    \f[
    x = x_b.
    \f]
    
    \note New search directions are created in processNewBestPoint().
    */
    void generateTrialPoints (const bool  bPrintDetails);

    //! Process a list of trial points that has been returned from the Conveyor.
    /*!
    First, check to see if the best of all these points is better than
    the current best point. If so, replace that best point (see
    processNewBestPoint()) and throw away the rest of the list.
    Otherwise, process each point and delete it.
    Finally, check if all the directions have converged and check if
    we have exhausted the maximum number of function evaluations.
    
    Returns true of a new best point has been found indicating the
    queue may be pruned.
    */
    bool processEvaluatedTrialPoints();
  

    //! State of the iterator.
    enum IteratorState {
      NOT_CONVERGED = 0,
      STEPLENGTH_CONVERGED,
      OBJECTIVE_REACHED,
      MAX_EVALS_FOR_CTZN,
      STOP_FROM_ERROR
    };

    //! Input parameters.
    ParameterList &  _cParams;

    //! Problem definition.
    const ProblemDef &  _cProbDef;

    //! Penalty function for nonlinear constraints.
    const NonlConstrPenalty &  _cPenalty;

    //! Pointer to the best trial point thus far.
    GssPoint *  _pBestPointPtr;

    //! True if the iterator has completed.
    bool  _bIsFinished;

    //! The search directions.
    GssDirections  _cDirections;

    //! List of trial points generated during an iteration.
    GssList  _cExchangeList;

    //! List of submitted points owned by this citizen instance.
    vector< int >  _cPointsOwned;

    //! Initial point provided by the user.  Could be a list someday.
    GssList  _cInitialList;

    //! The state of the iterator.
    IteratorState  _nState;

    //! Tolerance for saying whether or not we are on a boundary.
    double  _dBoundsTolerance;

    //! True if GSS should randomize the order in which points are generated.
    /*!
      When points are generated and added to the queue in a fixed
      order, it is possible (when running in asynchronous mode) to
      initially optimize over a strict subspace.  The drawback of this
      is that small steps must then be taken to the true solution to
      avoid undoing the gains made in the initial subspace.
    */
    bool  _bUseRandomOrder;

    //! Step length to use when no parent direction is known.
    double  _dInitialStep;

    //! A temporary vector, declared once to improve performance.
    Vector  _tmpReusableVector;

    //! Linear constraints.
    const LinConstr &  _cLinConstr;

    //! True if trial points are modified to snap to nearby linear constraints.
    bool  _bUseSnap;

    //! Scaled distance threshold for choosing constraints to snap onto.
    double  _dSnapDistance;

    //! Sufficient improvement factor.
    double  _dSuffImprovementFactorAlpha;

    //! Maximum number of evaluations this citizen can make.
    int  _nMaxGssEvals;

    //! Total number of evaluated points submitted and received by this citizen.
    int  _nNumGssEvals;
};

}

#endif     //-- HOPSPACK_GSSITERATOR_HPP
