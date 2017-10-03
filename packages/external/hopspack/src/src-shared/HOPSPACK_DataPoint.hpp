// $Id: HOPSPACK_DataPoint.hpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_DataPoint.hpp $ 

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
  @file HOPSPACK_DataPoint.hpp
  @brief Class declaration for HOPSPACK::DataPoint.
*/
#ifndef HOPSPACK_DATAPOINT_HPP
#define HOPSPACK_DATAPOINT_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_ProblemDef.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Contains a trial location point plus evaluation results.
/*!
  An instance contains a vector that specifies the location.  If the point
  has been evaluated, then the instance also contains objective function(s)
  and nonlinear constraint results.  Linear constraint results are not
  included.

  Each DataPoint is assigned a unique tag identifier during construction.
  Subclasses should rely on the base class constructor to assign a tag.

  Methods and states for evaluating just F (not FC) are provided for
  simplifying evaluation of problems without nonlinear constraints.
*/
//----------------------------------------------------------------------
class DataPoint 
{
  public:

    //! The evaluation state of the trial point.
    /*!
     *  May extend this later to include evaluation of gradients "GA".
     */
    enum State
    {
        //! Trial point has not yet been evaluated.
        UNEVALUATED = 0,

        //! Trial point functions evaluated by a worker.
        EVALUATED_FC,

        //! Trial point functions obtained from the cache.
        EVALUATED_FC_FROM_CACHE
    };


    //! Constructor.
    /*!
     *  @param[in] nObjGoal  Optimization goal: minimize, maximize, etc.
     *  @param[in] cX        Location of the trial point.
     */
    DataPoint (const ProblemDef::ObjectiveType  nObjGoal,
               const Vector &                   cX);

    //! Copy Constructor, performs deep copy.
    DataPoint (const DataPoint &  cArg);

    //! Assignment operator, performs deep copy.
    DataPoint &  operator= (const DataPoint &  cArg);

    //! Destructor.
    ~DataPoint (void);


    //@{ @name Accessors

    //! Return x.
    const Vector &  getX (void) const;

    //! Return a vector of evaluated objective function values.
    /*!
     *  Multiple objectives are allowed.  If an evaluation failed, the returned
     *  Vector will be empty.  If an individual function element failed,
     *  the returned element will be HOPSPACK::dne().
     */
    const Vector &  getVecF (void) const;

    //! Return the single best evaluated objective function value.
    /*!
     *  If getVecF() contains a single objective, then this method returns
     *  its value.
     *  If getVecF() contains multiple objectives, then this method returns
     *  the best single value; i.e., the largest value if the goal is to
     *  maximize, or the smallest value if the goal is to minimize.
     *  If all objective evaluations failed, the method returns HOPSPACK::dne().
     *
     *  Citizens should subclass and override this method to return a "best"
     *  objective value according to their needs.  For example, nonlinear
     *  constraints might be treated as penalty terms by using an instance of
     *  NonlConstrPenalty.  The penalty object is not part of this base class
     *  so that different citizens can run simultaneously with different
     *  penalty functions.
     */
    double  getBestF (void) const;

    //! Return a vector of evaluated nonlinear equality constraints.
    /*!
     *  Length always matches the number of nonlinear equalities.
     *  If an individual constraint failed to evaluate, the returned element
     *  will be HOPSPACK::dne().
     */
    const Vector &  getEqs (void) const;

    //! Return a vector of evaluated nonlinear inequality constraints.
    /*!
     *  Nonlinear inequalities have a negative value if in violation of the
     *  constraint.  A nonnegative value means there is no violation.
     *  Length always matches the number of nonlinear inequalities.
     *  If an individual constraint failed to evaluate, the returned element
     *  will be HOPSPACK::dne().
     */
    const Vector &  getIneqs (void) const;

    //! Return the L2 norm of nonlinear constraint violations at x.
    /*!
     *  Nonlinear equalities and inequalities are included.  The result is
     *  the L2 norm of the unscaled constraint violations.
     */
    double  getNonlConstrL2Norm (void) const;

    //! Return the infinity norm of nonlinear constraint violations at x.
    /*!
     *  Nonlinear equalities and inequalities are included.  The result is
     *  the largest unscaled violation of a constraint.
     */
    double  getNonlConstrLInfNorm (void) const;

    //! Return a number for multiplying a positive penalty term.
    /*!
     *  If the objective is to minimize, then return 1.
     *  If the objective is to maximize, then return -1.
     *  If the objective is to find a feasible point, then return 1.
     */
    double  getPenaltySign (void) const;

    //! Return the unique tag of the trial point.
    int  getTag (void) const;

    //! Return the evaluation state of the point.
    State  getState (void) const;

    //@}

    //@{ \name Manipulators

    //! Store evaluations results.
    /*!
     *  @param[in] cFns    Objective values.
     *                     Pass an empty vector if there are none.
     *  @param[in] cEqs    Nonlinear equality constraints.
     *                     Pass an empty vector if there are none.
     *  @param[in] cIneqs  Nonlinear inequality constraints.
     *                     Pass an empty vector if there are none.
     *  @param[in] sMsg    Message associated with the evaluation.
     */
    void  setEvalFC (const Vector &  cFns,
                     const Vector &  cEqs,
                     const Vector &  cIneqs,
                     const string &  sMsg);

    //! Store evaluations results obtained from the cache.
    /*!
     *  @param[in] cFns    Objective values.
     *  @param[in] cEqs    Nonlinear equality constraints.
     *  @param[in] cIneqs  Nonlinear inequality constraints.
     *  @param[in] sMsg    Message associated with the evaluation.
     */
    void  setCachedFC (const Vector &  cFns,
                       const Vector &  cEqs,
                       const Vector &  cIneqs,
                       const string &  sMsg);

    //@}

    //@{ \name Comparisons

    //! Return true if this point's objective is better than the argument's.
    /*!
     *  @param[in] cOther               Other point to compare against.
     *  @param[out] bAreObjsComparable  True if both objectives are
     *                                  comparable (see the cases below).
     *
     *  The method assumes both points are linearly feasible and ignores
     *  nonlinear constraints.
     *
     *  Generally, use getBestF() to combine multiple objectives and then
     *  make a scalar inequality comparison.  If the objective is equal, then
     *  use a tie-breaker rule based on the tag number.  Special rules apply
     *  if one or both points have not been evaluated, or the objective does
     *  not exist.  The logic, in order, is:
     *
     *  <ul>
     *  <li> The objective goal is to become feasible (objective is irrelevant).
     *       Return true (arbitrarily).
     *  <li> One point has been evaluated while the other has not.
     *       Return that the evaluated point as better.
     *  <li> Neither point has been evaluated.
     *       Return the tie-breaker comparison.
     *  <li> Tags are equal (the two evaluated points are identical).
     *       Return false.
     *  <li> One point's objective exists while the other does not.
     *       Return the objective that exists as better.
     *  <li> Neither point's objective exists.
     *       Return the tie-breaker comparison.
     *  <li> Both objectives exist and are not equal.
     *       Return a value showing which is better, based on the objective goal
     *       (MINIMIZE or MAXIMIZE), and report the objectives are comparable.
     *  <li> Both objectives exist and are equal.
     *       Return the tie-breaker comparison, and report the objectives
     *       are comparable.
     *  </ul>
     */
    bool  isBetterObjThan (const DataPoint &  cOther,
                                 bool      &  bAreObjsComparable) const;

    //! Return true if no component differs by more than the tolerance.
    /*!
     *  @param[in] cOther      Other point to compare against.
     *  @param[in] dTolerance  Comparison tolerance.  If the two vectors
     *                         are \f$ x \f$ and \f$ y \f$ and dTolerance
     *                         is denoted by \f$ \tau \f$, then the method
     *                         returns true if
     *                         \f[
     *                             \| x - y \|_{\infty} \leq \tau .
     *                         \f]
     */
    bool  isSamePoint (const DataPoint &  cOther,
                       const double       dTolerance) const;

    //@}

    //@{ \name Printing
 
    //! Print to the given stream.  No terminating newline is printed.
    /*!
     *  @param[in,out] stream   Print to this stream.
     *  @param[in] bIncludeMsg  If true then also print the evaluator message
     *                          associated with the point.
     *  @param[in] bPrintAllX   If true then print X no matter how long;
     *                          otherwise, print only if X has less than 10
     *                          elements.
     */
    void  leftshift (      ostream &  stream,
                     const bool       bIncludeMsg,
                     const bool       bPrintAllX = false) const;

    //! Debug function for analyzing memory leaks.
    /*!
     *  Add calls to this method at various places in the Mediator to help
     *  find memory leaks.
     */
    static void  debugPrintMemoryLists (void);

    //@}

  protected:

    //! Optimization goal, visible to subclasses.
    ProblemDef::ObjectiveType  _nObjGoal;

  private:

    //! Static counter used to generate unique tags.
    static int  _nTagGlobalCounter;

    //! If debugging memory leaks, track of all allocated and deleted points.
    static bool           _bDebuggingLeaks;
    static int            _nDebugTagCounter;
    int                   _nDebugTag;
    static vector< int >  _cDebugCreateList;
    static vector< int >  _cDebugDeleteList;

    //! The x vector location of the point.
    Vector  _cX;

    //! Evaluated objective function values.
    Vector  _cFns;

    //! Evaluated nonlinear equality constraints.
    Vector  _cEqs;

    //! Evaluated nonlinear inequality constraints.
    Vector  _cIneqs;

    //! Unique tag value.
    int  _nTag;

    //! Message string from evaluator.
    string  _sMsg;

    //! Trial point state.
    State  _nState;
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_DATAPOINT_HPP
