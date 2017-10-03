// $Id: HOPSPACK_GssPoint.hpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss/HOPSPACK_GssPoint.hpp $

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
  @file HOPSPACK_GssPoint.hpp
  @brief Class declaration for HOPSPACK::GssPoint, subclass of DataPoint.
*/
#ifndef HOPSPACK_GSSPOINT_HPP
#define HOPSPACK_GSSPOINT_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_DataPoint.hpp"
#include "HOPSPACK_NonlConstrPenalty.hpp"
#include "HOPSPACK_ProblemDef.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Contains a trial location, evaluation results, and GSS parent information.
/*!
  An instance contains the DataPoint base class members, plus GSS parent
  information and a method to compare evaluation results as a scalar.
*/
//----------------------------------------------------------------------
class GssPoint : public DataPoint
{
  public:

    //! Tag value if a point has no GSS parent.
    static const int  NO_PARENT_TAG = -1;

    //! Direction index value if a point has no parent direction.
    static const int  NO_DIR_INDEX = -1;


    //! Constructor for a point with no parent.
    /*!
     *  @param[in] nObjGoal  Optimization goal: minimize, maximize, etc.
     *  @param[in] cPenalty  Penalty function for nonlinear constraints.
     *  @param[in] cX        Location of the trial point.
     *  @param[in] dStep     Step length used to generate this point.
     */
    GssPoint (const ProblemDef::ObjectiveType  nObjGoal,
              const NonlConstrPenalty &        cPenalty,
              const Vector &                   cX,
              const double                     dStep);

    //! Constructor for a point generated from a parent point.
    /*!
     *  @param[in] nObjGoal       Optimization goal: minimize, maximize, etc.
     *  @param[in] cPenalty       Penalty function for nonlinear constraints.
     *  @param[in] cX             Location of the trial point.
     *  @param[in] dStep          Step length used to generate this point.
     *  @param[in] nParentTag     Tag of the parent of this point.
     *  @param[in] dParentObj     Objective of the parent of this point.
     *  @param[in] dPenaltyTerm   Nonlinear constraint penalty of the parent
     *                            at this point.
     *  @param[in] dSuffImprvAmt  Sufficient improvement amount (>= 0) that
     *                            the point must make over dParentObj.
     *  @param[in] nDirIndex      Direction index that generated this point.
     */
    GssPoint (const ProblemDef::ObjectiveType  nObjGoal,
              const NonlConstrPenalty &        cPenalty,
              const Vector &                   cX,
              const double                     dStep,
              const int                        nParentTag,
              const double                     dParentObj,
              const double                     dPenaltyTerm,
              const double                     dSuffImprvAmt,
              const int                        nDirIndex);

    //! Constructor that makes an instance from a DataPoint not from GSS.
    /*!
     *  The GSS citizen may use DataPoint instances generated from other
     *  citizens or the framework.  The new GssPoint has no parent, and is
     *  typically assigned the GSS default initial step length for purposes
     *  of computing more trial points.
     *
     *  @param[in] cArg      Existing point to be copied.
     *  @param[in] cPenalty  Penalty function for nonlinear constraints.
     *  @param[in] dStep     Step length to use.
     */
    GssPoint (const DataPoint         &  cArg,
              const NonlConstrPenalty &  cPenalty,
              const double               dStep);

    //! Copy Constructor that makes a deep copy.
    GssPoint (const GssPoint &  cArg);

    //! Destructor.
    ~GssPoint (void);


    //! Return the parent tag of the point.
    int  getParentTag (void) const;

    //! Return the step length used to generate this point.
    double  getStepLength (void) const;

    //! Return the index of the direction that generated this point.
    int  getDirIndex (void) const;

    //! Return the evaluated objective function value plus penalty.
    /*!
     *  This overrides the method in HOPSPACK_DataPoint to allow the citizen
     *  to add nonlinear constraints as penalty terms.
     *
     *  If getVecF() contains multiple objectives, then the "best" objective
     *  is determined by calling HOPSPACK_DataPoint::getBestF().
     *  The NonlConstrPenalty associated with this point is used to compute
     *  a penalty term for violated nonlinear constraints.  The term causes
     *  the objective to become worse if constraints are violated.
     *  For instance, if optimization is trying to minimize the objective,
     *  then a positive penalty term is added to the output of getBestF().
     *  If maximizing, then a positive penalty term is subtracted.
     */
    double  getBestF (void) const;

    //! Return true if this point's objective plus penalty is better.
    /*!
     *  The method assumes both points are linearly feasible.
     *  It calls DataPoint::getBestF() to combine multiple objectives and then
     *  make a scalar inequality comparison.  Special rules apply
     *  if the point has not been evaluated, or the objective does not exist.
     *  If nonlinear constraints are present, then the objective includes
     *  a penalty term.
     */
    bool  isBetterObjThan (const GssPoint &  other) const;

    //! Return true if the evaluated objective makes sufficient improvement.
    /*!
     *  The objective of the point must improve over its parent's by an
     *  amount that exceeds the "sufficient decrease" forcing function.
     *  The forcing function is defined in GssIterator where it computes
     *  the sufficient improvement amount for the GssPoint constructor.
     *  The function is
        \f[
          \rho (\Delta) = \alpha \Delta^2
        \f]
     *  where \f$ \Delta \f$ is step length, and \f$ \alpha \f$ is the
     *  user parameter "Sufficient Improvement Factor".
     *  If \f$ \alpha = 0 \f$, then the point need only make "simple decrease"
     *  over the parent point.
     *  If the point has no GSS parent, then the method returns true.
     */
    bool  hasSufficientImprovement (void) const;

    //! Print to the given stream.
    void  print (ostream &   stream,
                 const bool  bIncludeMsg = true) const;
    
  private:

    //! By design, there is no assignment operator.
    GssPoint & operator= (const GssPoint &);


    //! GSS parent point tag, or the value NO_PARENT_TAG if none.
    int  _nParentTag;

    //! Index of the direction that generated this point.
    int  _nDirIndex;

    //! Step length that was used to generate this point.
    double  _dStep;

    //! Parent's objective value.
    double  _dParentObjective;

    //! Sufficient improvement amount needed compared with parent's objective.
    double  _dSufficientImprovementAmount;

    //! Penalty function for nonlinear constraints.
    const NonlConstrPenalty &  _cPenalty;
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_GSSPOINT_HPP
