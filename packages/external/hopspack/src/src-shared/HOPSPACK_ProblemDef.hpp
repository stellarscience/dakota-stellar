// $Id: HOPSPACK_ProblemDef.hpp 166 2010-03-22 19:58:07Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_ProblemDef.hpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2010 Sandia Corporation
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
  @file HOPSPACK_ProblemDef.hpp
  @brief Class declaration for HOPSPACK::ProblemDef.
*/
#ifndef HOPSPACK_PROBLEMDEF_HPP
#define HOPSPACK_PROBLEMDEF_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_ParameterList.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Contains the optimization problem definition.
/*!
    The definition includes the number of objectives and variables,
    objective type, variable types, upper and lower variable bounds, etc.
    All members of the definition are determined from user input parameters.
    The number of general nonlinear constraints are part of the problem
    definition, but general linear constraints are not (see HOPSPACK_LinConstr).
 */
//----------------------------------------------------------------------
class ProblemDef 
{
  public:

    //! Objective type of the optimization problem.
    enum ObjectiveType
    {
        //! Minimize the single objective function.
        MINIMIZE = 0,

        //! Maximize the single objective function.
        MAXIMIZE,

        //! No objective function, just find feasible points.
        FIND_FEASIBLE_PT
    };

    //! Variable types.
    enum VariableType
    {
        //! Variable is continuous and real.
        CONTINUOUS = 0,

        //! Variable must be an integer value, but can temporarily be real.
        INTEGER,

        //! Variable can assume only integer values (aka "categorical").
        ORDINAL
    };


    //! Constructor.
    ProblemDef (void);

    //! Copy constructor.
    ProblemDef (const ProblemDef &);

    //! Destructor.
    ~ProblemDef (void);

    //! Initialize the definition from user input parameters.
    /*!
     *  If an initial point is defined but violates variable bounds,
     *  then it is modified to satisfy bound constraints.
     *
     *  @param[in] cProbDefParams  Sublist of "Problem Def" user parameters.
     *  @return                    false if fatal error encountered.
     */
    bool  initialize (const ParameterList &  cProbDefParams);

  //@{ @name Accessors

    //! Return the number of objective functions (0 or more).
    int  getNumObjs (void) const;

    //! Return the type of objective.
    ObjectiveType  getObjType (void) const;

    //! Return the objective target value, or HOPSPACK::dne() if none.
    /*!
     *  If minimizing, then the goal is to find a point whose objective
     *  is smaller than the target; if maximizing, then the goal is a point
     *  with larger objective value.
     */
    double  getObjTarget (void) const;

    //! Return the objective percent error of the target value.
    /*!
     *  The parameter can be defined only if an objective target is supplied.
     *  If undefined, then the method returns HOPSPACK::dne().
     */
    double  getObjPercentErrorThreshold (void) const;

    //! Return a vector of variable types.
    const vector< VariableType >  getVarTypes (void) const;

    //! Return true if the problem domain is continuous.
    /*!
     *  The domain is continuous if all variables are continuous, or if
     *  integer-valued variables are fixed at a single value by the bounds.
     */
    bool  isDomainContinuous (void) const;

    //! Return a vector of scaling factors.
    /*!
     *  Variable scaling can improve performance dramatically.  The user should
     *  provide scaling factors with parameter "Scaling".  If the parameter
     *  is not provided, then scaling is set from the variable bounds:
     *  \f[
     *    s_i = u_i - l_i
     *  \f]
     *  where \f$s\f$ is a vector of scaling factors, \f$u\f$ the vector of
     *  upper bounds, and \f$l\f$ the vector of lower bounds.
     *  If the upper or lower bound is not finite, then the user <b>must</b>
     *  provide a scaling factor using the "Scaling" parameter.
     */
    const Vector &  getVarScaling (void) const;

    //! Return true if scaling factors were set automatically, not by the user.
    bool  isAutoScaled (void) const;

    //! Return a vector of lower bounds.
    const Vector &  getLowerBnds (void) const;

    //! Return a vector of upper bounds.
    const Vector &  getUpperBnds (void) const;

    //! Return the user-defined initial point, or an empty vector if none.
    const Vector &  getInitialX (void) const;

    //! Return user-defined objective values at the initial point.
    /*!
     *  If no objective is defined or no initial point,
     *  then return an empty vector.
     */
    const Vector &  getInitialF (void) const;

    //! Return user-defined nonlinear equality constraint values.
    /*!
     *  If no equality constraints are defined or no initial point,
     *  then return an empty vector.
     */
    const Vector &  getInitialEqs (void) const;

    //! Return user-defined nonlinear inequality constraint values.
    /*!
     *  If no inequality constraints are defined or no initial point,
     *  then return an empty vector.
     */
    const Vector &  getInitialIneqs (void) const;

    //! Return true if the problem has at least one nonlinear constraint.
    bool  hasNonlinearConstr (void) const;

    //! Return the nonlinear active constraint tolerance.
    /*!
     *  The tolerance determines constraint tests made in the method
     *  isNonlinearlyFeasible().
     */
    double  getNonlinearActiveTol (void) const;

    //! Return the number of nonlinear equality constraints.
    int  getNumNonlinearEqs (void) const;

    //! Return the number of nonlinear inequality constraints.
    int  getNumNonlinearIneqs (void) const;

  //@}

    //! Return true if the point is feasible with respect to bound constraints.
    bool  isBndsFeasible (const Vector &  cX) const;

    //! Return true if the point is feasible for all nonlinear constraints.
    /*!
     *  @param[in] cEqs    Nonlinear equality constraint values at a point.
     *  @param[in] cIneqs  Nonlinear inequality constraint values at a point.
     *
     *  The feasibility test is made against the Nonlinear Active Tolerance.
     *  A point is feasible with respect to a given inequality if it satisfies
     *  the constraint or its unscaled distance is within the tolerance.
     *  A point is feasible with respect to a given equality if its unscaled
     *  distance is within the tolerance.
     *
     *  Note that nonlinear constraints are unscaled distances, versus scaled
     *  distances for linear constraints.  It is not possible to compute a
     *  scaled distance without knowing the normal vector to the constraint,
     *  which requires computing gradients for the nonlinear constraints.
     */
    bool  isNonlinearlyFeasible (const Vector &  cEqs,
                                 const Vector &  cIneqs) const;

    //! Project the argument onto the bounds to make it feasible.
    /*!
     *  @param[in]     dTol  Maximum distance any component can move.
     *                       If negative, then treat dTol as infinity.
     *  @param[in,out] cX    Point to be projected, unscaled.
     *  @return              true if feasibility was reached.
     *
     *  Each component of cX is checked for infeasibility with respect to
     *  variable bounds, and modified if infeasible by at most dTol.
     *  If the infeasibility exceeds dTol in any component, then return false.
     *  The tolerance and cX are not scaled.
     */
    bool  makeBndsFeasible (const double    dTol,
                                  Vector &  cX) const;

    //! Overwrite any user-defined initial point.
    /*!
     *  @param[in] newX  Position of the initial point, which must be feasible
     *                   with respect to variable bounds and any linear
     *                   constraints.  Pass an empty vector if the purpose
     *                   is to erase any initial point.
     *
     *  Any user-defined initial values for the objective or nonlinear
     *  constraints is erased.
     */
    void  resetInitialX (const Vector &  newX);

    //! Overwrite any user-defined initial point and initial values.
    /*!
     *  @param[in] newX      Position of the initial point, which must be
     *                       feasible with respect to variable bounds and any
     *                       linear constraints.
     *  @param[in] newX      Position of the initial point.
     *  @param[in] newF      Objective function values at newX.
     *  @param[in] newEqs    Nonlinear equality constraint values at newX.
     *  @param[in] newIneqs  Nonlinear inequality constraint values at newX.
     *  The number of initial values must match the problem size.
     */
    void  resetInitialX (const Vector &  newX,
                         const Vector &  newF,
                         const Vector &  newEqs,
                         const Vector &  newIneqs);

    //! Return true if dObjValue satisfies the user's target objective.
    /*!
     *  @param[in]  dObjValue  Scalar objective value to test.
     *  @param[out] dPercent   Difference from the target as a percent.
     *                         A value of zero means the target was reached
     *                         or exceeded (i.e., zero percent error).
     *                         Undefined if the method returns false.
     *  @return true           If dObjValue satisfies "Objective Target" or
     *                         is within a specified percent of it.
     *
     *  The point is assumed to be completely feasible.  This is especially
     *  important if dobjValue includes penalty terms for infeasibility.
     */
    bool  isObjTargetReached (const double    dObjValue,
                                    double &  dPercent) const;

    //! Print the problem definition.
    /*!
     *  @param [in]  bDisplayFull  If true, print all variable details
     *                             as directed by the "Display" parameter,
     *                             else print just a problem summary.
     */
    void  printDefinition (const bool  bDisplayFull) const;

  private:

    //! By design, there is no assignment operator.
    ProblemDef & operator= (const ProblemDef &);

    //! Helper methods for the constructor.  Most can throw a FATAL_ERROR.
    bool  setupObj_ (const ParameterList &  cParams);
    bool  setupVars_ (const ParameterList &  cParams);
    bool  setupVarBndsAndScaling_ (const ParameterList &  cParams);
    bool  setupInitialPoint_ (const ParameterList &  cParams);
    bool  setupMisc_ (const ParameterList &  cParams);

    void  printObjDefinition_ (void) const;
    void  printVarSummary_ (void) const;
    void  printVarName_ (const int  nVarNum) const;
    void  printInitPointSummary_ (void) const;


    int            _nNumObjs;          //-- NUMBER OF OBJECTIVE FNS
    ObjectiveType  _nObjType;          //-- TYPE OF OBJECTIVE
    double         _dObjTarget;        //-- STOP IF REACHED
    double         _dObjTgtPercent;    //-- STOP IF REACHED

    int            _nNumVars;          //-- NUMBER OF VARIABLES
    Vector         _cScaling;          //-- VARIABLE SCALING FACTORS
    bool           _bIsAutoScaled;     //-- TRUE IF SCALING SET AUTOMATICALLY
    Vector         _cLoBnds;           //-- VARIABLE LOWER BOUNDS
    Vector         _cUpBnds;           //-- VARIABLE UPPER BOUNDS

    Vector         _cInitialX;         //-- INITIAL POINT
    Vector         _cInitialF;         //-- OBJECTIVES AT INITIAL POINT
    Vector         _cInitialEqs;       //-- NONLINEAR EQS AT INITIAL POINT
    Vector         _cInitialIneqs;     //-- NONLINEAR INEQS AT INITIAL POINT

    int            _nNumNonlEqs;       //-- NONLINEAR EQ CONSTRAINTS
    int            _nNumNonlIneqs;     //-- NONLINEAR INEQ CONSTRAINTS
    double         _dNonlActTol;       //-- ACTIVE CONSTRAINT TOLERANCE

    int            _nDisplayFlag;      //-- VALUE OF 'Display'

    vector< VariableType >  _naVarTypes;
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_PROBLEMDEF_HPP
