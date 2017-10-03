// $Id: HOPSPACK_SolveLinConstrProj.hpp 164 2010-03-15 18:53:15Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_SolveLinConstrProj.hpp $

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
  @file HOPSPACK_SolveLinConstrProj.hpp
  @brief Class declaration for HOPSPACK::SolveLinConstrProj.
*/
#ifndef HOPSPACK_SOLVELINCONSTRPROJ_HPP
#define HOPSPACK_SOLVELINCONSTRPROJ_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_Matrix.hpp"
#include "HOPSPACK_ProblemDef.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Provides a method for solving the linear constraint projection problem.
/*!
 *  The goal is to find the nearest feasible point to a given point \f$ y \f$
 *  that satisfies all linear constraints.  Mathematically:
 *
 *  \f[
 *    \begin{array}{ccl}
 *      \mbox{minimize} & || y - x ||_2  \\
 *      \mbox{s.t.}     & A_{eq} x = b_{eq}
 *                          &  \mbox{(linear equalities)}  \\
 *                      & A_{ineq} x >= b_{ineq}
 *                          &  \mbox{(linear inequalities)}  \\
 *                      & b_{up} >= x >= b_{lo}
 *                          &  \mbox{(variable bounds)}
 *    \end{array}
 *  \f]
 *
 *  All variables are assumed continuous.
 *
 *  There are many possible ways to solve the problem.  The default
 *  implementation is a standard active set method that solves equality
 *  constrained subproblems using LAPACK dgglse.  Thus, computations use
 *  dense matrices and performance does not scale well for large numbers
 *  of variables.  Furthermore, the algorithm assumes linear independence
 *  of the active set of constraints at the solution and all intermediate
 *  iterates.
 */
//----------------------------------------------------------------------
class SolveLinConstrProj
{
  public:

    //! Constructor.
    SolveLinConstrProj (void);

    //! Destructor.
    ~SolveLinConstrProj (void);

    //! Solve to project a given point onto all linear constraints.
    /*!
     *  Linear equalities and inequalities are passed in cLinConstr, and
     *  variables bounds are passed in cProbDef.
     *  The test for constraint satisfaction is delegated to cLinConstr, which
     *  makes use of the "Active Tolerance" configuration parameter.
     *
     *  @param[in]  cProbDef     Problem definition, including variable bounds.
     *  @param[in]  cLinConstr   Linear constraints object.
     *  @param[in]  vX           The point to project, unscaled.
     *  @param[out] vProjection  The projected point, unscaled.
     *  @return                  True if successful.
     */
    bool  solve (const ProblemDef &  cProbDef,
                 const LinConstr  &  cLinConstr,
                 const Vector     &  vX,
                       Vector     &  vProjection);


  private:

     //! By design, there is no copy constructor.
    SolveLinConstrProj (const SolveLinConstrProj &);
    //! By design, there is no assignment operator.
    SolveLinConstrProj & operator= (const SolveLinConstrProj &);

    //! Find an inequality feasible point (ignore equalities).
    /*!
     *  Formulate a subproblem with violation variables \f$ v \f$ for all general
     *  inequalities.  Assume the initial point is already feasible with
     *  respect to variable bounds.  Equalities are ignored.
     *
     *  Solve the subproblem with an active set method to minimize
     *  \f$ || v ||_2 \f$ subject to inequalities.  A simpler linear objective
     *  \f$ \sum v_i \f$ would work, but LAPACK provides only a least squares
     *  solver.  Initialize violation variables to make the subproblem's
     *  initial point feasible.
     *
     *  The method fails if constraints are inconsistent or linearly dependent.
     *
     *  @param[in]     cLinConstr  Linear constraints for testing feasibility.
     *  @param[in]     mMatIneqs   Matrix of linear inequality constraints.
     *  @param[in]     vLoIneqs    Lower bounds on cMatIneqs.
     *  @param[in]     vUpIneqs    Upper bounds on cMatIneqs.
     *  @param[in,out] vX          On entry this is the initial point.
     *                             On exit this is a feasible point, assuming
     *                             the method returns true.
     *  @return                    True if successful.
     */
    bool findFeasibleIneqPoint_ (const LinConstr &  cLinConstr,
                                 const Matrix    &  mMatIneqs,
                                 const Vector    &  vLoIneqs,
                                 const Vector    &  vUpIneqs,
                                       Vector    &  vX) const;

    //! Starting from an inequality feasible point, find the closest point to the target.
    /*!
     *  Use an active set method to minimize \f$ || x_{target} - x ||_2 \f$
     *  subject to all constraints.  Note that \f$ x \f$ must be feasible
     *  with respect to all inequalities on entry.
     *  The method fails if constraints are linearly dependent.
     *
     *  @param[in]     mMatEqs    Matrix of linear equality constraints.
     *  @param[in]     vRhsEqs    Right-hand side of mMatEqs.
     *  @param[in]     mMatIneqs  Matrix of linear inequality constraints.
     *  @param[in]     vLoIneqs   Lower bounds on cMatIneqs.
     *  @param[in]     vUpIneqs   Upper bounds on cMatIneqs.
     *  @param[in]     vXtarget   Target to get close to.
     *  @param[in,out] vX         On entry this is the initial point.
     *                            On exit this is an optimal point, assuming
     *                            the method returns true.
     *  @return                   True if successful.
     */
    bool findClosestPoint_ (const Matrix &  mMatEqs,
                            const Vector &  vRhsEqs,
                            const Matrix &  mMatIneqs,
                            const Vector &  vLoIneqs,
                            const Vector &  vUpIneqs,
                            const Vector &  vXtarget,
                                  Vector &  vX) const;

    //! Minimize \f$ ||c - d^T x||_2 \f$ subject to equalities and inequalities.
    /*!
     *  The initial point must be feasible with respect to all inequalities.
     *  The method can fail if constraints are linearly dependent at an iterate.
     *
     *  @param[in]  vC         Fixed cost vector.
     *  @param[in]  vD         Scaling cost vector.
     *  @param[in]  vXinit     Feasible initial point.
     *  @param[in]  mMatEqs    Matrix of linear equality constraints.
     *  @param[in]  vRhsEqs    Right-hand side on mMatEqs.
     *  @param[in]  mMatIneqs  Matrix of linear inequality constraints.
     *  @param[in]  vLoIneqs   Lower bounds on cMatIneqs.
     *  @param[in]  vUpIneqs   Upper bounds on cMatIneqs.
     *  @param[out] vXsol      Solution point.
     *  @return                True if successful.
     */
    bool  computeActiveSetSolution_ (const Vector &  vC,
                                     const Vector &  vD,
                                     const Vector &  vXinit,
                                     const Matrix &  mMatEqs,
                                     const Vector &  vRhsEqs,
                                     const Matrix &  mMatIneqs,
                                     const Vector &  vLoIneqs,
                                     const Vector &  vUpIneqs,
                                           Vector &  vXsol) const;

    //! Calculate the unconstrained solution to min \f$ ||c - d^T x||_2 \f$.
    /*!
     *  If a component d_i is zero, then use c_i for that variable.
     *
     *  @param[in]  vC      Fixed cost vector.
     *  @param[in]  vD      Scaling cost vector.
     *  @param[out] vXsol   Solution point.
     */
    void  calcUnconstrainedSolution_ (const Vector &  vC,
                                      const Vector &  vD,
                                            Vector &  vXsol) const;

    //! Compute Lagrange multipliers for the active set.
    /*!
     *  This method is called by computeActiveSetSolution_ to help solve
     *  the constrained least squares problem.
     *
     *  Multipliers are found by solving the overdetermined least squares
     *  problem \f$ A^T \lambda = \nabla f \f$, where \f$ A \f$ is the matrix
     *  of active constraints.
     *  If an inequality multiplier is negative, it means the constraint
     *  should be removed so the objective can be reduced.  The method
     *  computes all multipliers and returns the index of the most negative
     *  multiplier for an inequality.
     *
     *  @param[in]  vC         Fixed cost vector.
     *  @param[in]  vD         Scaling cost vector.
     *  @param[in]  mMatCons   Matrix of active constraints, eqs before ineqs.
     *  @param[in]  nNumEqs    Number of equalities in mMatCons.
     *  @param[in]  vX         Feasible point at which to compute multipliers.
     *  @param[out] nConIndex  Index of most negative inequality multiplier,
     *                         or -1 if none are negative.
     *  @return                True if successful.
     */
    bool  computeMultipliers_ (const Vector &  vC,
                               const Vector &  vD,
                               const Matrix &  mMatCons,
                               const int       nNumEqs,
                               const Vector &  vX,
                                     int    &  nConIndex) const;

    //! Return true if the point roughly satisfies a set of inequalities.
    /*!
     *  The feasibility tolerance in this method is somewhat loose and not
     *  configurable.  The intent is to provide a sanity check while solving
     *  scaled subproblems.
     *
     *  @param[in]  vX         Point to check.
     *  @param[in]  mMatIneqs  Matrix of linear inequality constraints.
     *  @param[in]  vLoIneqs   Lower bounds on cMatIneqs.
     *  @param[in]  vUpIneqs   Upper bounds on cMatIneqs.
     *  @return                True if feasible.
     */
    bool  isIneqFeasible_ (const Vector &  vX,
                           const Matrix &  mMatIneqs,
                           const Vector &  vLoIneqs,
                           const Vector &  vUpIneqs) const;

    //! Tolerance for assessing bad steps in computeActiveSetSolution_.
    double  _dActiveTol;
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_SOLVELINCONSTRPROJ_HPP
