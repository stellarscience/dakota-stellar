// $Id: HOPSPACK_LinConstr.hpp 208 2012-08-01 22:59:33Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_LinConstr.hpp $

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
  \file HOPSPACK_LinConstr.hpp
  \brief Class declaration for HOPSPACK::LinConstr.
*/

#ifndef HOPSPACK_LINCONSTR_HPP
#define HOPSPACK_LINCONSTR_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Matrix.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_ProblemDef.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{

//! Contains general linear constraints and methods for using them.
/*!
  HOPSPACK can be configured at build time to exclude linking an LAPACK library,
  which means linear constraints cannot be used.  The configuration appears
  in the code as preprocessing symbol HAVE_LAPACK.  However, simple variable
  bounds are also a type of linear constraint, so instances of this
  class can be constructed even if general linear constraints are excluded.

  Consider the following linearly constrained problem
  \f[
    \begin{array}{rlcrcl}
    \displaystyle\mbox{minimize} & \multicolumn{5}{c}{f(x)} \\
    x \in R^n \\
    \mbox{subject to} & b^{\rm lower} & \le & x & \le & b^{\rm upper}\\
    & b^{\rm ineq\_lower} & \le & A^{\rm ineq} x & \le & b^{\rm ineq\_upper} \\
    & & & A^{\rm eq} x & = & b^{\rm eq}
    \end{array}
  \f]
  Note that inequality bounds may be infinite (i.e., unbounded).
  A variable may be fixed by setting identical lower and upper bounds, or by
  defining a linear equality.

  Variables are always scaled and shifted, and linear constraint computations
  are made in the scaled and shifted problem space.
  The scaling vector is denoted by \f$ s \f$, specified in ProblemDef, and
  defined so that all components are positive and finite.
  Variable shifting moves the lower bound to zero.

  Scaled and shifted quantities are:
  \f[
    \begin{array}{lcl}
    \tilde{x}                   &=& S^{-1}(x - \hat{\ell})  \\
    {\tilde f(\tilde x)}        &=& f(S \tilde x + \hat{\ell}) = f(x) \\
    \tilde{b}^{\rm lower}       &=& S^{-1}(b^{\rm lower} - \hat{\ell}) \\
    \tilde{b}^{\rm upper}       &=& S^{-1}(b^{\rm upper} - \hat{\ell}) \\
    \tilde{b}^{\rm ineq\_lower}
        &=& b^{\rm ineq\_lower} - A^{\rm ineq} \hat{\ell} \\
    \tilde{b}^{\rm ineq\_upper}
        &=& b^{\rm ineq\_upper} - A^{\rm ineq} \hat{\ell} \\
    \tilde{b}^{\rm eq}          &=& b^{\rm eq} - A^{\rm eq} \hat{\ell} \\
    \tilde{A}^{\rm ineq}        &=& A^{\rm ineq} S \\
    \tilde{A}^{\rm eq}          &=& A^{\rm eq} S
    \end{array}
  \f]
  where
  \f[
    S = \left(\begin{array}{ccc}
    s_1 \\
    & \ddots \\
    & &  s_n
    \end{array}\right)
  \f]
  \f[
    \hat{\ell}_ = \left\{
      \begin{array}{cl}
      b^{\rm lower}_i & {\rm if } \; b^{\rm lower}_i > -\infty \\
      0 & {\rm otherwise}.
      \end{array}
    \right.
  \f]

  The scaled problem is then given by
  \f[
    \begin{array}{rlcrcl}
    \displaystyle\mbox{minimize} & \multicolumn{5}{c}{\tilde f(\tilde x)} \\
    \tilde x \in R^n \\
    \mbox{subject to}
      & \tilde b^{\rm lower} & \le & \tilde x & \le & \tilde b^{\rm upper} \\
      & \tilde b^{\rm ineq\_lower} & \le & \tilde A^{\rm ineq} \tilde x
          & \le & \tilde b^{\rm ineq\_upper} \\
      & & & \tilde A^{\rm eq} \tilde x & = & \tilde b^{\rm eq}
    \end{array}
  \f]

  Internally, the code combines bound and linear inequality constraints into
  one matrix, giving
  \f[
    \begin{array}{rlcrcl}
    \displaystyle\mbox{minimize} & \multicolumn{5}{c}{\tilde f(\tilde x)} \\
    \tilde x \in R^n \\
    \mbox{subject to}
      & \hat b^{\rm lower} & \le & \hat A^{\rm ineq} \tilde x & \le
          & \hat b^{\rm upper} \\
      & & & \tilde A^{\rm eq} \tilde x & = & \tilde b^{\rm eq}
    \end{array}
  \f]
  where
  \f[
    \begin{array}{ccccc}
    \hat{b}^{\rm lower} = \left(
      \begin{array}{l}
      \tilde{b}^{\rm lower} \\
      \tilde{b}^{\rm ineq\_lower}
      \end{array}
    \right), 
    & &
    \hat{b}^{\rm upper} = \left(
      \begin{array}{l}
      \tilde{b}^{\rm upper} \\
      \tilde{b}^{\rm ineq\_upper}
      \end{array}
    \right),
    & &
    \hat{A} = \left(
      \begin{array}{l}
      I \\
      \tilde{A}^{\rm ineq}
      \end{array}
    \right).
    \end{array}
  \f]
*/
class LinConstr
{

public:

  /*! \brief State of one constraint with respect to both its upper and lower
    bounds and a given point. 

    A linear constraint is declared active at a point if the scaled distance
    to the constraint is within tolerance \f$ \epsilon \f$.
    The point does not have to be feasible.
    If \f$ \epsilon \f$ is sufficiently large, it is possible that both upper
    and lower bounds may be considered active.
  */
  enum ActiveType {
    //! Constraint is not active at the point.
    NEITHER_ACTIVE = 0,
    //! Lower bound on constraint is active at the point.
    LOWER_ACTIVE,
    //! Upper bound on constraint is active at the point.
    UPPER_ACTIVE,
    //! Both upper and lower bounds on constraint are active at the point.
    BOTH_ACTIVE
  };
  

  //! Constructor.
  /*!
   *  @param[in] probDef  Problem definition with scaling and bounds.
   */
  LinConstr (const ProblemDef &  probDef);

  //! Copy constructer that can exclude equalities.
  /*!
   *  @param[in] cOther    Linear constraint set that will be copied over.
   *  @param[in] bDropEqs  True means equalities are suppressed.
   *  @return              False if fatal error encountered.
   */
  LinConstr (const LinConstr &  cOther,
             const bool         bDropEqs);

  //! Destructor.
  ~LinConstr();

  //! Initialize the definition from user input parameters.
  /*!
   *  @param[in] cLinConstrParams  Sublist of "Linear Constraints" user
   *                               parameters (may be empty).
   *  @return                      False if fatal error encountered.
   */
  bool  initialize (const ParameterList &  cLinConstrParams);

  //@{ \name Accessors
  //! Return value that determines feasibility and activity of a constraint.
  /*!
   *  The default is a small multiple of machine epsilon, but this might
   *  be too tight depending on the condition of the active constraints.
   */
  double getActiveTol() const;

  //! Return true if there are linear constraints in addition to variable bounds.
  bool hasLinearConstraints() const;

  //! Return coefficient matrix for scaled and shifted inequality constraints.
  /*!
   *  The first n rows are the identity matrix, representing bound constraints.
   *  The next m rows are linear inequality constraints.
   *  Methods getBhatLower() and getBhatUpper() return vectors of length n+m
   *  that provide lower and upper bounds on Ahat.
   */
  const Matrix & getAhat() const;

  //! Return lower bounds for scaled and shifted inequality constraints.
  /*!
   *  The first n rows are lower bounds on constraints.
   *  The next m rows are lower bounds on linear inequality constraints.
   *  A value of DNE means there is no bound on the constraint.
   */
  const Vector & getBhatLower() const;

  //! Return upper bounds for scaled and shifted inequality constraints.
  /*!
   *  The first n rows are upper bounds on constraints.
   *  The next m rows are upper bounds on linear inequality constraints.
   *  A value of DNE means there is no bound on the constraint.
   */
  const Vector & getBhatUpper() const;

  //! Return coefficient matrix for scaled and shifted equality constraints.
  const Matrix & getAtildeEq() const;

  //! Return right-hand side for scaled and shifted equality constraints.
  const Vector & getBtildeEq() const;

  //@}
  
  //@{ \name  Feasibility verification methods.

  //! Return true if feasible, false otherwise.
  /*!  
    Simple bounds, linear inequalities, and linear equalities are checked.
    A point is feasible with respect to a given inequality if it satisfies
    the constraint or its scaled distance is within getActiveTol().
    A point is feasible with respect to a given equality if its scaled
    distance is within getActiveTol().

    \param[in] x                    The point to check, unscaled.
    \param[in] bPrintViolationInfo  If true, print the extent that a constraint
                                    is violated.
  */
  bool isFeasible(const Vector& x,
                  const bool bPrintViolationInfo = false) const;

  //! Return the L2 norm of the constraint violations at x.
  /*!
   *  Variable bounds, linear equalities, and linear inequalities are
   *  included.  The result is the L2 norm of the unscaled bound and constraint
   *  violations.
   */
  double  getL2Norm (const Vector &  x) const;

  //! Return the unscaled infinity norm of the constraint violations at x.
  /*!
   *  Variable bounds, linear equalities, and linear inequalities are
   *  included.  The result is the largest unscaled violation of a bound or
   *  general linear constraint.
   */
  double  getLInfNorm (const Vector &  x) const;
  
  //@}


  //@{ \name Transformation methods.
  //! Replace x with a scaled and shifted version.
  /*!
    The affine transformation is \f$x\f$ with \f$ S^{-1} (x - \hat{\ell}) \f$.
   */
  void scale(Vector& x) const;
  
  //! Replace w with an unscaled version.
  /*!
    The affine transformation is \f$w\f$ with \f$ Sw + \hat{\ell} \f$.
  */
  void unscale(Vector& w) const;

  //@}


  //@{ \name Constraint analysis.

  //! Returns maximum feasible step in the interval [0, step] along direction d.  
  /*!
    \param[in] x Current point, unscaled.
    \param[in] d Search direction, which must be in the null space of all
                 equality constraints.
    \param[in] maxLength Maximum step length, unscaled.

    On exit a nonnegative scalar \f$ \alpha \f$ is returned 
    that gives the maximum feasible step that can be made 
    along direction d.
    The method finds all \a blocking \a constraints and then the minimum
    distance to the first such constraint.
    An inequality is a blocking constraint if the search direction has a
    numerically significant component pointing towards its infeasible region.
    Let \f$ \epsilon \f$ be the value returned by getActiveTol().
    An inequality is blocking if it satisfies the following condition:

    \f[
      \begin{array}{lccl}
        \mbox{For finite variable lower bound $i$:}
          & d_i & < & - \epsilon \times scaling_i \\
        \mbox{For finite variable upper bound $i$:}
          & d_i & > &   \epsilon \times scaling_i \\
        \mbox{For linear inequality $a_i^T x \leq b_i$:}
          & a_i^T d & > & \epsilon
      \end{array}
    \f]

    If any blocking inequality constraint is deemed active according to
    getIneqState(), then a value of zero is returned.
    Also, if the search direction d is not in the null space of the equality
    constraints, then a value of zero is returned.
    Otherwise,
    \f[
      \alpha = {\rm min} \left( \{\Delta\} \cup  
      \left\{ \frac{b_i - a_i^T x}{a_i^T d} \Big |\;
          i \in \mathcal{B} \right\} \right)
    \f]
    is returned, where \f$ \mathcal{B} \f$ denotes the set of blocking
    inequality constraints, and \f$ \Delta \f$ is maxLength.

    \note
    Linear algebra roundoff error may cause the step to go slightly beyond
    the variable bounds (as well as slightly off the general linear
    constraints).  Variable bound violations are often considered serious;
    hence, new points computed from the step should be checked using
    isBndsFeasible() in HOPSPACK::ProblemDef.
  */
  double maxStep(const Vector& x, const Vector& d, double maxLength) const;

  //! Fill in the active/inactive state of all inequality constraints.
  /*!
    The ith inequality constraint is considered active at its lower bound if
    \f[
      {\rm  xdist}(i)  <  {\rm epsilon}.
    \f]

    \param xdist Contains the scaled distance to each inequality constraint from
                 a given point, as computed by formDistanceVector().
    \param epsilon Tolerance for determining if a constraint is active.
    \param index Vector with \f$ n+m \f$ entries, where \f$ n \f$ is the number
    of variables and \f$ m \f$ the number of linear inequality constraints.
    On exit each index[i] is set to one ActiveType value:
    - index[i] = NEITHER_ACTIVE if the \f$i_{ith}\f$ constraint is
    inactive at both its lower and upper bounds. 
    - index[i] = LOWER_ACTIVE if the \f$i_{ith}\f$ constraint is
    active at its lower bound, but not its upper bound. 
    - index[i] = UPPER_ACTIVE if the \f$i_{ith}\f$ constraint is
    active at its upper bound, but not its lower bound. 
    - index[i] = BOTH_ACTIVE if the \f$i_{ith}\f$ constraint is
    active at both its upper and lower bound.
  */
  void getActiveIneqIndices(const Vector& xdist,
                            double epsilon,
                            vector<ActiveType>& index) const;

  //! Return the scaled, projected distance from x to each inequality.
  /*!
    \param x (input) The given point, unscaled.
    \param xdist (output) Vector of length \f$ 2 (n+m) \f$, where \f$ n \f$ is
    the number of variables and \f$ m \f$ the number of inequality constraints.
    The first \f$ n+m \f$ elements are filled with the distance from x to each
    lower bound, and the second \f$ n+m \f$ elements filled with distances
    to upper bounds.

    The method computes distance constrained to movement within the
    nullspace of the equality constraints.
    The point is scaled, and the distance is computed for the scaled version
    of the constraints.

    The lower and upper inequality bounds are treated separately by stacking
    as follows:
    \f[
      \left( \begin{array}{r} -\hat A^{\rm ineq} \\
             \hat A^{\rm ineq} \end{array}
      \right)
      \tilde x \le 
      \left( \begin{array}{r} -\hat b^{\rm lower} \\
             \hat b^{\rm upper}\end{array}
      \right).
    \f]
    Let \f$ Z \f$ be a matrix whose columns form an orthonormal basis for
    the null space of \f$ \tilde A^{eq} \f$,
    let \f$ \epsilon \f$ be the value returned by getActiveTol(),
    and define
    \f[
      C = \left( \begin{array}{r} -\hat A^{\rm ineq} \\
                 \hat A^{\rm ineq} \end{array}
          \right),
      \; {\rm and } \quad \;
      d = \left( \begin{array}{r} -\hat b^{\rm lower} \\
                 \hat b^{\rm upper}\end{array}
          \right),
    \f]
    Then xdist is
    \f[
      {\rm xdist}(i) = \left\{
        \begin{array}{ll}
          | c_i^T \tilde x - d_i | / \|Z^T c_i\|_2
              & {\rm if } \; \|Z^T c_i\| > \epsilon, \\
          0   & {\rm if }  \; \|Z^T c_i \| < \epsilon,
                \;{\rm and }\; |c_i^T \tilde x - d_i| < \epsilon, \\
          \infty
              & {\rm if } \; \|Z^T c_i \| < \epsilon,
                \;{\rm and }\; |c_i^T \tilde x - d_i| > \epsilon.
        \end{array}
      \right.
    \f]
    If \f$ |d_i| = \infty \f$ then \f$ {\rm xdist}(i) = \infty \f$.  

    \note
    In exact arithmetic, the
    distance \f$ r(x,a,b) \f$ from a point \f$ x \f$ to the constraint
    \f$ a^T y \le b \f$ within the space spanned by orthonormal
    matrix \f$ Z \f$ is given by
    \f[
      r(x,a,b) = \left\{
        \begin{array}{ll}
        | a^T x - b | / \|Z^T a\|_2 & {\rm if }\;  Z^T a \ne 0, \\
        0 & {\rm if } \;  Z^T a = 0 \; {\rm and } \; |a^T x - b| = 0, \\
        \infty & {\rm if }\; Z^T a  = 0 \; {\rm and } \; |a^T x - b| > 0. \\
        \end{array}
      \right.
    \f]
  */
  void formDistanceVector(const Vector& x, Vector& xdist) const;

  //! Find the closest point that lies on all linear constraints within a given distance.
  /*! 
     \param[in,out] x Current point, replaced by snapped point.
                      In both cases the point is unscaled.
     \param[in] esnap Radius of a ball about x that defines "nearby" constraints.

     The primary purpose of this method is to quickly guess a solution point
     at the "corner" of the feasible region.
     Many derivative-free algorithms will "half-step" towards a boundary,
     asymptotically approaching it.  If the solution is on the boundary, then
     convergence to a specified tolerance may take a long time.  Snapping to
     the constraint is faster, and more accurate if the solution is at a corner.

     This method computes the closest point to the current x that satisfies
     all linear constraints lying within the distance esnap.
     If x violates any linear equality constraints, then the method only
     computes a point satisfying the equality constraints.

     The method begins by scaling the variables and determining all "nearby"
     inequality constraints using the argument esnap.
     A matrix \f$ A_{\rm esnap} \f$ is then formed that consists of the rows
     \f$ \tilde A^{Eq} \f$ followed by the rows of \f$ \hat A \f$ that
     correspond to "nearby" constraints.
     An LQ factorization is formed and \f$ A_{\rm esnap} \f$ redefined if
     necessary to have full row rank, based on the value of getActiveTol().
     The corresponding right-hand side vector \f$ b_{\rm esnap} \f$ is
     also formed.

     To find the point, the method solves an equality-constrained linear least
     squares problem and replaces x with the solution \f$ y^{*} \f$:
     \f[
       \begin{array}{cc}
         \mbox{minimize}    & \|y - x\| \\
         \mbox{subject to}  & A_{\rm esnap} y = b_{\rm esnap}.
       \end{array}
     \f]
  */
  void snapToBoundary(Vector& x, double esnap) const;

  //! Find the closest point to x that satisfies all constraints.
  /*! 
     \param[in,out] x Current point, replaced by feasible point.
                      In both cases the point is unscaled.
     \return          true if successful.

     The modified point will be feasible with respect to linear constraints
     and bound constraints.
  */
  bool  projectToFeasibility (Vector &  x) const;

  //@}

  //! Print the linear constraint definition.
  /*!
   *  @param [in]  bDisplayFull  If true, display all constraint details,
   *                             else just a summary.
   */
  void printDefinition (const bool  bDisplayFull) const;

private:

  /*! \brief State of either the upper or lower (not both) bound of a
    constraint with respect to a given point. 
  */
  enum StateType {
    //! Constraint does not exist, i.e., bound is \f$ \pm \infty \f$.
    DNE,
    //! Constraint is violated at the point.
    VIOLATED,
    //! Constraint is active at the point (within a tolerance).
    ACTIVE,
    //! Constraint is not active and the point is feasible.
    INACTIVE
  };

  //! Specify upper or lower bound type
  enum BoundType {
    //! Upper bound
    UPPER_BOUND,
    //! Lower bound
    LOWER_BOUND
  };

  //! By design, there is no copy constructor.
  LinConstr (const LinConstr &);
  //! By design, there is assignment operator.
  LinConstr & operator= (const LinConstr &);

  //@{ \name Constructor helper methods.
  
  //! Read the coefficient matrix for linear constraints.
  bool setupMatrix(const ParameterList& params);
  
  //! Read the right-hand side for linear constraints.
  bool setupRhs(const ParameterList& params);

  //! Form the scaled system of constraints.
  bool setupScaledSystem();  

  //@}

  //@{ \name Constraint categorizers.

  //! Return the state of inequality constraint i with respect to a given point.
  /*!
    Computes the state of inequality constraint i with respect to the scaled
    point \f$ \tilde x \f$.  The scaled distance to the constraint is calculated
    as

    \f[
      \frac{ | (\hat a_i^{\rm ineq})^T \tilde x - \hat b_i^{\rm lower} | }
           { \| \hat a_i^{\rm ineq} \| } \leq \epsilon
    \f]
    If the scaled distance is within plus or minus getActiveTol(), then the
    constraint is deemed ACTIVE.  If feasible and more than getActiveTol() away,
    then the constraint is INACTIVE.  If infeasible and more than getActiveTol()
    away, then it is VIOLATED.

    \param i (input) Constraint index.  Let \f$ n \f$ be the number
             of variables and \f$m\f$ the number of linear inequality
             constraints.  Constraints 0 thru \f$ n-1 \f$ correspond to
             variable bounds, and constraints \f$ n \f$ thru \f$ n+m-1 \f$
             correspond to linear inequalities.
    \param bType (input) Request UPPER or LOWER bound.
    \param xTilde (input) Scaled and shifted vector.
    \param bPrintViolationInfo (input) If true, then print violation info.
  */
  StateType getIneqState(const int       i,
                         const BoundType bType,
                         const Vector&   xTilde,
                         const bool      bPrintViolationInfo = false) const;

  //! Return the state of equality constraint i with respect to a given point.
  /*!
    Computes the state of equality constraint i with respect to the scaled
    point \f$ \tilde x \f$.  The scaled distance to the constraint is calculated
    as

    \f[
      \frac{ | (\tilde a_i^{\rm eq})^T \tilde x - \tilde b_i^{\rm eq} | }
           { \| \tilde a_i^{\rm eq} \| } \leq epsilon
    \f]
    If the scaled distance is within plus or minus getActiveTol(), then the
    constraint is deemed ACTIVE.
    
    \param i (input) Constraint index.  Let \f$ n \f$ be the number
             of variables and \f$ m \f$ the number of linear equality
             constraints  Constraints 0 thru \f$ n-1 \f$ correspond to
             variable bounds, and constraints \f$ n \f$ thru \f$ n+m-1 \f$
             correspond to linear equalities.
    \param xTilde (input) Scaled and shifted vector.
    \param bPrintViolationInfo (input) If true, then print violation info.
  */
  StateType getEqState(const int     i,
                       const Vector& xTilde,
                       const bool    bPrintViolationInfo = false) const;

  //! Return true if getEqState() finds all equality constraints are feasible.
  bool isEqualityFeasible(const Vector& xtilde,
                          const bool    bPrintViolationInfo = false) const;

  //! Return true if getIneqState() finds all inequality constraints are feasible.
  bool isInequalityFeasible(const Vector& xtilde,
                            const bool    bPrintViolationInfo = false) const;

  //@}

  //! Form quantities needed to snap a point to the boundary.
  /*!
   *  @param[in] xtilde  The point to be snapped to the boundary.
   *  @param[in] esnap  Radius of a ball that defines "nearby" constraints.
   *  @param[out] Asnap  Scaled matrix of equalities and nearby inequalities.
   *  @param[out] bsnap  Corresponding right-hand side of constraints.
  */
  void formSnapSystem(const Vector& xtilde, double esnap,
                      Matrix& Asnap, Vector& bsnap) const;

  //! Print an error message and throw an exception.
  void throwError(const string& fname, const string& msg) const;

  //! Print a count of all constraints.
  void printCounts_ (void) const;

  //! Print the name of an inequality constraint.
  void printIneqName_ (const int  nIneqNum) const;

  //! Print the name of an equality constraint.
  void printEqName_ (const int  nEqNum) const;


  //! Problem definition, including bound constraints and variable scaling.
  const ProblemDef &  probDef;

  //! Feasibility and active set tolerance.
  double _dActiveTol;

  //! Value of 'Display' parameter.
  int _nDisplayFlag;

  //! Variable scaling
  const Vector & scaling;

//@{ \name Unscaled-sytem members.
  //! The coefficient matrix of the unscaled linear inequality constraints.
  Matrix aIneq;

  //! The coefficient matrix of unscaled linear equality constraints.
  Matrix aEq;

  //! The unscaled lower bounds on linear inequality constraints.
  Vector bIneqLower;

  //! The unscaled upper bounds on linear inequality constraints.
  Vector bIneqUpper;

  //! The right-hand side of unscaled equality constraints.
  Vector bEq;
//@}

//@{ \name Scaled-sytem members.
  //! The scaled coefficient matrix of variable bounds and linear inequality constraints.
  Matrix aHat;

  //! Holds two norms of rows of aHat projected into the nullspace of the scaled equality constraints.
  Vector aHatZNorm;
 
  //! The scaled lower bounds on the variable bounds and linear inequality constraints.
  Vector bHatLower;

  //! The scaled upper bounds on the variable bounds and linear inequality constraints.
  Vector bHatUpper;

  //! The scaled coefficient matrix of linear equality constraints.
  Matrix aTildeEq;

  //! The scaled right-hand side of equality coefficient matrix.
  Vector bTildeEq;

  //! Used in transforming x to scaled and shifted space.
  /*!
    Scaled \f$ \tilde{x} \f$ is related to unscaled \f$ x \f$ via the
    affine transformation \f$ x = S \tilde{x} + \hat{l} \f$.
  */
  Vector lHat;
//@}

};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_LINCONSTR_HPP
