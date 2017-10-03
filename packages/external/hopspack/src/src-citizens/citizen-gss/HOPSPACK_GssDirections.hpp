// $Id: HOPSPACK_GssDirections.hpp 216 2013-11-13 23:34:51Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss/HOPSPACK_GssDirections.hpp $

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
  \file HOPSPACK_GssDirections.hpp
  \brief Class declaration for HOPSPACK::GssDirections.
*/

#ifndef HOPSPACK_GSSDIRECTIONS_HPP
#define HOPSPACK_GSSDIRECTIONS_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_GssPoint.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_Matrix.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_ProblemDef.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{

//! Generate and manage GSS search directions and associated information.
/*!  
  This object generates an appropriate set of search directions from a given
  "new best point", then tracks and updates the associated step lengths as
  the search from that progresses.  A new set of directions is initialized
  when computeNewDirections() is called.

  For each direction, track the desired step length, the true step length
  (if the desired step is outside the feasible region, then correct with a
  step length that maintains feasibility), and the trial point tag to the step.
  An internal cache of directions is maintained for greater efficiency.

  The true step and tag are actually calculated by
  GssIterator::generateTrialPoints().

  If the problem contains general linear constraints (as opposed to simple
  variable bounds), then the solution requires access to an LAPACK library.
  If constraints become nearly degenerate at a point, then the CDD library
  can be called to improve robustness and performance.

  <b>Parameters</b>

  These user parameters are parsed from the ParameterList passed to the
  constructor.

  <ul>
  <li>"Step Tolerance"
  <li>"Minimum Step"
  <li>"Contraction Factor"
  <li>"Epsilon Max"
  <li>"Add Projected Normals"
  <li>"Add Projected Compass"
  </ul>
*/
class GssDirections
{

public:

  //! Constructor 
  /*!
    Reads and validates relevant user parameters.
  */
  GssDirections(const ProblemDef    &  probDef,
                const LinConstr     &  linConstr,
                      ParameterList &  params);

  //! Destructor 
  ~GssDirections();
  
  //@{ \name Accessors

  //! Returns the indices of directions that are ready for new trial points.
  /*!  
    A direction is ready for a new trial point if its associated
    step length is greater than or equal to the convergence tolerance
    and it doesn't already have an active trial point.
  */
  void  getDirectionIndices(vector<int>& idvector) const;

  //! Returns the i-th search direction
  const Vector& getDirection(int i) const;

  //! Returns the i-th step
  double getStep(int i) const;

  //! Returns true if every step length is less than step convergence tolerance.
  /*!
    The steps are stored in the step vector.
    The convergence tolerance is stored in stepTolerance.
  */
  bool isStepConverged() const;

  //! Returns true if direction is empty, false otherwise.
  bool empty() const;
  //@}


  //@{ \name Manipulators

  /*! \brief Computes a new set of search directions based on a new best point.

    Search directions are computed by one of three private member functions:
    - buildWithNothing(),
    - buildWithLapack() (requires LAPACK library),
    - buildWithCddLib() (requires CDD library).

    If only bound constraints are present, search direction are computed
    exclusively by buildWithNothing().  
    If general linear constraints are present, a call is first made
    to buildNormalCone().  This method returns the following:
    - either a boolean value of true signifying that search directions for
      the current normal cone have already been generated and is stored
      in private member \e directionCache, or
    - a boolean value of false and the current normal cone.

    If the search directions are already in directionCache, computeNewDirections()
    returns these directions.
    If buildNormalCone() returns false, search directions must be computed.

   \b Computing \b new \b search \b directions:

   generateUnconstrained() is called whenever the problem looks locally
   unconstrained, i.e. the normal cone \f$ = \{ \emptyset \} \f$.
   Otherwise the search directions are computed with buildWithLapack()
   if possible.  If buildWithLapack() is unsuccessful (this can happen if
   nearby constraints are degenerate), then buildWithCddLib() is called.

   Finally, if user parameter "Add Projected Normals" is true, the projected
   normals are also added to the list of search directions.

   On exit the set of newly computed search direction are added to the direction  
   cache and the current search directions are updated.

    computeNewDirections() also resets the step information.
    The initial step (which is the same for all directions) is calculated as
    \f[
    \Delta_{\rm new} = \max \{ \Delta_{\rm old}, \Delta_{\min} \}.
    \f]
    Here, \f$\Delta_{\min}\f$ is minStep.
  */
  void computeNewDirections(const GssPoint& newPoint);

  //! Adds new directions when point remains the same but step is reduced.
  /*! The purpose of appendNewDirections is to handle the case when the step size
      is reduced.  Because HOPSPACK is asynchronous, the step size may be reduced
      for a given direction before the current point is rejected.  In a sense, we
      do not know yet if \f$ x_{k+1} = x_k \f$, signifying a failed iteration, or
      if a new best point will be found.  For failed iterations, the current
      point remains the same and the step size is reduced.  Reducing the step
      size may alter the relevant tangent cone, i.e. fewer constraints will be
      seen as blocking, increasing the number of relevant search directions.
      Thus, to remain theoretically convergent, trial-points must be appended to
      the queue corresponding to these new directions in the advent that the
      current iteration fails.  For this reason, appendNewDirections() will be
      called each time the reduced step
      \f$ \hat \Delta_{k}^i \f$ corresponding to the ith direction \f$d_i\f$
      satisfies
      \f[
         \hat \Delta_{k}^i < \min_{ 1 \le i \le n_d} \Delta_k^i
      \f]
      where \f$ n_d = \f$ equals the number of search directions.  
      \see computeNewDirections().
  */
  void appendNewDirections();

  //! Sets the step corresponding to the ith direction = stepTolerance/2.
  void setStepConverged(int i);

  //! Set the true step and tag for the trial point corresponding to direction i.
  void setTrueStepAndTag(int i, double trueStep_in, int tag_in);

  //! Reduce the step corresponding to direction i. 
  /*! 
    The new step is calculated as 
    \f[
    \Delta_{\rm new} = 0.5 \Delta_{\rm old}.
    \f]
    Also resets the corresponding trueStep and tag. 
  */
  void reduceStep(int i);

  //@}

  //! Print the directions, preceded with a label.
  void print(const string label) const;

  //! Return the smallest step size for directions that have not converged.
  double getSmallestStep() const;

private:

  //! By design, there is no copy constructor.
  GssDirections (const GssDirections &);
  //! By design, there is no assignment operator.
  GssDirections & operator= (const GssDirections &);

  //@{ \name Cone builders.

  /*! \brief This method builds the normal cone generators based on the current
    setting of the state vector constraintState.

    Let \f$V_p\f$ denote a matrix whose columns are obtained from the outward
    pointing normals of the active inequality constraints which are active
    on a single side.

    Let \f$ V_\ell \f$ denote a matrix
    whose columns are obtained from the normals of the inequality constraints
    which are active at both their upper and lower bounds and the equality
    constraints.
    Then the normal cone is given by
    \f[
    \mathcal{N}(x,\epsilon) = \{v |\; v = V_p \lambda + V_\ell \alpha, \lambda > 0 \}.
    \f]
    
    \param VpT (output) The transpose of \f$V_p\f$.
    \param VlT (output) The transpose of \f$V_\ell\f$.

    For more on how the active set is determined see updateConstraintState().
  */
  void buildNormalCone(Matrix& VpT, Matrix& VlT) const;

  //! This method builds the tangent cone base upon VpT and VlT. 
  /*! This method assumes the normal cone is given by
    \f[
    \mathcal{N}(x,\epsilon) = \{v |\; v = V_p \lambda + V_\ell \alpha,
                                      \lambda > 0 \}.
    \f]
    \param VpT (input) the transpose of \f$ V_p \f$
    \param VlT (input) the transpose of \f$ V_\ell \f$
    \param T (output) generator for the tangent cone.

    If VpT and VlT are empty, the method will simply call generateUnconstrained().    Otherwise buildWithLapack() is called.  If buildWithLapack() is unable to
    compute the tangent cone, buildWithCddLib() will next be called.  This can
    happen if the normal cone corresponds to degenerate constraints.
  */
  void buildTangentCone(const Matrix& VpT, const Matrix& VlT,
                        Matrix& T);

  //! Generate compass search directions for the variable bounds only case.
  /*!
    This method in unique in that it does not require LAPACK or CDDLIB.
    Thus GSS can reliably run, if only variable bounds are present,
    without either package.
  */
  void buildWithNothing(Matrix& D);

  //! Generates pattern search directions using CDDLIB.
  /*! A positive spanning set for the tangent to the normal cone
    \f[
    \mathcal{N} = \{v |\;  v = V_p \lambda + V_\ell \alpha, \lambda > 0 \}
    \f]
    is generated by a call to compute_cone_generators(),
    where \f$ V_p \f$ and \f$ V_\ell\f$ are inputed in transpose form.
    
    \note Method is called only if general linear constraints exist and a
    call to buildWithLapack() returns false.
    \param VpT (input) Transpose of \f$ V_p \f$.
    \param VlT (input) Transpose of \f$ V_\ell \f$.
    \param T (output) The tangent cone generators.
   */
#ifdef HAVE_CDDLIB
  bool buildWithCddLib(const Matrix& VpT, const Matrix& VlT,
                       Matrix& T);
#endif

  //! Generate search directions using LAPACK
  /*!
    The inputs define the \f$\epsilon\f$-normal cone
    \f[
    \mathcal{N} = \{v |\; v =  V_p \lambda + V_\ell \alpha, \lambda \geq 0 \}.
    \f]
    (Note that \f$V_p\f$ and \f$ V_\ell\f$ are provided in transpose
    form.)  The directions we need are the generators of the tangent
    cone \f$ \mathcal{T} = \mathcal{N}^\circ \f$.
   
    Let \f$ Z \f$ be a matrix whose columns form a basis for \f${\rm
    null}(V_\ell)\f$; let \f$ R \f$ be a right inverse for \f$ V_p^T Z
    \f$ (if it exists), i.e., \f$ V_p^T Z R = I\f$; and let \f$ N\f$
    be a matrix whose columns form a basis for \f${\rm null}(V_p^T
    Z)\f$.  

    Then: \f[ \mathcal{T} = \{w |\; w = -ZRu + ZN\xi, u \geq 0
    \}, \f]
    
    \param VpT (input) Transpose of \f$ V_p \f$. 
    \param VlT (input) Transpose of \f$ V_\ell \f$.
    \param T (output) The tangent cone generators.

    \return False if a right inverse for \f$ V_p^T Z \f$ does not
    exist; the search directions cannot be calculated with
    LAPACK. Otherwise, it returns true.
  */
  bool buildWithLapack(const Matrix& VpT, const Matrix& VlT,
                       Matrix& T);

  //@}
  
  //@{ \name Direction generating methods

  //! Generates compass search directions for the "looks locally unconstrained" case.
  void generateUnconstrained(Matrix& D);
 
  /*!\brief This method is called to generate search directions whenever
    general linear constraints are present. 
    
    The method first checks the direction cache, directionCache, directions for
    the current setting of constraintState have already been computed.  If the
    corresponding search directions have been found within the cache, input
    parameter D is updated accordingly.  If search directions are not found in
    the cache, the normal cone generators N are formed and the tangent cone T
    computed (stored by rows).

    If no tangent directions are found, i.e. T is empty, then
    - D = T.
    - D (an empty matrix) is cached with constraintState.
    - If T is empty, the projected normals are never add to D.

    If nontrivial tangent directions are found, then
    - D = T if withNormals = true.  
    - D = [T; PN] otherwise. Here PN is used to denote the projected normals, see
    addNormalDirections().
    - Prior to exiting, D is cached with constraintState.

    \note \li withNormals is a user-definable parameter via string
              "Add Projected Normals".
          \li By always returning D empty whenever T is empty we ensure that
              we never return a nontrivial set of infeasible directions.
   */
  void generateForLinear(Matrix& D);

  //! Add in the projected normals of the active linear inequality constraints.
  /*!
    This method assumes the normal cone is given by
    \f[
    \mathcal{N}(x,\epsilon) = \{v |\; v = V_p \lambda + V_\ell \alpha, \lambda > 0 \}.
    \f]
    \param VpT (input) the transpose of \f$ V_p \f$
    \param VlT (input) the transpose of \f$ V_\ell \f$
    \param D (input/output)  on exit, the rows of \f$ V_p^TZZ^T \f$ are added
                             to the bottom of matrix D after first being
                             normalized, then scaled.
    Here \f$ Z \f$ denotes an orthonormal basis matrix for the null space of
    \f$ V_\ell^T \f$.

    \warning addNormalDirections() requires LAPACK to perform the projection
    into the null-space of \f$ V_\ell^T \f$.  Thus if \f$ V_\ell\f$ is not empty
    and withNormals equals true, the method will throw an error and exit.
    A subtle point is that \f$ V_\ell \f$ can be non-empty despite the absence
    of equality constraints. See buildNormalCone().
  */
  void addNormalDirections(const Matrix& VpT, const Matrix& VlT,
                           Matrix& D);

  //! Method adds projected compass search direction to D.
  /*!
    \param VlT (input) the transpose of \f$ V_\ell \f$
    \param D (input/output)  on exit, the rows of \f$ \pm ZZ^T\f$ are added to
                             the bottom of matrix D after first being normalized,
                             then scaled.

    Here \f$ Z \f$ denotes an orthonormal basis matrix for the null space of
    \f$ V_\ell^T \f$.

    \note If \f$VlT\f$ is empty, then the classical compass search directions
    \f[
    \{\pm e_1, \pm e_2, \ldots, \pm e_n \}
    \f] are added after being scaled.

    \warning addNormalDirections() requires LAPACK to perform the projection
    into the null-space of \f$ V_\ell^T \f$.  Thus if \f$ V_\ell\f$ is not empty
    and withCompass equals true, the method will throw an error and exit.
    A subtle point is that \f$ V_\ell \f$ can be non-empty despite the absence
    of equality constraints.
  */
  void addCompassDirections(const Matrix& VlT, Matrix& D);

  //! Update the state according to distances stored within xDistance. 
  /*! 
    To determine which constraints are active (or near active), 
    an epsilon ball \f$\mathcal{B}(x,\epsilon)\f$ is defined about the current
    point, with x denoting the current point and \f$ \epsilon \f$ defined by
    \f[
    \epsilon = \min(\epsilon_{\rm max}, \Delta_k).
    \f]
    Here \f$\epsilon_{\rm max} = \f$epsilonMax, and \f$ \Delta_k \f$ is given by 
    input parameter newStep.
    All constraints which pass through this epsilon ball are considered active.
    Essentially, the \f$i_{th}\f$ constraint is considered active if 
    \f[
    \frac{|a_i^T x - b_i|}{\|a_i\|_2} \le \epsilon.
    \f]
    Activity of constraints is always determined in the scaled space.  
    This method updates constraintState via a call to
    LinConstr::getActiveIneqIndices(), passing in as parameters xDistance and
    epsilon.  The parameter xDistance provides a measure of the distance to each
    inequality constraint from the current point; it is assumed that
    xDistance has already been updated prior to function call.
    \returns False, if no changes were made to constraintState, i.e. the active
    index remained the same, signifying no updates to the current set of search
    directions are needed.  
    Otherwise the method returns true.

    \see For more on how active constraints are determined and described by
    constraintState see LinConstr::getActiveIneqIndices().
*/
  bool updateConstraintState(double newStep);
  
  //! Updates the step, trueStep, and tag information.
  /*!
    \param newStep Step value given to each new directions.
    \param isAppend If false, direction information will
    be updated for all rows in direction.  Otherwise, direction
    information will be updated only for recently added rows.

    \note When isAppend=true, the first nDirections rows of direction
    are considered old; the remaining rows are considered new additions.
  */
  void updateDirectionInfo(double newStep, bool isAppend=false);
  //@}


  //! The problem definition, containing bound constraints.
  const ProblemDef &  probDef;

  //! The constraints object
  const LinConstr &  linConstr;

  //! The number of dimensions
  int nDimensions;

  //! The vector of all zeros
  const Vector zero;

  //! The step tolerance used to test convergence, etc.
  double stepTolerance;

  //! Minimum size of the step after a new minimum is found
  double minStep;

  //! Contraction parameter
  double theta;

  //! Maximum radius about a current point, for determining active constraints.
  double epsilonMax;

  //! The current number of search directions
  int nDirections;

  //! The search directions are stored as rows of this matrix
  Matrix directionsMatrix;

  //! The steps associated with each search direction
  Vector step;

  //! The actual step associated with each search direction
  /*! The step represents the longest <em>feasible</em> step that is
    less than or equal to the corresponding delta.
  */
  Vector trueStep;

  //! The trial point tags corresponding to the search directions
  /*! A value of -1 means that there is no trial point associated with
    the given direction.
  */
  vector<int> tag;

  //! A temporary vector
  Vector tmpVector;

  //! A temporary integer vector
  mutable vector<int> idxVector;

  //! Stores a cache of previously computed search directions.
  /*!
    A vector of LinConstr::ActiveType values is paired with a corresponding
    set of search directions (stored as a Matrix) in a map to form a
    direction cache.  Each Matrix stored in \e directionCache consists of
    - generators for the tangent cone, and possibly,
    - projected generators for the normal cone (if parameter
    "Add Projected Normals" is true).

    \see LinConstr::getActiveIneqIndices()

    \note caching is only used if general linear constraints are present.
  */
  map< vector< LinConstr::ActiveType >, Matrix>  directionCache;

  //! An iterator to be used in conjunction with the cached directions.
  map< vector< LinConstr::ActiveType>, Matrix>::iterator  cacheIter;
  
  //! Stores a unique description of the tangent cone at the current point.
  /*!
    A call to LinConstr::getActiveIneqIndices()
    returns a vector of type \e enum LinConstr::ActiveType
    which uniquely identifies the tangent cone at a given point.
 */
  vector< LinConstr::ActiveType > constraintState;

  //! Records the last value of epsilonMin used to determine active set.
  double epsilonMin;
  
  //! If true, the projected constraints normals are added to direction. 
  bool withNormals;  

  //! If true, the projected compass search directions are added to direction. 
  bool withCompass;  

  //! Records number of time cached directions are used.
  int nCached;
  
  //! Records number of times directions are generated using LAPACK.
  int nLapack;
  
  //! Records number of times directions are generated using CDDLIB.
  int nCddLib;

  //! Records the maximum number of directions used per iteration.
  int nMaxDirections;

  //! Records the number of times appendNewDirections() is called and nontrivial directions added.
  int nAppend;

  //! Stores distance from current point to each boundary wrt LinConstr::aHat.
  Vector xDistance;

};

}

#endif     //-- HOPSPACK_GSSDIRECTIONS_HPP
