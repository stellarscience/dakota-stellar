/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LRDHierarchInterpPolyApproximation
//- Description:  Class for Nodal Interpolation Polynomial Approximation
//-               
//- Owner:        Chris Miller

#ifndef HIERARCH_INTERP_POLY_APPROXIMATION_HPP
#define HIERARCH_INTERP_POLY_APPROXIMATION_HPP

#include "InterpPolyApproximation.hpp"
#include "LocalRefinableDriver.hpp"

namespace Pecos {

/// Derived approximation class for piecewise linear and cubic hierarchical 
/// interpolation polynomials (local approximation interpolating function 
/// values and potentially gradients at collocation points).

/** The LRDHierarchInterpPolyApproximation class provides a local piecewise 
    polynomial approximation based on the hierarchical approach described in 
    X. Ma and N. Zabaras "An adaptive hierarchical sparse grid collocation 
    algorithm for the solution of stochastic differential equations", Journal 
    of Computational Physics, 228 (2009), 3084-3113.  Both piecewise linear 
    basis functions using function values at the collocation points and cubic 
    Hermite basis functions using both values and derivatives are available.  
    It is used primarily for stochastic collocation approaches to uncertainty 
    quantification. */

class LRDHierarchInterpPolyApproximation: public InterpPolyApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// Default constructor
  /** @param basis_type Has no effect on the derived class.
      @param num_vars The number of variables that define the problem.
      @param use_derivs If true the interpolant consists of piecewise cubic
        Hermite polynomials which interpolate the function values and gradients.
        If false the interpolant is piecewise linear and only interpolates
        function values. */ 
  LRDHierarchInterpPolyApproximation(short basis_type, size_t num_vars,
				     bool use_derivs);
  /// destructor
  ~LRDHierarchInterpPolyApproximation();

  //
  //- Heading: <out of date>
  //

  /// Update the coefficients using new data.
  void increment_coefficients();

  //
  //- Heading: Virtual function redefinitions
  //

  void allocate_expansion_coefficients();
  void compute_expansion_coefficients();
  void restore_expansion_coefficients();
  void compute_numerical_response_moments(size_t num_moments);
  void compute_numerical_expansion_moments(size_t num_moments);

  void compute_total_sobol_indices();
  void member_coefficients_weights(int set_value,
    const UShortArray& quad_order, const UShortArray& lev_index,
    const UShort2DArray& key, const SizetArray& colloc_index,
    RealVector& member_coeffs, RealVector& member_wts);

  /// retrieve the response expansion value for a given parameter vector.
  /** @param x A RealVector of size numVars.  The interpolant is evaluated at
        this point.
      @return The value of the interpolant at x.
  */
  Real value(const RealVector& x);
  
  /// retrieve the response expansion gradient for a given parameter vector
  /// and default DVV
  /** @param x A RealVector of size numVars.  The gradient of the interpolant
        is evaluated at this point.
     @return A const reference to approxGradient.  This will be a RealVector
       of length numVars containing the gradient of the interpolant at x.
  */
  const RealVector& gradient_basis_variables(const RealVector& x);
  
  /// retrieve the response expansion gradient for a given parameter vector
  /// and given DVV
  /** @param x A RealVector of size numVars.  The gradient of the interpolatnt
        with respect to the variables indicated in the array dvv is evaluated
        at this point.
      @param dvv A SizetArray of size at least 1 and at most numVars.  Each
        entry in dvv indicates a direction to take the gradient in.
      @return A const reference to approxGradient.  This will be a RealVector
        of length dvv.size() containing the gradient of the interpolant at x
        with respect to the variables indicated in the dvv input parameter.
  */
  const RealVector& gradient_basis_variables(const RealVector& x,
					     const SizetArray& dvv);

  const RealVector& gradient_nonbasis_variables(const RealVector& x);

  Real stored_value(const RealVector& x);
  const RealVector& stored_gradient_basis_variables(const RealVector& x);
  const RealVector& stored_gradient_nonbasis_variables(const RealVector& x);

  /// Returns the mean of the expansion, treating all variables as random.
  /** @return The mean of the expansion.
   */
  Real mean();

  /** @brief Returns the mean of the expansion, treating a subset of the 
      variables as random.
      @param x The non-random coordinates
      @return The mean of the expansion with respect to the random variables.
  */
  Real mean(const RealVector& x);
  
  const RealVector& mean_gradient();
  const RealVector& mean_gradient(const RealVector& x,
				  const SizetArray& dvv);

  Real variance();
  Real variance(const RealVector& x);
  const RealVector& variance_gradient();
  const RealVector& variance_gradient(const RealVector& x,
				      const SizetArray& dvv);

  Real covariance(PolynomialApproximation* poly_approx_2);
  Real covariance(const RealVector& x,
		      PolynomialApproximation* poly_approx_2);

protected:

  /// returns an int array containing the indices of the points whose
  /// support includes x
  const IntArray& in_support_of(const RealVector& x);

  /// compute the value at a point using a lower level than the full
  /// approximation.
  Real value(const RealVector& x, unsigned int max_level);

  /// compute the approximate gradient at a point using a lower level
  /// than the full approximation.
  const RealVector& gradient_basis_variables(const RealVector& x,
					     unsigned int max_level);

  /// The largest computed coefficient.
  unsigned int maxComputedCoeff;

private:
  
  ///Pecos:PIECEWISE_INTERP_POLYNOMIAL or Pecos:PIECEWISE_CUBIC_INTERP
  short polyType;

  //Array of ints indicating which basis functions have support
  //containing a given point.
  IntArray supportIndicator;

  /// the type1 coefficients of the expansion for interpolating values
  RealVector expansionType1Coeffs;
  /// the type2 coefficients of the expansion for interpolating gradients
  RealMatrix expansionType2Coeffs;
  /// the gradients of the type1 expansion coefficients
  /** may be interpreted as either the gradients of the expansion
      coefficients or the coefficients of expansions for the response
      gradients.  This array is used when sensitivities of moments are
      needed with respect to variables that do not appear in the
      expansion (e.g., with respect to design variables for an
      expansion only over the random variables). */
  RealMatrix expansionType1CoeffGrads;
};


} // namespace Pecos

#endif
