/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        InterpPolyApproximation
//- Description:  Class for Lagrange Interpolation Polynomial Approximation
//-               
//- Owner:        Mike Eldred

#ifndef INTERP_POLY_APPROXIMATION_HPP
#define INTERP_POLY_APPROXIMATION_HPP

#include "PolynomialApproximation.hpp"
#include "BasisPolynomial.hpp"

namespace Pecos {


/// Derived approximation class for interpolation polynomials (global
/// approximation).

/** The InterpPolyApproximation class provides a global approximation
    based on interpolation polynomials.  It is used primarily for
    stochastic collocation approaches to uncertainty quantification. */

class InterpPolyApproximation: public PolynomialApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  InterpPolyApproximation(const SharedBasisApproxData& shared_data);
  /// destructor
  ~InterpPolyApproximation();

  //
  //- Heading: member functions
  //

  //
  //- Heading: Virtual function redefinitions
  //

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  int min_coefficients() const;

  void allocate_arrays();

  /// computes component (main and interaction) Sobol' indices
  void compute_component_sobol();
  /// computes total Sobol' indices
  void compute_total_sobol();

  /// compute numerical moments to order 4
  void compute_moments(bool full_stats = true);
  /// compute numerical moments in all-variables mode to order 2
  void compute_moments(const RealVector& x, bool full_stats = true);

  //
  //- Heading: New virtual functions
  //

  /// compute moments of response using numerical integration
  virtual void integrate_response_moments(size_t num_moments) = 0;
  /// compute moments of expansion using numerical integration
  virtual void integrate_expansion_moments(size_t num_moments) = 0;

  virtual void compute_total_sobol_indices() = 0;
  virtual void compute_partial_variance(const BitArray& set_value);

  //
  //- Heading: Convenience functions
  //

  /// test accuracy of the interpolants
  void test_interpolation();

  //
  //- Heading: Data
  //

  /// the partialVariances of subset functions f_alpha
  RealVector partialVariance;

private:

  //
  //- Heading: Convenience functions
  //

  /// recursively identifies constituent subsets that are children of
  /// a parent set
  void proper_subsets(const BitArray& parent_set, BitArraySet& children);

  //
  //- Heading: Data
  //
};


inline InterpPolyApproximation::
InterpPolyApproximation(const SharedBasisApproxData& shared_data):
  PolynomialApproximation(shared_data)
{ }


inline InterpPolyApproximation::~InterpPolyApproximation()
{ }


inline void InterpPolyApproximation::compute_moments(bool full_stats)
{
  // standard variables mode supports four moments using the collocation rules
  // as integration rules
  integrate_response_moments(4);

  // do this second so that clearing any existing rules does not cause rework
  if (full_stats) //&& expConfigOptions.outputLevel >= VERBOSE_OUTPUT)
    integrate_expansion_moments(4);
  else
    expansionMoments.resize(0);
}


inline void InterpPolyApproximation::
compute_moments(const RealVector& x, bool full_stats)
{
  // all variables mode only supports first two moments
  mean(x); variance(x);
  //standardize_moments(numericalMoments);

  //if (full_stats) integrate_expansion_moments(4, x);
  //else            expansionMoments.resize(0);
  // Note: it would be feasible to implement an all_variables version of
  // integrate_expansion_moments() by evaluating the combined expansion at
  // {design/epistemic=initialPtU,aleatory=Gauss points}
  // --> can't do this for integrate_response_moments() (lacking reqd resp data)
  // --> would require generation of new TPQ/SSG grid only over aleatory vars
  // --> could possibly retire redundant all_vars functions
}

} // namespace Pecos

#endif
