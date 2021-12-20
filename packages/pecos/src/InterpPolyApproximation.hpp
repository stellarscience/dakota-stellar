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
#include "SharedInterpPolyApproxData.hpp"

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

  bool update_active_iterators(const ActiveKey& key);

  int min_coefficients() const;

  void allocate_arrays();

  void combined_to_active(bool clear_combined);

  /// computes component (main and interaction) Sobol' indices
  void compute_component_sobol();
  /// computes total Sobol' indices
  void compute_total_sobol();

  /// compute numerical moments (to order 4 if full_stats)
  void compute_moments(bool full_stats = true, bool combined_stats = false);
  // compute numerical moments in all-variables mode to order 2
  //void compute_moments(const RealVector& x, bool full_stats = true,
  //		         bool combined_stats = false);

  const RealVector& expansion_moments() const;
  const RealVector& numerical_integration_moments() const;

  //
  //- Heading: New virtual functions
  //

  /// update surr_data with synthetic data (from evaluating the interpolant)
  /// in order to ease downstream processing
  virtual void synthetic_surrogate_data(SurrogateData& surr_data) = 0;

  /// compute moments of response using numerical integration
  virtual void integrate_response_moments(size_t num_moments,
					  bool combined_stats) = 0;
  /// compute moments of expansion using numerical integration
  virtual void integrate_expansion_moments(size_t num_moments,
					   bool combined_stats) = 0;

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


inline bool InterpPolyApproximation::
update_active_iterators(const ActiveKey& key)
{
  surrData.active_key(key);
  PolynomialApproximation::update_active_iterators(key);
  return true;
}


/** TO DO: all current cases should use expansion moments (with the
    distinction drawn between interpolation and integration rules, we
    prefer Gauss integration rules on interpolant expansions). */
inline const RealVector& InterpPolyApproximation::expansion_moments() const
{ return secondaryMoments; }


inline const RealVector& InterpPolyApproximation::
numerical_integration_moments() const
{ return primaryMomIter->second; }


inline void InterpPolyApproximation::
compute_moments(bool full_stats, bool combined_stats)
{
  if (full_stats) {
    // std variables mode supports four moments using the collocation rules
    // as integration rules
    integrate_response_moments(4, combined_stats);
    // do this second so that clearing any existing rules does not cause rework
    //if (expConfigOptions.outputLevel >= VERBOSE_OUTPUT)
    integrate_expansion_moments(4, combined_stats);
  }
  else { // only two moments required for incremental metrics
    //integrate_response_moments(2, combined_stats);

    // this approach utilizes bit trackers for computed moments
    PolynomialApproximation::compute_moments(full_stats, combined_stats);
  }
}


/*
inline void InterpPolyApproximation::
compute_moments(const RealVector& x, bool full_stats, bool combined_stats)
{
  PolynomialApproximation::compute_moments(x, full_stats, combined_stats);

  //if (full_stats) integrate_expansion_moments(4, x, combined_stats);
  //else if (!secondaryMoments.empty()) secondaryMoments.resize(0);
  // Note: it would be feasible to implement an all_variables version of
  // integrate_expansion_moments() by evaluating the combined expansion at
  // {design/epistemic=initialPtU,aleatory=Gauss points}
  // --> can't do this for integrate_response_moments() (lacking reqd resp data)
  // --> would require generation of new TPQ/SSG grid only over aleatory vars
  // --> could possibly retire redundant all_vars functions
}
*/

} // namespace Pecos

#endif
