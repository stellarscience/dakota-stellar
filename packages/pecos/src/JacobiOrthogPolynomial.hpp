/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        JacobiOrthogPolynomial
//- Description:  Class for Jacobi Orthogonal Polynomial
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef JACOBI_ORTHOG_POLYNOMIAL_HPP
#define JACOBI_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"
#include "BetaRandomVariable.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Jacobi polynomials

/** The JacobiOrthogPolynomial class evaluates a univariate Jacobi
    polynomial P^(alpha,beta)_n of a particular order.  These
    polynomials are orthogonal with respect to the weight function
    (1-x)^alpha (1+x)^beta when integrated over the support range of
    [-1,+1].  This corresponds to the probability density function
    f(x) = (1-x)^alpha (1+x)^beta / (2^(alpha+beta+1) B(alpha+1,beta+1))
    for the beta distribution for [L,U]=[-1,1], where common
    statistical PDF notation conventions (see, e.g., the uncertain
    variables section in the DAKOTA Reference Manual) and the
    Abramowiz and Stegun orthogonal polynomial conventions are
    inverted and require conversion in this case (alpha_poly =
    beta_stat - 1; beta_poly = alpha_stat - 1 with the poly
    definitions used in both cases above).  It enables (mixed)
    multidimensional orthogonal polynomial basis functions within
    OrthogPolyApproximation.  A special case is the
    LegendreOrthogPolynomial (implemented separately), for which
    alpha_poly = beta_poly = 0. */

class JacobiOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  JacobiOrthogPolynomial();
  /// standard constructor
  JacobiOrthogPolynomial(Real alpha_stat, Real beta_stat);
  /// destructor
  ~JacobiOrthogPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

  /// calculate and return wtFactor based on alphaPoly and betaPoly
  Real weight_factor();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  Real type1_value(Real x, unsigned short order);
  Real type1_gradient(Real x, unsigned short order);
  Real type1_hessian(Real x, unsigned short order);
  Real norm_squared(unsigned short order);

  const RealArray& collocation_points(unsigned short order);
  const RealArray& type1_collocation_weights(unsigned short order);

  /// return alphaPoly
  Real alpha_polynomial() const;
  /// return betaPoly
  Real beta_polynomial() const;
  /// set betaPoly using the conversion betaPoly = alpha_stat - 1.
  void alpha_stat(Real alpha);
  /// set alphaPoly using the conversion alphaPoly = beta_stat - 1.
  void beta_stat(Real beta);

  /// override default definition (false) since Jacobi is parameterized
  bool parameterized() const;

  Real length_scale() const;

private:

  //
  //- Heading: Data
  //

  /// the alpha parameter for the Jacobi polynomial as defined by
  /// Abramowitz and Stegun (differs from statistical PDF notation)
  Real alphaPoly;
  /// the beta parameter for the Jacobi polynomial as defined by
  /// Abramowitz and Stegun (differs from statistical PDF notation)
  Real betaPoly;
};


inline JacobiOrthogPolynomial::JacobiOrthogPolynomial():
  alphaPoly(0.), betaPoly(0.)
{ collocRule = GAUSS_JACOBI; }


// TO DO
inline JacobiOrthogPolynomial::
JacobiOrthogPolynomial(Real alpha_stat, Real beta_stat):
  alphaPoly(beta_stat-1.), betaPoly(alpha_stat-1.) // inverted conventions
{ collocRule = GAUSS_JACOBI; }


inline JacobiOrthogPolynomial::~JacobiOrthogPolynomial()
{ }


inline Real JacobiOrthogPolynomial::alpha_polynomial() const
{ return alphaPoly; }


inline Real JacobiOrthogPolynomial::beta_polynomial() const
{ return betaPoly; }


inline void JacobiOrthogPolynomial::alpha_stat(Real alpha)
{
  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Therefore, set parametricUpdate to false unless an actual parameter change.
  // Logic for first pass included for completeness, but should not be needed.
  if (collocPoints.empty() || collocWeights.empty()) { // first pass
    parametricUpdate = true; // prevent false if default value assigned
    betaPoly = alpha - 1.;
  }
  else {
    parametricUpdate = false;
    Real bp = alpha - 1.;
    if (!real_compare(betaPoly, bp))
      { betaPoly = bp; parametricUpdate = true; reset_gauss(); }
  }
}


inline void JacobiOrthogPolynomial::beta_stat(Real beta)
{
  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Therefore, set parametricUpdate to false unless an actual parameter change.
  // Logic for first pass included for completeness, but should not be .
  if (collocPoints.empty() || collocWeights.empty()) { // first pass
    parametricUpdate = true; // prevent false if default value assigned
    alphaPoly = beta - 1.;
  }
  else {
    parametricUpdate = false;
    Real ap = beta - 1.;
    if (!real_compare(alphaPoly, ap))
      { alphaPoly = ap; parametricUpdate = true; reset_gauss(); }
  }
}


inline bool JacobiOrthogPolynomial::parameterized() const
{ return true; }


/** return max(mean, stdev) on [-1,1]. */
inline Real JacobiOrthogPolynomial::length_scale() const
{
  Real mean, stdev;
  // BetaRandomVariable uses alpha_stat, beta_stat:
  BetaRandomVariable::
    moments_from_params(betaPoly+1., alphaPoly+1., -1., 1., mean, stdev);
  return std::max(mean, stdev);
}

} // namespace Pecos

#endif
