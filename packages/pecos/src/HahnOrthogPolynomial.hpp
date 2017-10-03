/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HahnOrthogPolynomial
//- Description:  Class for Hahn Orthogonal Polynomial
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#ifndef HAHN_ORTHOG_POLYNOMIAL_HPP
#define HAHN_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Hahn polynomials

/** The HahnOrthogPolynomial class evaluates a univariate Hahn
    polynomial Q^(alpha,beta,N)_n of a particular order.  These polynomials
    are orthogonal with respect to the weight function 

    (K choose k)(N-K choose n-k)/( N choose n).

    This corresponds to the hypergeometric probability mass function, which
    describes the probability of k successes in n draws, without 
    replacement, from a finite population of size N that contains exactly 
    K successes.
    See appendix in Xiu & Karniadakis, Siam J. Sci. Comp., v24, n2,
    pp. 619-644, 2002 for more details.  */

class HahnOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  HahnOrthogPolynomial();
  /// destructor
  ~HahnOrthogPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

  //
  //- Heading: Noninherited memeber functions
  //

  //  These will need to be rolled into the API - RWH
  /// return gammaPoly
  Real gamma_polynomial() const;
  /// set gammaStat (total population)
  void gamma_stat(Real gamma);

protected:

  //
  //- Heading: Virtual function redefinitions
  //
  Real type1_value(Real x, unsigned short order);

  /// return alphaPoly
  Real alpha_polynomial() const;
  /// return betaPoly
  Real beta_polynomial() const;

  /// set alphaStat (selected population)
  void alpha_stat(Real alpha);
  /// set betaStat (num draws)
  void beta_stat(Real beta);

private:

  //
  //- Heading: Data
  //

  /// the hypergeometric alpha parameter
  Real alphaPoly;
  /// the hypergeometric beta parameter
  Real betaPoly;
  /// the number of discrete points on which to base the polynomial 
  Real gammaPoly;
};


inline HahnOrthogPolynomial::HahnOrthogPolynomial() :
  alphaPoly(-1.0), betaPoly(-1.0), gammaPoly(-1.0)
{ }

inline HahnOrthogPolynomial::~HahnOrthogPolynomial()
{ }

inline void HahnOrthogPolynomial::alpha_stat(Real alpha)
{
  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Therefore, set parametricUpdate to false unless an actual parameter change.
  // Logic for first pass included for completeness, but should not be needed.
  if (collocPoints.empty() || collocWeights.empty()) { // first pass
    parametricUpdate = true; // prevent false if default value assigned
    alphaPoly = alpha;
  }
  else {
    parametricUpdate = false;
    Real ap = alpha;
    if (!real_compare(alphaPoly, ap))
      { alphaPoly = ap; parametricUpdate = true; reset_gauss(); }
  }
}

inline void HahnOrthogPolynomial::beta_stat(Real beta)
{
  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Therefore, set parametricUpdate to false unless an actual parameter change.
  // Logic for first pass included for completeness, but should not be needed.
  if (collocPoints.empty() || collocWeights.empty()) { // first pass
    parametricUpdate = true; // prevent false if default value assigned
    betaPoly = beta;
  }
  else {
    parametricUpdate = false;
    Real bp = beta;
    if (!real_compare(betaPoly, bp))
      { betaPoly = bp; parametricUpdate = true; reset_gauss(); }
  }
}

inline void HahnOrthogPolynomial::gamma_stat(Real gamma)
{
  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Therefore, set parametricUpdate to false unless an actual parameter change.
  // Logic for first pass included for completeness, but should not be needed.
  if (collocPoints.empty() || collocWeights.empty()) { // first pass
    parametricUpdate = true; // prevent false if default value assigned
    gammaPoly = gamma;
  }
  else {
    parametricUpdate = false;
    Real bp = gamma;
    if (!real_compare(gammaPoly, bp))
      { gammaPoly = bp; parametricUpdate = true; reset_gauss(); }
  }
}

inline Real HahnOrthogPolynomial::alpha_polynomial() const
{ return alphaPoly; }

inline Real HahnOrthogPolynomial::beta_polynomial() const
{ return betaPoly; }

inline Real HahnOrthogPolynomial::gamma_polynomial() const
{ return gammaPoly; }

} // namespace Pecos

#endif
