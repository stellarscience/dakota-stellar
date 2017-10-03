/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        KrawtchoukOrthogPolynomial
//- Description:  Class for Krawtchouk Orthogonal Polynomial
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#ifndef KRAWTCHOUK_ORTHOG_POLYNOMIAL_HPP
#define KRAWTCHOUK_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Krawtchouk polynomials

/** The KrawtchoukOrthogPolynomial class evaluates a univariate Krawtchouk
    polynomial K^(p,N)_n of a particular order.  These polynomials
    are orthogonal with respect to the weight function 
    (N choose k)*p^k*(1-p)^(n-k).
    This corresponds to the binomial probability mass function,  
    which is the probability of k successes from N trials.
    See appendix in Xiu & Karniadakis, Siam J. Sci. Comp., v24, n2,
    pp. 619-644, 2002 for more details.  */

class KrawtchoukOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  KrawtchoukOrthogPolynomial();
  /// destructor
  ~KrawtchoukOrthogPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

  //
  //- Heading: Noninherited memeber functions
  //

protected:

  //
  //- Heading: Virtual function redefinitions
  //
  Real type1_value(Real x, unsigned short order);

  /// return alphaPoly
  Real alpha_polynomial() const;
  /// return betaPoly
  Real beta_polynomial() const;

  /// set alphaStat (probability per trial)
  void alpha_stat(Real alpha);
  /// set betaStat (num_trials)
  void beta_stat(Real beta);

private:

  //
  //- Heading: Data
  //

  /// the probability of a "success" for each experiment
  Real alphaPoly;
  /// the number of discrete points on which to base the polynomial 
  short betaPoly;
};


inline KrawtchoukOrthogPolynomial::KrawtchoukOrthogPolynomial() :
  alphaPoly(-1.0), betaPoly(-1.0)
{ }

inline KrawtchoukOrthogPolynomial::~KrawtchoukOrthogPolynomial()
{ }

inline void KrawtchoukOrthogPolynomial::alpha_stat(Real alpha)
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

inline void KrawtchoukOrthogPolynomial::beta_stat(Real beta)
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

inline Real KrawtchoukOrthogPolynomial::alpha_polynomial() const
{ return alphaPoly; }

inline Real KrawtchoukOrthogPolynomial::beta_polynomial() const
{ return betaPoly; }

} // namespace Pecos

#endif
