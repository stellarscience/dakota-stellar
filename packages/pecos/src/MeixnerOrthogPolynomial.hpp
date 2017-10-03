/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        MeixnerOrthogPolynomial
//- Description:  Class for Meixner Orthogonal Polynomial
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#ifndef MEIXNER_ORTHOG_POLYNOMIAL_HPP
#define MEIXNER_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Meixner polynomials

/** The MeixnerOrthogPolynomial class evaluates a univariate Meixner
    polynomial M^(c,Beta)_n of a particular order.  These polynomials
    are orthogonal with respect to the weight function 

    (k+r-1 choose k)*p^k*(1-p)^n, where p is the the probability of success
    of a trial, n is the number of total successes and k is the trial number.
    The geometric distribution (n==1) is a special case of the negative 
    binomial distribution.
    This corresponds to the negative binomial probability mass function.
    See appendix in Xiu & Karniadakis, Siam J. Sci. Comp., v24, n2,
    pp. 619-644, 2002 for more details.  */

class MeixnerOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  MeixnerOrthogPolynomial();
  /// destructor
  ~MeixnerOrthogPolynomial();

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
  /// set betaStat (num trials)
  void beta_stat(Real beta);

private:

  //
  //- Heading: Data
  //

  /// the probability of a "success" for each experiment
  Real alphaPoly;
  /// the number of failures allowed
  Real betaPoly;
};


inline MeixnerOrthogPolynomial::MeixnerOrthogPolynomial() :
  alphaPoly(-1.0), betaPoly(-1.0)
{ }

inline MeixnerOrthogPolynomial::~MeixnerOrthogPolynomial()
{ }

inline void MeixnerOrthogPolynomial::alpha_stat(Real alpha)
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

inline void MeixnerOrthogPolynomial::beta_stat(Real beta)
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

inline Real MeixnerOrthogPolynomial::alpha_polynomial() const
{ return alphaPoly; }

inline Real MeixnerOrthogPolynomial::beta_polynomial() const
{ return betaPoly; }

} // namespace Pecos

#endif
