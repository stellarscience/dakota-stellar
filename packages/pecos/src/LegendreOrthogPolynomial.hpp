/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LegendreOrthogPolynomial
//- Description:  Class for Legendre Orthogonal Polynomial
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef LEGENDRE_ORTHOG_POLYNOMIAL_HPP
#define LEGENDRE_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Legendre polynomials

/** The LegendreOrthogPolynomial class evaluates a univariate Legendre
    polynomial of a particular order.  These polynomials are
    orthogonal with respect to the weight function 1 when integrated
    over the support range of [-1,+1].  This corresponds to the
    probability density function f(x) = 1/(U-L) = 1/2 for the uniform
    distribution for [L,U]=[-1,1].  It enables (mixed)
    multidimensional orthogonal polynomial basis functions within
    OrthogPolyApproximation.  Legendre polynomials are a special case
    (alpha = beta = 0) of the more general Jacobi polynomials
    (implemented separately) which correspond to the beta distribution. */

class LegendreOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  LegendreOrthogPolynomial(short colloc_rule); ///< extended constructor
  LegendreOrthogPolynomial();                  ///< default constructor
  ~LegendreOrthogPolynomial();                 ///< destructor

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

  Real length_scale() const;

private:

  //
  //- Heading: Data
  //
};


// collocRule may be GAUSS_LEGENDRE (default), GAUSS_PATTERSON,
// CLENSHAW_CURTIS, or FEJER2
inline LegendreOrthogPolynomial::LegendreOrthogPolynomial(short colloc_rule)
{ collocRule = colloc_rule;    wtFactor = 0.5; }


inline LegendreOrthogPolynomial::LegendreOrthogPolynomial()
{ collocRule = GAUSS_LEGENDRE; wtFactor = 0.5; }


inline LegendreOrthogPolynomial::~LegendreOrthogPolynomial()
{ }


/** [-1,1]: mean is zero; return std deviation = 2/sqrt(12). */
inline Real LegendreOrthogPolynomial::length_scale() const
{ return std::pow(3., -0.5); }

} // namespace Pecos

#endif
