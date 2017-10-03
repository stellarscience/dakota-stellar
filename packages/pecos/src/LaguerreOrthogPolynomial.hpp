/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LaguerreOrthogPolynomial
//- Description:  Class for Laguerre Orthogonal Polynomial
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef LAGUERRE_ORTHOG_POLYNOMIAL_HPP
#define LAGUERRE_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Laguerre polynomials

/** The LaguerreOrthogPolynomial class evaluates a univariate Laguerre
    polynomial of a particular order.  These polynomials are
    orthogonal with respect to the weight function exp(-x) when
    integrated over the support range of [0,+infinity].  This
    corresponds to the probability density function for the standard
    exponential distribution.  It enables (mixed) multidimensional
    orthogonal polynomial basis functions within
    OrthogPolyApproximation.  Laguerre polynomials are a special case
    (alpha = 0) of the generalized Laguerre polynomials (implemented
    separately) which correspond to the standard gamma distribution. */

class LaguerreOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  LaguerreOrthogPolynomial();  ///< default constructor
  ~LaguerreOrthogPolynomial(); ///< destructor

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


inline LaguerreOrthogPolynomial::LaguerreOrthogPolynomial()
{ collocRule = GAUSS_LAGUERRE; }


inline LaguerreOrthogPolynomial::~LaguerreOrthogPolynomial()
{ }


/** return mean value */
inline Real LaguerreOrthogPolynomial::length_scale() const
{ return 1.; } // mean = std dev = 1

} // namespace Pecos

#endif
