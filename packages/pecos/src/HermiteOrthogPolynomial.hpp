/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HermiteOrthogPolynomial
//- Description:  Class for Hermite Orthogonal Polynomial
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef HERMITE_ORTHOG_POLYNOMIAL_HPP
#define HERMITE_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"
#include "pecos_global_defs.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for Hermite polynomials

/** The HermiteOrthogPolynomial class evaluates a univariate Hermite
    polynomial of a particular order.  It uses the "probabilist's"
    formulation for which the polynomials are orthogonal with respect
    to the weight function 1/std::sqrt(2*PI) exp(-x^2/2) when integrated
    over the support range of [-infinity,+infinity].  It enables
    (mixed) multidimensional orthogonal polynomial basis functions
    within OrthogPolyApproximation. */

class HermiteOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  HermiteOrthogPolynomial(short colloc_rule); ///< extended constructor
  HermiteOrthogPolynomial();                  ///< default constructor
  ~HermiteOrthogPolynomial();                 ///< destructor

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


inline HermiteOrthogPolynomial::HermiteOrthogPolynomial(short colloc_rule)
{
  collocRule = colloc_rule;
  ptFactor   = std::sqrt(2.);
  wtFactor   = 1./std::sqrt(PI);
}


inline HermiteOrthogPolynomial::HermiteOrthogPolynomial()
{
  collocRule = GAUSS_HERMITE;
  ptFactor   = std::sqrt(2.);
  wtFactor   = 1./std::sqrt(PI);
}


inline HermiteOrthogPolynomial::~HermiteOrthogPolynomial()
{ }


/** mean is zero; return std deviation = 1. */
inline Real HermiteOrthogPolynomial::length_scale() const
{ return 1.; }

} // namespace Pecos

#endif
