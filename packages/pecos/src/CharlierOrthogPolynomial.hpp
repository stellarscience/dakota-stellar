/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        CharlierOrthogPolynomial
//- Description:  Class for Charlier Orthogonal Polynomial
//-               
//- Owner:        John Jakeman, Sandia National Laboratories

#ifndef CHARLIER_ORTHOG_POLYNOMIAL_HPP
#define CHARLIER_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"

namespace Pecos {

/**
 * \class CharlierOrthogPolynomial
 * \brief One-dimensional Charlier polynomial
 */
class CharlierOrthogPolynomial : public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //
  
  /// default constructor
  CharlierOrthogPolynomial();
  /// destructor
  ~CharlierOrthogPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

protected:

  Real type1_value( Real x, unsigned short order );
  Real type1_gradient( Real x, unsigned short order );
  Real type1_hessian( Real x, unsigned short order );
  Real norm_squared( unsigned short order );

  /// return alphaPoly
  Real alpha_polynomial() const;
  /// set betaPoly
  void alpha_stat(Real alpha);

private: 
  
  /// Poisson distributioon is the probability that a realziations of 
  /// a random variable X with mean alpha_stat occurring k times in a 
  /// fixed interval.
  /// expected value of the random variable X
  Real alphaPoly;

};

inline CharlierOrthogPolynomial::CharlierOrthogPolynomial()
{};

inline CharlierOrthogPolynomial::~CharlierOrthogPolynomial()
{};

inline Real CharlierOrthogPolynomial::alpha_polynomial() const
{ return alphaPoly; }

inline void CharlierOrthogPolynomial::alpha_stat(Real alpha)
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

#endif // CHARLIER_ORTHOG_POLYNOMIAL_HPP

} // namespace Pecos
