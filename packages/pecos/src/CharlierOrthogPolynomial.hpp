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

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  Real type1_value( Real x, unsigned short order ) override;
  Real type1_gradient( Real x, unsigned short order ) override;
  Real type1_hessian( Real x, unsigned short order ) override;
  Real norm_squared( unsigned short order ) override;

  void pull_parameter(short dist_param, Real& param) const override;
  void push_parameter(short dist_param, Real  param) override;
  bool parameterized() const override;

private: 
  
  /// Poisson distribution is the probability that a realization of a random
  /// variable X with mean lambda occurring k times in a fixed interval.

  /// expected value of the random variable X
  Real lambdaStat;
};


inline CharlierOrthogPolynomial::CharlierOrthogPolynomial():
  lambdaStat(0.) // dummy value prior to update
{ }


inline CharlierOrthogPolynomial::~CharlierOrthogPolynomial()
{ }


inline void CharlierOrthogPolynomial::
pull_parameter(short dist_param, Real& param) const
{
  switch (dist_param) {
  case P_LAMBDA: param = lambdaStat; break;
  default:
    PCerr << "Error: unsupported distribution parameter in CharlierOrthog"
	  << "Polynomial::pull_parameter()." << std::endl;
    abort_handler(-1);
  }
}


inline void CharlierOrthogPolynomial::
push_parameter(short dist_param, Real param)
{
  if (dist_param != P_LAMBDA) {
    PCerr << "Error: unsupported distribution parameter in CharlierOrthog"
	  << "Polynomial::push_parameter()." << std::endl;
    abort_handler(-1);
  }

  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Logic for first pass included for completeness, but should not be needed.
  if (collocPointsMap.empty() || collocWeightsMap.empty()) // first pass
    lambdaStat = param;
  else if (!real_compare(lambdaStat, param))
    { lambdaStat = param;  reset_gauss(); }
}


inline bool CharlierOrthogPolynomial::parameterized() const
{ return true; }

} // namespace Pecos

#endif // CHARLIER_ORTHOG_POLYNOMIAL_HPP
