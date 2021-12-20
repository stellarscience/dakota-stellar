/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        GenLaguerreOrthogPolynomial
//- Description:  Class for GenLaguerre Orthogonal Polynomial
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#ifndef GEN_LAGUERRE_ORTHOG_POLYNOMIAL_HPP
#define GEN_LAGUERRE_ORTHOG_POLYNOMIAL_HPP

#include "OrthogonalPolynomial.hpp"
#include "GammaRandomVariable.hpp"


namespace Pecos {

/// Derived orthogonal polynomial class for generalized Laguerre polynomials

/** The GenLaguerreOrthogPolynomial class evaluates a univariate
    generalized/associated Laguerre polynomial L^(alpha)_n of a
    particular order.  These polynomials are orthogonal with respect
    to the weight function x^alpha exp(-x) when integrated over the
    support range of [0,+infinity].  This corresponds to the
    probability density function f(x) = x^alpha exp(-x) / Gamma(alpha+1)
    for the standard gamma distribution, although common statistical
    PDF parameter conventions (see, e.g., the uncertain variables
    section in the DAKOTA Reference Manual) and the Abramowitz and
    Stegun orthogonal polynomial parameter conventions require an
    offset conversion in this case (alpha_poly = alpha_stat - 1 with
    the poly definition used in both cases above).  It enables (mixed)
    multidimensional orthogonal polynomial basis functions within
    OrthogPolyApproximation.  A special case is the
    LaguerreOrthogPolynomial (implemented separately), for which
    alpha_poly = 0 and weight function = exp(-x) (the standard
    exponential distribution). */

class GenLaguerreOrthogPolynomial: public OrthogonalPolynomial
{
public:

  //
  //- Heading: Constructor and destructor
  //

  GenLaguerreOrthogPolynomial();                ///< default constructor
  GenLaguerreOrthogPolynomial(Real alpha_stat); ///< standard constructor
  ~GenLaguerreOrthogPolynomial();               ///< destructor

  //
  //- Heading: Virtual function redefinitions
  //

  /// calculate and return wtFactor based on alphaPoly
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

  void pull_parameter(short dist_param, Real& param) const;
  void push_parameter(short dist_param, Real  param);
  bool parameterized() const;

  Real length_scale() const;

private:

  //
  //- Heading: Data
  //

  /// the alpha parameter for the generalized Laguerre polynomial as defined
  /// by Abramowitz and Stegun (differs from statistical PDF notation)
  Real alphaPoly;
};


inline GenLaguerreOrthogPolynomial::GenLaguerreOrthogPolynomial(): alphaPoly(0.)
{ collocRule = GEN_GAUSS_LAGUERRE; }


inline GenLaguerreOrthogPolynomial::
GenLaguerreOrthogPolynomial(Real alpha_stat): alphaPoly(alpha_stat-1.)
{ collocRule = GEN_GAUSS_LAGUERRE; }


inline GenLaguerreOrthogPolynomial::~GenLaguerreOrthogPolynomial()
{ }


inline void GenLaguerreOrthogPolynomial::
pull_parameter(short dist_param, Real& param) const
{
  // return BE_ALPHA,BETA (unconverted)
  switch (dist_param) {
  case     GA_ALPHA: param = alphaPoly + 1.; break; // alpha poly to alpha stat
  case GENLAG_ALPHA: param = alphaPoly;      break;
  default:
    PCerr << "Error: unsupported distribution parameter in GenLaguerreOrthog"
	  << "Polynomial::parameter()." << std::endl;
    abort_handler(-1); break;
  }
}


inline void GenLaguerreOrthogPolynomial::
push_parameter(short dist_param, Real param)
{
  // *_stat() routines are called for each approximation build from
  // PolynomialApproximation::update_basis_distribution_parameters().
  // Logic for first pass included for completeness, but should not be needed.
  if (collocPointsMap.empty() || collocWeightsMap.empty()) { // first pass
    switch (dist_param) {
    case     GA_ALPHA: alphaPoly = param - 1.; break;// alpha stat to alpha poly
    case GENLAG_ALPHA: alphaPoly = param;      break;
    }
  }
  else
    switch (dist_param) {
    case GA_ALPHA: {
      Real ap = param - 1.; // alpha stat to beta poly
      if (!real_compare(alphaPoly, ap))
	{ alphaPoly = ap;  reset_gauss(); }
      break;
    }
    case GENLAG_ALPHA:
      if (!real_compare(alphaPoly, param))
	{ alphaPoly = param;  reset_gauss(); }
      break;
    }
}


inline bool GenLaguerreOrthogPolynomial::parameterized() const
{ return true; }


/** return max(mean,stdev) */
inline Real GenLaguerreOrthogPolynomial::length_scale() const
{
  Real mean, stdev;
  // GammaRandomVariable uses alpha_stat, beta_stat
  GammaRandomVariable::moments_from_params(alphaPoly + 1., 1., mean, stdev);
  return std::max(mean, stdev);

  //Real alpha_stat = alphaPoly + 1.;
  //return std::max(alpha_stat, std::sqrt(alpha_stat)); // too opaque
}

} // namespace Pecos

#endif
