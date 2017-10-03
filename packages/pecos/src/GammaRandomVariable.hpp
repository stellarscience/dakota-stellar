/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 GammaRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef GAMMA_RANDOM_VARIABLE_HPP
#define GAMMA_RANDOM_VARIABLE_HPP

#include "ExponentialRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for gamma random variables.

/** Manages alpha and inherits beta.  This follows the
    Gamma(alpha,beta) definition in Law & Kelton, which differs from
    the LHS definition (beta_LK = 1/beta_LHS). */

class GammaRandomVariable: public ExponentialRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  GammaRandomVariable();
  /// alternate constructor
  GammaRandomVariable(Real alpha, Real beta);
  /// destructor
  ~GammaRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  Real pdf_hessian(Real x) const;
  Real log_pdf(Real x) const;
  Real log_pdf_gradient(Real x) const;
  Real log_pdf_hessian(Real x) const;

  Real inverse_standard_cdf(Real p_cdf) const;

  Real standard_pdf(Real z) const;
  Real log_standard_pdf(Real z) const;
  Real log_standard_pdf_gradient(Real z) const;
  Real log_standard_pdf_hessian(Real z) const;

  // inherited from ExponentialRandomVariable
  //Real to_standard(Real x) const;
  //Real from_standard(Real z) const;

  Real parameter(short dist_param) const;
  void parameter(short dist_param, Real val);

  Real mean() const;
  Real median() const;
  Real mode() const;
  Real standard_deviation() const;
  Real variance() const;

  Real correlation_warping_factor(const RandomVariable& rv, Real corr) const;

  Real dx_ds(short dist_param, short u_type, Real x, Real z) const;
  Real dz_ds_factor(short u_type, Real x, Real z) const;

  //
  //- Heading: Member functions
  //

  void update(Real alpha, Real beta);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real alpha, Real beta);
  static Real cdf(Real x, Real alpha, Real beta);
  static Real inverse_cdf(Real cdf, Real alpha, Real beta);

  static void moments_from_params(Real alpha, Real beta,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Member functions
  //

  /// create a new gammaDist instance
  void update_boost();

  //
  //- Heading: Data
  //

  /// alpha shape parameter of gamma random variable (statistical PDF
  /// convention; differs from GenLaguerre polynomial convention)
  Real alphaShape;

  /// pointer to the Boost gamma_distribution instance
  gamma_dist* gammaDist;
};


inline GammaRandomVariable::GammaRandomVariable():
  ExponentialRandomVariable(), alphaShape(1.), gammaDist(NULL)
{ ranVarType = GAMMA; }


inline GammaRandomVariable::GammaRandomVariable(Real alpha, Real beta):
  ExponentialRandomVariable(beta), alphaShape(alpha),
  gammaDist(new gamma_dist(alphaShape, betaScale))
{ ranVarType = GAMMA; }


inline GammaRandomVariable::~GammaRandomVariable()
{ if (gammaDist) delete gammaDist; }


inline Real GammaRandomVariable::cdf(Real x) const
{ return bmth::cdf(*gammaDist, x); }


inline Real GammaRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*gammaDist, x)); }


inline Real GammaRandomVariable::inverse_cdf(Real p_cdf) const
{ return bmth::quantile(*gammaDist, p_cdf); }


inline Real GammaRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*gammaDist, p_ccdf)); }


//  F(x) = Boost
//  f(x) = x^(alpha-1) e^(-x/beta) / GammaFn(alpha) / beta^alpha
// f'(x) = [e^(-x/beta) (alpha-1) x^(alpha-2) - x^(alpha-1) e^(-x/beta)/beta]
//       / GammaFn(alpha) / beta^alpha
//       = f(x) [ (alpha-1)/x - 1/beta ]
inline Real GammaRandomVariable::pdf(Real x) const
{ return bmth::pdf(*gammaDist, x); }


inline Real GammaRandomVariable::pdf_gradient(Real x) const
{
  if (x <= 0.) {
    if      (alphaShape < 1.) return -std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return  std::numeric_limits<Real>::quiet_NaN();
    else return -ExponentialRandomVariable::pdf(x) / betaScale;
  }
  else
    return pdf(x) * ((alphaShape-1.)/x - 1./betaScale);
}


inline Real GammaRandomVariable::pdf_hessian(Real x) const
{
  if (x <= 0.) {
    if      (alphaShape < 1.) return std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return std::numeric_limits<Real>::quiet_NaN();
    else return ExponentialRandomVariable::pdf(x) / (betaScale*betaScale);
  }
  else {
    Real am1 = alphaShape - 1., term = am1 / x - 1. / betaScale;
    return pdf(x) * (term*term - am1 / (x*x));
  }
}


inline Real GammaRandomVariable::log_pdf(Real x) const
{
  if (x <= 0.) {
    if      (alphaShape < 1.) return  std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return -std::numeric_limits<Real>::infinity();
    else return ExponentialRandomVariable::log_pdf(x);
  }
  else
    return (alphaShape-1.)*std::log(x) - x/betaScale
      - bmth::lgamma(alphaShape) - alphaShape*std::log(betaScale);
}


inline Real GammaRandomVariable::log_pdf_gradient(Real x) const
{
  if (x <= 0.) {
    if      (alphaShape < 1.) return -std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return  std::numeric_limits<Real>::infinity();
    else return ExponentialRandomVariable::log_pdf_gradient(x);
  }
  else
    return (alphaShape-1.)/x - 1./betaScale;
}


inline Real GammaRandomVariable::log_pdf_hessian(Real x) const
{
  if (x <= 0.) {
    if      (alphaShape < 1.) return  std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return -std::numeric_limits<Real>::infinity();
    else return ExponentialRandomVariable::log_pdf_hessian(x);
  }
  else
    return (1.-alphaShape)/(x*x);
}


inline Real GammaRandomVariable::inverse_standard_cdf(Real p_cdf) const
{
  gamma_dist gamma1(alphaShape, 1.);
  return bmth::quantile(gamma1, p_cdf);
}


inline Real GammaRandomVariable::standard_pdf(Real z) const
{
  gamma_dist gamma1(alphaShape, 1.);
  return bmth::pdf(gamma1, z);
}


inline Real GammaRandomVariable::log_standard_pdf(Real z) const
{
  if (z <= 0.) {
    if      (alphaShape < 1.) return  std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return -std::numeric_limits<Real>::infinity();
    else return ExponentialRandomVariable::log_standard_pdf(z);
  }
  else
    return (alphaShape-1.)*std::log(z) - z - bmth::lgamma(alphaShape);
}


inline Real GammaRandomVariable::log_standard_pdf_gradient(Real z) const
{
  if (z <= 0.) {
    if      (alphaShape < 1.) return -std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return  std::numeric_limits<Real>::infinity();
    else return ExponentialRandomVariable::log_standard_pdf_gradient(z);
  }
  else
    return (alphaShape-1.)/z - 1.;
}


inline Real GammaRandomVariable::log_standard_pdf_hessian(Real z) const
{
  if (z <= 0.) {
    if      (alphaShape < 1.) return  std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return -std::numeric_limits<Real>::infinity();
    else return ExponentialRandomVariable::log_standard_pdf_hessian(z);
  }
  else
    return (1.-alphaShape)/(z*z);
}


inline Real GammaRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case GA_ALPHA: return alphaShape; break;
  case GA_BETA:  return betaScale;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in GammaRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void GammaRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case GA_ALPHA: alphaShape = val; break;
  case GA_BETA:  betaScale  = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in GammaRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
  update_boost(); // create a new gammaDist instance
}


inline Real GammaRandomVariable::mean() const
{ return bmth::mean(*gammaDist); }


inline Real GammaRandomVariable::median() const
{ return bmth::median(*gammaDist); }


inline Real GammaRandomVariable::mode() const
{ return bmth::mode(*gammaDist); }


inline Real GammaRandomVariable::standard_deviation() const
{ return bmth::standard_deviation(*gammaDist); }


inline Real GammaRandomVariable::variance() const
{ return bmth::variance(*gammaDist); }


inline Real GammaRandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  // correlation warping factor for transformations to STD_NORMAL space
  // Der Kiureghian and Liu, ASCE JEM 112:1, 1986
  Real COV = coefficient_of_variation(), COV_rv;
  switch (rv.type()) { // x-space types mapped to STD_NORMAL u-space

  // Der Kiureghian & Liu: Table 3 (quadratic approximations in COV)
  case GUMBEL: // Max Error 0.3%
    return 1.031 + (0.001 +  0.003*corr)*corr
      + (-0.007 + 0.131*COV - 0.132*corr)*COV; break;

  // Der Kiureghian & Liu: Table 6
  case GAMMA: // Max Error 4.0%
    COV_rv = rv.coefficient_of_variation();
    return 1.002 + 0.022*corr - 0.012*(COV + COV_rv) + 0.001*corr*corr
      + 0.125*(COV*COV + COV_rv*COV_rv) - 0.077*corr*(COV + COV_rv)
      + 0.014*COV*COV_rv; break;
  case FRECHET: // Max Error 4.2%
    COV_rv = rv.coefficient_of_variation();
    return 1.029 + (0.056 +  0.012*corr)*corr
      + (-0.03 + 0.174*COV - 0.313*corr)*COV
      + (0.225 + 0.379*COV_rv +  0.075*COV - 0.182*corr)*COV_rv; break;
  case WEIBULL: // Max Error 4.0%
    COV_rv = rv.coefficient_of_variation();
    return 1.032 + 0.034*corr
      + (-0.007 + 0.121*COV - 0.006*corr + 0.003*COV_rv)*COV
      + (-0.202 + 0.339*COV_rv - 0.111*corr)*COV_rv; break;

  // warping factors are defined once for lower triangle based on uv order
  case NORMAL: case LOGNORMAL: case UNIFORM: case EXPONENTIAL:
    return rv.correlation_warping_factor(*this, corr); break;

  default: // Unsupported warping (should be prevented upsteam)
    PCerr << "Error: unsupported correlation warping for GammaRV." << std::endl;
    abort_handler(-1); return 1.; break;
  }
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  dz/ds is zero if uncorrelated, 
    while dz_ds_factor() manages contributions in the correlated case. */
inline Real GammaRandomVariable::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  bool u_type_err = false, dist_err = false;
  switch (u_type) {
  case STD_GAMMA:
    switch (dist_param) { // x = z*beta
    // For distributions without simple closed-form CDFs (beta, gamma), dx/ds
    // is computed numerically in NatafTransformation::jacobian_dX_dS():
    //case GA_ALPHA:
    case GA_BETA: return z;   break;
    //case GA_LOCATION: - TO DO
    //case GA_SCALE:    - TO DO
    default: dist_err = true; break;
    }
    break;
  default: u_type_err = true; break;
  }

  if (u_type_err)
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in GammaRandomVariable::dx_ds()." << std::endl;
  if (dist_err)
    PCerr << "Error: mapping failure for distribution parameter " << dist_param
	  << " in GammaRandomVariable::dx_ds()." << std::endl;
  if (u_type_err || dist_err)
    abort_handler(-1);
  return 0.;
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  For the uncorrelated case,
    u and z are constants.  For the correlated case, u is a constant, but 
    z(s) = L(s) u due to Nataf dependence on s and dz/ds = dL/ds u. */
inline Real GammaRandomVariable::
dz_ds_factor(short u_type, Real x, Real z) const
{
  // x = z*beta --> add beta * dz/ds for nonzero dz/ds arising from correlation
  switch (u_type) {
  case STD_GAMMA: return betaScale; break;
  default:
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in GammaRandomVariable::dz_ds_factor()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void GammaRandomVariable::update_boost()
{
  if (gammaDist) delete gammaDist;
  gammaDist = new gamma_dist(alphaShape, betaScale);
}


inline void GammaRandomVariable::update(Real alpha, Real beta)
{
  if (!gammaDist || alphaShape != alpha || betaScale != beta)
    { alphaShape = alpha; betaScale = beta; update_boost(); }
}


inline Real GammaRandomVariable::pdf(Real x, Real alpha, Real beta)
{
  gamma_dist gamma1(alpha, beta);
  return bmth::pdf(gamma1, x);
}


inline Real GammaRandomVariable::cdf(Real x, Real alpha, Real beta)
{
  gamma_dist gamma1(alpha, beta);
  return bmth::cdf(gamma1, x);
}


inline Real GammaRandomVariable::inverse_cdf(Real cdf, Real alpha, Real beta)
{
  gamma_dist gamma1(alpha, beta);
  return bmth::quantile(gamma1, cdf);
}


inline void GammaRandomVariable::
moments_from_params(Real alpha, Real beta, Real& mean, Real& std_dev)
{ mean = alpha*beta; std_dev = std::sqrt(alpha)*beta; }

} // namespace Pecos

#endif
