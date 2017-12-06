/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 InvGammaRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef INV_GAMMA_RANDOM_VARIABLE_HPP
#define INV_GAMMA_RANDOM_VARIABLE_HPP

#include "ExponentialRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for gamma random variables.

/** Manages alphaShape and inherits betaScale.  This follows the definition 
    at https://en.wikipedia.org/wiki/Inverse-gamma_distribution.  This 
    implementation also supports a standard inverse-gamma with beta = 1. */

class InvGammaRandomVariable: public ExponentialRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  InvGammaRandomVariable();
  /// alternate constructor
  InvGammaRandomVariable(Real alpha, Real beta);
  /// destructor
  ~InvGammaRandomVariable();

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

  /// create a new invGammaDist instance
  void update_boost();

  //
  //- Heading: Data
  //

  /// alpha shape parameter of inverse gamma random variable (statistical
  /// PDF convention; differs from GenLaguerre polynomial convention)
  Real alphaShape;

  /// pointer to the Boost inv_gamma_dist instance
  inv_gamma_dist* invGammaDist;
};


inline InvGammaRandomVariable::InvGammaRandomVariable():
  ExponentialRandomVariable(), alphaShape(3.), invGammaDist(NULL)
{ ranVarType = INV_GAMMA; }


inline InvGammaRandomVariable::InvGammaRandomVariable(Real alpha, Real beta):
  ExponentialRandomVariable(beta), alphaShape(alpha),
  invGammaDist(new inv_gamma_dist(alphaShape, betaScale))
{ ranVarType = INV_GAMMA; }


inline InvGammaRandomVariable::~InvGammaRandomVariable()
{ if (invGammaDist) delete invGammaDist; }


inline Real InvGammaRandomVariable::cdf(Real x) const
{ return bmth::cdf(*invGammaDist, x); }


inline Real InvGammaRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*invGammaDist, x)); }


inline Real InvGammaRandomVariable::inverse_cdf(Real p_cdf) const
{ return bmth::quantile(*invGammaDist, p_cdf); }


inline Real InvGammaRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*invGammaDist, p_ccdf)); }


inline Real InvGammaRandomVariable::pdf(Real x) const
{ return bmth::pdf(*invGammaDist, x); }


inline Real InvGammaRandomVariable::pdf_gradient(Real x) const
{
  /*
  if (x <= 0.) {
    if      (alphaShape < 1.) return -std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return  std::numeric_limits<Real>::quiet_NaN();
    else return -ExponentialRandomVariable::pdf(x) / betaScale;
  }
  else
    return pdf(x) * ((alphaShape-1.)/x - 1./betaScale);
  */
  PCerr << "Error: InvGammaRandomVariable::pdf_gradient() not implemented."
	<< std::endl;
  abort_handler(-1);
  return std::numeric_limits<Real>::quiet_NaN();
}


inline Real InvGammaRandomVariable::pdf_hessian(Real x) const
{
  /*
  if (x <= 0.) {
    if      (alphaShape < 1.) return std::numeric_limits<Real>::infinity();
    else if (alphaShape > 1.) return std::numeric_limits<Real>::quiet_NaN();
    else return ExponentialRandomVariable::pdf(x) / (betaScale*betaScale);
  }
  else {
    Real am1 = alphaShape - 1., term = am1 / x - 1. / betaScale;
    return pdf(x) * (term*term - am1 / (x*x));
  }
  */
  PCerr << "Error: InvGammaRandomVariable::pdf_hessian() not implemented."
	<< std::endl;
  abort_handler(-1);
  return std::numeric_limits<Real>::quiet_NaN();
}


inline Real InvGammaRandomVariable::log_pdf(Real x) const
{
  // BMA TODO: This function should be reimplemented in terms of
  // gamma_p_derivative to avoid overflow.  For now, I made a small
  // improvement: log(tgamma(alpha)) -> lgamma(alpha)
  if (x <= 0.) // throw domain error?
    return std::numeric_limits<Real>::quiet_NaN();
  else
    return alphaShape*std::log(betaScale) - bmth::lgamma(alphaShape)
      - (alphaShape+1.)*std::log(x) - betaScale / x;
}


inline Real InvGammaRandomVariable::log_pdf_gradient(Real x) const
{
  if (x <= 0.) // throw domain error?
    return std::numeric_limits<Real>::quiet_NaN();
  else
    return (betaScale / x - alphaShape - 1.) / x;
}


inline Real InvGammaRandomVariable::log_pdf_hessian(Real x) const
{
  if (x <= 0.) // throw domain error?
    return std::numeric_limits<Real>::quiet_NaN();
  else
    return (alphaShape + 1. - 2.*betaScale / x) / (x*x);
}


inline Real InvGammaRandomVariable::inverse_standard_cdf(Real p_cdf) const
{
  inv_gamma_dist inv_gamma1(alphaShape, 1.);
  return bmth::quantile(inv_gamma1, p_cdf);
}


inline Real InvGammaRandomVariable::standard_pdf(Real z) const
{
  inv_gamma_dist inv_gamma1(alphaShape, 1.);
  return bmth::pdf(inv_gamma1, z);
}


inline Real InvGammaRandomVariable::log_standard_pdf(Real z) const
{
  if (z <= 0.) // throw domain error?
    return std::numeric_limits<Real>::quiet_NaN();
  else
    return -bmth::lgamma(alphaShape)  // log(gamma(alpha))
      - (alphaShape+1.)*std::log(z) - 1. / z;
}


inline Real InvGammaRandomVariable::log_standard_pdf_gradient(Real z) const
{
  if (z <= 0.) // throw domain error?
    return std::numeric_limits<Real>::quiet_NaN();
  else
    return (1. / z - alphaShape - 1.) / z;
}


inline Real InvGammaRandomVariable::log_standard_pdf_hessian(Real z) const
{
  if (z <= 0.) // throw domain error?
    return std::numeric_limits<Real>::quiet_NaN();
  else
    return (alphaShape + 1. - 2. / z) / (z*z);
}


inline Real InvGammaRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case IGA_ALPHA: return alphaShape; break;
  case IGA_BETA:  return betaScale;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in InvGammaRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void InvGammaRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case IGA_ALPHA: alphaShape = val; break;
  case IGA_BETA:  betaScale  = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in InvGammaRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
  update_boost(); // create a new invGammaDist instance
}


inline Real InvGammaRandomVariable::mean() const
{ return bmth::mean(*invGammaDist); }


inline Real InvGammaRandomVariable::median() const
{ return bmth::median(*invGammaDist); }


inline Real InvGammaRandomVariable::mode() const
{ return bmth::mode(*invGammaDist); }


inline Real InvGammaRandomVariable::standard_deviation() const
{ return bmth::standard_deviation(*invGammaDist); }


inline Real InvGammaRandomVariable::variance() const
{ return bmth::variance(*invGammaDist); }


inline Real InvGammaRandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  // hide ExponentialRandomVariable implementation
  PCerr << "Error: InvGammaRandomVariable::correlation_warping_factor() not "
	<< "implemented." << std::endl;
  abort_handler(-1);
  return std::numeric_limits<Real>::quiet_NaN();
}


inline Real InvGammaRandomVariable::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  // hide ExponentialRandomVariable implementation
  PCerr << "Error: InvGammaRandomVariable::dx_ds() not implemented."<<std::endl;
  abort_handler(-1);
  return std::numeric_limits<Real>::quiet_NaN();
}


inline Real InvGammaRandomVariable::
dz_ds_factor(short u_type, Real x, Real z) const
{
  // hide ExponentialRandomVariable implementation
  PCerr << "Error: InvGammaRandomVariable::dz_ds_factor() not implemented."
	<< std::endl;
  abort_handler(-1);
  return std::numeric_limits<Real>::quiet_NaN();
}


inline void InvGammaRandomVariable::update_boost()
{
  if (invGammaDist) delete invGammaDist;
  invGammaDist = new inv_gamma_dist(alphaShape, betaScale);
}


inline void InvGammaRandomVariable::update(Real alpha, Real beta)
{
  if (!invGammaDist || alphaShape != alpha || betaScale != beta)
    { alphaShape = alpha; betaScale = beta; update_boost(); }
}


inline Real InvGammaRandomVariable::pdf(Real x, Real alpha, Real beta)
{
  inv_gamma_dist inv_gamma1(alpha, beta);
  return bmth::pdf(inv_gamma1, x);
}


inline Real InvGammaRandomVariable::cdf(Real x, Real alpha, Real beta)
{
  inv_gamma_dist inv_gamma1(alpha, beta);
  return bmth::cdf(inv_gamma1, x);
}


inline Real InvGammaRandomVariable::inverse_cdf(Real cdf, Real alpha, Real beta)
{
  inv_gamma_dist inv_gamma1(alpha, beta);
  return bmth::quantile(inv_gamma1, cdf);
}


inline void InvGammaRandomVariable::
moments_from_params(Real alpha, Real beta, Real& mean, Real& std_dev)
{
  inv_gamma_dist inv_gamma1(alpha, beta);
  mean    = bmth::mean(inv_gamma1);                // domain_error if alpha <= 1
  std_dev = std::sqrt(bmth::variance(inv_gamma1)); // domain_error if alpha <= 2
}

} // namespace Pecos

#endif
