/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 FrechetRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef FRECHET_RANDOM_VARIABLE_HPP
#define FRECHET_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for frechet random variables.

/** Manages alpha and beta parameters. */

class FrechetRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  FrechetRandomVariable();
  /// alternate constructor
  FrechetRandomVariable(Real alpha, Real beta);
  /// destructor
  ~FrechetRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;

  Real log_pdf(Real x) const;

  Real parameter(short dist_param) const;
  void parameter(short dist_param, Real val);

  Real mean() const;
  Real median() const;
  Real mode() const;
  Real standard_deviation() const;
  Real variance() const;
  
  RealRealPair bounds() const;

  Real correlation_warping_factor(const RandomVariable& rv, Real corr) const;

  Real dx_ds(short dist_param, short u_type, Real x, Real z) const;
  Real dz_ds_factor(short u_type, Real x, Real z) const;

  //
  //- Heading: Member functions
  //

  void update(Real alpha, Real beta);

  /// inactive Z_to_X mapping option in NatafTransformation
  Real inverse_log_cdf(Real log_p) const;

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real alpha, Real beta);
  static Real cdf(Real x, Real alpha, Real beta);

  static void moments_from_params(Real alpha, Real beta,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// alpha parameter of frechet random variable (shape)
  Real alphaStat;
  /// beta parameter of frechet random variable (scale)
  Real betaStat;
};


inline FrechetRandomVariable::FrechetRandomVariable():
  RandomVariable(BaseConstructor()), alphaStat(1.), betaStat(1.)
{ ranVarType = FRECHET; }


inline FrechetRandomVariable::FrechetRandomVariable(Real alpha, Real beta):
  RandomVariable(BaseConstructor()), alphaStat(alpha), betaStat(beta)
{ ranVarType = FRECHET; }


inline FrechetRandomVariable::~FrechetRandomVariable()
{ }


inline Real FrechetRandomVariable::cdf(Real x) const
{ return std::exp(-std::pow(betaStat/x, alphaStat)); }


inline Real FrechetRandomVariable::ccdf(Real x) const
{ return -bmth::expm1(-std::pow(betaStat/x, alphaStat)); }


inline Real FrechetRandomVariable::inverse_cdf(Real p_cdf) const
{
  // p = std::exp(-std::pow(beta/x, alpha))
  // beta/x = std::pow( -log p, 1/alpha)
  return betaStat * std::pow(-std::log(p_cdf), -1./alphaStat);
}


inline Real FrechetRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return betaStat * std::pow(-bmth::log1p(-p_ccdf), -1./alphaStat); }


inline Real FrechetRandomVariable::inverse_log_cdf(Real log_p) const
{ return betaStat * std::pow(-log_p, -1./alphaStat); }


//  F(x) = e^(-(beta/x)^alpha)
//  f(x) = F(x) alpha (beta/x)^(alpha-1) beta/x^2
//       = F(x) alpha/beta (beta/x)^(alpha+1)
// f'(x) = alpha/beta ((beta/x)^(alpha+1) f(x) -
//                     F(x) (alpha+1)/beta (beta/x)^(alpha+2))
inline Real FrechetRandomVariable::pdf(Real x) const
{
  Real num = std::pow(betaStat/x, alphaStat);
  return alphaStat/x*num*std::exp(-num);
}


inline Real FrechetRandomVariable::pdf_gradient(Real x) const
{
  Real num = betaStat/x, cdf = std::exp(-std::pow(num, alphaStat)),
    ab_ratio = alphaStat/betaStat,
    pdf = ab_ratio * std::pow(num, alphaStat+1.) * cdf;
  return ab_ratio * (std::pow(num,alphaStat+1.) * pdf - 
		     cdf*(alphaStat+1.)/betaStat * std::pow(num,alphaStat+2.));
}


//inline Real FrechetRandomVariable::pdf_hessian(Real x) const
//{
//  return pdf(x, alphaStat, betaStat) * ...; // TO DO
//}


inline Real FrechetRandomVariable::log_pdf(Real x) const
{
  Real num = std::pow(betaStat/x, alphaStat);
  return std::log(alphaStat/x*num) - num; // fewer operations

  // more decomposed -> less likelihood of overflow?
  //Real num = betaStat/x;
  //return std::log(alphaStat/x) + alphaStat * std::log(num)
  //  - std::pow(num, alphaStat);
}


inline Real FrechetRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case F_ALPHA: return alphaStat; break;
  case F_BETA:  return betaStat;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in FrechetRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void FrechetRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case F_ALPHA: alphaStat = val; break;
  case F_BETA:  betaStat  = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in FrechetRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
}


inline Real FrechetRandomVariable::mean() const
{ return betaStat * bmth::tgamma(1.-1./alphaStat); }


inline Real FrechetRandomVariable::median() const
{ return betaStat / std::pow(std::log(2.), 1./alphaStat); } //inverse_cdf(.5)


inline Real FrechetRandomVariable::mode() const
{ return betaStat * std::pow(alphaStat/(1.+alphaStat), 1./alphaStat); }


inline Real FrechetRandomVariable::standard_deviation() const
{
  Real gam = bmth::tgamma(1.-1./alphaStat);
  return betaStat*std::sqrt(bmth::tgamma(1.-2./alphaStat)-gam*gam);
}


inline Real FrechetRandomVariable::variance() const
{
  Real gam = bmth::tgamma(1.-1./alphaStat);
  return betaStat*betaStat*(bmth::tgamma(1.-2./alphaStat)-gam*gam);
}


inline RealRealPair FrechetRandomVariable::bounds() const
{ return RealRealPair(0., std::numeric_limits<Real>::infinity()); }


inline Real FrechetRandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  // correlation warping factor for transformations to STD_NORMAL space
  // Der Kiureghian and Liu, ASCE JEM 112:1, 1986
  Real COV = coefficient_of_variation(), COV_rv;
  switch (rv.type()) { // x-space types mapped to STD_NORMAL u-space

  // Der Kiureghian & Liu: Table 6
  case FRECHET: { // Max Error 4.3%
    COV_rv = rv.coefficient_of_variation();
    Real COV2 = COV*COV, COV_rv2 = COV_rv*COV_rv, corr2 = corr*corr;
    return 1.086 + 0.054*corr +  0.104*(COV + COV_rv) - 0.055*corr2
      + 0.662*(COV2 + COV_rv2) - 0.57*corr*(COV + COV_rv) +  0.203*COV*COV_rv
      - 0.02*corr2*corr - 0.218*(COV2*COV+COV_rv2*COV_rv)
      - 0.371*corr*(COV2 + COV_rv2) +  0.257*corr2*(COV + COV_rv)
      + 0.141*COV*COV_rv*(COV + COV_rv); break;
  }
  case WEIBULL: // Max Error 3.8%
    COV_rv = rv.coefficient_of_variation();
    return 1.065 + (0.146 + 0.013*corr)*corr
      + (-0.259 + 0.435*COV_rv + 0.034*COV - 0.481*corr)*COV_rv
      + ( 0.241 + 0.372*COV + 0.005*corr)*COV; break;

  // warping factors are defined once for lower triangle based on uv order
  case NORMAL: case LOGNORMAL: case UNIFORM: case EXPONENTIAL: case GAMMA:
  case GUMBEL:
    return rv.correlation_warping_factor(*this, corr); break;

  default: // Unsupported warping (should be prevented upsteam)
    PCerr << "Error: unsupported correlation warping for FrechetRV."<<std::endl;
    abort_handler(-1); return 1.; break;
  }
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  dz/ds is zero if uncorrelated, 
    while dz_ds_factor() manages contributions in the correlated case. */
inline Real FrechetRandomVariable::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  // to STD_NORMAL: x = beta (-ln(Phi(z)))^(-1/alpha)
  bool u_type_err = false, dist_err = false;
  switch (u_type) {
  case STD_NORMAL: {
    switch (dist_param) {
    case F_ALPHA:
      return x * std::log(-NormalRandomVariable::log_std_cdf(z)) /
	(alphaStat*alphaStat);         break;
    case F_BETA:  return x / betaStat; break;
    // Frechet Mean          - TO DO
    // Frechet Std Deviation - TO DO
    default: dist_err = true;          break;
    }
    break;
  }
  //case FRECHET:  TO DO;              break;
  default:         u_type_err = true;  break;
  }

  if (u_type_err)
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in FrechetRandomVariable::dx_ds()." << std::endl;
  if (dist_err)
    PCerr << "Error: mapping failure for distribution parameter " << dist_param
	  << " in FrechetRandomVariable::dx_ds()." << std::endl;
  if (u_type_err || dist_err)
    abort_handler(-1);
  return 0.;
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  For the uncorrelated case,
    u and z are constants.  For the correlated case, u is a constant, but 
    z(s) = L(s) u due to Nataf dependence on s and dz/ds = dL/ds u. */
inline Real FrechetRandomVariable::
dz_ds_factor(short u_type, Real x, Real z) const
{
  switch (u_type) {
  case STD_NORMAL:
    return -x * NormalRandomVariable::std_pdf(z) / (alphaStat *
      NormalRandomVariable::std_cdf(z) * NormalRandomVariable::log_std_cdf(z));
    break;
  //case FRECHET: TO DO; break;
  default:
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in FrechetRandomVariable::dz_ds_factor()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void FrechetRandomVariable::update(Real alpha, Real beta)
{ alphaStat = alpha; betaStat = beta; }

// static functions:

inline Real FrechetRandomVariable::pdf(Real x, Real alpha, Real beta)
{
  Real num = std::pow(beta/x, alpha);
  return alpha/x*num*std::exp(-num);
}


inline Real FrechetRandomVariable::cdf(Real x, Real alpha, Real beta)
{ return std::exp(-std::pow(beta/x, alpha)); }


inline void FrechetRandomVariable::
moments_from_params(Real alpha, Real beta, Real& mean, Real& std_dev)
{
  // See Haldar and Mahadevan, p. 91-92
  Real gam = bmth::tgamma(1.-1./alpha);
  mean     = beta*gam;
  std_dev  = beta*std::sqrt(bmth::tgamma(1.-2./alpha)-gam*gam);
}

} // namespace Pecos

#endif
