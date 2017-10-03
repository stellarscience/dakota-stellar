/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 WeibullRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef WEIBULL_RANDOM_VARIABLE_HPP
#define WEIBULL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for weibull random variables.

/** Manages alpha and beta parameters. */

class WeibullRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  WeibullRandomVariable();
  /// alternate constructor
  WeibullRandomVariable(Real alpha, Real beta);
  /// destructor
  ~WeibullRandomVariable();

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
  Real inverse_log_ccdf(Real log_p_ccdf) const;

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real alpha, Real beta);
  static Real cdf(Real x, Real alpha, Real beta);

  static void moments_from_params(Real alpha, Real beta,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Member functions
  //

  /// create a new weibullDist instance
  void update_boost();

  //
  //- Heading: Data
  //

  /// alpha parameter of weibull random variable (shape)
  Real alphaStat;
  /// beta parameter of weibull random variable (scale)
  Real betaStat;

  /// pointer to the Boost weibull_distribution instance
  weibull_dist* weibullDist;
};


inline WeibullRandomVariable::WeibullRandomVariable():
  RandomVariable(BaseConstructor()), alphaStat(1.), betaStat(1.),
  weibullDist(NULL)
{ ranVarType = WEIBULL; }


inline WeibullRandomVariable::WeibullRandomVariable(Real alpha, Real beta):
  RandomVariable(BaseConstructor()), alphaStat(alpha), betaStat(beta),
  weibullDist(new weibull_dist(alpha, beta))
{ ranVarType = WEIBULL; }


inline WeibullRandomVariable::~WeibullRandomVariable()
{ if (weibullDist) delete weibullDist; }


inline Real WeibullRandomVariable::cdf(Real x) const
{ return bmth::cdf(*weibullDist, x); }


inline Real WeibullRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*weibullDist, x)); }


inline Real WeibullRandomVariable::inverse_cdf(Real p_cdf) const
{
  return bmth::quantile(*weibullDist, p_cdf);
  //return betaStat * std::pow(-bmth::log1p(-p), 1./alphaStat);
}


inline Real WeibullRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  return bmth::quantile(complement(*weibullDist, p_ccdf));
  //return betaStat * std::pow(-std::log(p_ccdf), 1./alphaStat);
}


inline Real WeibullRandomVariable::inverse_log_ccdf(Real log_p_ccdf) const
{ return betaStat * std::pow(-log_p_ccdf, 1./alphaStat); }


//  F(x) = 1.-e^(-(x/beta)^alpha)
//  f(x) = alpha/beta e^(-(x/beta)^alpha) (x/beta)^(alpha-1)
// f'(x) = alpha/beta (e^(-(x/beta)^alpha) (alpha-1)/beta
//                     (x/beta)^(alpha-2) - (x/beta)^(alpha-1) f(x))
inline Real WeibullRandomVariable::pdf(Real x) const
{
  return bmth::pdf(*weibullDist, x);
  //return alpha/beta * std::pow(x/beta,alpha-1.) *
  //  std::exp(-std::pow(x/beta,alpha));
}


inline Real WeibullRandomVariable::pdf_gradient(Real x) const
{
  Real num = x / betaStat, num2 = std::exp(-std::pow(num, alphaStat)),
    ab_ratio = alphaStat/betaStat,
    pdf = ab_ratio * num2 * std::pow(num, alphaStat - 1.);
  return ab_ratio * (num2 * (alphaStat-1.) / betaStat * 
		     std::pow(num, alphaStat - 2.) -
		     std::pow(num, alphaStat - 1.) * pdf);
}


//inline Real WeibullRandomVariable::pdf_hessian(Real x) const
//{
//  return pdf(x) * ...; // TO DO
//}


inline Real WeibullRandomVariable::log_pdf(Real x) const
{
  Real num = x/betaStat;
  return std::log(alphaStat/betaStat) + (alphaStat-1.) * std::log(num)
    - std::pow(num,alphaStat);
}


inline Real WeibullRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case W_ALPHA: return alphaStat; break;
  case W_BETA:  return betaStat;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in WeibullRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void WeibullRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case W_ALPHA: alphaStat = val; break;
  case W_BETA:  betaStat  = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in WeibullRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
  update_boost(); // create a new weibullDist instance
}


inline Real WeibullRandomVariable::mean() const
{ return bmth::mean(*weibullDist); }


inline Real WeibullRandomVariable::median() const
{ return bmth::median(*weibullDist); }


inline Real WeibullRandomVariable::mode() const
{ return bmth::mode(*weibullDist); }


inline Real WeibullRandomVariable::standard_deviation() const
{ return bmth::standard_deviation(*weibullDist); }


inline Real WeibullRandomVariable::variance() const
{ return bmth::variance(*weibullDist); }


inline RealRealPair WeibullRandomVariable::bounds() const
{ return RealRealPair(0., std::numeric_limits<Real>::infinity()); }


inline Real WeibullRandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  // correlation warping factor for transformations to STD_NORMAL space
  // Der Kiureghian and Liu, ASCE JEM 112:1, 1986
  switch (rv.type()) { // x-space types mapped to STD_NORMAL u-space

  // Der Kiureghian & Liu: Table 6 (quadratic approximations in COV)
  case WEIBULL: {  // Max Error 2.6%
    Real COV    =    coefficient_of_variation(),
         COV_rv = rv.coefficient_of_variation();
    return 1.063 + (-0.004 - 0.001*corr)*corr - 0.007*COV*COV_rv
      + (COV + COV_rv) * (0.007*corr - 0.2) + 0.337*(COV*COV + COV_rv*COV_rv);
    break;
  }

  // warping factors are defined once for lower triangle based on uv order
  case NORMAL: case LOGNORMAL: case UNIFORM: case EXPONENTIAL: case GAMMA:
  case GUMBEL: case FRECHET:
    return rv.correlation_warping_factor(*this, corr); break;

  default: // Unsupported warping (should be prevented upsteam)
    PCerr << "Error: unsupported correlation warping for WeibullRV."<<std::endl;
    abort_handler(-1); return 1.; break;
  }
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  dz/ds is zero if uncorrelated, 
    while dz_ds_factor() manages contributions in the correlated case. */
inline Real WeibullRandomVariable::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  // to STD_NORMAL: x = beta (-ln(1-Phi(z)))^(1/alpha)
  bool u_type_err = false, dist_err = false;
  switch (u_type) {
  case STD_NORMAL:
    switch (dist_param) {
    case W_ALPHA:
      return -x * std::log(-NormalRandomVariable::log_std_ccdf(z)) /
	(alphaStat*alphaStat);        break;
    case W_BETA: return x / betaStat; break;
    // Weibull Mean          - TO DO
    // Weibull Std Deviation - TO DO
    default:      dist_err = true;    break;
    }
    break;
  //case WEIBULL: TO DO;              break;
  default:        u_type_err = true;  break;
  }

  if (u_type_err)
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in WeibullRandomVariable::dx_ds()." << std::endl;
  if (dist_err)
    PCerr << "Error: mapping failure for distribution parameter " << dist_param
	  << " in WeibullRandomVariable::dx_ds()." << std::endl;
  if (u_type_err || dist_err)
    abort_handler(-1);
  return 0.;
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  For the uncorrelated case,
    u and z are constants.  For the correlated case, u is a constant, but 
    z(s) = L(s) u due to Nataf dependence on s and dz/ds = dL/ds u. */
inline Real WeibullRandomVariable::
dz_ds_factor(short u_type, Real x, Real z) const
{
  switch (u_type) {
  case STD_NORMAL:
    return -x * NormalRandomVariable::std_pdf(z) / (alphaStat *
      NormalRandomVariable::std_ccdf(z) *
      NormalRandomVariable::log_std_ccdf(z)); break;
  //case WEIBULL:   TO DO;                    break;
  default:
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in WeibullRandomVariable::dz_ds_factor()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void WeibullRandomVariable::update_boost()
{
  if (weibullDist) delete weibullDist;
  weibullDist = new weibull_dist(alphaStat, betaStat);
}


inline void WeibullRandomVariable::update(Real alpha, Real beta)
{
  if (!weibullDist || alphaStat != alpha || betaStat != beta)
    { alphaStat = alpha; betaStat = beta; update_boost(); }
}

// Static member functions:

inline Real WeibullRandomVariable::pdf(Real x, Real alpha, Real beta)
{
  weibull_dist weibull1(alpha, beta);
  return bmth::pdf(weibull1, x);
  //return alpha/beta * std::pow(x/beta,alpha-1.) *
  //  std::exp(-std::pow(x/beta,alpha));
}


inline Real WeibullRandomVariable::cdf(Real x, Real alpha, Real beta)
{
  weibull_dist weibull1(alpha, beta);
  return bmth::cdf(weibull1, x);
  // avoid numerical probs when exp()~1
  //return -std::expm1(-std::pow(x/beta, alpha));
}


inline void WeibullRandomVariable::
moments_from_params(Real alpha, Real beta, Real& mean, Real& std_dev)
{
  // See Haldar and Mahadevan, p. 97
  Real gam = bmth::tgamma(1.+1./alpha),
       COV = std::sqrt(bmth::tgamma(1.+2./alpha)/gam/gam - 1.);
  mean     = beta*gam;
  std_dev  = COV*mean;
}

} // namespace Pecos

#endif
