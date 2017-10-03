/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 BoundedLognormalRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef BOUNDED_LOGNORMAL_RANDOM_VARIABLE_HPP
#define BOUNDED_LOGNORMAL_RANDOM_VARIABLE_HPP

#include "LognormalRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for bounded lognormal random variables.

/** Manages lower and upper bounds. */

class BoundedLognormalRandomVariable: public LognormalRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  BoundedLognormalRandomVariable();
  /// alternate constructor
  BoundedLognormalRandomVariable(Real lambda, Real zeta, Real lwr, Real upr);
  /// destructor
  ~BoundedLognormalRandomVariable();

  //
  //- Heading: Member functions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  //Real pdf_gradient(Real x) const;
  //Real pdf_hessian(Real x) const;
  Real log_pdf(Real x) const;
  Real log_pdf_gradient(Real x) const;
  Real log_pdf_hessian(Real x) const;

  Real parameter(short dist_param) const;
  void parameter(short dist_param, Real val);

  Real mean() const;
  Real median() const;
  Real mode() const;
  Real standard_deviation() const;
  Real variance() const;
  
  RealRealPair moments() const;
  RealRealPair bounds() const;

  Real coefficient_of_variation() const;

  Real dx_ds(short dist_param, short u_type, Real x, Real z) const;
  Real dz_ds_factor(short u_type, Real x, Real z) const;

  void update(Real lambda, Real zeta, Real lwr, Real upr);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real lambda, Real zeta, Real lwr, Real upr);
  static Real cdf(Real x, Real lambda, Real zeta, Real lwr, Real upr);

protected:

  //
  //- Heading: Data
  //

  /// lower bound of bounded_lognormal random variable
  Real lowerBnd;
  /// upper bound of bounded_lognormal random variable
  Real upperBnd;
};


inline BoundedLognormalRandomVariable::BoundedLognormalRandomVariable():
  LognormalRandomVariable(), lowerBnd(0.),
  upperBnd(std::numeric_limits<Real>::infinity())
{ ranVarType = BOUNDED_LOGNORMAL; }


inline BoundedLognormalRandomVariable::
BoundedLognormalRandomVariable(Real lambda, Real zeta, Real lwr, Real upr):
  LognormalRandomVariable(lambda, zeta), lowerBnd(lwr), upperBnd(upr)
{ ranVarType = BOUNDED_LOGNORMAL; }


inline BoundedLognormalRandomVariable::~BoundedLognormalRandomVariable()
{ }


inline Real BoundedLognormalRandomVariable::
pdf(Real x, Real lambda, Real zeta, Real lwr, Real upr)
{
  if (x < lwr || x > upr) return 0.;
  else {
    Real Phi_lms = (lwr > 0.) ?
      NormalRandomVariable::std_cdf((std::log(lwr)-lambda)/zeta) : 0.;
    Real Phi_ums = (upr < std::numeric_limits<Real>::infinity()) ?
      NormalRandomVariable::std_cdf((std::log(upr)-lambda)/zeta) : 1.;
    return NormalRandomVariable::std_pdf((std::log(x)-lambda)/zeta) /
      (Phi_ums-Phi_lms)/x/zeta;
  }
}


inline Real BoundedLognormalRandomVariable::
cdf(Real x, Real lambda, Real zeta, Real lwr, Real upr)
{
  if      (x < lwr) return 0.;
  else if (x > upr) return 1.;
  else {
    Real Phi_lms = (lwr > 0.) ?
      NormalRandomVariable::std_cdf((std::log(lwr)-lambda)/zeta) : 0.;
    Real Phi_ums = (upr < std::numeric_limits<Real>::infinity()) ?
      NormalRandomVariable::std_cdf((std::log(upr)-lambda)/zeta) : 1.;
    return (NormalRandomVariable::std_cdf((std::log(x)-lambda)/zeta) - Phi_lms)
      / (Phi_ums - Phi_lms);
  }
}


inline Real BoundedLognormalRandomVariable::cdf(Real x) const
{ return cdf(x, lnLambda, lnZeta, lowerBnd, upperBnd); }


inline Real BoundedLognormalRandomVariable::ccdf(Real x) const
{
  if      (x < lowerBnd) return 1.;
  else if (x > upperBnd) return 0.;
  else {
    Real Phi_lms = (lowerBnd > 0.) ?
      NormalRandomVariable::std_cdf((std::log(lowerBnd)-lnLambda)/lnZeta) : 0.;
    Real Phi_ums = (upperBnd < std::numeric_limits<Real>::infinity()) ?
      NormalRandomVariable::std_cdf((std::log(upperBnd)-lnLambda)/lnZeta) : 1.;
    return (Phi_ums - NormalRandomVariable::
	    std_cdf((std::log(x)-lnLambda)/lnZeta)) / (Phi_ums - Phi_lms);
  }
}


inline Real BoundedLognormalRandomVariable::inverse_cdf(Real p_cdf) const
{
  if      (p_cdf <= 0.) return lowerBnd;
  else if (p_cdf >= 1.) return upperBnd;
  else {
    // p = (Phi((log(x)-lambda)/zeta) - Phi_lms)/(Phi_ums - Phi_lms)
    // log(x) = Phi_inverse[p * (Phi_ums - Phi_lms) + Phi_lms] * zeta + lambda
    Real Phi_lms = (lowerBnd > 0.) ?
      NormalRandomVariable::std_cdf((std::log(lowerBnd)-lnLambda)/lnZeta) : 0.;
    Real Phi_ums = (upperBnd < std::numeric_limits<Real>::infinity()) ?
      NormalRandomVariable::std_cdf((std::log(upperBnd)-lnLambda)/lnZeta) : 1.;
    return std::exp(lnLambda + lnZeta * NormalRandomVariable::
		    inverse_std_cdf(p_cdf * (Phi_ums - Phi_lms) + Phi_lms));
  }
}


inline Real BoundedLognormalRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  if      (p_ccdf >= 1.) return lowerBnd;
  else if (p_ccdf <= 0.) return upperBnd;
  else {
    // p = (Phi_ums - Phi((log(x)-lambda)/zeta))/(Phi_ums - Phi_lms)
    // log(x) = lambda + zeta * Phi_inverse[Phi_ums - p * (Phi_ums - Phi_lms)]
    Real Phi_lms = (lowerBnd > 0.) ?
      NormalRandomVariable::std_cdf((std::log(lowerBnd)-lnLambda)/lnZeta) : 0.;
    Real Phi_ums = (upperBnd < std::numeric_limits<Real>::infinity()) ?
      NormalRandomVariable::std_cdf((std::log(upperBnd)-lnLambda)/lnZeta) : 1.;
    return std::exp(lnLambda + lnZeta * NormalRandomVariable::
		    inverse_std_cdf(Phi_ums - p_ccdf * (Phi_ums - Phi_lms)));
  }
}


inline Real BoundedLognormalRandomVariable::pdf(Real x) const
{ return pdf(x, lnLambda, lnZeta, lowerBnd, upperBnd); }


// Same form as base definition for different virtual pdf(x):
//inline Real BoundedLognormalRandomVariable::pdf_gradient(Real x) const
//{ return -pdf(x) / x * (1. + (std::log(x) - lnLambda) / (lnZeta*lnZeta) ); }


// Same form as base definition for different virtual pdf(x):
//inline Real BoundedLognormalRandomVariable::pdf_hessian(Real x) const
//{
//  Real zeta_sq = lnZeta*lnZeta, num = (std::log(x) - lnLambda) / zeta_sq;
//  return pdf(x) / (x*x) * (num * (1. + num) - 1. / zeta_sq);
//}


inline Real BoundedLognormalRandomVariable::log_pdf(Real x) const
{
  Real dbl_inf = std::numeric_limits<Real>::infinity();
  if (x < lowerBnd || x > upperBnd) return -dbl_inf;
  else {
    Real Phi_lms = (lowerBnd > 0.) ?
      NormalRandomVariable::std_cdf((std::log(lowerBnd)-lnLambda)/lnZeta) : 0.;
    Real Phi_ums = (upperBnd < dbl_inf) ?
      NormalRandomVariable::std_cdf((std::log(upperBnd)-lnLambda)/lnZeta) : 1.;
    return LognormalRandomVariable::log_pdf(x) - std::log(Phi_ums-Phi_lms);
  }
}


inline Real BoundedLognormalRandomVariable::log_pdf_gradient(Real x) const
{
  if (x < lowerBnd || x > upperBnd) return 0.;
  else // same as base definition since std::log(Phi_ums-Phi_lms) term drops
    return LognormalRandomVariable::log_pdf_gradient(x);
}


inline Real BoundedLognormalRandomVariable::log_pdf_hessian(Real x) const
{
  if (x < lowerBnd || x > upperBnd) return 0.;
  else // same as base definition since std::log(Phi_ums-Phi_lms) term drops
    return LognormalRandomVariable::log_pdf_hessian(x);
}


inline Real BoundedLognormalRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case LN_MEAN: case LN_STD_DEV: case LN_LAMBDA: case LN_ZETA:
    return LognormalRandomVariable::parameter(dist_param); break;
  case LN_LWR_BND: return lowerBnd; break;
  case LN_UPR_BND: return upperBnd; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in BoundedLognormalRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
}


inline void BoundedLognormalRandomVariable::
parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case LN_MEAN: case LN_STD_DEV: case LN_LAMBDA: case LN_ZETA:
    LognormalRandomVariable::parameter(dist_param, val); break;
  case LN_LWR_BND: lowerBnd = val; break;
  case LN_UPR_BND: upperBnd = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in BoundedLognormalRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
}


inline Real BoundedLognormalRandomVariable::mean() const
{
  Real Phi_lms = 0., Phi_ums = 1., term = 0.;
  if (lowerBnd > 0.) {
    Real lms = (std::log(lowerBnd)-lnLambda)/lnZeta;
    Phi_lms  = NormalRandomVariable::std_cdf(lms);
    term    += NormalRandomVariable::std_cdf(lnZeta - lms);
  }
  if (upperBnd < std::numeric_limits<Real>::infinity()) {
    Real ums = (std::log(upperBnd)-lnLambda)/lnZeta;
    Phi_ums  = NormalRandomVariable::std_cdf(ums);
    term    -= NormalRandomVariable::std_cdf(lnZeta - ums);
  }
  return LognormalRandomVariable::mean() * term / (Phi_ums - Phi_lms);
}


inline Real BoundedLognormalRandomVariable::median() const
{ return inverse_cdf(.5); }


inline Real BoundedLognormalRandomVariable::mode() const
{
  Real ln_mode = LognormalRandomVariable::mode();
  if      (ln_mode < lowerBnd) return lowerBnd;
  else if (ln_mode > upperBnd) return upperBnd;
  else                         return ln_mode;
}


inline Real BoundedLognormalRandomVariable::standard_deviation() const
{ return std::sqrt(variance()); }


inline Real BoundedLognormalRandomVariable::variance() const
{ return moments().second; } // variance computation requires mean computation


inline RealRealPair BoundedLognormalRandomVariable::moments() const
{
  Real Phi_lms = 0., Phi_ums = 1., term1 = 0., term2 = 0.;
  if (lowerBnd > 0.) {
    Real lms = (std::log(lowerBnd)-lnLambda)/lnZeta;
    Phi_lms  = NormalRandomVariable::std_cdf(lms);
    term1   += NormalRandomVariable::std_cdf(lnZeta - lms);
    term2   += NormalRandomVariable::std_cdf(2.*lnZeta - lms);
  }
  if (upperBnd < std::numeric_limits<Real>::infinity()) {
    Real ums = (std::log(upperBnd)-lnLambda)/lnZeta;
    Phi_ums  = NormalRandomVariable::std_cdf(ums);
    term1   -= NormalRandomVariable::std_cdf(lnZeta - ums);
    term2   -= NormalRandomVariable::std_cdf(2.*lnZeta - ums);
  }
  Real Phi_range  = Phi_ums - Phi_lms,
       mean       = LognormalRandomVariable::mean() * term1 / Phi_range,
       ln_raw_mom = std::exp(2.*(lnLambda+lnZeta*lnZeta));
  return RealRealPair(mean, ln_raw_mom * term2 / Phi_range - mean * mean);
}


inline RealRealPair BoundedLognormalRandomVariable::bounds() const
{ return RealRealPair(lowerBnd, upperBnd); }


inline Real BoundedLognormalRandomVariable::coefficient_of_variation() const
{ RealRealPair mom = moments(); return mom.second / mom.first; }


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  dz/ds is zero if uncorrelated, 
    while dz_ds_factor() manages contributions in the correlated case. */
inline Real BoundedLognormalRandomVariable::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  bool u_type_err = false, dist_err = false;
  switch (u_type) {
  case STD_NORMAL: {
    Real dlambda_ds = 0., dzeta_ds = 0., dlwr_ds = 0., dupr_ds = 0.;
    switch (dist_param) {
    case LN_MEAN: {
      Real mean, stdev, mean_sq, var;
      moments_from_params(lnLambda, lnZeta, mean, stdev);
      mean_sq = mean*mean; var = stdev*stdev;
      dlambda_ds = (1.+var/(mean_sq+var))/mean;
      dzeta_ds   = -var/lnZeta/mean/(mean_sq+var);                      break;
    }
    case LN_STD_DEV: {
      Real mean, stdev, mean_sq, var;
      moments_from_params(lnLambda, lnZeta, mean, stdev);
      mean_sq = mean*mean; var = stdev*stdev;
      dlambda_ds = -stdev/(mean_sq+var);
      dzeta_ds   =  stdev/lnZeta/(mean_sq+var);                         break;
    }
    case LN_LAMBDA: dlambda_ds = 1.;                                    break;
    case LN_ZETA:   dzeta_ds   = 1.;                                    break;
    case LN_ERR_FACT: {
      Real inv_95 = NormalRandomVariable::inverse_std_cdf(0.95),
	err_fact = std::exp(inv_95 * lnZeta);
      dzeta_ds   = 1./inv_95/err_fact;
      dlambda_ds = -lnZeta*dzeta_ds;                                    break;
    }
    case LN_LWR_BND: dlwr_ds = 1.; break;
    case LN_UPR_BND: dupr_ds = 1.; break;
    default: dist_err = true;      break;
    }

    Real lms, phi_lms, ums, phi_ums, xms = (std::log(x)-lnLambda)/lnZeta,
      phi_xms = NormalRandomVariable::std_pdf(xms), dlms_ds, dums_ds;
    if (lowerBnd > 0.) {
      lms     = (std::log(lowerBnd)-lnLambda)/lnZeta;
      phi_lms = NormalRandomVariable::std_pdf(lms);
      dlms_ds = (dlwr_ds/lowerBnd - dlambda_ds - lms*dzeta_ds)/lnZeta;
    }
    else phi_lms = dlms_ds = 0.;
    if (upperBnd < std::numeric_limits<Real>::infinity()) { 
      ums     = (std::log(upperBnd)-lnLambda)/lnZeta;
      phi_ums = NormalRandomVariable::std_pdf(ums);
      dums_ds = (dupr_ds/upperBnd - dlambda_ds - ums*dzeta_ds)/lnZeta;
    }
    else phi_ums = dums_ds = 0.;
    Real dxms_ds = NormalRandomVariable::std_cdf(z)/phi_xms*
      (phi_ums*dums_ds - phi_lms*dlms_ds) + phi_lms/phi_xms*dlms_ds;
    return x*(lnZeta*dxms_ds + dlambda_ds + xms*dzeta_ds);            break;
  }
  //case BOUNDED_LOGNORMAL: TO DO;                                    break;
  default: u_type_err = true;                                         break;
  }

  if (u_type_err)
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in BoundedLognormalRandomVariable::dx_ds()." << std::endl;
  if (dist_err)
    PCerr << "Error: mapping failure for distribution parameter " << dist_param
	  << " in BoundedLognormalRandomVariable::dx_ds()." << std::endl;
  if (u_type_err || dist_err)
    abort_handler(-1);
  return 0.;
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  For the uncorrelated case,
    u and z are constants.  For the correlated case, u is a constant, but 
    z(s) = L(s) u due to Nataf dependence on s and dz/ds = dL/ds u. */
inline Real BoundedLognormalRandomVariable::
dz_ds_factor(short u_type, Real x, Real z) const
{
  Real num = 0., lms, ums, xms = (std::log(x)-lnLambda)/lnZeta;
  switch (u_type) {
  case STD_NORMAL:
    if (upperBnd < std::numeric_limits<Real>::infinity()) {
      ums = (std::log(upperBnd)-lnLambda)/lnZeta;
      num += NormalRandomVariable::std_cdf(ums);
    }
    else num += 1.;
    if (lowerBnd > 0.) {
      lms = (std::log(lowerBnd)-lnLambda)/lnZeta;
      num -= NormalRandomVariable::std_cdf(lms);
    }
    return num * NormalRandomVariable::std_pdf(z) /
      NormalRandomVariable::std_pdf(xms);
    break;
  //case BOUNDED_LOGNORMAL: TO DO; break;
  default:
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in BoundedLognormalRandomVariable::dz_ds_factor()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void BoundedLognormalRandomVariable::
update(Real lambda, Real zeta, Real lwr, Real upr)
{ lnLambda = lambda; lnZeta = zeta; lowerBnd = lwr; upperBnd = upr; }

} // namespace Pecos

#endif
