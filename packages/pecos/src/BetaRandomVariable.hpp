/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 BetaRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef BETA_RANDOM_VARIABLE_HPP
#define BETA_RANDOM_VARIABLE_HPP

#include "UniformRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for beta random variables.

/** Manages alpha and beta and inherits bounds.  See Haldar and
    Mahadevan, p. 72. */

class BetaRandomVariable: public UniformRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  BetaRandomVariable();
  /// alternate constructor
  BetaRandomVariable(Real alpha, Real beta, Real lwr, Real upr);
  /// destructor
  ~BetaRandomVariable();

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

  // inherited from UniformRandomVariable
  //Real to_standard(Real x) const;
  //Real from_standard(Real z) const;

  Real parameter(short dist_param) const;
  void parameter(short dist_param, Real val);

  Real mean() const;
  Real median() const;
  Real mode() const;
  Real standard_deviation() const;
  Real variance() const;

  Real dx_ds(short dist_param, short u_type, Real x, Real z) const;
  Real dz_ds_factor(short u_type, Real x, Real z) const;

  //
  //- Heading: Member functions
  //

  void update(Real alpha, Real beta, Real lwr, Real upr);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real std_pdf(Real x, Real alpha, Real beta);
  static Real std_cdf(Real x, Real alpha, Real beta);
  static Real inverse_std_cdf(Real cdf, Real alpha, Real beta);

  static Real pdf(Real x, Real alpha, Real beta, Real lwr, Real upr);
  static Real cdf(Real x, Real alpha, Real beta, Real lwr, Real upr);

  static void moments_from_params(Real alpha, Real beta, Real lwr, Real upr,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Member functions
  //

  /// create a new betaDist instance
  void update_boost();

  //
  //- Heading: Data
  //

  /// alpha parameter of beta random variable (statistical PDF
  /// convention; differs from Jacobi polynomial convention)
  Real alphaStat;
  /// beta parameter of beta random variable (statistical PDF
  /// convention; differs from Jacobi polynomial convention)
  Real betaStat;

  /// pointer to the Boost beta_distribution instance
  beta_dist* betaDist;
};


inline BetaRandomVariable::BetaRandomVariable():
  UniformRandomVariable(), alphaStat(1.), betaStat(1.), betaDist(NULL)
{ ranVarType = BETA; }


inline BetaRandomVariable::
BetaRandomVariable(Real alpha, Real beta, Real lwr, Real upr):
  UniformRandomVariable(lwr, upr), alphaStat(alpha), betaStat(beta),
  betaDist(new beta_dist(alphaStat, betaStat))
{ ranVarType = BETA; }


inline BetaRandomVariable::~BetaRandomVariable()
{ if (betaDist) delete betaDist; }


inline Real BetaRandomVariable::cdf(Real x) const
{
  //return cdf(x, alphaStat, betaStat, lowerBnd, upperBnd);

  Real std_x = (x - lowerBnd)/(upperBnd - lowerBnd);// scale from [L,U] to [0,1]
  return bmth::cdf(*betaDist, std_x);
}


inline Real BetaRandomVariable::ccdf(Real x) const
{
  //return cdf(x, alphaStat, betaStat, lowerBnd, upperBnd);

  Real std_x = (x - lowerBnd)/(upperBnd - lowerBnd);// scale from [L,U] to [0,1]
  return bmth::cdf(complement(*betaDist, std_x));
}


inline Real BetaRandomVariable::inverse_cdf(Real p_cdf) const
{
  //Real std_x = inverse_std_cdf(p, alphaStat, betaStat); // [0,1]

  Real std_x = bmth::quantile(*betaDist, p_cdf);
  return lowerBnd + (upperBnd - lowerBnd) * std_x;  // scale from [0,1] to [L,U]
}


inline Real BetaRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  //Real std_x = inverse_std_cdf(p, alphaStat, betaStat); // [0,1]

  Real std_x = bmth::quantile(complement(*betaDist, p_ccdf));
  return lowerBnd + (upperBnd - lowerBnd) * std_x;  // scale from [0,1] to [L,U]
}


//  F(x) = Boost
//  f(x) = (x-L)^{alpha-1} (U-x)^{beta-1}
//       / Beta(alpha,beta) / (U-L)^{alpha+beta-1}
// f'(x) = f(x) ((alpha-1)/(x-lwr) - (beta-1)/(upr-x))
inline Real BetaRandomVariable::pdf(Real x) const
{
  //return pdf(x, alphaStat, betaStat, lowerBnd, upperBnd);

  Real range = upperBnd - lowerBnd,
    scaled_x = (x - lowerBnd)/range; // from [L,U] to [0,1]
  return bmth::pdf(*betaDist, scaled_x) / range;
}


inline Real BetaRandomVariable::pdf_gradient(Real x) const
{
  // beta PDF is defined on an open interval (l,u)
  // (see also boost/math/special_functions/beta.hpp:ibeta_derivative_imp())
  if (x <= lowerBnd) { // eval @lowerBnd as approached from above
    if      (alphaStat > 1.) return  std::numeric_limits<Real>::quiet_NaN();
    else if (alphaStat < 1.) return -std::numeric_limits<Real>::infinity();
    else // x=L, alphaStat=1: drop (x-L)^{alpha-1}
      return pdf(x) * (1.-betaStat) / (upperBnd-x);
  }
  else if (x >= upperBnd) { // eval @upperBnd as approached from below
    if      (betaStat > 1.) return  std::numeric_limits<Real>::quiet_NaN();
    else if (betaStat < 1.) return  std::numeric_limits<Real>::infinity();
    else // x=U, betaStat=1: drop (U-x)^{beta-1}
      return pdf(x) * (alphaStat-1.) / (x-lowerBnd);
  }
  else
    return pdf(x) *
      ( (alphaStat-1.) / (x-lowerBnd) - (betaStat-1.) / (upperBnd-x) );
}


inline Real BetaRandomVariable::pdf_hessian(Real x) const
{
  // beta PDF is defined on an open interval (l,u)
  // (see also boost/math/special_functions/beta.hpp:ibeta_derivative_imp())
  if (x <= lowerBnd) { // eval @lowerBnd as approached from above
    if      (alphaStat > 1.) return std::numeric_limits<Real>::quiet_NaN();
    else if (alphaStat < 1.) return std::numeric_limits<Real>::infinity();
    else { // x=L, alphaStat=1: drop (x-L)^{alpha-1}
      Real umx = upperBnd - x,  bm1 = betaStat - 1., term =  -bm1 / umx;
      return pdf(x) * (term * term - bm1/(umx*umx));
    }
  }
  else if (x >= upperBnd) { // eval @upperBnd as approached from below
    if      (betaStat > 1.) return std::numeric_limits<Real>::quiet_NaN();
    else if (betaStat < 1.) return std::numeric_limits<Real>::infinity();
    else { // x=U, betaStat=1: drop (U-x)^{beta-1}
      Real xml = x - lowerBnd, am1 = alphaStat - 1., term = am1 / xml;
      return pdf(x) * (term * term - am1/(xml*xml));
    }
  }
  else {
    Real umx = upperBnd - x,  xml = x - lowerBnd, am1 = alphaStat - 1.,
         bm1 = betaStat - 1., term = am1 / xml - bm1 / umx;
    return pdf(x) * (term * term - bm1/(umx*umx) - am1/(xml*xml));
  }
}


inline Real BetaRandomVariable::log_pdf(Real x) const
{
  // beta PDF is defined on an open interval (l,u)
  // (see also boost/math/special_functions/beta.hpp:ibeta_derivative_imp())
  if (x <= lowerBnd) { // eval @lowerBnd as approached from above
    if      (alphaStat > 1.) return -std::numeric_limits<Real>::infinity();
    else if (alphaStat < 1.) return  std::numeric_limits<Real>::infinity();
    else // x=L, alphaStat=1: drop (x-L)^{alpha-1} and combine log(range)
      return -std::log(upperBnd-lowerBnd)
	- std::log(bmth::beta(alphaStat,betaStat));
  }
  else if (x >= upperBnd) { // eval @upperBnd as approached from below
    if      (betaStat > 1.) return -std::numeric_limits<Real>::infinity();
    else if (betaStat < 1.) return  std::numeric_limits<Real>::infinity();
    else // x=U, betaStat=1: drop (U-x)^{beta-1} and combine log(range)
      return -std::log(upperBnd-lowerBnd)
	- std::log(bmth::beta(alphaStat,betaStat));
  }
  else
    return (alphaStat-1.)*std::log(x-lowerBnd)
      + (betaStat-1.)*std::log(upperBnd-x)
      - (alphaStat+betaStat-1.)*std::log(upperBnd-lowerBnd)
      - std::log(bmth::beta(alphaStat,betaStat));
}


inline Real BetaRandomVariable::log_pdf_gradient(Real x) const
{
  // beta PDF is defined on an open interval (l,u).
  // (see also boost/math/special_functions/beta.hpp:ibeta_derivative_imp())
  if (x <= lowerBnd) { // eval @lowerBnd as approached from above
    if      (alphaStat > 1.) return  std::numeric_limits<Real>::infinity();
    else if (alphaStat < 1.) return -std::numeric_limits<Real>::infinity();
    else // x=L, alphaStat=1: drop (x-L)^{alpha-1}
      return (1.-betaStat)/(upperBnd - x);
  }
  else if (x >= upperBnd) { // eval @upperBnd as approached from below
    if      (betaStat > 1.) return -std::numeric_limits<Real>::infinity();
    else if (betaStat < 1.) return  std::numeric_limits<Real>::infinity();
    else // x=U, betaStat=1: drop (U-x)^{beta-1}
      return (alphaStat-1.)/(x - lowerBnd);
  }
  else
    return (alphaStat-1.)/(x - lowerBnd) + (1.-betaStat)/(upperBnd - x);
}


inline Real BetaRandomVariable::log_pdf_hessian(Real x) const
{
  // beta PDF is defined on an open interval (l,u).
  // (see also boost/math/special_functions/beta.hpp:ibeta_derivative_imp())
  if (x <= lowerBnd) { // eval @lowerBnd as approached from above
    if      (alphaStat > 1.) return -std::numeric_limits<Real>::infinity();
    else if (alphaStat < 1.) return  std::numeric_limits<Real>::infinity();
    else // x=L, alphaStat=1: drop (x-L)^{alpha-1}
      { Real umx = upperBnd-x; return (1.-betaStat)/(umx*umx); }
  }
  else if (x >= upperBnd) { // eval @upperBnd as approached from below
    if      (betaStat > 1.) return -std::numeric_limits<Real>::infinity();
    else if (betaStat < 1.) return  std::numeric_limits<Real>::infinity();
    else // x=U, betaStat=1: drop (U-x)^{beta-1}
      { Real xml = x-lowerBnd; return (1.-alphaStat)/(xml*xml); }
  }
  else {
    Real umx = upperBnd - x,  xml = x - lowerBnd;
    return (1.-alphaStat)/(xml*xml) + (1.-betaStat)/(umx*umx);
  }
}


inline Real BetaRandomVariable::inverse_standard_cdf(Real p_cdf) const
{
  Real scaled_z = bmth::quantile(*betaDist, p_cdf); // [0,1]
  return 2.*scaled_z - 1.; // [0,1] to [-1,1]
}


inline Real BetaRandomVariable::standard_pdf(Real z) const
{
  Real scaled_z = (z + 1.)/2.; // [-1,1] to [0,1]
  return bmth::pdf(*betaDist, scaled_z) / 2.;
}


inline Real BetaRandomVariable::log_standard_pdf(Real z) const
{
  // beta PDF is defined on an open interval (-1,1)
  // (see also boost/math/special_functions/beta.hpp:ibeta_derivative_imp())
  if (z <= -1.) { // eval @lowerBnd as approached from above
    if      (alphaStat > 1.) return -std::numeric_limits<Real>::infinity();
    else if (alphaStat < 1.) return  std::numeric_limits<Real>::infinity();
    else // z=L, alphaStat=1: drop (z-L)^{alpha-1} and combine log(range)
      return -std::log(2.) - std::log(bmth::beta(alphaStat,betaStat));
  }
  else if (z >= 1.) { // eval @upperBnd as approached from below
    if      (betaStat > 1.) return -std::numeric_limits<Real>::infinity();
    else if (betaStat < 1.) return  std::numeric_limits<Real>::infinity();
    else // z=U, betaStat=1: drop (U-z)^{beta-1} and combine log(range)
      return -std::log(2.) - std::log(bmth::beta(alphaStat,betaStat));
  }
  else
    return (alphaStat-1.)*bmth::log1p(z) + (betaStat-1.)*bmth::log1p(-z)
      - (alphaStat+betaStat-1.)*std::log(2.)
      - std::log(bmth::beta(alphaStat,betaStat));
}


inline Real BetaRandomVariable::log_standard_pdf_gradient(Real z) const
{
  // beta PDF is defined on an open interval (l,u).
  // (see also boost/math/special_functions/beta.hpp:ibeta_derivative_imp())
  if (z <= -1.) { // eval @lowerBnd as approached from above
    if      (alphaStat > 1.) return  std::numeric_limits<Real>::infinity();
    else if (alphaStat < 1.) return -std::numeric_limits<Real>::infinity();
    else // z=L, alphaStat=1: drop (z-L)^{alpha-1}
      return (1.-betaStat)/(1. - z);
  }
  else if (z >= 1.) { // eval @upperBnd as approached from below
    if      (betaStat > 1.) return -std::numeric_limits<Real>::infinity();
    else if (betaStat < 1.) return  std::numeric_limits<Real>::infinity();
    else // z=U, betaStat=1: drop (U-z)^{beta-1}
      return (alphaStat-1.)/(z + 1.);
  }
  else
    return (alphaStat-1.)/(z + 1.) + (1.-betaStat)/(1. - z);
}


inline Real BetaRandomVariable::log_standard_pdf_hessian(Real z) const
{
  // beta PDF is defined on an open interval (l,u).
  // (see also boost/math/special_functions/beta.hpp:ibeta_derivative_imp())
  if (z <= -1.) { // eval @lowerBnd as approached from above
    if      (alphaStat > 1.) return -std::numeric_limits<Real>::infinity();
    else if (alphaStat < 1.) return  std::numeric_limits<Real>::infinity();
    else // x=L, alphaStat=1: drop (x-L)^{alpha-1}
      { Real umz = 1. - z; return (1.-betaStat)/(umz*umz); }
  }
  else if (z >= 1.) { // eval @upperBnd as approached from below
    if      (betaStat > 1.) return -std::numeric_limits<Real>::infinity();
    else if (betaStat < 1.) return  std::numeric_limits<Real>::infinity();
    else // x=U, betaStat=1: drop (U-x)^{beta-1}
      { Real zml = z + 1.; return (1.-alphaStat)/(zml*zml); }
  }
  else {
    Real umz = 1. - z, zml = z + 1.;
    return (1.-alphaStat)/(zml*zml) + (1.-betaStat)/(umz*umz);
  }
}


inline Real BetaRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case BE_ALPHA:   return alphaStat; break;
  case BE_BETA:    return betaStat;  break;
  case BE_LWR_BND: return lowerBnd;  break;
  case BE_UPR_BND: return upperBnd;  break;
  //case BE_LOCATION: - TO DO
  //case BE_SCALE:    - TO DO
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in BetaRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void BetaRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case BE_ALPHA:   alphaStat = val; update_boost(); break;
  case BE_BETA:    betaStat  = val; update_boost(); break;
  case BE_LWR_BND: lowerBnd  = val; break;
  case BE_UPR_BND: upperBnd  = val; break;
  //case BE_LOCATION: - TO DO
  //case BE_SCALE:    - TO DO
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in BetaRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
}


inline Real BetaRandomVariable::mean() const
{ return lowerBnd + bmth::mean(*betaDist) * (upperBnd - lowerBnd); }


inline Real BetaRandomVariable::median() const
{ return lowerBnd + bmth::median(*betaDist) * (upperBnd - lowerBnd); }


inline Real BetaRandomVariable::mode() const
{ return lowerBnd + bmth::mode(*betaDist) * (upperBnd - lowerBnd); }


inline Real BetaRandomVariable::standard_deviation() const
{ return (upperBnd - lowerBnd) * bmth::standard_deviation(*betaDist); }


inline Real BetaRandomVariable::variance() const
{
  Real range = upperBnd - lowerBnd;
  return range * range * bmth::variance(*betaDist);
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  dz/ds is zero if uncorrelated, 
    while dz_ds_factor() manages contributions in the correlated case. */
inline Real BetaRandomVariable::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  bool u_type_err = false, dist_err = false;
  switch (u_type) {
  case STD_BETA:
    switch (dist_param) { // x = lwr + (upr - lwr)*(z+1.)/2.
    // For distributions without simple closed-form CDFs (beta, gamma), dx/ds
    // is computed numerically in NatafTransformation::jacobian_dX_dS():
    //case BE_ALPHA: case BE_BETA:
    case BE_LWR_BND: return (1. - z)/2.; break;
    case BE_UPR_BND: return (z + 1.)/2.; break;
    //case BE_LOCATION: - TO DO
    //case BE_SCALE:    - TO DO
    default: dist_err = true;                    break;
    }
    break;
  default:
    u_type_err = true; break;
  }

  if (u_type_err)
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in BetaRandomVariable::dx_ds()." << std::endl;
  if (dist_err)
    PCerr << "Error: mapping failure for distribution parameter " << dist_param
	  << " in BetaRandomVariable::dx_ds()." << std::endl;
  if (u_type_err || dist_err)
    abort_handler(-1);
  return 0.;
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  For the uncorrelated case,
    u and z are constants.  For the correlated case, u is a constant, but 
    z(s) = L(s) u due to Nataf dependence on s and dz/ds = dL/ds u. */
inline Real BetaRandomVariable::dz_ds_factor(short u_type, Real x, Real z) const
{
  switch (u_type) {
  //case STD_NORMAL:
  //case STD_UNIFORM:
  case STD_BETA: // x = lwr + (upr - lwr)*(z+1.)/2.
    // --> add (upr - lwr)/2. * dz/ds for nonzero dz/ds arising from correlation
    return (upperBnd-lowerBnd)/2.; break;
  default:
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in BetaRandomVariable::dz_ds_factor()." << std::endl;
    abort_handler(-1);
  }
}


inline void BetaRandomVariable::update_boost()
{
  if (betaDist) delete betaDist;
  betaDist = new beta_dist(alphaStat, betaStat);
}


inline void BetaRandomVariable::
update(Real alpha, Real beta, Real lwr, Real upr)
{
  lowerBnd = lwr; upperBnd = upr; // don't affect betaDist
  if (!betaDist || alphaStat != alpha || betaStat != beta)
    { alphaStat = alpha; betaStat = beta; update_boost(); }
}


inline Real BetaRandomVariable::std_pdf(Real z, Real alpha, Real beta)
{
  beta_dist beta1(alpha, beta);
  Real scaled_z = (z + 1.)/2.; // [-1,1] to [0,1]
  return bmth::pdf(beta1, scaled_z) / 2.;
}


inline Real BetaRandomVariable::std_cdf(Real z, Real alpha, Real beta)
{
  beta_dist beta1(alpha, beta);
  Real scaled_z = (z + 1.)/2.; // [-1,1] to [0,1]
  return bmth::cdf(beta1, scaled_z);
}


inline Real BetaRandomVariable::inverse_std_cdf(Real cdf, Real alpha, Real beta)
{
  beta_dist beta1(alpha, beta);
  return 2.*bmth::quantile(beta1, cdf) - 1.; // [0,1] to [-1,1]
}


inline Real BetaRandomVariable::
pdf(Real x, Real alpha, Real beta, Real lwr, Real upr)
{
  beta_dist beta1(alpha, beta);
  Real range = upr - lwr, scaled_x = (x - lwr)/range; // from [L,U] to [0,1]
  return bmth::pdf(beta1, scaled_x) / range;
}


inline Real BetaRandomVariable::
cdf(Real x, Real alpha, Real beta, Real lwr, Real upr)
{
  beta_dist beta1(alpha, beta);
  Real scaled_x = (x - lwr)/(upr - lwr); // from [L,U] to [0,1]
  return bmth::cdf(beta1, scaled_x);
}


inline void BetaRandomVariable::
moments_from_params(Real alpha, Real beta, Real lwr, Real upr,
		    Real& mean, Real& std_dev)
{
  Real range = upr - lwr;
  mean       = lwr + alpha/(alpha+beta)*range;
  std_dev    = std::sqrt(alpha*beta/(alpha+beta+1.))/(alpha+beta)*range;
}

} // namespace Pecos

#endif
