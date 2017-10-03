/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 UniformRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef UNIFORM_RANDOM_VARIABLE_HPP
#define UNIFORM_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"
#include "NormalRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for uniform random variables.

/** Manages lower and upper bounds. */

class UniformRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  UniformRandomVariable();                   ///< default constructor
  UniformRandomVariable(Real lwr, Real upr); ///< alternate constructor
  ~UniformRandomVariable();                  ///< destructor

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
  Real log_pdf_gradient(Real x) const;
  Real log_pdf_hessian(Real x) const;

  Real inverse_standard_cdf(Real p_cdf) const;

  Real standard_pdf(Real z) const;
  Real log_standard_pdf_gradient(Real z) const;
  Real log_standard_pdf_hessian(Real x) const;

  Real to_standard(Real x) const;
  Real from_standard(Real z) const;

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

  void update(Real lwr, Real upr);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real lwr, Real upr);
  static Real cdf(Real x, Real lwr, Real upr);
  static Real inverse_cdf(Real p_cdf, Real lwr, Real upr);

  static Real std_pdf();
  static Real log_std_pdf();
  static Real log_std_pdf_gradient();
  static Real log_std_pdf_hessian();

  static Real std_cdf(Real beta);
  static Real std_ccdf(Real beta);
  static Real inverse_std_cdf(Real p_cdf);
  static Real inverse_std_ccdf(Real p_ccdf);

  template <typename Engine> 
  static Real draw_std_sample(Engine& rng);

  static void moments_from_params(Real lwr, Real upr, Real& mean,
				  Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// lower bound of uniform random variable
  Real lowerBnd;
  /// upper bound of uniform random variable
  Real upperBnd;
};


inline UniformRandomVariable::UniformRandomVariable():
  RandomVariable(BaseConstructor()), lowerBnd(-1.), upperBnd(1.)
{ ranVarType = UNIFORM; }


inline UniformRandomVariable::UniformRandomVariable(Real lwr, Real upr):
  RandomVariable(BaseConstructor()), lowerBnd(lwr), upperBnd(upr)
{ ranVarType = UNIFORM; }


inline UniformRandomVariable::~UniformRandomVariable()
{ }


inline Real UniformRandomVariable::std_pdf()
{ return 0.5; } // equal probability on [-1,1]


inline Real UniformRandomVariable::log_std_pdf()
{ return std::log(std_pdf()); }


inline Real UniformRandomVariable::log_std_pdf_gradient()
{ return 0.; }


inline Real UniformRandomVariable::log_std_pdf_hessian()
{ return 0.; }


inline Real UniformRandomVariable::std_cdf(Real x)
{
  if      (x >=  1.) return 1.;
  else if (x <= -1.) return 0.;
  else               return (x + 1.)/2.; // linear x \in [-1,1] -> p \in [0,1]
}


inline Real UniformRandomVariable::std_ccdf(Real x)
{
  if      (x >=  1.) return 0.;
  else if (x <= -1.) return 1.;
  else               return (1. - x)/2.; // linear x \in [-1,1] -> p \in [1,0]
}


inline Real UniformRandomVariable::inverse_std_cdf(Real p_cdf)
{
  if      (p_cdf >= 1.) return  1.;
  else if (p_cdf <= 0.) return -1.;
  else return 2.*p_cdf - 1.; // linear p \in [0,1] -> x \in [-1,1]
}


inline Real UniformRandomVariable::inverse_std_ccdf(Real p_ccdf)
{
  if      (p_ccdf >= 1.) return -1.;
  else if (p_ccdf <= 0.) return  1.;
  else return 1. - 2.*p_ccdf; // linear p \in [1,0] -> x \in [-1,1]
}


inline Real UniformRandomVariable::pdf(Real lwr, Real upr)
{ return 1./(upr - lwr); } // equal probability on [lwr,upr]


inline Real UniformRandomVariable::cdf(Real x, Real lwr, Real upr)
{
  if      (x >= upr) return 1.;
  else if (x <= lwr) return 0.;
  else               return (x - lwr)/(upr - lwr); // linear [l,u] -> [0,1]
}


inline Real UniformRandomVariable::inverse_cdf(Real p_cdf, Real lwr, Real upr)
{
  if      (p_cdf >= 1.) return upr;
  else if (p_cdf <= 0.) return lwr;
  else                  return lwr + (upr - lwr) * p_cdf;
}


inline Real UniformRandomVariable::cdf(Real x) const
{ return cdf(x, lowerBnd, upperBnd); }


inline Real UniformRandomVariable::ccdf(Real x) const
{
  if      (x >= upperBnd) return 0.;
  else if (x <= lowerBnd) return 1.;
  else                    return (upperBnd - x)/(upperBnd - lowerBnd);
}


inline Real UniformRandomVariable::inverse_cdf(Real p_cdf) const
{ return inverse_cdf(p_cdf, lowerBnd, upperBnd); }


inline Real UniformRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  if      (p_ccdf >= 1.) return lowerBnd;
  else if (p_ccdf <= 0.) return upperBnd;
  else                   return upperBnd - (upperBnd - lowerBnd) * p_ccdf;
}


inline Real UniformRandomVariable::pdf(Real x) const
{ return pdf(lowerBnd, upperBnd); }


inline Real UniformRandomVariable::pdf_gradient(Real x) const
{ return 0.; }


inline Real UniformRandomVariable::pdf_hessian(Real x) const
{ return 0.; }


inline Real UniformRandomVariable::log_pdf_gradient(Real x) const
{ return 0.; }


inline Real UniformRandomVariable::log_pdf_hessian(Real x) const
{ return 0.; }


inline Real UniformRandomVariable::inverse_standard_cdf(Real p_cdf) const
{ return inverse_std_cdf(p_cdf); }


inline Real UniformRandomVariable::standard_pdf(Real z) const
{ return std_pdf(); }


inline Real UniformRandomVariable::log_standard_pdf_gradient(Real x) const
{ return 0.; }


inline Real UniformRandomVariable::log_standard_pdf_hessian(Real x) const
{ return 0.; }


inline Real UniformRandomVariable::to_standard(Real x) const
{
  // [L,U] -> [-1,1]
  if      (x >= upperBnd) return  1.;
  else if (x <= lowerBnd) return -1.;
  else                    return 2.*(x - lowerBnd)/(upperBnd - lowerBnd) - 1.;
}


inline Real UniformRandomVariable::from_standard(Real z) const
{
  // [-1,1] -> [L,U]
  if      (z >=  1.) return upperBnd;
  else if (z <= -1.) return lowerBnd;
  else               return lowerBnd + (upperBnd - lowerBnd) * (z + 1.) / 2.;
}


inline Real UniformRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case U_LWR_BND: case CDV_LWR_BND: case CSV_LWR_BND: return lowerBnd; break;
  case U_UPR_BND: case CDV_UPR_BND: case CSV_UPR_BND: return upperBnd; break;
  //case U_LOCATION: - TO DO
  //case U_SCALE:    - TO DO
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in UniformRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void UniformRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case U_LWR_BND: case CDV_LWR_BND: case CSV_LWR_BND: lowerBnd = val; break;
  case U_UPR_BND: case CDV_UPR_BND: case CSV_UPR_BND: upperBnd = val; break;
  //case U_LOCATION: - TO DO
  //case U_SCALE:    - TO DO
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in UniformRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
}


inline Real UniformRandomVariable::mean() const
{ return (lowerBnd + upperBnd)/2.; }


inline Real UniformRandomVariable::median() const
{ return mean(); }


inline Real UniformRandomVariable::mode() const
{ return mean(); } // not well-defined: any value in [L,U] is valid


inline Real UniformRandomVariable::standard_deviation() const
{ return (upperBnd - lowerBnd)/std::sqrt(12.); }


inline Real UniformRandomVariable::variance() const
{
  Real range = upperBnd - lowerBnd;
  return range*range/12.;
}


inline RealRealPair UniformRandomVariable::bounds() const
{ return RealRealPair(lowerBnd, upperBnd); }


inline Real UniformRandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  // correlation warping factor for transformations to STD_NORMAL space
  // Der Kiureghian and Liu, ASCE JEM 112:1, 1986
  Real COV;
  switch (rv.type()) { // x-space types mapped to STD_NORMAL u-space

  // Der Kiureghian & Liu: Table 4
  case UNIFORM:     return 1.047 - 0.047*corr*corr; break; // Max Error 0.0%
  case EXPONENTIAL: return 1.133 + 0.029*corr*corr; break; // Max Error 0.0%
  case GUMBEL:      return 1.055 + 0.015*corr*corr; break; // Max Error 0.0%

  // Der Kiureghian & Liu: Table 5 (quadratic approximations in COV,corr)
  case GAMMA:     // Max Error 0.1%
    COV = rv.coefficient_of_variation();
    return 1.023 + (-0.007 + 0.127*COV)*COV + 0.002*corr*corr; break;
  case FRECHET:   // Max Error 2.1%
    COV = rv.coefficient_of_variation();
    return 1.033 + ( 0.305 + 0.405*COV)*COV + 0.074*corr*corr; break;
  case WEIBULL:   // Max Error 0.5%
    COV = rv.coefficient_of_variation();
    return 1.061 + (-0.237 + 0.379*COV)*COV - 0.005*corr*corr; break;

  // warping factors are defined once for lower triangle based on uv order
  case NORMAL: case LOGNORMAL:
    return rv.correlation_warping_factor(*this, corr); break;

  default: // Unsupported warping (should be prevented upsteam)
    PCerr << "Error: unsupported correlation warping for UniformRV."<<std::endl;
    abort_handler(-1); return 1.; break;
  }
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  dz/ds is zero if uncorrelated, 
    while dz_ds_factor() manages contributions in the correlated case. */
inline Real UniformRandomVariable::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  // to STD_NORMAL:  x = L + Phi(z) (U - L)
  // to STD_UNIFORM: x = L + (z + 1)*(U - L)/2
  bool u_type_err = false, dist_err = false;
  switch (dist_param) {
  case U_LWR_BND: case CDV_LWR_BND: case CSV_LWR_BND:
    // Deriv of Uniform w.r.t. its Lower Bound
    switch (u_type) {
    //case STD_NORMAL:  return NormalRandomVariable::std_ccdf(z); break;
    //case STD_UNIFORM: return std_ccdf(z); break;
    // The following is equivalent with fewer operations:
    case STD_NORMAL: case STD_UNIFORM: return ccdf(x); break;
    default: u_type_err = true;                        break;
    }
    break;
  case U_UPR_BND: case CDV_UPR_BND: case CSV_UPR_BND:
    // Deriv of Uniform w.r.t. its Upper Bound
    switch (u_type) {
    //case STD_NORMAL:  return NormalRandomVariable::std_cdf(z); break;
    //case STD_UNIFORM: return std_cdf(z);                       break;
    // The following is equivalent with fewer operations:
    case STD_NORMAL: case STD_UNIFORM: return cdf(x); break;
    default:          u_type_err = true;                       break;
    }
    break;
  //case U_LOCATION: - TO DO
  //case U_SCALE:    - TO DO
  case NO_TARGET: // can occur for all_variables Jacobians
    if (ranVarType == CONTINUOUS_DESIGN || ranVarType == CONTINUOUS_INTERVAL ||
	ranVarType == CONTINUOUS_STATE)
      return 0.;
    else dist_err = true;
    break;
  default:
    dist_err = true; break;
  }

  if (u_type_err)
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in UniformRandomVariable::dx_ds()." << std::endl;
  if (dist_err)
    PCerr << "Error: mapping failure for distribution parameter " << dist_param
	  << " in UniformRandomVariable::dx_ds()." << std::endl;
  if (u_type_err || dist_err)
    abort_handler(-1);
  return 0.;
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  For the uncorrelated case,
    u and z are constants.  For the correlated case, u is a constant, but 
    z(s) = L(s) u due to Nataf dependence on s and dz/ds = dL/ds u. */
inline Real UniformRandomVariable::
dz_ds_factor(short u_type, Real x, Real z) const
{
  // to STD_NORMAL:  x = L + Phi(z) (U - L)
  // to STD_UNIFORM: x = L + (z + 1)*(U - L)/2
  switch (u_type) {
  case STD_NORMAL:
    return (upperBnd-lowerBnd) * NormalRandomVariable::std_pdf(z); break;
  case STD_UNIFORM:
    return (upperBnd-lowerBnd) * std_pdf();                        break;
  default:
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in UniformRandomVariable::dz_ds_factor()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void UniformRandomVariable::update(Real lwr, Real upr)
{ lowerBnd = lwr; upperBnd = upr; }


template <typename Engine> 
Real UniformRandomVariable::draw_std_sample(Engine& rng)
{
  // draw random number on [0,1] from a persistent RNG sequence
  boost::uniform_real<Real> uniform_sampler;
  Real u01 = uniform_sampler(rng);
  return inverse_std_cdf(u01);
}


inline void UniformRandomVariable::
moments_from_params(Real lwr, Real upr, Real& mean, Real& std_dev)
{ mean = (lwr + upr)/2.; std_dev = (upr - lwr)/std::sqrt(12.); }

} // namespace Pecos

#endif
