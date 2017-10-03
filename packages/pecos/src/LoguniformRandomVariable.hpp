/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 LoguniformRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef LOGUNIFORM_RANDOM_VARIABLE_HPP
#define LOGUNIFORM_RANDOM_VARIABLE_HPP

#include "UniformRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for loguniform random variables.

/** Manages lower and upper bounds.  See SAND98-0210 LHS manual, pp. 43-44. */

class LoguniformRandomVariable: public UniformRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  LoguniformRandomVariable();                   ///< default constructor
  LoguniformRandomVariable(Real lwr, Real upr); ///< alternate constructor
  ~LoguniformRandomVariable();                  ///< destructor

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

  //void update(Real lwr, Real upr); // inherits from UniformRV

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real lwr, Real upr);
  static Real cdf(Real x, Real lwr, Real upr);

  static void moments_from_params(Real lwr, Real upr, Real& mean,
				  Real& std_dev);

protected:

  //
  //- Heading: Data
  //
};


inline LoguniformRandomVariable::LoguniformRandomVariable():
  UniformRandomVariable()
{ ranVarType = LOGUNIFORM; }


inline LoguniformRandomVariable::LoguniformRandomVariable(Real lwr, Real upr):
  UniformRandomVariable(lwr, upr)
{ ranVarType = LOGUNIFORM; }


inline LoguniformRandomVariable::~LoguniformRandomVariable()
{ }


inline Real LoguniformRandomVariable::cdf(Real x) const
{
  return (std::log(x)        - std::log(lowerBnd)) /
         (std::log(upperBnd) - std::log(lowerBnd));
}


inline Real LoguniformRandomVariable::ccdf(Real x) const
{
  return (std::log(upperBnd) - std::log(x)) /
         (std::log(upperBnd) - std::log(lowerBnd));
}


inline Real LoguniformRandomVariable::inverse_cdf(Real p_cdf) const
{
  // p = (ln x - ln L)/(ln U - ln L)
  return lowerBnd * std::exp(p_cdf * (std::log(upperBnd) - std::log(lowerBnd)));
}


inline Real LoguniformRandomVariable::inverse_ccdf(Real p_ccdf) const
{
  // p = (ln U - ln x)/(ln U - ln L)
  return upperBnd /
    std::exp(p_ccdf * (std::log(upperBnd) - std::log(lowerBnd)));
}


//  F(x) = (ln x - ln L)/(ln U - ln L)
//  f(x) =  1/(ln U - ln L)/x
// f'(x) = -1/(ln U - ln L)/x^2
// f'(x) =  2/(ln U - ln L)/x^3
inline Real LoguniformRandomVariable::pdf(Real x) const
{ return  1./((std::log(upperBnd) - std::log(lowerBnd)) * x); }


inline Real LoguniformRandomVariable::pdf_gradient(Real x) const
{ return -1./((std::log(upperBnd) - std::log(lowerBnd)) * x * x); }


inline Real LoguniformRandomVariable::pdf_hessian(Real x) const
{ return  2./((std::log(upperBnd) - std::log(lowerBnd)) * x * x * x); }


inline Real LoguniformRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case LU_LWR_BND: return lowerBnd; break;
  case LU_UPR_BND: return upperBnd; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in LoguniformRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void LoguniformRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case LU_LWR_BND: lowerBnd = val; break;
  case LU_UPR_BND: upperBnd = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in LoguniformRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
}


inline Real LoguniformRandomVariable::mean() const
{ return (upperBnd - lowerBnd)/(std::log(upperBnd) - std::log(lowerBnd)); }


inline Real LoguniformRandomVariable::median() const
{ return inverse_cdf(.5); }


inline Real LoguniformRandomVariable::mode() const
{ return lowerBnd; }


inline Real LoguniformRandomVariable::standard_deviation() const
{
  Real  range = upperBnd - lowerBnd,
    log_range = std::log(upperBnd) - std::log(lowerBnd);
  return std::sqrt(range*(log_range*(upperBnd+lowerBnd)/2.-range))/log_range;
}


inline Real LoguniformRandomVariable::variance() const
{
  Real  range = upperBnd - lowerBnd,
    log_range = std::log(upperBnd) - std::log(lowerBnd);
  return range*(log_range*(upperBnd+lowerBnd)/2.-range)/(log_range*log_range);
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  dz/ds is zero if uncorrelated, 
    while dz_ds_factor() manages contributions in the correlated case. */
inline Real LoguniformRandomVariable::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  // to STD_NORMAL:  ln x = ln L +  Phi(z) (ln U - ln L)
  // to STD_UNIFORM: ln x = ln L + (z+1)/2 (ln U - ln L)
  bool u_type_err = false, dist_err = false;
  switch (dist_param) {
  case LU_LWR_BND: // Deriv of Loguniform w.r.t. its Lower Bound
    switch (u_type) {
    case STD_UNIFORM:
      return x*UniformRandomVariable::std_ccdf(z)/lowerBnd; break;
    case STD_NORMAL:
      return x* NormalRandomVariable::std_ccdf(z)/lowerBnd; break;
    //case LOGUNIFORM:  TO DO; break;
    default: u_type_err = true;                             break;
    }
    break;
  case LU_UPR_BND: // Deriv of Loguniform w.r.t. its Upper Bound
    switch (u_type) {
    case STD_NORMAL:
      return x* NormalRandomVariable::std_cdf(z)/upperBnd;  break;
    case STD_UNIFORM:
      return x*UniformRandomVariable::std_cdf(z)/upperBnd;  break;
    //case LOGUNIFORM:  TO DO; break;
    default: u_type_err = true;                             break;
    }
    break;
  default:   dist_err = true;                               break;
  }

  if (u_type_err)
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in LoguniformRandomVariable::dx_ds()." << std::endl;
  if (dist_err)
    PCerr << "Error: mapping failure for distribution parameter " << dist_param
	  << " in LoguniformRandomVariable::dx_ds()." << std::endl;
  if (u_type_err || dist_err)
    abort_handler(-1);
  return 0.;
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  For the uncorrelated case,
    u and z are constants.  For the correlated case, u is a constant, but 
    z(s) = L(s) u due to Nataf dependence on s and dz/ds = dL/ds u. */
inline Real LoguniformRandomVariable::
dz_ds_factor(short u_type, Real x, Real z) const
{
  // to STD_NORMAL:  ln x = ln L +  Phi(z) (ln U - ln L)
  // to STD_UNIFORM: ln x = ln L + (z+1)/2 (ln U - ln L)
  Real log_range = std::log(upperBnd)-std::log(lowerBnd);
  switch (u_type) {
  case STD_NORMAL:  return x*log_range*NormalRandomVariable::std_pdf(z); break;
  case STD_UNIFORM: return x*log_range*UniformRandomVariable::std_pdf(); break;
  //case LOGUNIFORM:  TO DO; break;
  default:
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in LoguniformRandomVariable::dz_ds_factor()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}

// static functions:

inline Real LoguniformRandomVariable::pdf(Real x, Real lwr, Real upr)
{ return 1./(std::log(upr) - std::log(lwr))/x; }


inline Real LoguniformRandomVariable::cdf(Real x, Real lwr, Real upr)
{ return (std::log(x) - std::log(lwr))/(std::log(upr) - std::log(lwr)); }


inline void LoguniformRandomVariable::
moments_from_params(Real lwr, Real upr, Real& mean, Real& std_dev)
{
  Real range = upr - lwr, log_range = std::log(upr) - std::log(lwr);
  mean       = range/log_range;
  std_dev    = std::sqrt(range*(log_range*(upr+lwr)/2.-range))/log_range;
}

} // namespace Pecos

#endif
