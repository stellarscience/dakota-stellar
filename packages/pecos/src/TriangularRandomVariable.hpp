/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 TriangularRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef TRIANGULAR_RANDOM_VARIABLE_HPP
#define TRIANGULAR_RANDOM_VARIABLE_HPP

#include "UniformRandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for triangular random variables.

/** Manages mode and inherits bounds.  See Haldar and Mahadevan, p. 99. */

class TriangularRandomVariable: public UniformRandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  TriangularRandomVariable();
  /// alternate constructor
  TriangularRandomVariable(Real lwr, Real mode, Real upr);
  /// destructor
  ~TriangularRandomVariable();

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

  void update(Real lwr, Real mode, Real upr);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real lwr, Real mode, Real upr);
  static Real cdf(Real x, Real lwr, Real mode, Real upr);

  static void moments_from_params(Real lwr, Real mode, Real upr,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Member functions
  //

  /// create a new triangDist instance
  void update_boost();

  //
  //- Heading: Data
  //

  /// mode of triangular random variable
  Real triangularMode;

  /// pointer to the Boost gamma_distribution instance
  triangular_dist* triangDist;
};


inline TriangularRandomVariable::TriangularRandomVariable():
  UniformRandomVariable(), triangularMode(0), triangDist(NULL)
{ ranVarType = TRIANGULAR; }


inline TriangularRandomVariable::
TriangularRandomVariable(Real lwr, Real mode, Real upr):
  UniformRandomVariable(lwr, upr), triangularMode(mode),
  triangDist(new triangular_dist(lwr, mode, upr))
{ ranVarType = TRIANGULAR; }


inline TriangularRandomVariable::~TriangularRandomVariable()
{ if (triangDist) delete triangDist; }


inline Real TriangularRandomVariable::pdf(Real x, Real lwr, Real mode, Real upr)
{
  triangular_dist tri1(lwr, mode, upr);
  return bmth::pdf(tri1, x);

  //return (x < mode) ? 2.*(x-lwr)/(upr-lwr)/(mode-lwr) :
  //                    2.*(upr-x)/(upr-lwr)/(upr-mode);
}


inline Real TriangularRandomVariable::cdf(Real x, Real lwr, Real mode, Real upr)
{
  triangular_dist tri1(lwr, mode, upr);
  return bmth::cdf(tri1, x);

  //return (x < mode) ? std::pow(x-lwr,2.)/(upr-lwr)/(mode-lwr) :
  //  ((mode-lwr) - (x+mode-2*upr)*(x-mode)/(upr-mode))/(upr-lwr);
}


inline Real TriangularRandomVariable::cdf(Real x) const
{ return bmth::cdf(*triangDist, x); }


inline Real TriangularRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*triangDist, x)); }


inline Real TriangularRandomVariable::inverse_cdf(Real p_cdf) const
{
  return bmth::quantile(*triangDist, p_cdf);

  /*
  // assume x < mode and then check
  Real range = upperBnd - lowerBnd,
       x = lowerBnd + std::sqrt(p*range*(triangularMode-lowerBnd));
  Real x_pdf = 2.*(x-lowerBnd)/range/(triangularMode-lowerBnd),
       m_pdf = 2./range;
  // check pdf value to ensure that proper equation used
  if ( x_pdf > m_pdf )
    x = upperBnd - std::sqrt((1.-p)*range*(upperBnd-triangularMode));
  return x;
  */
}


inline Real TriangularRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*triangDist, p_ccdf)); }


//             x < M                        x > M
//  F(x): (x-L)^2/(U-L)/(M-L)    (M-L)/(U-L) - (x+M-2U)(x-M)/(U-L)/(U-M)
//  f(x): 2(x-L)/(U-L)/(M-L)     2(U-x)/(U-L)/(U-M)
// f'(x): 2/(U-L)/(M-L)          -2/(U-L)/(U-M)
// Note: at x=M, F(x) and f(x) are continuous but f'(x) is not
inline Real TriangularRandomVariable::pdf(Real x) const
{ return bmth::pdf(*triangDist, x); }


inline Real TriangularRandomVariable::pdf_gradient(Real x) const
{
  Real range = upperBnd - lowerBnd;
  if (x < triangularMode)
    return  2. / ( range * (triangularMode - lowerBnd) );
  else if (x > triangularMode)
    return -2. / ( range * (upperBnd - triangularMode) );
  else // x == triangularMode
    return  0.; // f'(x) is undefined: use 0.
}


inline Real TriangularRandomVariable::pdf_hessian(Real x) const
{ return 0.; }


inline Real TriangularRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case T_MODE:    return triangularMode; break;
  case T_LWR_BND: return lowerBnd;       break;
  case T_UPR_BND: return upperBnd;       break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in TriangularRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void TriangularRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case T_MODE:    triangularMode = val; break;
  case T_LWR_BND: lowerBnd       = val; break;
  case T_UPR_BND: upperBnd       = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in TriangularRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
  update_boost(); // create a new triangDist instance
}


inline Real TriangularRandomVariable::mean() const
{ return bmth::mean(*triangDist); }


inline Real TriangularRandomVariable::median() const
{ return bmth::median(*triangDist); }


inline Real TriangularRandomVariable::mode() const
{ return bmth::mode(*triangDist); }


inline Real TriangularRandomVariable::standard_deviation() const
{ return bmth::standard_deviation(*triangDist); }


inline Real TriangularRandomVariable::variance() const
{ return bmth::variance(*triangDist); }


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  dz/ds is zero if uncorrelated, 
    while dz_ds_factor() manages contributions in the correlated case. */
inline Real TriangularRandomVariable::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  bool u_type_err = false, dist_error = false; Real sdf;
  if (x < triangularMode)
    switch (u_type) {
    case STD_NORMAL:  sdf =  NormalRandomVariable::std_cdf(z);  break;
    case STD_UNIFORM: sdf = UniformRandomVariable::std_cdf(z);  break;
    //case TRIANGULAR: break;
    default: u_type_err = true; break;
    }
  else
    switch (u_type) {
    case STD_NORMAL:  sdf =  NormalRandomVariable::std_ccdf(z); break;
    case STD_UNIFORM: sdf = UniformRandomVariable::std_ccdf(z); break;
    //case TRIANGULAR: break;
    default: u_type_err = true; break;
    }
  if (u_type_err) {
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in TriangularRandomVariable::dx_ds()." << std::endl;
    abort_handler(-1);
  }

  if (x < triangularMode) {
    Real denom = 2.*(x-lowerBnd);
    switch (dist_param) {
    case T_MODE:    return sdf*(upperBnd-lowerBnd)/denom;        break;
    case T_LWR_BND:
      return 1.+sdf*(2.*lowerBnd-upperBnd-triangularMode)/denom; break;
    case T_UPR_BND: return sdf*(triangularMode-lowerBnd)/denom;  break;
    // Triangular Mean          - TO DO
    // Triangular Std Deviation - TO DO
    default:        dist_error = true;                           break;
    }
  }
  else {
    Real denom = 2.*(upperBnd-x);
    switch (dist_param) {
    case T_MODE:    return sdf*(upperBnd-lowerBnd)/denom;        break;
    case T_LWR_BND: return (upperBnd-triangularMode)*sdf/denom;  break;
    case T_UPR_BND:
      return 1.-sdf*(2.*upperBnd-lowerBnd-triangularMode)/denom; break;
    // Triangular Mean          - TO DO
    // Triangular Std Deviation - TO DO
    default:        dist_error = true;                           break;
    }
  }
  if (dist_error) {
    PCerr << "Error: mapping failure for distribution parameter " << dist_param
	  << " in TriangularRandomVariable::dx_ds()." << std::endl;
    abort_handler(-1); return 0.;
  }
}


/** dx/ds is derived by differentiating NatafTransformation::trans_Z_to_X()
    with respect to distribution parameter s.  For the uncorrelated case,
    u and z are constants.  For the correlated case, u is a constant, but 
    z(s) = L(s) u due to Nataf dependence on s and dz/ds = dL/ds u. */
inline Real TriangularRandomVariable::
dz_ds_factor(short u_type, Real x, Real z) const
{
  Real pdf;
  switch (u_type) {
  case STD_NORMAL:  pdf = NormalRandomVariable::std_pdf(z); break;
  case STD_UNIFORM: pdf = UniformRandomVariable::std_pdf(); break;
  //case TRIANGULAR: break;
  default:
    PCerr << "Error: unsupported u-space type " << u_type
	  << " in TriangularRandomVariable::dz_ds_factor()." << std::endl;
    abort_handler(-1); break;
  }

  return (x < triangularMode) ?
    (upperBnd-lowerBnd)*(triangularMode-lowerBnd)*pdf/(2.*(x-lowerBnd)) :
    (upperBnd-triangularMode)*(upperBnd-lowerBnd)*pdf/(2.*(upperBnd-x));
}


inline void TriangularRandomVariable::update_boost()
{
  if (triangDist) delete triangDist;
  triangDist = new triangular_dist(lowerBnd, triangularMode, upperBnd);
}


inline void TriangularRandomVariable::update(Real lwr, Real mode, Real upr)
{
  if (!triangDist ||
      lowerBnd != lwr || triangularMode != mode || upperBnd != upr)
    { lowerBnd = lwr; triangularMode = mode; upperBnd = upr; update_boost(); }
}


inline void TriangularRandomVariable::
moments_from_params(Real lwr, Real mode, Real upr, Real& mean, Real& std_dev)
{
  mean = (lwr + mode + upr)/3.;
  std_dev
    = std::sqrt((lwr*(lwr - mode) + mode*(mode - upr) + upr*(upr - lwr))/18.);
}

} // namespace Pecos

#endif
