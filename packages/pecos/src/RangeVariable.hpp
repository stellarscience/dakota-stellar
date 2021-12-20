/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 RangeVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef RANGE_VARIABLE_HPP
#define RANGE_VARIABLE_HPP

#include "RandomVariable.hpp"
#include "UniformRandomVariable.hpp"

namespace Pecos {


/// Derived RandomVariable class for range variables.

/** This is distinct from UniformRandomVariable in that it is used for
    non-random types without associated probability densities.  As
    such, some statistical operations are suppressed. */

template <typename T>
class RangeVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  RangeVariable();
  /// alternate constructor
  RangeVariable(T lwr, T upr);
  /// destructor
  ~RangeVariable();

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

  Real mean() const;
  //Real median() const;
  Real mode() const;
  Real standard_deviation() const;
  Real variance() const;
  
  RealRealPair moments() const;
  RealRealPair distribution_bounds() const;

  void lower_bound(T l_bnd);
  void upper_bound(T l_bnd);

  //Real coefficient_of_variation() const;

  void pull_parameter(short dist_param, T& val) const;
  void push_parameter(short dist_param, T  val);

  void copy_parameters(const RandomVariable& rv);

  Real correlation_warping_factor(const RandomVariable& rv, Real corr) const;
  Real dx_ds(short dist_param, short u_type, Real x, Real z) const;
  Real dz_ds_factor(short u_type, Real x, Real z) const;

  //
  //- Heading: Member functions
  //

  void update(T lwr, T upr);

  //
  //- Heading: Static member functions (global utilities)
  //

  //static void moments_from_params(T lwr, T upr, Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Member functions
  //

  void no_template_specialization(String fn) const;

  //
  //- Heading: Data
  //

  /// lower bound of range variable
  T lowerBnd;
  /// upper bound of range variable
  T upperBnd;
};


//// GENERIC ////


template <typename T>
RangeVariable<T>::RangeVariable():
  RandomVariable(BaseConstructor())
{ }


template <typename T>
RangeVariable<T>::RangeVariable(T lwr, T upr):
  RandomVariable(BaseConstructor())
{ update(lwr, upr); }


template <typename T>
RangeVariable<T>::~RangeVariable()
{ }


template <typename T>
void RangeVariable<T>::update(T lwr, T upr)
{ lowerBnd = lwr; upperBnd = upr; }
// specializations used for assigning ranVarType, but could also employ
// std::is_same for type identification


template <typename T>
void RangeVariable<T>::pull_parameter(short dist_param, T& val) const
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case CR_LWR_BND: case DR_LWR_BND:  val = lowerBnd;  break;
  case CR_UPR_BND: case DR_UPR_BND:  val = upperBnd;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in RangeVariable::pull_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
void RangeVariable<T>::push_parameter(short dist_param, T val)
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case CR_LWR_BND: case DR_LWR_BND:  lowerBnd = val;  break;
  case CR_UPR_BND: case DR_UPR_BND:  upperBnd = val;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in RangeVariable::push_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
void RangeVariable<T>::copy_parameters(const RandomVariable& rv)
{
  switch (ranVarType) {
  case CONTINUOUS_RANGE:
    rv.pull_parameter(CR_LWR_BND, lowerBnd);
    rv.pull_parameter(CR_UPR_BND, upperBnd); break;
  case DISCRETE_RANGE:
    rv.pull_parameter(DR_LWR_BND, lowerBnd);
    rv.pull_parameter(DR_UPR_BND, upperBnd); break;
  }
}


template <typename T>
void RangeVariable<T>::lower_bound(T l_bnd)
{ lowerBnd = l_bnd; }


template <typename T>
void RangeVariable<T>::upper_bound(T u_bnd)
{ upperBnd = u_bnd; }


template <typename T>
void RangeVariable<T>::no_template_specialization(String fn) const
{
  PCerr << "Error: no template specialization of " << fn << "() for "
	<< "RangeVariable<T>." << std::endl;
  abort_handler(-1);
}


template <typename T>
Real RangeVariable<T>::pdf(Real x) const
{ no_template_specialization("pdf"); return 0.; }


template <typename T>
Real RangeVariable<T>::pdf_gradient(Real x) const
{ return 0.; }


template <typename T>
Real RangeVariable<T>::pdf_hessian(Real x) const
{ return 0.; }


template <typename T>
Real RangeVariable<T>::log_pdf_gradient(Real x) const
{ return 0.; }


template <typename T>
Real RangeVariable<T>::log_pdf_hessian(Real x) const
{ return 0.; }


template <typename T>
Real RangeVariable<T>::cdf(Real x) const
{ no_template_specialization("cdf"); return 0.; }


template <typename T>
Real RangeVariable<T>::ccdf(Real x) const
{ no_template_specialization("ccdf"); return 0.; }


template <typename T>
Real RangeVariable<T>::inverse_cdf(Real p_cdf) const
{ no_template_specialization("inverse_cdf"); return 0.; }


template <typename T>
Real RangeVariable<T>::inverse_ccdf(Real p_ccdf) const
{ no_template_specialization("inverse_ccdf"); return 0.; }


template <typename T>
RealRealPair RangeVariable<T>::moments() const
{ no_template_specialization("moments"); return RealRealPair(); }


template <typename T>
Real RangeVariable<T>::mean() const
{ return moments().first; }


//template <typename T>
//Real RangeVariable<T>::median() const
//{ return inverse_cdf(.5); }


template <typename T>
Real RangeVariable<T>::standard_deviation() const
{ return moments().second; }


template <typename T>
Real RangeVariable<T>::variance() const
{ Real std_dev = moments().second; return std_dev * std_dev; }


//template <typename T>
//Real RangeVariable<T>::coefficient_of_variation() const
//{ RealRealPair mom = moments(); return mom.second / mom.first; }


template <typename T>
Real RangeVariable<T>::mode() const
{ return moments().first; } // not well-defined: any value in [L,U] is valid


/*
/// for T-valued histogram, return a real-valued mean and std dev
template <typename T>
void RangeVariable<T>::
moments_from_params(const std::set<T>& vals, Real& mean, Real& std_dev)
{
  mean = 0.;
  Real val, raw2 = 0.;  size_t num_vals = vals.size();
  typename std::set<T>::const_iterator cit;
  for (cit = vals.begin(); cit != vals.end(); ++cit) {
    val   = (Real)(*cit);
    mean += val;
    raw2 += val * val;
  }
  mean /= num_vals;  raw2 /= num_vals;
  std_dev = std::sqrt(raw2 - mean * mean);
}
*/


template <typename T>
RealRealPair RangeVariable<T>::distribution_bounds() const
{ return RealRealPair((Real)lowerBnd, (Real)upperBnd); }


template <typename T>
Real RangeVariable<T>::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{ no_template_specialization("correlation_warping_factor"); return 0.; }


template <typename T>
Real RangeVariable<T>::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{ no_template_specialization("dx_ds"); return 0.; }


template <typename T>
Real RangeVariable<T>::dz_ds_factor(short u_type, Real x, Real z) const
{ no_template_specialization("dz_ds_factor"); return 0.; }


//// SPECIALIZATIONS ////


template <>
inline Real RangeVariable<Real>::pdf(Real x) const
{ return UniformRandomVariable::pdf(x, lowerBnd, upperBnd); }


template <>
inline Real RangeVariable<Real>::cdf(Real x) const
{ return UniformRandomVariable::cdf(x, lowerBnd, upperBnd); }


template <>
inline Real RangeVariable<Real>::ccdf(Real x) const
{ return UniformRandomVariable::ccdf(x, lowerBnd, upperBnd); }


template <>
inline Real RangeVariable<Real>::inverse_cdf(Real p_cdf) const
{ return UniformRandomVariable::inverse_cdf(p_cdf, lowerBnd, upperBnd); }


template <>
inline Real RangeVariable<Real>::inverse_ccdf(Real p_ccdf) const
{ return UniformRandomVariable::inverse_ccdf(p_ccdf, lowerBnd, upperBnd); }


template <>
inline RealRealPair RangeVariable<Real>::moments() const
{
  RealRealPair moms;
  UniformRandomVariable::
    moments_from_params(lowerBnd, upperBnd, moms.first, moms.second);
  return moms;
}


template <>
inline Real RangeVariable<Real>::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{ return UniformRandomVariable::corr_warp_fact(rv, corr); }


template <>
inline Real RangeVariable<Real>::
dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  return UniformRandomVariable::
    dx_ds_fact(dist_param, u_type, ranVarType, x, z);
}


template <>
inline Real RangeVariable<Real>::
dz_ds_factor(short u_type, Real x, Real z) const
{ return UniformRandomVariable::dz_ds_fact(u_type, upperBnd - lowerBnd, x, z); }

} // namespace Pecos

#endif
