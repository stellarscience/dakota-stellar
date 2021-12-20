/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 IntervalRandomVariable
//- Description: A random variable described by discrete values + probabilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef INTERVAL_RANDOM_VARIABLE_HPP
#define INTERVAL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"
#include "HistogramBinRandomVariable.hpp"
#include "DiscreteSetRandomVariable.hpp"
#include "pecos_stat_util.hpp"

namespace Pecos {


/// Derived random variable class for interval random variables

/** Manages basic probability assignments (BPAs: pairings of
    intervals with probabilities) for types int and real. */

template <typename T>
class IntervalRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  IntervalRandomVariable();
  /// alternate constructor
  IntervalRandomVariable(const std::map<std::pair<T, T>, Real>& bpa);
  /// destructor
  ~IntervalRandomVariable();

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

  Real mean() const;
  //Real median() const;
  Real mode() const;
  Real standard_deviation() const;
  Real variance() const;
  
  RealRealPair moments() const;
  RealRealPair distribution_bounds() const;

  //Real coefficient_of_variation() const;

  void pull_parameter(short dist_param,
		      std::map<std::pair<T, T>, Real>& bpa) const;
  void push_parameter(short dist_param,
		      const std::map<std::pair<T, T>, Real>& bpa);

  void copy_parameters(const RandomVariable& rv);

  //
  //- Heading: Member functions
  //

  /// update intervalBPA
  void update(const std::map<std::pair<T, T>, Real>& bpa);

  /// verify that valueProbPairs is defined
  void check_vpp() const;
  /// activate valueProbPairs tracking by copying intervalBPA
  void activate_vpp();
  /// update valueProbPairs if currently in use
  void update_vpp_if_active();

  //
  //- Heading: Static member functions (global utilities)
  //

  static void moments_from_params(const std::map<std::pair<T, T>, Real>& bpa,
				  Real& mean, Real& std_dev);

  //static Real pdf(Real x, const RealVector& dist_params);

protected:

  //
  //- Heading: Data
  //

  /// interval-probability pairs for basic probability assignment (BPA)
  std::map<std::pair<T, T>, Real> intervalBPA;

  /// collapse (overlapping, disjoint) intervals into a histogram-like format
  /// (defined if needed for use in moments/PDF/CDF)
  std::map<T, Real> valueProbPairs;
};


//// GENERIC ////


template <typename T>
IntervalRandomVariable<T>::IntervalRandomVariable():
  RandomVariable(BaseConstructor())
{ }


template <typename T>
IntervalRandomVariable<T>::
IntervalRandomVariable(const std::map<std::pair<T, T>, Real>& bpa):
  RandomVariable(BaseConstructor()), intervalBPA(bpa)
{ }


template <typename T>
IntervalRandomVariable<T>::~IntervalRandomVariable()
{ }


template <typename T>
void IntervalRandomVariable<T>::check_vpp() const
{
  if (valueProbPairs.empty()) {
    PCerr << "Error: valueProbPairs not activated in IntervalRandomVariable<T>."
	  << std::endl;
    abort_handler(-1);
  }
}


template <typename T>
void IntervalRandomVariable<T>::activate_vpp()
{
  if (valueProbPairs.empty()) // activate once to avoid copying repeatedly
    intervals_to_xy_pdf(intervalBPA, valueProbPairs);
}


template <typename T>
void IntervalRandomVariable<T>::update_vpp_if_active()
{
  if (!valueProbPairs.empty()) // update valueProbPairs if already activated
    intervals_to_xy_pdf(intervalBPA, valueProbPairs);
}


template <typename T>
void IntervalRandomVariable<T>::
update(const std::map<std::pair<T, T>, Real>& bpa)
{
  intervalBPA = bpa;
  update_vpp_if_active(); // update valueProbPairs if in use, else leave empty
}
// specializations could be used for assigning ranVarType, or could employ
// std::is_same for type identification.  Simplest: ranVarType assigned at
// bottom of RandomVariable::get_random_variable().


template <typename T>
void IntervalRandomVariable<T>::
pull_parameter(short dist_param, std::map<std::pair<T, T>, Real>& bpa) const
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case CIU_BPA: case DIU_BPA:  bpa = intervalBPA;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in IntervalRandomVariable::pull_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
void IntervalRandomVariable<T>::
push_parameter(short dist_param, const std::map<std::pair<T, T>, Real>& bpa)
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case CIU_BPA: case DIU_BPA:
    intervalBPA = bpa;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in IntervalRandomVariable::push_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
  update_vpp_if_active(); // update valueProbPairs if in use, else leave empty
}


template <typename T>
void IntervalRandomVariable<T>::copy_parameters(const RandomVariable& rv)
{
  switch (ranVarType) {
  case CONTINUOUS_INTERVAL_UNCERTAIN:
    rv.pull_parameter(CIU_BPA, intervalBPA); break;
  case DISCRETE_INTERVAL_UNCERTAIN:
    rv.pull_parameter(DIU_BPA, intervalBPA); break;
  default:
    PCerr << "Error: update failure for RandomVariable type " << rv.type()
	  << " in IntervalRandomVariable::copy_parameters(T)." << std::endl;
    abort_handler(-1); break;
  }
  update_vpp_if_active(); // update valueProbPairs if in use, else leave empty
}


template <typename T>
Real IntervalRandomVariable<T>::cdf(Real x) const
{
  //activate_vpp(); // can't do this since non-const

  if (valueProbPairs.empty()) {
    std::map<T, Real> xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    return DiscreteSetRandomVariable<T>::cdf(x, xy_map);
  }
  else
    return DiscreteSetRandomVariable<T>::cdf(x, valueProbPairs);
}


template <typename T>
Real IntervalRandomVariable<T>::ccdf(Real x) const
{
  //activate_vpp(); // can't do this since non-const

  if (valueProbPairs.empty()) {
    std::map<T, Real> xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    return DiscreteSetRandomVariable<T>::ccdf(x, xy_map);
  }
  else
    return DiscreteSetRandomVariable<T>::ccdf(x, valueProbPairs);
}


template <typename T>
Real IntervalRandomVariable<T>::inverse_cdf(Real p_cdf) const
{
  //activate_vpp(); // can't do this since non-const

  if (valueProbPairs.empty()) {
    std::map<T, Real> xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    return DiscreteSetRandomVariable<T>::inverse_cdf(p_cdf, xy_map);
  }
  else
    return DiscreteSetRandomVariable<T>::inverse_cdf(p_cdf, valueProbPairs);
}


template <typename T>
Real IntervalRandomVariable<T>::inverse_ccdf(Real p_ccdf) const
{
  //activate_vpp(); // can't do this since non-const

  if (valueProbPairs.empty()) {
    std::map<T, Real> xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    return DiscreteSetRandomVariable<T>::inverse_ccdf(p_ccdf, xy_map);
  }
  else
    return DiscreteSetRandomVariable<T>::inverse_ccdf(p_ccdf, valueProbPairs);
}


template <typename T>
Real IntervalRandomVariable<T>::pdf(Real x) const
{
  //activate_vpp(); // can't do this since non-const

  if (valueProbPairs.empty()) {
    std::map<T, Real> xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    return DiscreteSetRandomVariable<T>::pdf(x, xy_map);
  }
  else
    return DiscreteSetRandomVariable<T>::pdf(x, valueProbPairs);
}


template <typename T>
Real IntervalRandomVariable<T>::pdf_gradient(Real x) const
{ return 0.; }


template <typename T>
Real IntervalRandomVariable<T>::pdf_hessian(Real x) const
{ return 0.; }


template <typename T>
RealRealPair IntervalRandomVariable<T>::moments() const
{
  //activate_vpp(); // can't do this since non-const

  RealRealPair moms;
  if (valueProbPairs.empty()) {
    std::map<T, Real> xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    DiscreteSetRandomVariable<T>::
      moments_from_params(xy_map, moms.first, moms.second);
  }
  else
    DiscreteSetRandomVariable<T>::
      moments_from_params(valueProbPairs, moms.first, moms.second);
  return moms;
}


template <typename T>
Real IntervalRandomVariable<T>::mean() const
{ return moments().first; }


template <typename T>
Real IntervalRandomVariable<T>::standard_deviation() const
{ return moments().second; }


template <typename T>
Real IntervalRandomVariable<T>::variance() const
{ Real std_dev = moments().second; return std_dev * std_dev; }


template <typename T>
Real IntervalRandomVariable<T>::mode() const
{
  //activate_vpp(); // can't do this since non-const

  if (valueProbPairs.empty()) {
    std::map<T, Real> xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    return DiscreteSetRandomVariable<T>::mode(xy_map);
  }
  else
    return DiscreteSetRandomVariable<T>::mode(valueProbPairs);
}


template <typename T>
RealRealPair IntervalRandomVariable<T>::distribution_bounds() const
{
  T l_bnd, u_bnd;
  if (valueProbPairs.empty()) {
    typename std::map<std::pair<T, T>, Real>::const_iterator
      cit = intervalBPA.begin();
    const std::pair<T, T>& bnds0 = cit->first;
    l_bnd = bnds0.first;  u_bnd = bnds0.second;  ++cit;
    for (; cit!=intervalBPA.end(); ++cit) {
      const std::pair<T, T>& bnds = cit->first;
      if (bnds.first  < l_bnd) l_bnd = bnds.first;
      if (bnds.second > u_bnd) u_bnd = bnds.second;
    }
  }
  else {
    l_bnd =   valueProbPairs.begin()->first;
    u_bnd = (--valueProbPairs.end())->first;
  }

  return RealRealPair((Real)l_bnd, (Real)u_bnd);
}


/// for T-valued histogram, return a real-valued mean and std dev
template <typename T>
void IntervalRandomVariable<T>::
moments_from_params(const std::map<std::pair<T, T>, Real>& bpa,
		    Real& mean, Real& std_dev)
{
  typename std::map<T, Real> xy_map;
  intervals_to_xy_pdf(bpa, xy_map);
  DiscreteSetRandomVariable<T>::moments_from_params(xy_map, mean, std_dev);
}


//// SPECIALIZATIONS ////


template <>
inline Real IntervalRandomVariable<Real>::cdf(Real x) const
{
  //activate_vpp(); // can't do this since non-const

  if (valueProbPairs.empty()) {
    RealRealMap xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    return HistogramBinRandomVariable::cdf(x, xy_map);
  }
  else
    return HistogramBinRandomVariable::cdf(x, valueProbPairs);
}


template <>
inline Real IntervalRandomVariable<Real>::ccdf(Real x) const
{
  //activate_vpp(); // can't do this since non-const

  if (valueProbPairs.empty()) {
    RealRealMap xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    return HistogramBinRandomVariable::ccdf(x, xy_map);
  }
  else
    return HistogramBinRandomVariable::ccdf(x, valueProbPairs);
}


template <>
inline Real IntervalRandomVariable<Real>::inverse_cdf(Real x) const
{
  //activate_vpp(); // can't do this since non-const

  if (valueProbPairs.empty()) {
    RealRealMap xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    return HistogramBinRandomVariable::inverse_cdf(x, xy_map);
  }
  else
    return HistogramBinRandomVariable::inverse_cdf(x, valueProbPairs);
}


template <>
inline Real IntervalRandomVariable<Real>::inverse_ccdf(Real x) const
{
  //activate_vpp(); // can't do this since non-const

  if (valueProbPairs.empty()) {
    RealRealMap xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    return HistogramBinRandomVariable::inverse_ccdf(x, xy_map);
  }
  else
    return HistogramBinRandomVariable::inverse_ccdf(x, valueProbPairs);
}


template <>
inline Real IntervalRandomVariable<Real>::pdf(Real x) const
{
  //activate_vpp(); // can't do this since non-const

  if (valueProbPairs.empty()) {
    RealRealMap xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    return HistogramBinRandomVariable::pdf(x, xy_map);
  }
  else
    return HistogramBinRandomVariable::pdf(x, valueProbPairs);
}


template <>
inline RealRealPair IntervalRandomVariable<Real>::moments() const
{
  //activate_vpp(); // can't do this since non-const

  RealRealPair moms;
  if (valueProbPairs.empty()) {
    RealRealMap xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    HistogramBinRandomVariable::
      moments_from_params(xy_map, moms.first, moms.second);
  }
  else
    HistogramBinRandomVariable::
      moments_from_params(valueProbPairs, moms.first, moms.second);
  return moms;
}


template <>
inline Real IntervalRandomVariable<Real>::mode() const
{
  //activate_vpp(); // can't do this since non-const

  if (valueProbPairs.empty()) {
    RealRealMap xy_map;  intervals_to_xy_pdf(intervalBPA, xy_map);
    return HistogramBinRandomVariable::mode(xy_map);
  }
  else
    return HistogramBinRandomVariable::mode(valueProbPairs);
}


template <>
inline void IntervalRandomVariable<Real>::
moments_from_params(const RealRealPairRealMap& bpa, Real& mean, Real& std_dev)
{
  RealRealMap xy_map;  intervals_to_xy_pdf(bpa, xy_map);
  HistogramBinRandomVariable::moments_from_params(xy_map, mean, std_dev);
}

} // namespace Pecos

#endif
