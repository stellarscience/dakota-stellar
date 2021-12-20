/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 DiscreteSetRandomVariable
//- Description: A random variable described by discrete values + probabilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef DISCRETE_SET_RANDOM_VARIABLE_HPP
#define DISCRETE_SET_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for discrete set random variables.

/** Manages value-probability pairings for types int, string, and real.
    String values are managed by index rather than value, requiring
    template specializations. */

template <typename T>
class DiscreteSetRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  DiscreteSetRandomVariable();
  /// alternate constructor
  DiscreteSetRandomVariable(const std::map<T, Real>& vals_probs);
  /// destructor
  ~DiscreteSetRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;

  Real mean() const;
  //Real median() const;
  Real mode() const;
  Real standard_deviation() const;
  Real variance() const;
  
  RealRealPair moments() const;
  RealRealPair distribution_bounds() const;

  Real coefficient_of_variation() const;

  void pull_parameter(short dist_param, std::map<T, Real>& val) const;
  void push_parameter(short dist_param, const std::map<T, Real>& val);

  void copy_parameters(const RandomVariable& rv);

  //
  //- Heading: Member functions
  //

  void update(const std::map<T, Real>& vals_probs);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real cdf(Real x, const std::map<T, Real>& vals_probs);
  static Real ccdf(Real x, const std::map<T, Real>& vals_probs);
  static Real inverse_cdf(Real p_cdf, const std::map<T, Real>& vals_probs);
  static Real inverse_ccdf(Real p_ccdf, const std::map<T, Real>& vals_probs);

  static Real pdf(Real x, const std::map<T, Real>& vals_probs);

  static Real mode(const std::map<T, Real>& vals_probs);

  static RealRealPair distribution_bounds(const std::map<T, Real>& vals_probs);

  static void moments_from_params(const std::map<T, Real>& vals_probs,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// value-prob pairs for int values within a set
  std::map<T, Real> valueProbPairs;
};


//// GENERIC ////


template <typename T>
DiscreteSetRandomVariable<T>::DiscreteSetRandomVariable():
  RandomVariable(BaseConstructor())
{ }


template <typename T>
DiscreteSetRandomVariable<T>::
DiscreteSetRandomVariable(const std::map<T, Real>& vals_probs):
  RandomVariable(BaseConstructor()), valueProbPairs(vals_probs)
{ }


template <typename T>
DiscreteSetRandomVariable<T>::~DiscreteSetRandomVariable()
{ }


template <typename T>
void DiscreteSetRandomVariable<T>::update(const std::map<T, Real>& vals_probs)
{ valueProbPairs = vals_probs; }
// specializations could be used for assigning ranVarType, or could employ
// std::is_same for type identification.  Simplest: ranVarType assigned at
// bottom of RandomVariable::get_random_variable().


template <typename T>
void DiscreteSetRandomVariable<T>::
pull_parameter(short dist_param, std::map<T, Real>& val) const
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case H_PT_INT_PAIRS:    case H_PT_STR_PAIRS:    case H_PT_REAL_PAIRS:
  case DUSI_VALUES_PROBS: case DUSS_VALUES_PROBS: case DUSR_VALUES_PROBS:
    val = valueProbPairs; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in DiscreteSetRandomVariable::pull_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
void DiscreteSetRandomVariable<T>::
push_parameter(short dist_param, const std::map<T, Real>& val)
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case H_PT_INT_PAIRS:    case H_PT_STR_PAIRS:    case H_PT_REAL_PAIRS:
  case DUSI_VALUES_PROBS: case DUSS_VALUES_PROBS: case DUSR_VALUES_PROBS:
    valueProbPairs = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in DiscreteSetRandomVariable::push_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
void DiscreteSetRandomVariable<T>::copy_parameters(const RandomVariable& rv)
{
  switch (ranVarType) {
  case HISTOGRAM_PT_INT:
    rv.pull_parameter(H_PT_INT_PAIRS,    valueProbPairs); break;
  case DISCRETE_UNCERTAIN_SET_INT:
    rv.pull_parameter(DUSI_VALUES_PROBS, valueProbPairs); break;
  case HISTOGRAM_PT_STRING:
    rv.pull_parameter(H_PT_STR_PAIRS,    valueProbPairs); break;
  case DISCRETE_UNCERTAIN_SET_STRING:
    rv.pull_parameter(DUSS_VALUES_PROBS, valueProbPairs); break;
  case HISTOGRAM_PT_REAL:
    rv.pull_parameter(H_PT_REAL_PAIRS,   valueProbPairs); break;
  case DISCRETE_UNCERTAIN_SET_REAL:
    rv.pull_parameter(DUSR_VALUES_PROBS, valueProbPairs); break;
  default:
    PCerr << "Error: update failure for RandomVariable type " << rv.type()
	  << " in DiscreteSetRandomVariable::copy_parameters(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
Real DiscreteSetRandomVariable<T>::
cdf(Real x, const std::map<T, Real>& vals_probs)
{
  // increment CDF probability until x is surpassed
  typename std::map<T, Real>::const_iterator cit;
  Real p_cdf = 0.;
  for (cit=vals_probs.begin(); cit!=vals_probs.end(); ++cit) {
    if (x <= (Real)(cit->first)) // has not reached start of this CDF increment
      return p_cdf;              // return previous value
    else
      p_cdf += cit->second;
  }
  return 1.;
}


template <typename T>
Real DiscreteSetRandomVariable<T>::
ccdf(Real x, const std::map<T, Real>& vals_probs)
{
  // decrement CCDF probability until x is surpassed
  typename std::map<T, Real>::const_iterator cit;
  Real p_ccdf = 1.;
  for (cit=vals_probs.begin(); cit!=vals_probs.end(); ++cit) {
    if (x < (Real)(cit->first))
      return p_ccdf;
    else
      p_ccdf -= cit->second;
  }
  return 0.;
}


template <typename T>
Real DiscreteSetRandomVariable<T>::
inverse_cdf(Real p_cdf, const std::map<T, Real>& vals_probs)
{
  // increment CDF probability until x is surpassed
  typename std::map<T, Real>::const_iterator cit;
  Real accumulated_p = 0., x = 0.; // return of this x value should not happen
  for (cit=vals_probs.begin(); cit!=vals_probs.end(); ++cit) {
    if (p_cdf <= accumulated_p) // CDF prob = P(z <= z-bar)
      return x;
    else {
      accumulated_p += cit->second;
      x = (Real)(cit->first);
    }
  }
  return (--vals_probs.end())->first;
}


template <typename T>
Real DiscreteSetRandomVariable<T>::
inverse_ccdf(Real p_ccdf, const std::map<T, Real>& vals_probs)
{
  // decrement CCDF probability until x is surpassed
  typename std::map<T, Real>::const_iterator cit;
  Real decremented_p = 1., x = 0.; // return of this x value should not happen
  for (cit=vals_probs.begin(); cit!=vals_probs.end(); ++cit) {
    if (p_ccdf > decremented_p) // CCDF prob = P(z > z-bar)
      return x;
    else {
      decremented_p -= cit->second;
      x = (Real)(cit->first);
    }
  }
  return (--vals_probs.end())->first;
}


template <typename T>
Real DiscreteSetRandomVariable<T>::
pdf(Real x, const std::map<T, Real>& vals_probs)
{
  // check to see if x can be converted to type T without loss of Real data
  T abscissa = (T)x; // cast Real to type T
  if (!real_compare(x, (Real)abscissa)) return 0.; // did casting change value?

  // x is a valid discrete value, now search for presence in vals_probs
  typename std::map<T, Real>::const_iterator cit = vals_probs.find(abscissa);
  return (cit == vals_probs.end()) ? 0. : cit->second;
}


template <typename T>
Real DiscreteSetRandomVariable<T>::mode(const std::map<T, Real>& vals_probs)
{
  Real mode, mode_prob;
  typename std::map<T, Real>::const_iterator cit = vals_probs.begin();
  mode = (Real)cit->first;  mode_prob = cit->second;  ++cit;
  for (; cit != vals_probs.end(); ++cit)
    if (cit->second > mode_prob)
      { mode = (Real)cit->first;  mode_prob = cit->second; }
  return mode;
}


template <typename T>
RealRealPair DiscreteSetRandomVariable<T>::
distribution_bounds(const std::map<T, Real>& vals_probs)
{
  RealRealPair bnds;
  bnds.first  = (Real)vals_probs.begin()->first;   // lower bound
  bnds.second = (Real)(--vals_probs.end())->first; // upper bound
  return bnds;
}


template <typename T>
Real DiscreteSetRandomVariable<T>::cdf(Real x) const
{ return cdf(x, valueProbPairs); }


template <typename T>
Real DiscreteSetRandomVariable<T>::ccdf(Real x) const
{ return ccdf(x, valueProbPairs); }


template <typename T>
Real DiscreteSetRandomVariable<T>::inverse_cdf(Real p_cdf) const
{ return inverse_cdf(p_cdf, valueProbPairs); }


template <typename T>
Real DiscreteSetRandomVariable<T>::inverse_ccdf(Real p_ccdf) const
{ return inverse_ccdf(p_ccdf, valueProbPairs); }


template <typename T>
Real DiscreteSetRandomVariable<T>::pdf(Real x) const
{ return pdf(x, valueProbPairs); }


template <typename T>
Real DiscreteSetRandomVariable<T>::mode() const
{ return mode(valueProbPairs); }


template <typename T>
RealRealPair DiscreteSetRandomVariable<T>::moments() const
{
  RealRealPair moms;
  moments_from_params(valueProbPairs, moms.first, moms.second);
  return moms;
}


template <typename T>
Real DiscreteSetRandomVariable<T>::mean() const
{ return moments().first; }


//template <typename T>
//Real DiscreteSetRandomVariable<T>::median() const
//{ return inverse_cdf(.5); } // default


template <typename T>
Real DiscreteSetRandomVariable<T>::standard_deviation() const
{ return moments().second; }


template <typename T>
Real DiscreteSetRandomVariable<T>::variance() const
{ Real std_dev = moments().second; return std_dev * std_dev; }


template <typename T>
Real DiscreteSetRandomVariable<T>::coefficient_of_variation() const
{ RealRealPair mom = moments(); return mom.second / mom.first; }


template <typename T>
RealRealPair DiscreteSetRandomVariable<T>::distribution_bounds() const
{ return distribution_bounds(valueProbPairs); }


/// for T-valued histogram, return a real-valued mean and std dev
template <typename T>
void DiscreteSetRandomVariable<T>::
moments_from_params(const std::map<T, Real>& vals_probs,
		    Real& mean, Real& std_dev)
{
  // In point histogram case, (x,y) and (x,c) are equivalent since bins
  // have zero-width.  Assume normalization (prob values sum to 1.).
  mean = 0.;
  Real val, prod, raw2 = 0.;
  typename std::map<T, Real>::const_iterator cit;
  for (cit = vals_probs.begin(); cit != vals_probs.end(); ++cit) {
    val   = (Real)cit->first;
    prod  = cit->second * val; // prob * val
    mean += prod;
    raw2 += prod * val;        // prob * val^2
  }
  std_dev = std::sqrt(raw2 - mean * mean);
}


//// SPECIALIZATIONS ////


// for string vars, moments/bounds are based on weighting of set indices


template <>
inline Real DiscreteSetRandomVariable<String>::
cdf(Real x, const StringRealMap& vals_probs)
{
  // increment CDF probability until x is surpassed
  SRMCIter cit;  size_t index = 0;  Real p_cdf = 0.;
  for (cit=vals_probs.begin(); cit!=vals_probs.end(); ++cit, ++index) {
    if (x <= (Real)index) // has not reached start of this index increment
      return p_cdf;       // return previous value
    else
      p_cdf += cit->second;
  }
  return 1.;
}


template <>
inline Real DiscreteSetRandomVariable<String>::
ccdf(Real x, const StringRealMap& vals_probs)
{
  // decrement CCDF probability until x is surpassed
  SRMCIter cit;  size_t index = 0;  Real p_ccdf = 1.;
  for (cit=vals_probs.begin(); cit!=vals_probs.end(); ++cit, ++index) {
    if (x < (Real)index)
      return p_ccdf;
    else
      p_ccdf -= cit->second;
  }
  return 0.;
}


template <>
inline Real DiscreteSetRandomVariable<String>::
inverse_cdf(Real p_cdf, const StringRealMap& vals_probs)
{
  // increment CDF probability until x is surpassed
  SRMCIter cit;  size_t index = 0;
  Real accumulated_p = 0., x = 0.; // return of this x value should not happen
  for (cit=vals_probs.begin(); cit!=vals_probs.end(); ++cit, ++index) {
    if (p_cdf <= accumulated_p) // CDF prob = P(z <= z-bar)
      return x;
    else {
      accumulated_p += cit->second;
      x = (Real)index;
    }
  }
  return x;
}


template <>
inline Real DiscreteSetRandomVariable<String>::
inverse_ccdf(Real p_ccdf, const StringRealMap& vals_probs)
{
  // decrement CCDF probability until x is surpassed
  SRMCIter cit;  size_t index = 0;
  Real decremented_p = 1., x = 0.; // return of this x value should not happen
  for (cit=vals_probs.begin(); cit!=vals_probs.end(); ++cit, ++index) {
    if (p_ccdf > decremented_p) // CCDF prob = P(z > z-bar)
      return x;
    else {
      decremented_p -= cit->second;
      x = (Real)index;
    }
  }
  return x;
}


template <>
inline Real DiscreteSetRandomVariable<String>::
pdf(Real x, const StringRealMap& vals_probs)
{
  // check to see if x can be converted to type T without loss of Real data
  size_t index = (size_t)x; // cast Real to size_t
  if ( !real_compare(x, (Real)index) || // did casting change value?
       index >= vals_probs.size() )      // index out of range
    return 0.;

  // x is a valid discrete value, now return ith probability from vals_probs
  SRMCIter cit = vals_probs.begin();  std::advance(cit, index);
  return cit->second;
}


template <>
inline Real DiscreteSetRandomVariable<String>::
mode(const StringRealMap& vals_probs)
{
  Real mode, mode_prob;
  SRMCIter cit = vals_probs.begin();
  mode = 0.;  mode_prob = cit->second;  ++cit;
  for (size_t index=1; cit!=vals_probs.end(); ++cit, ++index)
    if (cit->second > mode_prob)
      { mode = (Real)index;  mode_prob = cit->second; }
  return mode;
}


template <>
inline RealRealPair DiscreteSetRandomVariable<String>::
distribution_bounds(const StringRealMap& vals_probs)
{
  RealRealPair bnds;
  bnds.first  = 0.;                            // index lower bound
  bnds.second = (Real)(vals_probs.size() - 1); // index upper bound
  return bnds;
}


template <>
inline void DiscreteSetRandomVariable<String>::
moments_from_params(const StringRealMap& vals_probs, Real& mean, Real& std_dev)
{
  // in point case, (x,y) and (x,c) are equivalent since bins have zero-width.
  // assume normalization (probs sum to 1.).
  mean = 0.;
  Real val, prod, raw2 = 0.;  size_t index = 0;  SRMCIter cit;
  for (cit = vals_probs.begin(); cit != vals_probs.end(); ++cit, ++index) {
    val   = (Real)index;
    prod  = cit->second * val; // normalized prob * val
    mean += prod;
    raw2 += prod * val;        // normalized prob * val^2
  }
  std_dev = std::sqrt(raw2 - mean * mean);
}

} // namespace Pecos

#endif
