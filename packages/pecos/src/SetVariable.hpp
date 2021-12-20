/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 SetVariable
//- Description: A random variable described by discrete values without
//-              associated probabilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef SET_VARIABLE_HPP
#define SET_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for discrete set random variables.

/** Manages a set of discrete values without associated probability pairings
    (refer to DiscreteSetRandomVariable for pairings) for types int, string,
    and real.  String values are managed by index rather than value,
    requiring template specializations. */

template <typename T>
class SetVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  SetVariable();
  /// alternate constructor
  SetVariable(const std::set<T>& vals);
  /// destructor
  ~SetVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  /*
  Real mean() const;
  //Real median() const;
  Real mode() const;
  Real standard_deviation() const;
  Real variance() const;
  
  RealRealPair moments() const;
  Real coefficient_of_variation() const;
  */
  RealRealPair distribution_bounds() const;

  void pull_parameter(short dist_param, std::set<T>& vals) const;
  void push_parameter(short dist_param, const std::set<T>& vals);

  void copy_parameters(const RandomVariable& rv);

  //
  //- Heading: Member functions
  //

  void update(const std::set<T>& vals);

  //
  //- Heading: Static member functions (global utilities)
  //

  static void moments_from_params(const std::set<T>& vals,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// value-prob pairs for int values within a set
  std::set<T> setValues;
};


//// GENERIC ////


template <typename T>
SetVariable<T>::SetVariable():
  RandomVariable(BaseConstructor())
{ }


template <typename T>
SetVariable<T>::SetVariable(const std::set<T>& vals):
  RandomVariable(BaseConstructor()), setValues(vals)
{ }


template <typename T>
SetVariable<T>::~SetVariable()
{ }


template <typename T>
void SetVariable<T>::update(const std::set<T>& vals)
{ setValues = vals; }
// specializations could be used for assigning ranVarType, or could employ
// std::is_same for type identification.  Simplest: ranVarType assigned at
// bottom of RandomVariable::get_random_variable().


template <typename T>
void SetVariable<T>::pull_parameter(short dist_param, std::set<T>& vals) const
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case DSI_VALUES: case DSS_VALUES: case DSR_VALUES:
    vals = setValues; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in SetVariable::pull_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
void SetVariable<T>::push_parameter(short dist_param, const std::set<T>& vals)
{
  // could specialize template, but case aggregation seems adequate

  switch (dist_param) {
  case DSI_VALUES: case DSS_VALUES: case DSR_VALUES:
    setValues = vals; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in SetVariable::push_parameter(T)." << std::endl;
    abort_handler(-1); break;
  }
}


template <typename T>
void SetVariable<T>::copy_parameters(const RandomVariable& rv)
{
  switch (ranVarType) {
  case DISCRETE_SET_INT:     rv.pull_parameter(DSI_VALUES, setValues);  break;
  case DISCRETE_SET_STRING:  rv.pull_parameter(DSS_VALUES, setValues);  break;
  case DISCRETE_SET_REAL:    rv.pull_parameter(DSR_VALUES, setValues);  break;
  }
}


/*
template <typename T>
RealRealPair SetVariable<T>::moments() const
{
  RealRealPair moms;
  moments_from_params(setValues, moms.first, moms.second);
  return moms;
}


template <typename T>
Real SetVariable<T>::mean() const
{ return moments().first; }


//template <typename T>
//Real SetVariable<T>::median() const
//{ return inverse_cdf(.5); } // default


template <typename T>
Real SetVariable<T>::standard_deviation() const
{ return moments().second; }


template <typename T>
Real SetVariable<T>::variance() const
{ Real std_dev = moments().second; return std_dev * std_dev; }


template <typename T>
Real SetVariable<T>::coefficient_of_variation() const
{ RealRealPair mom = moments(); return mom.second / mom.first; }


template <typename T>
Real SetVariable<T>::mode() const
{
  Real mode, mode_prob;
  typename std::set<T>::const_iterator cit = setValues.begin();
  mode = (Real)cit->first;  mode_prob = cit->second;  ++cit;
  for (; cit != setValues.end(); ++cit)
    if (cit->second > mode_prob)
      { mode = (Real)cit->first;  mode_prob = cit->second; }
  return mode;
}
*/


template <typename T>
RealRealPair SetVariable<T>::distribution_bounds() const
{
  // set values are sorted
  T l_bnd = *setValues.begin(), u_bnd = *(--setValues.end());
  return RealRealPair((Real)l_bnd, (Real)u_bnd);
}


/*
/// for T-valued histogram, return a real-valued mean and std dev
template <typename T>
void SetVariable<T>::
moments_from_params(const std::set<T>& vals, Real& mean, Real& std_dev)
{
  // In point histogram case, (x,y) and (x,c) are equivalent since bins
  // have zero-width.  Assume normalization (prob values sum to 1.).
  mean = 0.;
  Real val, prod, raw2 = 0.;
  typename std::set<T>::const_iterator cit;
  for (cit = vals.begin(); cit != vals.end(); ++cit) {
    val   = (Real)cit->first;
    prod  = cit->second * val; // prob * val
    mean += prod;
    raw2 += prod * val;        // prob * val^2
  }
  std_dev = std::sqrt(raw2 - mean * mean);
}
*/


//// SPECIALIZATIONS ////
// for string vars, moments/bounds are based on weighting of set indices


/*
template <>
inline Real SetVariable<String>::mode() const
{
  Real mode, mode_prob;
  SRMCIter cit = setValues.begin();
  mode = 0.;  mode_prob = cit->second;  ++cit;
  for (size_t index=1; cit!=setValues.end(); ++cit, ++index)
    if (cit->second > mode_prob)
      { mode = (Real)index;  mode_prob = cit->second; }
  return mode;
}
*/


template <>
inline RealRealPair SetVariable<String>::distribution_bounds() const
{
  size_t last_index = setValues.size() - 1;
  return RealRealPair(0., (Real)last_index);
}


/*
template <>
void SetVariable<String>::
inline moments_from_params(const StringRealMap& s_prs,
                           Real& mean, Real& std_dev)
{
  // in point case, (x,y) and (x,c) are equivalent since bins have zero-width.
  // assume normalization (probs sum to 1.).
  mean = 0.;
  Real val, prod, raw2 = 0.;  size_t index = 0;  SRMCIter cit;
  for (cit = s_prs.begin(); cit != s_prs.end(); ++cit, ++index) {
    val   = (Real)index;
    prod  = cit->second * val; // normalized prob * val
    mean += prod;
    raw2 += prod * val;        // normalized prob * val^2
  }
  std_dev = std::sqrt(raw2 - mean * mean);
}
*/

} // namespace Pecos

#endif
