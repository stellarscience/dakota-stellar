/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 GeometricRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef GEOMETRIC_RANDOM_VARIABLE_HPP
#define GEOMETRIC_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for geometric random variables.

/** Manages the probPerTrial parameter.  The geometric distribution is a
    special case of the negative binomial distribution for numTrials = 1. */

class GeometricRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  GeometricRandomVariable();
  /// alternate constructor
  GeometricRandomVariable(Real prob_per_trial);
  /// destructor
  ~GeometricRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;

  Real parameter(short dist_param) const;
  void parameter(short dist_param, Real val);

  Real mean() const;
  Real median() const;
  Real mode() const;
  Real standard_deviation() const;
  Real variance() const;
  RealRealPair bounds() const;

  //
  //- Heading: Member functions
  //

  void update(Real prob_per_trial);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real prob_per_trial);
  static Real cdf(Real x, Real prob_per_trial);

  static void moments_from_params(Real prob_per_trial,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Member functions
  //

  /// create a new geometricDist instance
  void update_boost();

  //
  //- Heading: Data
  //

  /// p parameter of geometric random variable
  Real probPerTrial;

  /// pointer to the Boost geometric_distribution instance
  geometric_dist* geometricDist;
};


inline GeometricRandomVariable::GeometricRandomVariable():
  RandomVariable(BaseConstructor()), probPerTrial(1.), geometricDist(NULL)
{ ranVarType = GEOMETRIC; }


inline GeometricRandomVariable::
GeometricRandomVariable(Real prob_per_trial):
  RandomVariable(BaseConstructor()), probPerTrial(prob_per_trial),
  geometricDist(new geometric_dist(prob_per_trial))
{ ranVarType = GEOMETRIC; }


inline GeometricRandomVariable::~GeometricRandomVariable()
{ if (geometricDist) delete geometricDist; }


inline Real GeometricRandomVariable::cdf(Real x) const
{ return bmth::cdf(*geometricDist, x); }


inline Real GeometricRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*geometricDist, x)); }


inline Real GeometricRandomVariable::inverse_cdf(Real p_cdf) const
{ return bmth::quantile(*geometricDist, p_cdf); }


inline Real GeometricRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*geometricDist, p_ccdf)); }


inline Real GeometricRandomVariable::pdf(Real x) const
{ return bmth::pdf(*geometricDist, x); }


inline Real GeometricRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case GE_P_PER_TRIAL: return probPerTrial; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in GeometricRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void GeometricRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case GE_P_PER_TRIAL: probPerTrial = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in GeometricRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
  update_boost(); // create a new geometricDist instance
}


inline Real GeometricRandomVariable::mean() const
{ return bmth::mean(*geometricDist); }


inline Real GeometricRandomVariable::median() const
{ return bmth::median(*geometricDist); }


inline Real GeometricRandomVariable::mode() const
{ return bmth::mode(*geometricDist); }


inline Real GeometricRandomVariable::standard_deviation() const
{ return bmth::standard_deviation(*geometricDist); }


inline Real GeometricRandomVariable::variance() const
{ return bmth::variance(*geometricDist); }


inline RealRealPair GeometricRandomVariable::bounds() const
{ return RealRealPair(0., std::numeric_limits<Real>::infinity()); }


inline void GeometricRandomVariable::update_boost()
{
  if (geometricDist) delete geometricDist;
  geometricDist = new geometric_dist(probPerTrial);
}


inline void GeometricRandomVariable::update(Real prob_per_trial)
{
  if (!geometricDist || probPerTrial != prob_per_trial)
    { probPerTrial = prob_per_trial; update_boost(); }
}

// Static member functions:

inline Real GeometricRandomVariable::pdf(Real x, Real prob_per_trial)
{
  geometric_dist geometric1(prob_per_trial);
  return bmth::pdf(geometric1, x);
}


inline Real GeometricRandomVariable::cdf(Real x, Real prob_per_trial)
{
  geometric_dist geometric1(prob_per_trial);
  return bmth::cdf(geometric1, x);
}


inline void GeometricRandomVariable::
moments_from_params(Real prob_per_trial, Real& mean, Real& std_dev)
{
  Real comp_p = 1. - prob_per_trial;
  mean = comp_p / prob_per_trial; std_dev = std::sqrt(comp_p) / prob_per_trial;
}

} // namespace Pecos

#endif
