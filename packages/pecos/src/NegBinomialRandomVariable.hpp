/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 NegBinomialRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef NEG_BINOMIAL_RANDOM_VARIABLE_HPP
#define NEG_BINOMIAL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for negative binomial random variables.

/** Manages numTrials and probPerTrial parameters.  Note that the geometric
    distribution is a special case of the negative binomial distribution
    for numTrials = 1; however, there is currently little benefit to
    deriving NegBinomialRandomVariable from GeometricRandomVariable (e.g.,
    the Boost distribution pointers are distinct). */

class NegBinomialRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  NegBinomialRandomVariable();
  /// alternate constructor
  NegBinomialRandomVariable(unsigned int num_trials, Real prob_per_trial);
  /// destructor
  ~NegBinomialRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;

  void pull_parameter(short dist_param, Real& val) const;
  void push_parameter(short dist_param, Real  val);
  void pull_parameter(short dist_param, unsigned int& val) const;
  void push_parameter(short dist_param, unsigned int  val);

  void copy_parameters(const RandomVariable& rv);

  Real mean() const;
  Real median() const;
  Real mode() const;
  Real standard_deviation() const;
  Real variance() const;
  RealRealPair distribution_bounds() const;

  //
  //- Heading: Member functions
  //

  void update(unsigned int num_trials, Real prob_per_trial);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, unsigned int num_trials, Real prob_per_trial);
  static Real cdf(Real x, unsigned int num_trials, Real prob_per_trial);

  static void moments_from_params(unsigned int num_trials, Real prob_per_trial,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Member functions
  //

  /// create a new negBinomialDist instance
  void update_boost();

  //
  //- Heading: Data
  //

  /// p parameter of negative binomial random variable
  Real probPerTrial;
  /// r parameter of negative binomial random variable
  unsigned int numTrials;

  /// pointer to the Boost negative_binomial_distribution instance
  std::unique_ptr<negative_binomial_dist> negBinomialDist;
};


inline NegBinomialRandomVariable::NegBinomialRandomVariable():
  RandomVariable(BaseConstructor()), probPerTrial(1.), numTrials(1),
  negBinomialDist(new negative_binomial_dist((Real)numTrials, probPerTrial))
{ ranVarType = NEGATIVE_BINOMIAL; }


inline NegBinomialRandomVariable::
NegBinomialRandomVariable(unsigned int num_trials, Real prob_per_trial):
  RandomVariable(BaseConstructor()),
  probPerTrial(prob_per_trial), numTrials(num_trials), 
  negBinomialDist(new negative_binomial_dist((Real)num_trials, prob_per_trial))
{ ranVarType = NEGATIVE_BINOMIAL; }


inline NegBinomialRandomVariable::~NegBinomialRandomVariable()
{ }


inline Real NegBinomialRandomVariable::cdf(Real x) const
{ return bmth::cdf(*negBinomialDist, x); }


inline Real NegBinomialRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*negBinomialDist, x)); }


inline Real NegBinomialRandomVariable::inverse_cdf(Real p_cdf) const
{ return bmth::quantile(*negBinomialDist, p_cdf); }


inline Real NegBinomialRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*negBinomialDist, p_ccdf)); }


inline Real NegBinomialRandomVariable::pdf(Real x) const
{ return bmth::pdf(*negBinomialDist, x); }


inline void NegBinomialRandomVariable::
pull_parameter(short dist_param, Real& val) const
{
  switch (dist_param) {
  case NBI_P_PER_TRIAL:  val = probPerTrial;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in NegBinomialRandomVariable::pull_parameter(Real)."<< std::endl;
    abort_handler(-1); break;
  }
}


inline void NegBinomialRandomVariable::
push_parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case NBI_P_PER_TRIAL:  probPerTrial = val;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in NegBinomialRandomVariable::push_parameter(Real)."<< std::endl;
    abort_handler(-1); break;
  }
  update_boost(); // create a new binomialDist instance
}


inline void NegBinomialRandomVariable::
pull_parameter(short dist_param, unsigned int& val) const
{
  switch (dist_param) {
  case NBI_TRIALS:  val = numTrials;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in NegBinomialRandomVariable::pull_parameter(unsigned int)."
	  << std::endl;
    abort_handler(-1); break;
  }
}


inline void NegBinomialRandomVariable::
push_parameter(short dist_param, unsigned int val)
{
  switch (dist_param) {
  case NBI_TRIALS:  numTrials = val;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in NegBinomialRandomVariable::push_parameter(unsigned int)."
	  << std::endl;
    abort_handler(-1); break;
  }
  update_boost(); // create a new negBinomialDist instance
}


inline void NegBinomialRandomVariable::copy_parameters(const RandomVariable& rv)
{
  rv.pull_parameter(NBI_P_PER_TRIAL, probPerTrial);
  rv.pull_parameter(NBI_TRIALS,      numTrials);
  update_boost(); // create a new negBinomialDist instance
}


inline Real NegBinomialRandomVariable::mean() const
{ return bmth::mean(*negBinomialDist); }


inline Real NegBinomialRandomVariable::median() const
{ return bmth::median(*negBinomialDist); }


inline Real NegBinomialRandomVariable::mode() const
{ return bmth::mode(*negBinomialDist); }


inline Real NegBinomialRandomVariable::standard_deviation() const
{ return bmth::standard_deviation(*negBinomialDist); }


inline Real NegBinomialRandomVariable::variance() const
{ return bmth::variance(*negBinomialDist); }


inline RealRealPair NegBinomialRandomVariable::distribution_bounds() const
{ return RealRealPair(0., std::numeric_limits<Real>::infinity()); }


inline void NegBinomialRandomVariable::update_boost()
{
  negBinomialDist.reset
    (new negative_binomial_dist((Real)numTrials, probPerTrial))
;
}


inline void NegBinomialRandomVariable::
update(unsigned int num_trials, Real prob_per_trial)
{
  if (!negBinomialDist ||
      numTrials != num_trials || probPerTrial != prob_per_trial)
    { numTrials = num_trials; probPerTrial = prob_per_trial; update_boost(); }
}

// Static member functions:

inline Real NegBinomialRandomVariable::
pdf(Real x, unsigned int num_trials, Real prob_per_trial)
{
  negative_binomial_dist neg_binomial1((Real)num_trials, prob_per_trial);
  return bmth::pdf(neg_binomial1, x);
}


inline Real NegBinomialRandomVariable::
cdf(Real x, unsigned int num_trials, Real prob_per_trial)
{
  negative_binomial_dist neg_binomial1((Real)num_trials, prob_per_trial);
  return bmth::cdf(neg_binomial1, x);
}


inline void NegBinomialRandomVariable::
moments_from_params(unsigned int num_trials, Real prob_per_trial,
		    Real& mean, Real& std_dev)
{
  Real n1mp = (Real)num_trials * (1. - prob_per_trial);
  mean = n1mp / prob_per_trial; std_dev = std::sqrt(n1mp) / prob_per_trial;
}

} // namespace Pecos

#endif
