/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 BinomialRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef BINOMIAL_RANDOM_VARIABLE_HPP
#define BINOMIAL_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for binomial random variables.

/** Manages numTrials and probPerTrial parameters. */

class BinomialRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  BinomialRandomVariable();
  /// alternate constructor
  BinomialRandomVariable(unsigned int num_trials, Real prob_per_trial);
  /// destructor
  ~BinomialRandomVariable();

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

  /// create a new binomialDist instance
  void update_boost();

  //
  //- Heading: Data
  //

  /// p parameter of binomial random variable
  Real probPerTrial;
  /// n parameter of binomial random variable
  unsigned int numTrials;

  /// pointer to the Boost binomial_distribution instance
  std::unique_ptr<binomial_dist> binomialDist;
};


inline BinomialRandomVariable::BinomialRandomVariable():
  RandomVariable(BaseConstructor()), probPerTrial(1.), numTrials(1), 
  binomialDist(new binomial_dist((Real)numTrials, probPerTrial))
{ ranVarType = BINOMIAL; }


inline BinomialRandomVariable::
BinomialRandomVariable(unsigned int num_trials, Real prob_per_trial):
  RandomVariable(BaseConstructor()),
  probPerTrial(prob_per_trial), numTrials(num_trials), 
  binomialDist(new binomial_dist((Real)num_trials, prob_per_trial))
{ ranVarType = BINOMIAL; }


inline BinomialRandomVariable::~BinomialRandomVariable()
{ }


inline Real BinomialRandomVariable::cdf(Real x) const
{ return bmth::cdf(*binomialDist, x); }


inline Real BinomialRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*binomialDist, x)); }


inline Real BinomialRandomVariable::inverse_cdf(Real p_cdf) const
{ return bmth::quantile(*binomialDist, p_cdf); }


inline Real BinomialRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*binomialDist, p_ccdf)); }


inline Real BinomialRandomVariable::pdf(Real x) const
{ return bmth::pdf(*binomialDist, x); }


inline void BinomialRandomVariable::
pull_parameter(short dist_param, Real& val) const
{
  switch (dist_param) {
  case BI_P_PER_TRIAL:  val = probPerTrial;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in BinomialRandomVariable::pull_parameter(Real)." << std::endl;
    abort_handler(-1); break;
  }
}


inline void BinomialRandomVariable::push_parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case BI_P_PER_TRIAL:  probPerTrial = val;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in BinomialRandomVariable::push_parameter(Real)." << std::endl;
    abort_handler(-1); break;
  }
  update_boost(); // create a new binomialDist instance
}


inline void BinomialRandomVariable::
pull_parameter(short dist_param, unsigned int& val) const
{
  switch (dist_param) {
  case BI_TRIALS:  val = numTrials;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in BinomialRandomVariable::pull_parameter(unsigned int)."
	  << std::endl;
    abort_handler(-1); break;
  }
}


inline void BinomialRandomVariable::
push_parameter(short dist_param, unsigned int val)
{
  switch (dist_param) {
  case BI_TRIALS:  numTrials = val;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in BinomialRandomVariable::push_parameter(unsigned int)."
	  << std::endl;
    abort_handler(-1); break;
  }
  update_boost(); // create a new binomialDist instance
}


inline void BinomialRandomVariable::copy_parameters(const RandomVariable& rv)
{
  rv.pull_parameter(BI_P_PER_TRIAL, probPerTrial);
  rv.pull_parameter(BI_TRIALS,      numTrials);
  update_boost(); // create a new binomialDist instance
}


inline Real BinomialRandomVariable::mean() const
{ return bmth::mean(*binomialDist); }


inline Real BinomialRandomVariable::median() const
{ return bmth::median(*binomialDist); }


inline Real BinomialRandomVariable::mode() const
{ return bmth::mode(*binomialDist); }


inline Real BinomialRandomVariable::standard_deviation() const
{ return bmth::standard_deviation(*binomialDist); }


inline Real BinomialRandomVariable::variance() const
{ return bmth::variance(*binomialDist); }


inline RealRealPair BinomialRandomVariable::distribution_bounds() const
{ return RealRealPair(0., std::numeric_limits<Real>::infinity()); }


inline void BinomialRandomVariable::update_boost()
{
  binomialDist.reset(new binomial_dist((Real)numTrials, probPerTrial));
}


inline void BinomialRandomVariable::
update(unsigned int num_trials, Real prob_per_trial)
{
  if (!binomialDist ||
      numTrials != num_trials || probPerTrial != prob_per_trial)
    { numTrials = num_trials; probPerTrial = prob_per_trial; update_boost(); }
}

// Static member functions:

inline Real BinomialRandomVariable::
pdf(Real x, unsigned int num_trials, Real prob_per_trial)
{
  binomial_dist binomial1((Real)num_trials, prob_per_trial);
  return bmth::pdf(binomial1, x);
}


inline Real BinomialRandomVariable::
cdf(Real x, unsigned int num_trials, Real prob_per_trial)
{
  binomial_dist binomial1((Real)num_trials, prob_per_trial);
  return bmth::cdf(binomial1, x);
}


inline void BinomialRandomVariable::
moments_from_params(unsigned int num_trials, Real prob_per_trial,
		    Real& mean, Real& std_dev)
{
  mean    = (Real)num_trials * prob_per_trial;
  std_dev = std::sqrt(mean * (1.-prob_per_trial));
}

} // namespace Pecos

#endif
