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

  /// n parameter of binomial random variable
  unsigned int numTrials;
  /// p parameter of binomial random variable
  Real probPerTrial;

  /// pointer to the Boost binomial_distribution instance
  binomial_dist* binomialDist;
};


inline BinomialRandomVariable::BinomialRandomVariable():
  RandomVariable(BaseConstructor()), numTrials(1), probPerTrial(1.),
  binomialDist(NULL)
{ ranVarType = BINOMIAL; }


inline BinomialRandomVariable::
BinomialRandomVariable(unsigned int num_trials, Real prob_per_trial):
  RandomVariable(BaseConstructor()),
  numTrials(num_trials), probPerTrial(prob_per_trial),
  binomialDist(new binomial_dist((Real)num_trials, prob_per_trial))
{ ranVarType = BINOMIAL; }


inline BinomialRandomVariable::~BinomialRandomVariable()
{ if (binomialDist) delete binomialDist; }


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


inline Real BinomialRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case BI_TRIALS:      return (Real)numTrials; break;
  case BI_P_PER_TRIAL: return probPerTrial;    break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in BinomialRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void BinomialRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case BI_TRIALS:      numTrials = (unsigned int)std::floor(val+.5); break;
  case BI_P_PER_TRIAL: probPerTrial = val;                           break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in BinomialRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
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


inline RealRealPair BinomialRandomVariable::bounds() const
{ return RealRealPair(0., std::numeric_limits<Real>::infinity()); }


inline void BinomialRandomVariable::update_boost()
{
  if (binomialDist) delete binomialDist;
  binomialDist = new binomial_dist((Real)numTrials, probPerTrial);
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
