/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 HypergeometricRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef HYPERGEOMETRIC_RANDOM_VARIABLE_HPP
#define HYPERGEOMETRIC_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for hypergeometric random variables.

/** Manages numTotalPop, numSelectPop, and numDrawn parameters. */

class HypergeometricRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  HypergeometricRandomVariable();
  /// alternate constructor
  HypergeometricRandomVariable(unsigned int num_total_pop,
			       unsigned int num_sel_pop,
			       unsigned int num_drawn);
  /// destructor
  ~HypergeometricRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;

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

  void update(unsigned int num_total_pop, unsigned int num_sel_pop,
	      unsigned int num_drawn);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, unsigned int num_total_pop, unsigned int num_sel_pop,
		  unsigned int num_drawn);
  static Real cdf(Real x, unsigned int num_total_pop, unsigned int num_sel_pop,
		  unsigned int num_drawn);

  static void moments_from_params(unsigned int num_total_pop,
				  unsigned int num_sel_pop,
				  unsigned int num_drawn,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Member functions
  //

  /// create a new hypergeomDist instance
  void update_boost();
  /// create a new hypergeomDist instance if parameter set is valid
  void update_boost_conditionally();

  //
  //- Heading: Data
  //

  /// N parameter of hypergeometric random variable
  unsigned int numTotalPop;
  /// n parameter of hypergeometric random variable
  unsigned int numSelectPop;
  /// r parameter of hypergeometric random variable
  unsigned int numDrawn;

  /// pointer to the Boost hypergeometric_distribution instance
  std::unique_ptr<hypergeometric_dist> hypergeomDist;
};


inline HypergeometricRandomVariable::HypergeometricRandomVariable():
  RandomVariable(BaseConstructor()),
  numTotalPop(1), numSelectPop(1), numDrawn(1),
  hypergeomDist(new hypergeometric_dist(numDrawn, numSelectPop, numTotalPop))
{ ranVarType = HYPERGEOMETRIC; }


inline HypergeometricRandomVariable::
HypergeometricRandomVariable(unsigned int num_total_pop,
			     unsigned int num_sel_pop, unsigned int num_drawn):
  RandomVariable(BaseConstructor()), numTotalPop(num_total_pop),
  numSelectPop(num_sel_pop), numDrawn(num_drawn),
  hypergeomDist(new hypergeometric_dist(num_drawn, num_sel_pop, num_total_pop))
{ ranVarType = HYPERGEOMETRIC; }


inline HypergeometricRandomVariable::~HypergeometricRandomVariable()
{ }


inline Real HypergeometricRandomVariable::cdf(Real x) const
{ return bmth::cdf(*hypergeomDist, x); }


inline Real HypergeometricRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*hypergeomDist, x)); }


inline Real HypergeometricRandomVariable::inverse_cdf(Real p_cdf) const
{ return bmth::quantile(*hypergeomDist, p_cdf); }


inline Real HypergeometricRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*hypergeomDist, p_ccdf)); }


inline Real HypergeometricRandomVariable::pdf(Real x) const
{ return bmth::pdf(*hypergeomDist, x); }


inline void HypergeometricRandomVariable::
pull_parameter(short dist_param, unsigned int& val) const
{
  switch (dist_param) {
  case HGE_TOT_POP: val = numTotalPop;  break;
  case HGE_SEL_POP: val = numSelectPop; break;
  case HGE_DRAWN:   val = numDrawn;     break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in HypergeometricRandomVariable::pull_parameter(unsigned int)."
	  << std::endl;
    abort_handler(-1); break;
  }
}


inline void HypergeometricRandomVariable::
push_parameter(short dist_param, unsigned int val)
{
  switch (dist_param) {
  case HGE_TOT_POP: numTotalPop  = val; break;
  case HGE_SEL_POP: numSelectPop = val; break;
  case HGE_DRAWN:   numDrawn     = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in HypergeometricRandomVariable::push_parameter(unsigned int)."
	  << std::endl;
    abort_handler(-1); break;
  }
  update_boost_conditionally(); // create a new hypergeomDist instance
}


inline void HypergeometricRandomVariable::
copy_parameters(const RandomVariable& rv)
{
  rv.pull_parameter(HGE_TOT_POP, numTotalPop);
  rv.pull_parameter(HGE_SEL_POP, numSelectPop);
  rv.pull_parameter(HGE_DRAWN,   numDrawn);
  update_boost(); // create a new hypergeomDist instance
}


inline Real HypergeometricRandomVariable::mean() const
{ return bmth::mean(*hypergeomDist); }


inline Real HypergeometricRandomVariable::median() const
{ return bmth::median(*hypergeomDist); }


inline Real HypergeometricRandomVariable::mode() const
{ return bmth::mode(*hypergeomDist); }


inline Real HypergeometricRandomVariable::standard_deviation() const
{ return bmth::standard_deviation(*hypergeomDist); }


inline Real HypergeometricRandomVariable::variance() const
{ return bmth::variance(*hypergeomDist); }


inline RealRealPair HypergeometricRandomVariable::distribution_bounds() const
{
  unsigned int npr = numSelectPop + numDrawn,
    l_bnd = (npr > numTotalPop) ? npr - numTotalPop : 0;// care w/ unsigned math
  unsigned int u_bnd = std::min(numSelectPop, numDrawn);
  return RealRealPair((Real)l_bnd, (Real)u_bnd);
}


inline void HypergeometricRandomVariable::update_boost()
{
  hypergeomDist.reset
    (new hypergeometric_dist(numDrawn, numSelectPop, numTotalPop));
}


inline void HypergeometricRandomVariable::update_boost_conditionally()
{
  // old is now invalid
  if (hypergeomDist) hypergeomDist.reset();
  // new may not be valid as of yet
  if (numDrawn <= numTotalPop && numSelectPop <= numTotalPop)
    hypergeomDist.reset
      (new hypergeometric_dist(numDrawn, numSelectPop, numTotalPop));
  // else wait for pending param updates
}


inline void HypergeometricRandomVariable::
update(unsigned int num_total_pop, unsigned int num_sel_pop,
       unsigned int num_drawn)
{
  if (!hypergeomDist || numTotalPop != num_total_pop ||
      numSelectPop != num_sel_pop || numDrawn != num_drawn) {
    numTotalPop = num_total_pop;  numSelectPop = num_sel_pop;
    numDrawn    = num_drawn;      update_boost();
  }
}

// Static member functions:

inline Real HypergeometricRandomVariable::
pdf(Real x, unsigned int num_total_pop, unsigned int num_sel_pop,
    unsigned int num_drawn)
{
  hypergeometric_dist hypergeometric1(num_drawn, num_sel_pop, num_total_pop);
  return bmth::pdf(hypergeometric1, x);
}


inline Real HypergeometricRandomVariable::
cdf(Real x, unsigned int num_total_pop, unsigned int num_sel_pop,
    unsigned int num_drawn)
{
  hypergeometric_dist hypergeometric1(num_drawn, num_sel_pop, num_total_pop);
  return bmth::cdf(hypergeometric1, x);
}


inline void HypergeometricRandomVariable::
moments_from_params(unsigned int num_total_pop, unsigned int num_sel_pop,
		    unsigned int num_drawn, Real& mean, Real& std_dev)
{
  mean    = (Real)(num_drawn*num_sel_pop)/(Real)num_total_pop;
  std_dev = std::sqrt(mean
	  * (Real)((num_total_pop-num_drawn)*(num_total_pop-num_sel_pop))
	  / (Real)(num_total_pop*(num_total_pop-1)));
}

} // namespace Pecos

#endif
