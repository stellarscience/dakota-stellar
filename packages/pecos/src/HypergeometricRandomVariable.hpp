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

/** Manages numTotalPop, numSelectPop, and numFail parameters. */

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
			       unsigned int num_sel_pop, unsigned int num_fail);
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

  void update(unsigned int num_total_pop, unsigned int num_sel_pop,
	      unsigned int num_fail);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, unsigned int num_total_pop, unsigned int num_sel_pop,
		  unsigned int num_fail);
  static Real cdf(Real x, unsigned int num_total_pop, unsigned int num_sel_pop,
		  unsigned int num_fail);

  static void moments_from_params(unsigned int num_total_pop,
				  unsigned int num_sel_pop,
				  unsigned int num_fail,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Member functions
  //

  /// create a new hypergeomDist instance
  void update_boost();

  //
  //- Heading: Data
  //

  /// N parameter of hypergeometric random variable
  unsigned int numTotalPop;
  /// n parameter of hypergeometric random variable
  unsigned int numSelectPop;
  /// r parameter of hypergeometric random variable
  unsigned int numFail;

  /// pointer to the Boost hypergeometric_distribution instance
  hypergeometric_dist* hypergeomDist;
};


inline HypergeometricRandomVariable::HypergeometricRandomVariable():
  RandomVariable(BaseConstructor()), numTotalPop(1), numSelectPop(1),
  numFail(1), hypergeomDist(NULL)
{ ranVarType = HYPERGEOMETRIC; }


inline HypergeometricRandomVariable::
HypergeometricRandomVariable(unsigned int num_total_pop,
			     unsigned int num_sel_pop, unsigned int num_fail):
  RandomVariable(BaseConstructor()), numTotalPop(num_total_pop),
  numSelectPop(num_sel_pop), numFail(num_fail),
  hypergeomDist(new hypergeometric_dist(numFail, numSelectPop, numTotalPop))
{ ranVarType = HYPERGEOMETRIC; }


inline HypergeometricRandomVariable::~HypergeometricRandomVariable()
{ if (hypergeomDist) delete hypergeomDist; }


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


inline Real HypergeometricRandomVariable::parameter(short dist_param) const
{
  switch (dist_param) {
  case HGE_TOT_POP: return (Real)numTotalPop;  break;
  case HGE_SEL_POP: return (Real)numSelectPop; break;
  case HGE_FAILED:  return (Real)numFail;      break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in HypergeometricRandomVariable::parameter()." << std::endl;
    abort_handler(-1); return 0.; break;
  }
}


inline void HypergeometricRandomVariable::parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case HGE_TOT_POP: numTotalPop  = (unsigned int)std::floor(val+.5); break;
  case HGE_SEL_POP: numSelectPop = (unsigned int)std::floor(val+.5); break;
  case HGE_FAILED:  numFail      = (unsigned int)std::floor(val+.5); break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in HypergeometricRandomVariable::parameter()." << std::endl;
    abort_handler(-1); break;
  }
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


inline RealRealPair HypergeometricRandomVariable::bounds() const
{
  unsigned int npr = numSelectPop + numFail,
    l_bnd = (npr > numTotalPop) ? npr - numTotalPop : 0;// care w/ unsigned math
  unsigned int u_bnd = std::min(numSelectPop, numFail);
  return RealRealPair((Real)l_bnd, (Real)u_bnd);
}


inline void HypergeometricRandomVariable::update_boost()
{
  if (hypergeomDist) delete hypergeomDist;
  hypergeomDist = new hypergeometric_dist(numFail, numSelectPop, numTotalPop);
}


inline void HypergeometricRandomVariable::
update(unsigned int num_total_pop, unsigned int num_sel_pop,
       unsigned int num_fail)
{
  if (!hypergeomDist || numTotalPop != num_total_pop ||
      numSelectPop != num_sel_pop || numFail != num_fail) {
    numTotalPop = num_total_pop; numSelectPop = num_sel_pop; numFail = num_fail;
    update_boost();
  }
}

// Static member functions:

inline Real HypergeometricRandomVariable::
pdf(Real x, unsigned int num_total_pop, unsigned int num_sel_pop,
    unsigned int num_fail)
{
  hypergeometric_dist hypergeometric1(num_fail, num_sel_pop, num_total_pop);
  return bmth::pdf(hypergeometric1, x);
}


inline Real HypergeometricRandomVariable::
cdf(Real x, unsigned int num_total_pop, unsigned int num_sel_pop,
    unsigned int num_fail)
{
  hypergeometric_dist hypergeometric1(num_fail, num_sel_pop, num_total_pop);
  return bmth::cdf(hypergeometric1, x);
}


inline void HypergeometricRandomVariable::
moments_from_params(unsigned int num_total_pop, unsigned int num_sel_pop,
		    unsigned int num_fail, Real& mean, Real& std_dev)
{
  mean    = (Real)(num_fail*num_sel_pop)/(Real)num_total_pop;
  std_dev = std::sqrt(mean
	  * (Real)((num_total_pop-num_fail)*(num_total_pop-num_sel_pop))
	  / (Real)(num_total_pop*(num_total_pop-1)));
}

} // namespace Pecos

#endif
