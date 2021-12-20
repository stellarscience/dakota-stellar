/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 PoissonRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef POISSON_RANDOM_VARIABLE_HPP
#define POISSON_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for poisson random variables.

/** Manages the lambda parameter. */

class PoissonRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  PoissonRandomVariable();
  /// alternate constructor
  PoissonRandomVariable(Real lambda);
  /// destructor
  ~PoissonRandomVariable();

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

  void update(Real lambda);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, Real lambda);
  static Real cdf(Real x, Real lambda);

  static void moments_from_params(Real lambda, Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Member functions
  //

  /// create a new poissonDist instance
  void update_boost();

  //
  //- Heading: Data
  //

  /// lambda parameter of poisson random variable
  Real poissonLambda;

  /// pointer to the Boost poisson_distribution instance
  std::unique_ptr<poisson_dist> poissonDist;
};


inline PoissonRandomVariable::PoissonRandomVariable():
  RandomVariable(BaseConstructor()), poissonLambda(1.),
  poissonDist(new poisson_dist(poissonLambda))
{ ranVarType = POISSON; }


inline PoissonRandomVariable::PoissonRandomVariable(Real lambda):
  RandomVariable(BaseConstructor()), poissonLambda(lambda),
  poissonDist(new poisson_dist(lambda))
{ ranVarType = POISSON; }


inline PoissonRandomVariable::~PoissonRandomVariable()
{ }


inline Real PoissonRandomVariable::cdf(Real x) const
{ return bmth::cdf(*poissonDist, x); }


inline Real PoissonRandomVariable::ccdf(Real x) const
{ return bmth::cdf(complement(*poissonDist, x)); }


inline Real PoissonRandomVariable::inverse_cdf(Real p_cdf) const
{ return bmth::quantile(*poissonDist, p_cdf); }


inline Real PoissonRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return bmth::quantile(complement(*poissonDist, p_ccdf)); }


inline Real PoissonRandomVariable::pdf(Real x) const
{ return bmth::pdf(*poissonDist, x); }


inline void PoissonRandomVariable::
pull_parameter(short dist_param, Real& val) const
{
  switch (dist_param) {
  case P_LAMBDA: val = poissonLambda; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in PoissonRandomVariable::pull_parameter(Real)." << std::endl;
    abort_handler(-1); break;
  }
}


inline void PoissonRandomVariable::push_parameter(short dist_param, Real val)
{
  switch (dist_param) {
  case P_LAMBDA: poissonLambda = val; break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in PoissonRandomVariable::push_parameter(Real)." << std::endl;
    abort_handler(-1); break;
  }
  update_boost(); // create a new poissonDist instance
}


inline void PoissonRandomVariable::copy_parameters(const RandomVariable& rv)
{
  rv.pull_parameter(P_LAMBDA, poissonLambda);
  update_boost(); // create a new poissonDist instance
}


inline Real PoissonRandomVariable::mean() const
{ return bmth::mean(*poissonDist); }


inline Real PoissonRandomVariable::median() const
{ return bmth::median(*poissonDist); }


inline Real PoissonRandomVariable::mode() const
{ return bmth::mode(*poissonDist); }


inline Real PoissonRandomVariable::standard_deviation() const
{ return bmth::standard_deviation(*poissonDist); }


inline Real PoissonRandomVariable::variance() const
{ return bmth::variance(*poissonDist); }


inline RealRealPair PoissonRandomVariable::distribution_bounds() const
{ return RealRealPair(0., std::numeric_limits<Real>::infinity()); }


inline void PoissonRandomVariable::update_boost()
{
  poissonDist.reset(new poisson_dist(poissonLambda));
}


inline void PoissonRandomVariable::update(Real lambda)
{
  if (!poissonDist || poissonLambda != lambda)
    { poissonLambda = lambda; update_boost(); }
}

// Static member functions:

inline Real PoissonRandomVariable::pdf(Real x, Real lambda)
{ poisson_dist poisson1(lambda); return bmth::pdf(poisson1, x); }


inline Real PoissonRandomVariable::cdf(Real x, Real lambda)
{ poisson_dist poisson1(lambda); return bmth::cdf(poisson1, x); }


inline void PoissonRandomVariable::
moments_from_params(Real lambda, Real& mean, Real& std_dev)
{ mean = lambda; std_dev = std::sqrt(lambda); }

} // namespace Pecos

#endif
