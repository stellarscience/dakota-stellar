/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 RandomVariable
//- Description: 
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef RANDOM_VARIABLE_HPP
#define RANDOM_VARIABLE_HPP

#include "pecos_stat_util.hpp"

namespace Pecos {


/// base class for random variable hierarchy

/** This class enables cdf(), ccdf(), inverse_cdf(), inverse_ccdf(),
    pdf(), pdf_gradient(), pdf_hessian(), and related random variable
    utilities from contained distribution parameters. */

// ExtremeValueBase: Gumbel, Frechet, Weibull -> alpha, beta, no bounds
// EmpiricalBase: histogram bin, KDE, ...

class RandomVariable
{
public:

  //
  //- Heading: Constructors, destructor, and operator=
  //

  /// default constructor
  RandomVariable();
  /// standard constructor for envelope
  RandomVariable(short ran_var_type);
  /// copy constructor
  RandomVariable(const RandomVariable& ran_var);

  /// destructor
  virtual ~RandomVariable();

  /// assignment operator
  RandomVariable operator=(const RandomVariable& ran_var);

  //
  //- Heading: Virtual functions
  //

  /// return the cumulative distribution function value of the random
  /// variable at x
  virtual Real cdf(Real x) const;
  /// return the x value corresponding to a cumulative probability
  virtual Real inverse_cdf(Real p_cdf) const;

  /// return the complementary cumulative distribution function value
  /// of the random variable at x
  virtual Real ccdf(Real x) const;
  /// return the x value corresponding to a complementary cumulative probability
  virtual Real inverse_ccdf(Real p_ccdf) const;

  /// return the value of the random variable's probability density
  /// function at x
  virtual Real pdf(Real x) const;
  /// return the gradient of the random variable's probability density
  /// function at x
  virtual Real pdf_gradient(Real x) const;
  /// return the hessian of the random variable's probability density
  /// function at x
  virtual Real pdf_hessian(Real x) const;
  /// return the value of the natural log of the random variable's probability
  /// density function at x (useful for calculations of log density in Bayesian
  /// methods)
  virtual Real log_pdf(Real x) const;
  /// return the gradient of the natural log of the random variable's
  /// probability density function at x (useful for defining MCMC proposal
  /// distributions in Bayesian methods)
  virtual Real log_pdf_gradient(Real x) const;
  /// return the Hessian of the natural log of the random variable's probability
  /// density function at x (useful for defining MCMC proposal distributions in
  /// Bayesian methods)
  virtual Real log_pdf_hessian(Real x) const;

  /// return the x value for a standardized probability distribution
  /// corresponding to a cumulative probability
  virtual Real inverse_standard_cdf(Real p_cdf) const;

  /// return the value of a standardized random variable's probability density
  /// function at x
  virtual Real standard_pdf(Real z) const;
  /// return the natural log of a standardized random variable's probability
  /// density function at x (useful for calculations of log density in
  /// Bayesian methods)
  virtual Real log_standard_pdf(Real z) const;
  /// return the gradient of the natural log of a standardized random
  /// variable's probability density function at x (useful for
  /// calculations of log density in Bayesian methods)
  virtual Real log_standard_pdf_gradient(Real z) const;
  /// return the Hessian of the natural log of a standardized random
  /// variable's probability density function at x (useful for
  /// calculations of log density in Bayesian methods)
  virtual Real log_standard_pdf_hessian(Real z) const;

  /// scale variable value x from current to standardized distribution
  virtual Real to_standard(Real x) const;
  /// scale variable value z from standardized to current distribution
  virtual Real from_standard(Real z) const;

  /// return the value of the named distribution parameter
  virtual void pull_parameter(short dist_param, Real& val) const;
  /// return the value of the named distribution parameter
  virtual void pull_parameter(short dist_param, int& val) const;
  /// return the value of the named distribution parameter
  virtual void pull_parameter(short dist_param, unsigned int& val) const;
  /// return the value of the named distribution parameter
  virtual void pull_parameter(short dist_param, IntSet& val) const;
  /// return the value of the named distribution parameter
  virtual void pull_parameter(short dist_param, StringSet& val) const;
  /// return the value of the named distribution parameter
  virtual void pull_parameter(short dist_param, RealSet& val) const;
  /// return the value of the named distribution parameter
  virtual void pull_parameter(short dist_param, IntRealMap& val) const;
  /// return the value of the named distribution parameter
  virtual void pull_parameter(short dist_param, StringRealMap& val) const;
  /// return the value of the named distribution parameter
  virtual void pull_parameter(short dist_param, RealRealMap& val) const;
  /// return the value of the named distribution parameter
  virtual void pull_parameter(short dist_param, IntIntPairRealMap& val) const;
  /// return the value of the named distribution parameter
  virtual void pull_parameter(short dist_param, RealRealPairRealMap& val) const;

  /// update the value of the named distribution parameter
  virtual void push_parameter(short dist_param, Real val);
  /// update the value of the named distribution parameter
  virtual void push_parameter(short dist_param, int val);
  /// update the value of the named distribution parameter
  virtual void push_parameter(short dist_param, unsigned int val);
  /// update the value of the named distribution parameter
  virtual void push_parameter(short dist_param, const IntSet& val);
  /// update the value of the named distribution parameter
  virtual void push_parameter(short dist_param, const StringSet& val);
  /// update the value of the named distribution parameter
  virtual void push_parameter(short dist_param, const RealSet& val);
  /// update the value of the named distribution parameter
  virtual void push_parameter(short dist_param, const IntRealMap& val);
  /// update the value of the named distribution parameter
  virtual void push_parameter(short dist_param, const StringRealMap& val);
  /// update the value of the named distribution parameter
  virtual void push_parameter(short dist_param, const RealRealMap& val);
  /// update the value of the named distribution parameter
  virtual void push_parameter(short dist_param, const IntIntPairRealMap& val);
  /// update the value of the named distribution parameter
  virtual void push_parameter(short dist_param, const RealRealPairRealMap& val);

  /// copy the distribution parameter values from rv into *this
  virtual void copy_parameters(const RandomVariable& rv);

  /// return the distribution mean
  virtual Real mean() const;
  /// return the distribution mode
  virtual Real median() const;
  /// return the distribution mode
  virtual Real mode() const;
  /// return the distribution variance
  virtual Real standard_deviation() const;
  /// return the distribution variance
  virtual Real variance() const;

  /// return the distribution mean and standard deviation as a pair
  /** default is only overridden when more efficient to compute together */
  virtual RealRealPair moments() const;

  /// return the distribution lower and upper bounds as a pair
  virtual RealRealPair distribution_bounds() const;

  //virtual Real pull_lower_bound() const;
  //virtual int  pull_lower_bound() const;
  //virtual void push_lower_bound(Real l_bnd);
  //virtual void push_lower_bound(int  l_bnd);
  virtual void lower_bound(Real l_bnd);
  virtual void lower_bound(int  l_bnd);
  virtual void upper_bound(Real l_bnd);
  virtual void upper_bound(int  l_bnd);

  /// compute the coefficient of variation (used to compute selected
  /// correlation warping factors); defined for semi-infinite distributions
  /// with nonzero mean (lognormal, exponential, gamma, frechet, weibull)
  /** default is only overridden when more efficient to compute together */
  virtual Real coefficient_of_variation() const;
  /// compute the warping factor for correlation between the current
  /// variable and the one passed in (used in NatafTransformation)
  virtual Real correlation_warping_factor(const RandomVariable& rv,
					  Real corr) const;

  /// compute the design Jacobian from differentiating the X->Z mapping with
  /// respect to the distibution parameter s
  virtual Real dx_ds(short dist_param, short u_type, Real x, Real z) const;
  /// compute the mapping-specific factor that is multiplied by dz/ds for
  /// contributions to the dx/ds design Jacobian in the case of correlated
  /// random variables (dz/ds is evaluated numerically and multiplied by
  /// this analytic factor)
  virtual Real dz_ds_factor(short u_type, Real x, Real z) const;

  //
  //- Heading: Member functions
  //

  // Invoke virtual pull_parameter(short, T) and return result
  //template <typename T>
  //T return_parameter(short dist_param) const;

  /// Draw a sample from the distribution using inverse_cdf on uniform[0,1]
  template <typename Engine>
  Real draw_sample(Engine& rng) const;
  /// Draw a sample from the distribution using inverse_cdf on uniform[0,1]
  template <typename Engine>
  Real draw_standard_sample(Engine& rng) const;
  
  /// set ranVarType
  void type(short ran_var_type);
  /// get ranVarType
  short type() const;

  /// returns ranVarRep for access to derived class member functions
  /// that are not mapped to the base level
  std::shared_ptr<RandomVariable> random_variable_rep() const;
  /// function to check modelRep (does this envelope contain a letter)
  bool is_null() const;

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  RandomVariable(BaseConstructor);

  //
  //- Heading: Member functions
  //

  //
  //- Heading: Data
  //

  /// enumeration value indicating type of random variable
  short ranVarType;

private:

  //
  //- Heading: Member functions
  //

  /// Used only by the standard envelope constructor to initialize
  /// ranVarRep to the appropriate derived type.
  std::shared_ptr<RandomVariable> get_random_variable(short ran_var_type);

  //
  //- Heading: Data members
  //

  /// draws real samples on [0,1]
  boost::uniform_real<Real> uniformSampler;
  
  /// pointer to the letter (initialized only for the envelope)
  std::shared_ptr<RandomVariable> ranVarRep;
};


/*
template <typename T>
T RandomVariable::return_parameter(short dist_param) const
{
  T val;
  pull_parameter(dist_param, val);
  return val;
}
*/


template <typename Engine> 
Real RandomVariable::draw_sample(Engine& rng) const
{
  if (ranVarRep)
    return ranVarRep->draw_sample(rng);
  else {
    // draw random number on [0,1] from a persistent RNG sequence
    Real u01 = uniformSampler(rng);
    return inverse_cdf(u01);
  }
}


template <typename Engine> 
Real RandomVariable::draw_standard_sample(Engine& rng) const
{
  if (ranVarRep)
    return ranVarRep->draw_standard_sample(rng);
  else {
    // draw random number on [0,1] from a persistent RNG sequence
    Real u01 = uniformSampler(rng);
    return inverse_standard_cdf(u01);
  }
}


inline void RandomVariable::type(short ran_var_type)
{
  if (ranVarRep) ranVarRep->ranVarType = ran_var_type;
  else           ranVarType = ran_var_type;
}


inline short RandomVariable::type() const
{ return (ranVarRep) ? ranVarRep->ranVarType : ranVarType; }


inline std::shared_ptr<RandomVariable>
RandomVariable::random_variable_rep() const
{ return ranVarRep; }


inline bool RandomVariable::is_null() const
{ return (ranVarRep) ? false : true; }

} // namespace Pecos

#endif
