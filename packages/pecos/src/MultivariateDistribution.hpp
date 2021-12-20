/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef MULTIVARIATE_DISTRIBUTION_HPP
#define MULTIVARIATE_DISTRIBUTION_HPP

#include "pecos_data_types.hpp"

namespace Pecos {

class RandomVariable;


/// Base class for all multivariate distributions

/** The base class for multivariate distributions, including joint
    distributions and distributions based on marginals plus correlation. */

class MultivariateDistribution
{
public:

  /// default constructor
  MultivariateDistribution();
  /// standard constructor for envelope
  MultivariateDistribution(short mv_dist_type);
  /// copy constructor
  MultivariateDistribution(const MultivariateDistribution& mv_dist);

  /// destructor
  virtual ~MultivariateDistribution();

  /// assignment operator
  MultivariateDistribution operator=(const MultivariateDistribution& mv_dist);

  //
  //- Heading: Virtual functions
  //

  /// return randomVars (marginals, when present)
  virtual const std::vector<RandomVariable>& random_variables() const;
  /// return randomVars (marginals, when present)
  virtual std::vector<RandomVariable>& random_variables();

  /// return randomVars[i] (marginal, when present)
  virtual const RandomVariable& random_variable(size_t i) const;
  /// return randomVars[i] (marginal, when present)
  virtual RandomVariable& random_variable(size_t i);

  /// return ranVarTypes (marginals, when present)
  virtual const ShortArray& random_variable_types() const;
  /// set ranVarTypes (marginals, when present)
  virtual void random_variable_types(const ShortArray& rv_types);

  /// return ranVarTypes[i] (marginal, when present)
  virtual short random_variable_type(size_t i) const;
  /// set ranVarTypes[i] (marginal, when present)
  virtual void random_variable_type(short rv_type, size_t i);

  /// return active subset of variables
  virtual const BitArray& active_variables() const;
  /// return subset of variables to which correlation matrix applies
  virtual const BitArray& active_correlations() const;

  /// return corrMatrix
  virtual const RealSymMatrix& correlation_matrix() const;
  /// set corrMatrix
  virtual void correlation_matrix(const RealSymMatrix& corr);

  /// pull non-standardized distribution parameters from mv_dist to this,
  /// assuming a consistent variable ordering
  virtual void
    pull_distribution_parameters(const MultivariateDistribution& mv_dist);
  /// pull non-standardized distribution parameters from mv_dist to this,
  /// using label lookups to support potentially inconsistent ordering
  virtual void
    pull_distribution_parameters(const MultivariateDistribution& mv_dist,
      const StringArray& pull_labels, const StringArray& push_labels);

  /// return marginal means,standard deviations for multivariate distribution
  virtual RealRealPairArray moments() const;
  /// return marginal means from multivariate distribution
  virtual RealVector means() const;
  /// return marginal standard deviations for multivariate distribution
  virtual RealVector std_deviations() const;
  /// return marginal variances for multivariate distribution
  virtual RealVector variances() const;

  /// return lower and upper bounds for multivariate distribution
  virtual RealRealPairArray distribution_bounds() const;
  /// return lower bounds for multivariate distribution
  virtual RealVector distribution_lower_bounds() const;
  /// return upper bounds for multivariate distribution
  virtual RealVector distribution_upper_bounds() const;

  /// indicate presence of data beyond distribution bounds, which define the
  /// "global" parameter space
  virtual bool global_bounds() const;

  /// set lower bounds for multivariate distribution
  virtual void lower_bounds(const RealVector& l_bnds,
			    const BitArray& mask = BitArray());
  /// set lower bounds for multivariate distribution
  virtual void lower_bounds(const IntVector& l_bnds,
			    const BitArray& mask = BitArray());
  /// set lower bound for a variable within multivariate distribution
  virtual void lower_bound(Real l_bnd, size_t rv_index);
  /// set lower bound for a variable within multivariate distribution
  virtual void lower_bound(int l_bnd, size_t rv_index);

  /// set upper bounds for multivariate distribution
  virtual void upper_bounds(const RealVector& u_bnds,
			    const BitArray& mask = BitArray());
  /// set upper bounds for multivariate distribution
  virtual void upper_bounds(const IntVector& u_bnds,
			    const BitArray& mask = BitArray());
  /// set upper bound for a variable within multivariate distribution
  virtual void upper_bound(Real u_bnd, size_t rv_index);
  /// set upper bound for a variable within multivariate distribution
  virtual void upper_bound(int u_bnd, size_t rv_index);

  /// return the multivariate PDF value for the random variables
  virtual Real pdf(const RealVector& pt) const;
  /// return the multivariate log PDF value for the random variables
  virtual Real log_pdf(const RealVector& pt) const;

  /// return the univariate PDF value for a random variable
  virtual Real pdf(Real val, size_t i) const;
  /// return the gradient of the univariate PDF value for a random variable
  virtual Real pdf_gradient(Real val, size_t i) const;
  /// return the the Hessian of the univariate PDF value for a random variable
  virtual Real pdf_hessian(Real val, size_t i) const;
  /// return the univariate log PDF value for a random variable
  virtual Real log_pdf(Real val, size_t i) const;
  /// return the gradient of the univariate log PDF for a random variable
  virtual Real log_pdf_gradient(Real val, size_t i) const;
  /// return the Hessian of the univariate log PDF for a random variable
  virtual Real log_pdf_hessian(Real val, size_t i) const;

  //
  //- Heading: Member functions
  //

  /// return index of i-th active variable
  size_t active_variable_index(size_t i) const;
  /// finds index of i-th active variable and returns randomVars[index]
  const RandomVariable& active_random_variable(size_t i) const;
  /// finds index of i-th active variable and returns randomVars[index]
  RandomVariable& active_random_variable(size_t i);
  /// finds index of i-th active variable and returns ranVarTypes[index]
  short active_random_variable_type(size_t i) const;

  /// return correlationFlag
  bool correlation() const;

  /// perform a deep copy of this object and return by value
  MultivariateDistribution copy() const;

  /// returns ranVarRep for access to derived class member functions
  /// that are not mapped to the base level
  MultivariateDistribution* multivar_dist_rep() const;
  /// function to check modelRep (does this envelope contain a letter)
  bool is_null() const;

  /// return mvDistType
  short type() const;

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  MultivariateDistribution(BaseConstructor);
  
  //
  //- Heading: Virtual functions
  //

  /// copy the data from mvd_rep to the current representation
  virtual void copy_rep(MultivariateDistribution* mvd_rep);

  //
  //- Heading: Data members
  //

  /// derived type of MultivariateDistribution
  short mvDistType;

  /// flag for indicating if correlation exists among the random variables
  bool correlationFlag;

private:

  //
  //- Heading: Member functions
  //

  /// Used only by the standard envelope constructor to initialize
  /// mvDistRep to the appropriate derived type.
  MultivariateDistribution* get_distribution(short mv_dist_type) const;

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  MultivariateDistribution* mvDistRep;
  /// number of objects sharing mvDistRep
  int referenceCount;
};


inline short MultivariateDistribution::type() const
{ return mvDistType; }


inline bool MultivariateDistribution::correlation() const
{ return (mvDistRep) ? mvDistRep->correlationFlag : correlationFlag; }


inline MultivariateDistribution* MultivariateDistribution::
multivar_dist_rep() const
{ return mvDistRep; }


inline bool MultivariateDistribution::is_null() const
{ return (mvDistRep) ? false : true; }

} // namespace Pecos

#endif
