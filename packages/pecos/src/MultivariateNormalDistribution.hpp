/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef MULTIVARIATE_NORMAL_DISTRIBUTION_HPP
#define MULTIVARIATE_NORMAL_DISTRIBUTION_HPP

#include "MultivariateDistribution.hpp"
#include "NormalRandomVariable.hpp"

namespace Pecos {


/// Derived multivariate distribution class for a multivariate normal (MVN)

/** An MVN can be compactly defined by a mean vector and a symmetric
    covariance matric. */

class MultivariateNormalDistribution: public MultivariateDistribution
{
public:

  //
  //- Heading: Constructors and destructor
  //

  MultivariateNormalDistribution();                         ///< constructor
  MultivariateNormalDistribution(const RealVector& means,
				 const RealSymMatrix& cov); ///< alt constructor
  ~MultivariateNormalDistribution();                        ///< destructor

  //
  //- Heading: Member functions
  //

  void update(const RealVector& means, const RealSymMatrix& cov);

  void initialize_correlations();

  /// update vector values
  void push_parameters(short dist_param, const RealVector& values);
  /// update symmetric matrix values
  void push_parameters(short dist_param, const RealSymMatrix& values);

  // migrated from NormalRandomVariable static fns:
  static Real std_pdf(Real beta, size_t n);
  static Real std_pdf(const RealVector& u);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// return the multivariate PDF value for full set of random variables
  Real pdf(const RealVector& pt) const;
  /// return the multivariate log PDF value for full set of random variables
  Real log_pdf(const RealVector& pt) const;

  //
  //- Heading: Data
  //

  /// vector of mean values for multivariate Normal distribution
  RealVector mvnMeans;
  /// symmetric covariance matrix for multivariate Normal distribution
  RealSymMatrix mvnCovariance;
};


inline MultivariateNormalDistribution::MultivariateNormalDistribution():
  MultivariateDistribution(BaseConstructor())
{ }


inline MultivariateNormalDistribution::
MultivariateNormalDistribution(const RealVector& means,
			       const RealSymMatrix& cov):
  MultivariateDistribution(BaseConstructor()),
  mvnMeans(means), mvnCovariance(cov)
{ initialize_correlations(); }


inline MultivariateNormalDistribution::~MultivariateNormalDistribution()
{ }


inline void MultivariateNormalDistribution::
update(const RealVector& means, const RealSymMatrix& cov)
{ mvnMeans = means; mvnCovariance = cov; initialize_correlations(); }


inline void MultivariateNormalDistribution::initialize_correlations()
{
  size_t i, j, num_rv = mvnCovariance.numRows();
  correlationFlag = false;
  for (i=1; i<num_rv; i++)
    for (j=0; j<i; j++)
      if (std::abs(mvnCovariance(i,j)) > SMALL_NUMBER)
	{ correlationFlag = true; break; }
}


inline void MultivariateNormalDistribution::
push_parameters(short dist_param, const RealVector& values)
{
  switch (dist_param) {
  case N_MEAN: case N_LOCATION:  mvnMeans = values;  break;
  default:
    PCerr << "Error: lookup failure for distribution parameter " << dist_param
	  << " in MultivariateNormalDistribution::parameters(RealVector)."
	  << std::endl;
    abort_handler(-1); break;
  }
}


inline void MultivariateNormalDistribution::
push_parameters(short dist_param, const RealSymMatrix& values)
{
  bool err_flag = false;
  switch (dist_param) {
  case N_VARIANCE: //case N_STD_DEV: case N_SCALE:
    mvnCovariance = values;  initialize_correlations();  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in MultivariateNormalDistribution::parameters(RealSymMatrix)."
	  << std::endl;
    abort_handler(-1); break;
  }
}


inline Real MultivariateNormalDistribution::pdf(const RealVector& x) const
{
  // *** TO DO:
  //   std::exp(- (x-mvnMeans)^T inverse_mvnCovar * (x-mvnMeans) / 2.)
  // / std::sqrt( determinant(mvnCovar) * (2.*Pi)^d )
  return 0.; //
}


// Multivariate standard normal density function from vector.
inline Real MultivariateNormalDistribution::log_pdf(const RealVector& x) const
{
  // *** TO DO:
  //   std::log(term) - (x-mvnMeans)^T inverse_mvnCovar * (x-mvnMeans) / 2.
  return 0.;
}


// *** The following are fns migrated from NormalRandomVariable::mvn_*()

// Multivariate standard normal density function with aggregate distance.
inline Real MultivariateNormalDistribution::std_pdf(Real beta, size_t n)
{
  // need n instances of 1/sqrt(2Pi), but 1D pdf only includes 1:
  return (n > 1) ? // correct the 1D pdf for n dimensions
    NormalRandomVariable::std_pdf(beta) * std::pow(2.*PI, -((Real)(n-1))/2.) :
    NormalRandomVariable::std_pdf(beta);
}


// Multivariate standard normal density function from vector.
inline Real MultivariateNormalDistribution::std_pdf(const RealVector& u)
{
  return std_pdf(u.normFrobenius(), u.length());

  // Alternate implementation invokes exp() repeatedly:
  //normal_dist norm(0., 1.);
  //size_t i, n = u.length(); Real pdf = 1.;
  //for (i=0; i<n; ++i)
  //  pdf *= bmth::pdf(norm, u[i]);
}

} // namespace Pecos

#endif
