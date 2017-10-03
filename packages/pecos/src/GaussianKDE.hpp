/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef GAUSSIAN_KDE_HPP_
#define GAUSSIAN_KDE_HPP_

#include "DensityEstimator.hpp"
#include "RosenblattTransformation.hpp"
#include "pecos_data_types.hpp"

namespace Pecos {

#ifndef M_SQRT2PI
#define M_SQRT2PI       2.5066282746310002             /* sqrt(2*pi) */
#endif

  /// Class for kernel density estimation with gaussian kernels.

  /** The kernel density estimation method estimates an unknown density based
   * on realizations of it. The estimated probability density function is
   * defined as
   *
   *        f(x) = 1/n \sum_{i=1}^n K((x - \mu_i) / \sigma)
   *
   * where \mu_i are the realizations and \sigma is the bandwidth of the gaussian
   * kernels K.
   **/

  class GaussianKDE: public DensityEstimator {
  public:

    //
    //- Heading: Constructors and destructor
    //

    /// default constructor
    GaussianKDE();

    /// destructor
    ~GaussianKDE();

    /// initialize the density estimator
    void initialize(RealMatrix& samples,
		    Teuchos::ETransp trans = Teuchos::NO_TRANS );
    void initialize(RealVectorArray& samples);

    /// get the dimensionality
    size_t getDim();

    /// computes the sample mean
    Real mean();
    /// computed according to [Goodman, 1962]
    Real variance();
    /// computes the standard deviation
    Real std_deviation();

    /// computes the covariance matrix
    void cov(RealMatrix& cov);

    /// operations for single samples
    Real pdf(const RealVector& x) const;

    /// operations for a set of samples
    void pdf(const RealMatrix& data, RealVector& res,
	     Teuchos::ETransp trans = Teuchos::NO_TRANS) const;

    /// marginalization operations
    virtual void marginalize(size_t dim, DensityEstimator& estimator);
    virtual void margToDimXs(const IntVector& dims,
			     DensityEstimator& estimator);
    virtual void margToDimX(size_t dim, DensityEstimator& estimator);

    /// conditionalization operations
    virtual void conditionalize(const RealVector& x, const IntVector& dims,
				DensityEstimator& estimator);
    virtual void condToDimX(const RealVector& x, size_t dim,
			    DensityEstimator& estimator);

    /// getter and setter functions
    void getConditionalizationFactor(RealVector& pcond);
    void setConditionalizationFactor(const RealVector& pcond);
    void getBandwidths(RealVector& bandwidths);

    /// get samples
    const RealVectorArray& getSamples() const;

  private:
    void computeOptKDEbdwth();

    Real getSampleMean(RealVector& data);
    Real getSampleVariance(RealVector& data);
    Real getSampleStd(RealVector& data);

    void updateConditionalizationFactors(const RealVector& x,
					 const IntVector& dims, 
					 RealVector& cond);

  protected:
    /// samples
    RealVectorArray samplesVec;

    size_t nsamples;
    size_t ndim;

    /// standard deviations for the kernels in 1d
    RealVector bandwidths;
    /// normalization factor for 1d kernels
    RealVector norm;
    /// conditionalization factors
    RealVector cond;
    Real sumCond;
  }
    ;

} /* namespace Pecos */

#endif /* GAUSSIAN_KDE_HPP_ */
