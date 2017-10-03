/*
 * DensityEstimator.hpp
 *
 *  Created on: Apr 6, 2015
 *      Author: franzefn
 */

#ifndef DENSITY_ESTIMATOR_HPP_
#define DENSITY_ESTIMATOR_HPP_

#include "pecos_data_types.hpp"

namespace Pecos {

  /// Base class for density estimation strategies

  class DensityEstimator {
  public:
    /// default constructor
    DensityEstimator();
    /// standard constructor for envelope
    DensityEstimator(const String& density_estimator_type);
    /// copy constructor
    DensityEstimator(const DensityEstimator& density_estimator);

    /// destructor
    virtual ~DensityEstimator();

    /// assignment operator
    DensityEstimator operator=(const DensityEstimator& density_estimator);

    /// initialize the densities
    virtual void initialize(RealMatrix& samples,  
			    Teuchos::ETransp trans = Teuchos::NO_TRANS );
    virtual void initialize(RealVectorArray& samples);

    /// get the dimensionality
    virtual size_t getDim();

    /// mean
    virtual Real mean();
    /// variance
    virtual Real variance();
    /// standard deviation
    virtual Real std_deviation();
    /// computes the covariance matrix
    virtual void cov(RealMatrix& cov);
    /// computes the correlation matrix
    virtual void corrcoeff(RealMatrix& corr);

    /// operations for one sample
    virtual Real pdf(const RealVector& x) const;
    /// operations for a set of samples
    virtual void pdf(const RealMatrix& data, RealVector& res,
		     Teuchos::ETransp trans = Teuchos::NO_TRANS ) const;

    /// marginalization operations

    /// marginalizes over dim
    virtual void marginalize(size_t dim, DensityEstimator& estimator);
    /// creates a maginalized density with remaining dims
    virtual void margToDimXs(const IntVector& dims,
			     DensityEstimator& estimator);
    /// creates a maginalized density with one remaining dimension
    virtual void margToDimX(size_t dim, DensityEstimator& estimator);

    /// conditionalization operations
    virtual void conditionalize(const RealVector& x, const IntVector& dims,
				DensityEstimator& estimator);
    virtual void condToDimX(const RealVector& x, size_t dim,
			    DensityEstimator& estimator);

    String getType();
    DensityEstimator* getEnvelope();
    bool is_null();

  protected:
    String density_estimator_type;

  private:
    static DensityEstimator* get_density_estimator(
						   const String& density_estimator_type);

    DensityEstimator* densityEstimator;

    /// number of objects sharing samples
    int referenceCount;
  };

} /* namespace Pecos */

#endif /* DENSITY_ESTIMATOR_HPP_ */
