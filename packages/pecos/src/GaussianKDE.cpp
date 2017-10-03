/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Fabian Franzelin (fabian.franzelin@ipvs.uni-stuttgart.de)
#include "GaussianKDE.hpp"
#include "RosenblattTransformation.hpp"
#include <cmath>
#include <limits>
#include <iostream>

using namespace std;

// #define REFCOUNT_DEBUG

namespace Pecos {

// -------------------- constructors and desctructors --------------------
GaussianKDE::GaussianKDE() :
        nsamples(0), ndim(0), sumCond(1.0) {
    // set the density estimator type
    density_estimator_type = "gaussian_kde";
}

GaussianKDE::~GaussianKDE() {
#ifdef REFCOUNT_DEBUG
    std::cout << "deleting GaussianKDE" << std::endl;
#endif
}
// ----------------------------------------------------------------------

void GaussianKDE::initialize(RealVectorArray& samples) {
    ndim = samples.size();
    if (ndim > 0) {
        nsamples = samples[0].length();
        if (nsamples > 1) {
            // copy 1d samples to vector
            samplesVec.resize(ndim);
            for (size_t idim = 0; idim < ndim; idim++) {
                samplesVec[idim] = samples[idim]; // copy
            }

            // init the bandwidths
            bandwidths.resize(ndim);
            computeOptKDEbdwth();

            // initialize normalization factors
            norm.resize(ndim);
            for (size_t d = 0; d < ndim; d++) {
                norm[d] = 1. / (bandwidths[d] * M_SQRT2PI);
            }

            // initialize conditionalization factor
            cond.resize(nsamples);
            cond.putScalar(1.0);
            sumCond = nsamples;
        } else {
            PCerr<< "Error: KDE needs at least two samples to estimate the bandwidth\n"
            << std::endl;
            abort_handler(-1);
        }
    } else {
        PCerr << "Error: KDE needs at least one dimensional data\n"
        << std::endl;
        abort_handler(-1);
    }
}

  void GaussianKDE::initialize(RealMatrix& samples, Teuchos::ETransp trans ) {
    if ( trans == Teuchos::NO_TRANS ){
      nsamples = samples.numRows();
      ndim = samples.numCols();
    }else{
      nsamples = samples.numCols();
      ndim = samples.numRows();
    }

    if (ndim > 0) {
        if (nsamples > 1) {
	  
            // copy 1d samples to vector
            samplesVec.resize(ndim);

	    for (size_t idim = 0; idim < ndim; idim++) {
	      samplesVec[idim].resize(nsamples);
	      for (size_t isample = 0; isample < nsamples; isample++) {
		if ( trans == Teuchos::NO_TRANS )
		  samplesVec[idim][isample] = samples(isample, idim);
		else
		  samplesVec[idim][isample] = samples(idim, isample);
	      }
	    }

            // init the bandwidths
            bandwidths.resize(ndim);
            computeOptKDEbdwth();

            // initialize normalization factors
            norm.resize(ndim);
            for (size_t d = 0; d < ndim; d++) {
                norm[d] = 1. / (bandwidths[d] * M_SQRT2PI);
            }

            // initialize conditionalization factor
            cond.resize(nsamples);
            cond.putScalar(1.0);
            sumCond = nsamples;
        } else {
            PCerr<< "Error: KDE needs at least two samples to estimate the bandwidth\n"
            << std::endl;
            abort_handler(-1);
        }
    } else {
        PCerr << "Error: KDE needs at least one dimensional data\n"
        << std::endl;
        abort_handler(-1);
    }
}

size_t GaussianKDE::getDim() {
    return ndim;
}

void GaussianKDE::getBandwidths(RealVector& bandwidths) {
    // copy
    bandwidths.resize(this->bandwidths.length());
    for (size_t i = 0; i < this->bandwidths.length(); i++) {
        bandwidths[i] = this->bandwidths[i];
    }
}

  void GaussianKDE::pdf(const RealMatrix& data,RealVector& res,
			Teuchos::ETransp trans ) const {
    // init variables
    RealVector x(ndim);

    int num_data = ( trans==Teuchos::NO_TRANS ) ? data.numRows():data.numCols();

    // resize result vector
    res.resize(num_data);
    res.putScalar(0.0);

    // run over all data points
    for (size_t idata = 0; idata < num_data; idata++) {
        // copy samples
      for (size_t idim = 0; idim < ndim; idim++) {
	if ( trans==Teuchos::NO_TRANS )
	  x[idim] = data(idata, idim);
	else
	  x[idim] = data(idim, idata);
      }
      res[idata] = pdf(x);
    }
}

Real GaussianKDE::pdf(const RealVector& x) const {
    // init variables
    Real res = 0.0;
    Real kern = 0, y = 0.0;

    // run over all data points
    for (size_t isample = 0; isample < nsamples; isample++) {
        kern = 1.;
        for (size_t idim = 0; idim < ndim; idim++) {
            // normalize x
            y = (x[idim] - samplesVec[idim][isample]) / bandwidths[idim];
            // evaluate kernel
            kern *= norm[idim] * std::exp(-(y * y) / 2.);
        }
        res += cond[isample] * kern;
    }
    return res / sumCond;
}

void GaussianKDE::cov(RealMatrix& cov) {
    if ((cov.numRows() != ndim) || (cov.numCols() != ndim)) {
        // throw error -> covariance matrix has wrong size
        cout << "covariance matrix has the wrong size" << endl;
        exit(-1);
    }

    // prepare covariance marix
    cov.putScalar(0.0);

    // generate 1d densities and compute means and variances
    vector<Real> means(ndim);
    vector<Real> variances(ndim);

    DensityEstimator kdeMarginalized("gaussian_kde");
    for (size_t idim = 0; idim < ndim; idim++) {
        margToDimX(idim, kdeMarginalized);
        // store moments
        means[idim] = kdeMarginalized.mean();
        variances[idim] = kdeMarginalized.variance();
    }

    // helper variables
    IntVector mdims(2);
    double covij = 0.0;

    DensityEstimator kdeijdim("gaussian_kde");
    for (size_t idim = 0; idim < ndim; idim++) {
        // diagonal is equal to the variance of the marginalized densities
        cov(idim, idim) = variances[idim];
        for (size_t jdim = idim + 1; jdim < ndim; jdim++) {
            // marginalize the density
            mdims[0] = idim;
            mdims[1] = jdim;
            margToDimXs(mdims, kdeijdim);
            // -----------------------------------------------------
            // compute the covariance of Cov(X_i, X_j)
            covij = kdeijdim.mean() - means[idim] * means[jdim];
            cov(idim, jdim) = covij;
            cov(jdim, idim) = covij;
            // -----------------------------------------------------
        }
    }
}

Real GaussianKDE::mean() {
    Real res = 0, kernelMean = 1.;
    for (size_t isample = 0; isample < nsamples; isample++) {
        kernelMean = 1.;
        for (size_t idim = 0; idim < ndim; idim++) {
            kernelMean *= samplesVec[idim][isample];
        }
        res += kernelMean;
    }
    return res / static_cast<Real>(nsamples);
}

Real GaussianKDE::variance() {
    Real meansquared = 0, kernelVariance = 1., x = 0.0, sigma = 0.0;
    for (size_t isample = 0; isample < nsamples; isample++) {
        kernelVariance = 1.;
        for (size_t idim = 0; idim < ndim; idim++) {
            x = samplesVec[idim][isample];
            sigma = bandwidths[idim];
            kernelVariance *= sigma * sigma + x * x;
        }
        meansquared += kernelVariance;
    }
    meansquared /= static_cast<double>(nsamples);

    Real mu = mean();
    Real var = meansquared - mu * mu;

    return var;
}

Real GaussianKDE::std_deviation() {
    return sqrt(variance());
}

void GaussianKDE::computeOptKDEbdwth() {
    if (ndim != bandwidths.length()) {
        std::cerr << "KDEBdwth dimension error" << std::endl;
        exit(-1);
    }

    RealVector flag(ndim);
    flag.putScalar(1.);

    // get min and max in each direction
    RealVector datamin(ndim);
    datamin.putScalar(std::numeric_limits<Real>::max());
    RealVector datamax(ndim);
    datamax.putScalar(std::numeric_limits<Real>::min());

    Real stdd;
    for (size_t idim = 0; idim < ndim; idim++) {
        size_t numBorder = 0;
        // search for maximum in current dimension
        for (size_t isample = 0; isample < nsamples; isample++) {
            if (samplesVec[idim][isample] < datamin[idim]) {
                datamin[idim] = samplesVec[idim][isample];
            }
            if (samplesVec[idim][isample] > datamax[idim]) {
                datamax[idim] = samplesVec[idim][isample];
            }
        }
        Real nearBorder = (datamax[idim] - datamin[idim]) / 20.;

        // count how many values are close to the border
        for (size_t isample = 0; isample < nsamples; isample++) {
            if (samplesVec[idim][isample] - datamin[idim] < nearBorder
                    || datamax[idim] - samplesVec[idim][isample] < nearBorder) {
                numBorder++;
            }
        }
        if (numBorder > static_cast<Real>(nsamples) / 20.) {
            flag[idim] = 0.5;
        }

        // compute the standard deviation
        stdd = getSampleStd(samplesVec[idim]);

        // compute the bandwidth in dimension idim
        bandwidths[idim] = flag[idim]
                * std::pow(4. / (static_cast<Real>(ndim) + 2),
                        1. / (static_cast<Real>(ndim) + 4.)) * stdd
                * std::pow(static_cast<Real>(nsamples),
                        -1. / (static_cast<Real>(ndim) + 4.));
    }
//
//    // compute the bandwidths
//    for (size_t idim = 0; idim < ndim; idim++) {
//        samples1d = samplesVec[idim];
//        stdd = getSampleStd(*samples1d);
//
//        // compute the bandwidth in dimension idim
//        bandwidths[idim] = flag[idim]
//                * std::pow(4. / (static_cast<Real>(ndim) + 2),
//                        1. / (static_cast<Real>(ndim) + 4.)) * stdd
//                * std::pow(static_cast<Real>(nsamples),
//                        -1. / (static_cast<Real>(ndim) + 4.));
//    }

    return;
}

Real GaussianKDE::getSampleMean(RealVector& data) {
    Real res = 0.;
    size_t n = data.length();

    for (size_t i = 0; i < n; i++) {
        res += data[i];
    }
    return res / static_cast<Real>(n);
}

Real GaussianKDE::getSampleVariance(RealVector& data) {
    Real mean = getSampleMean(data);
    Real diff1 = 0.0;
    Real diff2 = 0.0;

    size_t n = data.length();
    for (size_t i = 0; i < n; i++) {
        diff1 += (data[i] - mean) * (data[i] - mean);
        diff2 += (data[i] - mean);
    }
    return 1. / (static_cast<Real>(n) - 1.)
            * (diff1 - 1. / static_cast<Real>(n) * diff2 * diff2);
}

Real GaussianKDE::getSampleStd(RealVector& data) {
    return sqrt(getSampleVariance(data));
}

// ------------------------- additional operations ---------------------------

void GaussianKDE::getConditionalizationFactor(RealVector& pcond) {
    pcond.resize(nsamples);
    for (size_t isample = 0; isample < nsamples; isample++) {
        pcond[isample] = cond[isample];
    }
}

void GaussianKDE::setConditionalizationFactor(const RealVector& pcond) {
    sumCond = 0.0;
    for (size_t isample = 0; isample < nsamples; isample++) {
        cond[isample] = pcond[isample];
        sumCond += cond[isample];
    }
}

void GaussianKDE::marginalize(size_t dim, DensityEstimator& estimator) {
    // dimensionality of new set
    RealVectorArray newSamplesVec(ndim - 1);

    // copy all values but the ones in dim
    for (size_t idim = 0; idim < ndim; idim++) {
        if (idim != dim) {
            newSamplesVec[idim] = samplesVec[idim];
        }
    }

    // initialize kde with new samples
    estimator.initialize(newSamplesVec);
}

void GaussianKDE::margToDimXs(const IntVector& dims,
        DensityEstimator& estimator) {
    // dimensionality of new set
    size_t ndimsNew = dims.length();
    RealVectorArray newSamplesVec(ndimsNew);

    // get the subset of the data for marginalized density
    size_t jdim = 0;
    for (size_t idim = 0; idim < ndimsNew; idim++) {
        newSamplesVec[idim] = samplesVec[dims[idim]];
    }

    // initialize kde with new samples
    estimator.initialize(newSamplesVec);
}

void GaussianKDE::margToDimX(size_t dim, DensityEstimator& estimator) {
    if (dim < 0 || dim >= ndim) {
        PCerr<< "Error: can not marginalize to dim " << dim << "\n"
        << std::endl;
        abort_handler(-1);
    }

    // dimensionality of new set
    RealVectorArray newSamplesVec(1);
    newSamplesVec[0] = samplesVec[dim];// copy

    // initialize marginalized kde
    estimator.initialize(newSamplesVec);
}

/// conditionalization operations
void GaussianKDE::conditionalize(const RealVector& x, const IntVector& dims,
        DensityEstimator& estimator) {
    // compute the dimensions to conditionalize
    IntVector condDims(ndim - dims.length());

    size_t idim = 0;
    size_t jdim = 0;
    size_t kdim = 0;

    while (idim < ndim) {
        jdim = 0;
        while (jdim < dims.length() && idim != dims[jdim]) {
            jdim++;
        }

        // check if current dim has been found in dims
        if (jdim == dims.length()) {
            condDims[kdim] = idim;
            kdim++;
        }

        idim++;
    }

    // update the conditionalization factors
    RealVector pcond(cond.length());
    getConditionalizationFactor(pcond);
    updateConditionalizationFactors(x, condDims, pcond);

    // marginalize the kde to the desired dimensions
    margToDimXs(dims, estimator);
    // set the conditionalization coefficients
    static_cast<GaussianKDE*>(estimator.getEnvelope())->setConditionalizationFactor(
            pcond);
}

void GaussianKDE::condToDimX(const RealVector& x, size_t dim,
        DensityEstimator& estimator) {
    // compute the dimensions to conditionalize over
    IntVector condDims(ndim - 1);
    size_t jdim = 0;
    for (size_t idim = 0; idim < ndim; idim++) {
        if (idim != dim) {
            condDims[jdim] = idim;
            jdim++;
        }
    }

    // update the conditionalization factors
    RealVector pcond(cond.length());
    getConditionalizationFactor(pcond);
    updateConditionalizationFactors(x, condDims, pcond);

    // marginalize the kde to the desired dimensions
    margToDimX(dim, estimator);
    // set the conditionalization coefficients
    static_cast<GaussianKDE*>(estimator.getEnvelope())->setConditionalizationFactor(
            pcond);
}

void GaussianKDE::updateConditionalizationFactors(const RealVector& x,
        const IntVector& dims, RealVector& pcond) {
    // run over all samples and evaluate the kernels in each dimension
    // that should be conditionalized
    size_t idim = 0;
    Real xi = 0.0;
    for (size_t i = 0; i < dims.length(); i++) {
        idim = dims[i];
        if ((idim >= 0) && (idim < ndim)) {
            for (size_t isample = 0; isample < nsamples; isample++) {
                xi = (x[idim] - samplesVec[idim][isample]) / bandwidths[idim];
                pcond[isample] *= norm[idim] * std::exp(-(xi * xi) / 2.);
            }
        } else {
            PCerr<< "Error: can not conditionalize in non existing dimension\n"
            << std::endl;
            abort_handler(-1);
        }
    }
}

const RealVectorArray& GaussianKDE::getSamples() const {
    return samplesVec;
}

}
/* namespace Pecos */
