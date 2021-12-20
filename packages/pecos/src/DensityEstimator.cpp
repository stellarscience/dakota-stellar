/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "DensityEstimator.hpp"
#include "GaussianKDE.hpp"
//#include "NatafDensity.hpp" // Deactivate until Fabian resolves errors with
// Nataf transformation

// #define REFCOUNT_DEBUG

namespace Pecos {

// -------------------- constructors and desctructors --------------------

/** The default constructor: samples is NULL in this case. */
DensityEstimator::DensityEstimator()
{ /* empty ctor */ }

/** Envelope constructor only needs to extract enough data to properly
 execute get_density_estimator, since DensityEstimator(BaseConstructor)
 builds the actual base class data for the derived transformations. */
DensityEstimator::DensityEstimator(const String& density_estimator):
  // Set the rep pointer to the appropriate derived type
  density_estimator_type(density_estimator),
  densityEstimator(get_density_estimator(density_estimator))
{
  if (!densityEstimator) // bad type or insufficient memory
    abort_handler(-1);
}

/** Copy constructor manages sharing of densityEstimator. */
DensityEstimator::DensityEstimator(const DensityEstimator& density_estimator):
  densityEstimator(density_estimator.densityEstimator)
{ /* empty ctor */ }

/** Assignment operator decrements referenceCount for old dataTransRep,
 assigns new dataTransRep, and increments referenceCount for new
 dataTransRep. */
DensityEstimator 
DensityEstimator::operator=(const DensityEstimator& density_estimator)
{
  densityEstimator = density_estimator.densityEstimator;
  return *this; // calls copy constructor since returned by value
}


DensityEstimator::~DensityEstimator()
{ /* empty dtor */}


// ----------------------------------------------------------------------

/** Used only by the envelope constructor to initialize densityEstimator to the
 appropriate derived type. */
std::shared_ptr<DensityEstimator>
DensityEstimator::get_density_estimator(const String& density_estimator_type)
{
  if (density_estimator_type == "gaussian_kde") {
    return std::make_shared<GaussianKDE>();
    // Deactivate until Fabian resolves errors with
    // Nataf transformation 
    //} else if (density_estimator_type == "nataf") {
    //return new NatafDensity();
  } else {
    PCerr << "Error: DensityEstimator type '" << density_estimator_type
	  << "' not available." << std::endl;
    return std::make_shared<DensityEstimator>();
  }
}

String DensityEstimator::getType() {
    if (densityEstimator) { // envelope fwd to letter
        return densityEstimator->getType();
    } else { // empty envelope
        // if envelope is letter
        return density_estimator_type;
    }
}

DensityEstimator* DensityEstimator::getEnvelope() {
  if (densityEstimator) {
    // return envelope
    return densityEstimator.get();
  } else {
    // if envelope is letter -> return this
    return this;
  }
}

//bool DensityEstimator::is_null() {
//    return densityEstimator == NULL;
//}

  void DensityEstimator::initialize(RealMatrix& samples,Teuchos::ETransp trans){
    if (densityEstimator) { // envelope fwd to letter
      densityEstimator->initialize(samples, trans);
    } else { // letter lacking redefinition of virtual fn
        PCerr
                << "Error: derived class does not redefine initialize(RealMatrix& samples) virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

void DensityEstimator::initialize(RealVectorArray& samples) {
    if (densityEstimator) { // envelope fwd to letter
        densityEstimator->initialize(samples);
    } else { // letter lacking redefinition of virtual fn
        PCerr
                << "Error: derived class does not redefine initialize(RealVectorArray& samples) virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

// ----------------------------------------------------------------------

size_t DensityEstimator::getDim() {
    if (densityEstimator) { // envelope fwd to letter
        return densityEstimator->getDim();
    } else { // letter lacking redefinition of virtual fn
        PCerr << "Error: derived class does not redefine getDim() virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

Real DensityEstimator::mean() {
    if (densityEstimator) { // envelope fwd to letter
        return densityEstimator->mean();
    } else { // letter lacking redefinition of virtual fn
        PCerr << "Error: derived class does not redefine mean() virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

Real DensityEstimator::variance() {
    if (densityEstimator) { // envelope fwd to letter
        return densityEstimator->variance();
    } else { // letter lacking redefinition of virtual fn
        PCerr
                << "Error: derived class does not redefine variance() virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

Real DensityEstimator::std_deviation() {
    if (densityEstimator) { // envelope fwd to letter
        return densityEstimator->std_deviation();
    } else { // letter lacking redefinition of virtual fn
        PCerr
                << "Error: derived class does not redefine std_deviation() virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

/// computes the covariance matrix
void DensityEstimator::cov(RealMatrix& cov) {
    if (densityEstimator) { // envelope fwd to letter
        return densityEstimator->cov(cov);
    } else { // letter lacking redefinition of virtual fn
        PCerr << "Error: derived class does not redefine cov() virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}
/// computes the correlation matrix
void DensityEstimator::corrcoeff(RealMatrix& corr) {
    if (densityEstimator) { // envelope fwd to letter
        // get covariance matrix from the letter
        densityEstimator->cov(corr);
    } else {
        // envelope is letter
        cov(corr);
    }

    // normalize it
    double corrij = 0.0;
    double sigmai = 0.0, sigmaj = 0.0;
    size_t ndim = corr.numCols();

    for (size_t idim = 0; idim < ndim; idim++) {
        // normalize row idim but the diagonal element
        sigmai = sqrt(corr(idim, idim));
        for (size_t jdim = idim + 1; jdim < ndim; jdim++) {
            sigmaj = sqrt(corr(jdim, jdim));
            corrij = corr(idim, jdim) / (sigmai * sigmaj);
            corr(idim, jdim) = corrij;
            corr(jdim, idim) = corrij;
        }
        // set the diagonal element
        corr(idim, idim) = 1.0;
    }
}

/// operations for one sample
Real DensityEstimator::pdf(const RealVector& x) const{
    if (densityEstimator) { // envelope fwd to letter
        return densityEstimator->pdf(x);
    } else { // letter lacking redefinition of virtual fn
        PCerr << "Error: derived class does not redefine pdf() virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

/// operations for a set of samples
void DensityEstimator::pdf(const RealMatrix& data,RealVector& res,
			   Teuchos::ETransp trans) const{
    if (densityEstimator) { // envelope fwd to letter
      densityEstimator->pdf(data, res, trans);
    } else { // letter lacking redefinition of virtual fn
        PCerr << "Error: derived class does not redefine pdf() virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

/// marginalization operations
void DensityEstimator::marginalize(size_t dim, DensityEstimator& estimator) {
    if (densityEstimator) { // envelope fwd to letter
        densityEstimator->marginalize(dim, estimator);
    } else { // letter lacking redefinition of virtual fn
        PCerr
                << "Error: derived class does not redefine marginalize() virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

void DensityEstimator::margToDimXs(const IntVector& dims,
        DensityEstimator& estimator) {
    if (densityEstimator) { // envelope fwd to letter
        densityEstimator->margToDimXs(dims, estimator);
    } else { // letter lacking redefinition of virtual fn
        PCerr
                << "Error: derived class does not redefine marginalize() virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

void DensityEstimator::margToDimX(size_t dim, DensityEstimator& estimator) {
    if (densityEstimator) { // envelope fwd to letter
        return densityEstimator->margToDimX(dim, estimator);
    } else { // letter lacking redefinition of virtual fn
        PCerr
                << "Error: derived class does not redefine margToDimX() virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

/// conditionalization operations
void DensityEstimator::conditionalize(const RealVector& x,
        const IntVector& dims, DensityEstimator& estimator) {
    if (densityEstimator) { // envelope fwd to letter
        return densityEstimator->conditionalize(x, dims, estimator);
    } else { // letter lacking redefinition of virtual fn
        PCerr
                << "Error: derived class does not redefine conditionalize() virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

void DensityEstimator::condToDimX(const RealVector& x, size_t dim,
        DensityEstimator& estimator) {
    if (densityEstimator) { // envelope fwd to letter
        densityEstimator->condToDimX(x, dim, estimator);
    } else { // letter lacking redefinition of virtual fn
        PCerr
                << "Error: derived class does not redefine margToDimX() virtual fn.\n"
                << "       No default defined at DensityEstimator base class.\n"
                << std::endl;
        abort_handler(-1);
    }
}

}
/* namespace Pecos */
