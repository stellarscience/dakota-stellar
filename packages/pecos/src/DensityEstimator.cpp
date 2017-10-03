#include "DensityEstimator.hpp"
#include "GaussianKDE.hpp"
//#include "NatafDensity.hpp" // Deactivate until Fabian resolves errors with
// Nataf transformation

// #define REFCOUNT_DEBUG

namespace Pecos {

// -------------------- constructors and desctructors --------------------

/** The default constructor: samples is NULL in this case.  This
 makes it necessary to check for NULL in the copy constructor,
 assignment operator, and destructor. */
DensityEstimator::DensityEstimator() :
        densityEstimator(NULL), referenceCount(1) {
#ifdef REFCOUNT_DEBUG
    PCout << "DensityEstimator::DensityEstimator(BaseConstructor) called to "
    << "build empty envelope." << std::endl;
#endif
}

/** Envelope constructor only needs to extract enough data to properly
 execute get_density_estimator, since DensityEstimator(BaseConstructor)
 builds the actual base class data for the derived transformations. */
DensityEstimator::DensityEstimator(const String& density_estimator) :
        densityEstimator(NULL), referenceCount(1) {
#ifdef REFCOUNT_DEBUG
    PCout << "DensityEstimator::DensityEstimator(string&) called to "
    << "instantiate envelope." << std::endl;
#endif

    // Set the rep pointer to the appropriate derived type
    density_estimator_type = density_estimator;
    densityEstimator = get_density_estimator(density_estimator);
    if (!densityEstimator) // bad type or insufficient memory
        abort_handler(-1);
}

/** Copy constructor manages sharing of densityEstimatorand incrementing
 of referenceCount. */
DensityEstimator::DensityEstimator(const DensityEstimator& density_estimator) {
    // Increment new (no old to decrement)
    densityEstimator = density_estimator.densityEstimator;

    if (densityEstimator) // Check for an assignment of NULL
        densityEstimator->referenceCount++;

#ifdef REFCOUNT_DEBUG
    PCout << "DensityEstimator::DensityEstimator(DensityEstimator&)"
    << std::endl;
    if (densityEstimator)
    PCout << "densityEstimator referenceCount = " << densityEstimator->referenceCount
    << std::endl;
#endif
}

/** Assignment operator decrements referenceCount for old dataTransRep,
 assigns new dataTransRep, and increments referenceCount for new
 dataTransRep. */
DensityEstimator DensityEstimator::operator=(
        const DensityEstimator& density_estimator) {
    if (densityEstimator != density_estimator.densityEstimator) { // normal case: old != new
        // Decrement old
        if (densityEstimator) {            // Check for null pointer
            if (--densityEstimator->referenceCount == 0) {
                delete densityEstimator;
            }
        }
        // Assign and increment new
        densityEstimator = density_estimator.densityEstimator;
        if (densityEstimator) // Check for an assignment of NULL
            densityEstimator->referenceCount++;
    }
    // else if assigning same rep, then do nothing since referenceCount
    // should already be correct

#ifdef REFCOUNT_DEBUG
    PCout << "DensityEstimator::operator=(DensityEstimator&)" << std::endl;
    if (densityEstimator)
    PCout << "densityEstimator referenceCount = " << densityEstimator->referenceCount
    << std::endl;
#endif

    return *this; // calls copy constructor since returned by value
}

/** Destructor decrements referenceCount and only deletes dataTransRep
 when referenceCount reaches zero. */
DensityEstimator::~DensityEstimator() {
    // Check for NULL pointer
    if (densityEstimator) {
        --densityEstimator->referenceCount;
#ifdef REFCOUNT_DEBUG
        PCout << "densityEstimator referenceCount decremented to "
        << densityEstimator->referenceCount << std::endl;
#endif
        if (densityEstimator->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
            PCout << "deleting densityEstimator" << std::endl;
#endif
            delete densityEstimator;
        }
    }
}

// ----------------------------------------------------------------------

/** Used only by the envelope constructor to initialize densityEstimator to the
 appropriate derived type. */
DensityEstimator* DensityEstimator::get_density_estimator(
        const String& density_estimator_type) {
#ifdef REFCOUNT_DEBUG
    PCout << "Envelope instantiating letter in get_density_estimator(string&)."
    << std::endl;
#endif

    if (density_estimator_type == "gaussian_kde") {
        return new GaussianKDE();
	// Deactivate until Fabian resolves errors with
	// Nataf transformation 
	//} else if (density_estimator_type == "nataf") {
        //return new NatafDensity();
    } else {
        PCerr << "Error: DensityEstimator type '" << density_estimator_type
                << "' not available." << std::endl;
        return NULL;
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
        return densityEstimator;
    } else {
        // if envelope is letter -> return this
        return this;
    }
}

bool DensityEstimator::is_null() {
    return densityEstimator == NULL;
}

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
