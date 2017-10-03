/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "BasisApproximation.hpp"
#include "NodalInterpPolyApproximation.hpp"
#include "HierarchInterpPolyApproximation.hpp"
#include "RegressOrthogPolyApproximation.hpp"
#include "ProjectOrthogPolyApproximation.hpp"

static const char rcsId[]="@(#) $Id: BasisApproximation.cpp 4768 2007-12-17 17:49:32Z mseldre $";

namespace Pecos {


/** This constructor is the one which must build the base class data
    for all derived classes.  get_basis_approx() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_basis_approx() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~BasisApproximation). */
BasisApproximation::
BasisApproximation(BaseConstructor, const SharedBasisApproxData& shared_data):
  sharedDataRep(shared_data.data_rep()), basisApproxRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "BasisApproximation::BasisApproximation(BaseConstructor) called to "
        << "build base class for letter." << std::endl;
#endif
}


/** The default constructor: basisApproxRep is NULL in this case.  This
    makes it necessary to check for NULL in the copy constructor,
    assignment operator, and destructor. */
BasisApproximation::BasisApproximation():
  basisApproxRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "BasisApproximation::BasisApproximation() called to build empty "
        << "envelope." << std::endl;
#endif
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_basis_approx, since BasisApproximation(BaseConstructor)
    builds the actual base class data for the derived basis functions. */
BasisApproximation::
BasisApproximation(const SharedBasisApproxData& shared_data):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "BasisApproximation::BasisApproximation(string&) called to "
        << "instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  basisApproxRep = get_basis_approx(shared_data);
  if ( !basisApproxRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize basisApproxRep to the 
    appropriate derived type. */
BasisApproximation* BasisApproximation::
get_basis_approx(const SharedBasisApproxData& shared_data)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_basis_approx(string&)."
        << std::endl;
#endif

  switch (shared_data.data_rep()->basisType) {
  case GLOBAL_NODAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL:
    return new NodalInterpPolyApproximation(shared_data);    break;
  case GLOBAL_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
    return new HierarchInterpPolyApproximation(shared_data); break;
  case GLOBAL_REGRESSION_ORTHOGONAL_POLYNOMIAL:
  //case PIECEWISE_REGRESSION_ORTHOGONAL_POLYNOMIAL:
    // L1 or L2 regression
    return new RegressOrthogPolyApproximation(shared_data);  break;
  case GLOBAL_PROJECTION_ORTHOGONAL_POLYNOMIAL:
  //case PIECEWISE_PROJECTION_ORTHOGONAL_POLYNOMIAL:
    // projection via numerical integration of inner products
    return new ProjectOrthogPolyApproximation(shared_data);  break;
  case GLOBAL_ORTHOGONAL_POLYNOMIAL: //case PIECEWISE_ORTHOGONAL_POLYNOMIAL:
    // coefficient import -- no coefficient computation required
    return new OrthogPolyApproximation(shared_data);         break;
  //case FOURIER_BASIS:
  //  return new FourierBasisApproximation();                break;
  //case EIGEN_BASIS:
  //  return new SVDLeftEigenBasisApproximation();           break;
  default:
    PCerr << "Error: BasisApproximation type "
	  << shared_data.data_rep()->basisType << " not available."<< std::endl;
    return NULL; break;
  }
}


/** Copy constructor manages sharing of basisApproxRep and incrementing
    of referenceCount. */
BasisApproximation::BasisApproximation(const BasisApproximation& basis_approx)
{
  // Increment new (no old to decrement)
  basisApproxRep = basis_approx.basisApproxRep;
  if (basisApproxRep) // Check for an assignment of NULL
    basisApproxRep->referenceCount++;

#ifdef REFCOUNT_DEBUG
  PCout << "BasisApproximation::BasisApproximation(BasisApproximation&)"
	<< std::endl;
  if (basisApproxRep)
    PCout << "basisApproxRep referenceCount = "
	  << basisApproxRep->referenceCount << std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old basisApproxRep,
    assigns new basisApproxRep, and increments referenceCount for new
    basisApproxRep. */
BasisApproximation BasisApproximation::
operator=(const BasisApproximation& basis_approx)
{
  if (basisApproxRep != basis_approx.basisApproxRep) { // std case: old != new
    // Decrement old
    if (basisApproxRep) // Check for null pointer
      if (--basisApproxRep->referenceCount == 0) 
	delete basisApproxRep;
    // Assign and increment new
    basisApproxRep = basis_approx.basisApproxRep;
    if (basisApproxRep) // Check for an assignment of NULL
      basisApproxRep->referenceCount++;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  PCout << "BasisApproximation::operator=(BasisApproximation&)" << std::endl;
  if (basisApproxRep)
    PCout << "basisApproxRep referenceCount = "
	  << basisApproxRep->referenceCount << std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes basisApproxRep
    when referenceCount reaches zero. */
BasisApproximation::~BasisApproximation()
{ 
  // Check for NULL pointer 
  if (basisApproxRep) {
    --basisApproxRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    PCout << "basisApproxRep referenceCount decremented to "
	  << basisApproxRep->referenceCount << std::endl;
#endif
    if (basisApproxRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      PCout << "deleting basisApproxRep" << std::endl;
#endif
      delete basisApproxRep;
    }
  }
}


void BasisApproximation::
assign_rep(BasisApproximation* approx_rep, bool ref_count_incr)
{
  if (basisApproxRep == approx_rep) {
    // if ref_count_incr = true (rep from another envelope), do nothing as
    // referenceCount should already be correct (see also operator= logic).
    // if ref_count_incr = false (rep from on the fly), then this is an error.
    if (!ref_count_incr) {
      PCerr << "Error: duplicated approx_rep pointer assignment without "
	    << "reference count increment in BasisApproximation::assign_rep()."
	    << std::endl;
      abort_handler(-1);
    }
  }
  else { // normal case: old != new
    // Decrement old
    if (basisApproxRep) // Check for NULL
      if ( --basisApproxRep->referenceCount == 0 ) 
	delete basisApproxRep;
    // Assign new
    basisApproxRep = approx_rep;
    // Increment new
    if (basisApproxRep && ref_count_incr)// Check for NULL; honor ref_count_incr
      basisApproxRep->referenceCount++;
  }

#ifdef REFCOUNT_DEBUG
  PCout << "BasisApproximation::assign_rep(BasisApproximation*)" << std::endl;
  if (basisApproxRep)
    PCout << "basisApproxRep referenceCount = "
	  << basisApproxRep->referenceCount << std::endl;
#endif
}


Real BasisApproximation::value(const RealVector& x)
{
  if (!basisApproxRep) {
    PCerr << "Error: value() not available for this basis approximation "
	  << "type." << std::endl;
    abort_handler(-1);
  }

  return basisApproxRep->value(x);
}


const RealVector& BasisApproximation::gradient(const RealVector& x)
{
  if (!basisApproxRep) {
    PCerr << "Error: gradient() not available for this basis approximation "
	  << "type." << std::endl;
    abort_handler(-1);
  }

  return basisApproxRep->gradient(x);
}


const RealSymMatrix& BasisApproximation::hessian(const RealVector& x)
{
  if (!basisApproxRep) {
    PCerr << "Error: hessian() not available for this basis approximation "
	  << "type." << std::endl;
    abort_handler(-1);
  }
    
  return basisApproxRep->hessian(x);
}


void BasisApproximation::surrogate_data(const SurrogateData& data)
{
  if (basisApproxRep)
    basisApproxRep->surrogate_data(data);
  else {
    PCerr << "Error: surrogate_data(SurrogateData&) not available for this "
	  << "basis approximation type." << std::endl;
    abort_handler(-1);
  }
}


const SurrogateData& BasisApproximation::surrogate_data() const
{
  if (!basisApproxRep) {
    PCerr << "Error: surrogate_data() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
    
  return basisApproxRep->surrogate_data();
}


int BasisApproximation::min_coefficients() const
{
  if (!basisApproxRep) { // no default implementation
    PCerr << "Error: min_coefficients() not defined for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }

  return basisApproxRep->min_coefficients(); // fwd to letter
}


void BasisApproximation::compute_coefficients()
{
  if (basisApproxRep)
    basisApproxRep->compute_coefficients(); 
  else {
    PCerr << "Error: compute_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::increment_coefficients()
{
  if (basisApproxRep)
    basisApproxRep->increment_coefficients(); 
  else {
    PCerr << "Error: increment_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::decrement_coefficients()
{
  if (basisApproxRep)
    basisApproxRep->decrement_coefficients(); 
  else {
    PCerr << "Error: decrement_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::push_coefficients()
{
  if (basisApproxRep)
    basisApproxRep->push_coefficients(); 
  else {
    PCerr << "Error: push_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::finalize_coefficients()
{
  if (basisApproxRep)
    basisApproxRep->finalize_coefficients(); 
  else {
    PCerr << "Error: finalize_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::store_coefficients(size_t index)
{
  if (basisApproxRep)
    basisApproxRep->store_coefficients(index); 
  else {
    PCerr << "Error: store_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::restore_coefficients(size_t index)
{
  if (basisApproxRep)
    basisApproxRep->restore_coefficients(index); 
  else {
    PCerr << "Error: restore_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::swap_coefficients(size_t index)
{
  if (basisApproxRep)
    basisApproxRep->swap_coefficients(index);
  else {
    PCerr << "Error: swap_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::remove_stored_coefficients(size_t index)
{
  if (basisApproxRep)
    basisApproxRep->remove_stored_coefficients(index); 
  else {
    PCerr << "Error: remove_stored_coefficients() not available for this "
	  << "basis approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::
combine_coefficients(short combine_type, size_t swap_index)
{
  if (basisApproxRep)
    basisApproxRep->combine_coefficients(combine_type, swap_index);
  else {
    PCerr << "Error: combine_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::print_coefficients(std::ostream& s, bool normalized)
{
  if (basisApproxRep)
    basisApproxRep->print_coefficients(s, normalized);
  else {
    PCerr << "Error: print_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


RealVector BasisApproximation::approximation_coefficients(bool normalized) const
{
  if (!basisApproxRep) {
    PCerr << "Error: approximation_coefficients() not available for this "
	  << "basis approximation type." << std::endl;
    abort_handler(-1);
  }
   
  return basisApproxRep->approximation_coefficients(normalized);// fwd to letter
}


void BasisApproximation::
approximation_coefficients(const RealVector& approx_coeffs, bool normalized)
{
  if (basisApproxRep) // fwd to letter
    basisApproxRep->approximation_coefficients(approx_coeffs, normalized);
  else {
    PCerr << "Error: approximation_coefficients() not available for this "
	  << "basis approximation type." << std::endl;
    abort_handler(-1);
  }
}

void BasisApproximation::
coefficient_labels(std::vector<std::string>& coeff_labels) const
{
  if (basisApproxRep)
    basisApproxRep->coefficient_labels(coeff_labels);
  else {
    PCerr << "Error: coefficient_labels() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}

} // namespace Pecos
