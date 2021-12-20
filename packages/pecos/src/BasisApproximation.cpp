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
    letter IS the representation, its rep pointer is set to NULL. */
BasisApproximation::
BasisApproximation(BaseConstructor, const SharedBasisApproxData& shared_data):
  sharedDataRep(shared_data.data_rep())
{ /* empty ctor */ }


/** The default constructor: basisApproxRep is NULL in this case. */
BasisApproximation::BasisApproximation()
{ /* empty ctor */ }


/** Envelope constructor only needs to extract enough data to properly
    execute get_basis_approx, since BasisApproximation(BaseConstructor)
    builds the actual base class data for the derived basis functions. */
BasisApproximation::
BasisApproximation(const SharedBasisApproxData& shared_data):
  // Set the rep pointer to the appropriate derived type
  basisApproxRep(get_basis_approx(shared_data))
{
  if ( !basisApproxRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize basisApproxRep to the 
    appropriate derived type. */
std::shared_ptr<BasisApproximation> BasisApproximation::
get_basis_approx(const SharedBasisApproxData& shared_data)
{
  switch (shared_data.data_rep()->basisType) {
  case GLOBAL_NODAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL:
    return std::make_shared<NodalInterpPolyApproximation>(shared_data); break;
  case GLOBAL_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
    return std::make_shared<HierarchInterpPolyApproximation>(shared_data); break;
  case GLOBAL_REGRESSION_ORTHOGONAL_POLYNOMIAL:
  //case PIECEWISE_REGRESSION_ORTHOGONAL_POLYNOMIAL:
    // L1 or L2 regression
    return std::make_shared<RegressOrthogPolyApproximation>(shared_data); break;
  case GLOBAL_PROJECTION_ORTHOGONAL_POLYNOMIAL:
  //case PIECEWISE_PROJECTION_ORTHOGONAL_POLYNOMIAL:
    // projection via numerical integration of inner products
    return std::make_shared<ProjectOrthogPolyApproximation>(shared_data); break;
  case GLOBAL_ORTHOGONAL_POLYNOMIAL: //case PIECEWISE_ORTHOGONAL_POLYNOMIAL:
    // coefficient import -- no coefficient computation required
    return std::make_shared<OrthogPolyApproximation>(shared_data); break;
  //case FOURIER_BASIS:
  //  return std::make_shared<FourierBasisApproximation>(); break;
  //case EIGEN_BASIS:
  //  return std::make_shared<SVDLeftEigenBasisApproximation>(); break;
  default:
    PCerr << "Error: BasisApproximation type "
	  << shared_data.data_rep()->basisType << " not available."<< std::endl;
    return std::shared_ptr<BasisApproximation>(); break;
  }
}


/** Copy constructor manages sharing of basisApproxRep. */
BasisApproximation::BasisApproximation(const BasisApproximation& basis_approx):
  basisApproxRep(basis_approx.basisApproxRep)
{ /* empty ctor */ }


BasisApproximation BasisApproximation::
operator=(const BasisApproximation& basis_approx)
{
  basisApproxRep = basis_approx.basisApproxRep;
  return *this; // calls copy constructor since returned by value
}


BasisApproximation::~BasisApproximation()
{ /* empty dtor */ }


void BasisApproximation::
assign_rep(std::shared_ptr<BasisApproximation> approx_rep)
{
  basisApproxRep = approx_rep;
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
    PCerr << "Error: surrogate_data(SurrogateData&) not available "
	  << "for this basis approximation type." << std::endl;
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


SurrogateData& BasisApproximation::surrogate_data()
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


void BasisApproximation::pop_coefficients(bool save_data)
{
  if (basisApproxRep)
    basisApproxRep->pop_coefficients(save_data); 
  else {
    PCerr << "Error: pop_coefficients() not available for this basis "
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


bool BasisApproximation::advancement_available()
{
  // default is no saturation in refinement candidates
  return (basisApproxRep) ? basisApproxRep->advancement_available() : true;
}


/*
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
*/


void BasisApproximation::clear_inactive()
{
  if (basisApproxRep)
    basisApproxRep->clear_inactive(); 
  //else
  //  default: no stored approx data to clear
}


void BasisApproximation::combine_coefficients()
{
  if (basisApproxRep)
    basisApproxRep->combine_coefficients();
  else {
    PCerr << "Error: combine_coefficients() not available for this basis "
	  << "approximation type." << std::endl;
    abort_handler(-1);
  }
}


void BasisApproximation::combined_to_active(bool clear_combined)
{
  if (basisApproxRep)
    basisApproxRep->combined_to_active(clear_combined);
  else {
    PCerr << "Error: combined_to_active() not available for this basis "
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
