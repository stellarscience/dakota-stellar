/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "SharedBasisApproxData.hpp"
#include "SharedNodalInterpPolyApproxData.hpp"
#include "SharedHierarchInterpPolyApproxData.hpp"
#include "SharedRegressOrthogPolyApproxData.hpp"
#include "SharedProjectOrthogPolyApproxData.hpp"

static const char rcsId[]="@(#) $Id: SharedBasisApproxData.cpp 4768 2007-12-17 17:49:32Z mseldre $";

namespace Pecos {


/** This constructor is the one which must build the base class data
    for all derived classes.  get_shared_data() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_shared_data() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~SharedBasisApproxData). */
SharedBasisApproxData::
SharedBasisApproxData(BaseConstructor, short basis_type, size_t num_vars):
  basisType(basis_type), numVars(num_vars), dataRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "SharedBasisApproxData::SharedBasisApproxData(BaseConstructor) "
        << "called to build base class for letter." << std::endl;
#endif
}


/** The default constructor: dataRep is NULL in this case.  This
    makes it necessary to check for NULL in the copy constructor,
    assignment operator, and destructor. */
SharedBasisApproxData::SharedBasisApproxData(): dataRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "SharedBasisApproxData::SharedBasisApproxData() called to build "
        << "empty envelope." << std::endl;
#endif
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_shared_data, since SharedBasisApproxData(BaseConstructor)
    builds the actual base class data for the derived basis functions. */
SharedBasisApproxData::
SharedBasisApproxData(short basis_type, const UShortArray& approx_order,
		      size_t num_vars):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "SharedBasisApproxData::SharedBasisApproxData(short) called to "
        << "instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  dataRep = get_shared_data(basis_type, approx_order, num_vars);
  if ( !dataRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize dataRep to the 
    appropriate derived type. */
SharedBasisApproxData* SharedBasisApproxData::
get_shared_data(short basis_type, const UShortArray& approx_order,
		size_t num_vars)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_shared_data()." << std::endl;
#endif

  switch (basis_type) {
  case GLOBAL_NODAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL:
    return new SharedNodalInterpPolyApproxData(basis_type, num_vars);
    break;
  case GLOBAL_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
    return new SharedHierarchInterpPolyApproxData(basis_type, num_vars);
    break;
  case GLOBAL_REGRESSION_ORTHOGONAL_POLYNOMIAL:
  //case PIECEWISE_REGRESSION_ORTHOGONAL_POLYNOMIAL:
    // L1 or L2 regression
    return new
      SharedRegressOrthogPolyApproxData(basis_type, approx_order, num_vars);
    break;
  case GLOBAL_PROJECTION_ORTHOGONAL_POLYNOMIAL:
  //case PIECEWISE_PROJECTION_ORTHOGONAL_POLYNOMIAL:
    // projection via numerical integration of inner products
    return new
      SharedProjectOrthogPolyApproxData(basis_type, approx_order, num_vars);
    break;
  case GLOBAL_ORTHOGONAL_POLYNOMIAL: //case PIECEWISE_ORTHOGONAL_POLYNOMIAL:
    // coefficient import -- no coefficient computation required
    return new SharedOrthogPolyApproxData(basis_type, approx_order, num_vars);
    break;
  //case FOURIER_BASIS:
  //  return new SharedFourierBasisApproxData(num_vars); break;
  //case EIGEN_BASIS:
  //  return new SharedEigenBasisApproxData(num_vars);   break;
  default:
    PCerr << "Error: SharedBasisApproxData type " << basis_type
	  << " not available." << std::endl;
    return NULL; break;
  }
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_shared_data, since SharedBasisApproxData(BaseConstructor)
    builds the actual base class data for the derived basis functions. */
SharedBasisApproxData::
SharedBasisApproxData(short basis_type, const UShortArray& approx_order,
		      size_t num_vars, const ExpansionConfigOptions& ec_options,
		      const BasisConfigOptions& bc_options,
		      const RegressionConfigOptions& rc_options):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "SharedBasisApproxData::SharedBasisApproxData(short) called to "
        << "instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  dataRep = get_shared_data(basis_type, approx_order, num_vars,
			    ec_options, bc_options, rc_options);
  if ( !dataRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize dataRep to the 
    appropriate derived type. */
SharedBasisApproxData* SharedBasisApproxData::
get_shared_data(short basis_type, const UShortArray& approx_order,
		size_t num_vars, const ExpansionConfigOptions& ec_options,
		const BasisConfigOptions& bc_options,
		const RegressionConfigOptions& rc_options)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_shared_data()." << std::endl;
#endif

  switch (basis_type) {
  case GLOBAL_NODAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL:
    return new SharedNodalInterpPolyApproxData(basis_type, num_vars,
					       ec_options, bc_options);
    break;
  case GLOBAL_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
    return new SharedHierarchInterpPolyApproxData(basis_type, num_vars,
						  ec_options, bc_options);
    break;
  case GLOBAL_REGRESSION_ORTHOGONAL_POLYNOMIAL:
  //case PIECEWISE_REGRESSION_ORTHOGONAL_POLYNOMIAL:
    // L1 or L2 regression
    return new SharedRegressOrthogPolyApproxData(basis_type, approx_order,
						 num_vars, ec_options,
						 bc_options, rc_options);
    break;
  case GLOBAL_PROJECTION_ORTHOGONAL_POLYNOMIAL:
  //case PIECEWISE_PROJECTION_ORTHOGONAL_POLYNOMIAL:
    // projection via numerical integration of inner products
    return new SharedProjectOrthogPolyApproxData(basis_type, approx_order,
						 num_vars, ec_options,
						 bc_options);
    break;
  case GLOBAL_ORTHOGONAL_POLYNOMIAL: //case PIECEWISE_ORTHOGONAL_POLYNOMIAL:
    // coefficient import -- no coefficient computation required
    return new SharedOrthogPolyApproxData(basis_type, approx_order, num_vars,
					  ec_options, bc_options);
    break;
  //case FOURIER_BASIS:
  //  return new SharedFourierBasisApproxData(num_vars); break;
  //case EIGEN_BASIS:
  //  return new SharedEigenBasisApproxData(num_vars);   break;
  default:
    PCerr << "Error: SharedBasisApproxData type " << basis_type
	  << " not available." << std::endl;
    return NULL; break;
  }
}


/** Copy constructor manages sharing of dataRep and incrementing
    of referenceCount. */
SharedBasisApproxData::
SharedBasisApproxData(const SharedBasisApproxData& shared_data)
{
  // Increment new (no old to decrement)
  dataRep = shared_data.dataRep;
  if (dataRep) // Check for an assignment of NULL
    dataRep->referenceCount++;

#ifdef REFCOUNT_DEBUG
  PCout << "SharedBasisApproxData::SharedBasisApproxData"
	<< "(SharedBasisApproxData&)" << std::endl;
  if (dataRep)
    PCout << "dataRep referenceCount = " << dataRep->referenceCount<< std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old dataRep,
    assigns new dataRep, and increments referenceCount for new
    dataRep. */
SharedBasisApproxData SharedBasisApproxData::
operator=(const SharedBasisApproxData& shared_data)
{
  if (dataRep != shared_data.dataRep) { // std case: old != new
    // Decrement old
    if (dataRep) // Check for null pointer
      if (--dataRep->referenceCount == 0) 
	delete dataRep;
    // Assign and increment new
    dataRep = shared_data.dataRep;
    if (dataRep) // Check for an assignment of NULL
      dataRep->referenceCount++;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  PCout << "SharedBasisApproxData::operator=(SharedBasisApproxData&)"
	<< std::endl;
  if (dataRep)
    PCout << "dataRep referenceCount = " << dataRep->referenceCount<< std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes dataRep
    when referenceCount reaches zero. */
SharedBasisApproxData::~SharedBasisApproxData()
{
  // Check for NULL pointer 
  if (dataRep) {
    --dataRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    PCout << "dataRep referenceCount decremented to "
	  << dataRep->referenceCount << std::endl;
#endif
    if (dataRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      PCout << "deleting dataRep" << std::endl;
#endif
      delete dataRep;
    }
  }
}


void SharedBasisApproxData::
assign_rep(SharedBasisApproxData* data_rep, bool ref_count_incr)
{
  if (dataRep == data_rep) {
    // if ref_count_incr = true (rep from another envelope), do nothing as
    // referenceCount should already be correct (see also operator= logic).
    // if ref_count_incr = false (rep from on the fly), then this is an error.
    if (!ref_count_incr) {
      PCerr << "Error: duplicated data_rep pointer assignment without "
	    << "reference count increment in SharedBasisApproxData::"
	    << "assign_rep()." << std::endl;
      abort_handler(-1);
    }
  }
  else { // normal case: old != new
    // Decrement old
    if (dataRep) // Check for NULL
      if ( --dataRep->referenceCount == 0 ) 
	delete dataRep;
    // Assign new
    dataRep = data_rep;
    // Increment new
    if (dataRep && ref_count_incr) // Check for NULL & honor ref_count_incr
      dataRep->referenceCount++;
  }

#ifdef REFCOUNT_DEBUG
  PCout << "SharedBasisApproxData::assign_rep(BasisApproximation*)" <<std::endl;
  if (dataRep)
    PCout << "dataRep referenceCount = " << dataRep->referenceCount <<std::endl;
#endif
}

} // namespace Pecos
