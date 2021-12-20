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
    letter IS the representation, its rep pointer is set to NULL. */
SharedBasisApproxData::
SharedBasisApproxData(BaseConstructor, short basis_type, size_t num_vars):
  basisType(basis_type), numVars(num_vars)
{ /* empty ctor */ }


/** The default constructor: dataRep is NULL in this case. */
SharedBasisApproxData::SharedBasisApproxData()
{ /* empty ctor */ }


/** Envelope constructor only needs to extract enough data to properly
    execute get_shared_data, since SharedBasisApproxData(BaseConstructor)
    builds the actual base class data for the derived basis functions. */
SharedBasisApproxData::
SharedBasisApproxData(short basis_type, const UShortArray& approx_order,
		      size_t num_vars):
  // Set the rep pointer to the appropriate derived type
  dataRep(get_shared_data(basis_type, approx_order, num_vars))
{
  if ( !dataRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize dataRep to the 
    appropriate derived type. */
std::shared_ptr<SharedBasisApproxData> SharedBasisApproxData::
get_shared_data(short basis_type, const UShortArray& approx_order,
		size_t num_vars)
{
  switch (basis_type) {
  case GLOBAL_NODAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL:
    return std::make_shared<SharedNodalInterpPolyApproxData>
      (basis_type, num_vars);
    break;
  case GLOBAL_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
    return std::make_shared<SharedHierarchInterpPolyApproxData>
      (basis_type, num_vars);
    break;
  case GLOBAL_REGRESSION_ORTHOGONAL_POLYNOMIAL:
  //case PIECEWISE_REGRESSION_ORTHOGONAL_POLYNOMIAL:
    // L1 or L2 regression
    return std::make_shared<SharedRegressOrthogPolyApproxData>
      (basis_type, approx_order, num_vars);
    break;
  case GLOBAL_PROJECTION_ORTHOGONAL_POLYNOMIAL:
  //case PIECEWISE_PROJECTION_ORTHOGONAL_POLYNOMIAL:
    // projection via numerical integration of inner products
    return std::make_shared<SharedProjectOrthogPolyApproxData>
      (basis_type, approx_order, num_vars);
    break;
  case GLOBAL_ORTHOGONAL_POLYNOMIAL: //case PIECEWISE_ORTHOGONAL_POLYNOMIAL:
    // coefficient import -- no coefficient computation required
    return std::make_shared<SharedOrthogPolyApproxData>
      (basis_type, approx_order, num_vars);
    break;
  //case FOURIER_BASIS:
  //  return std::make_shared<SharedFourierBasisApproxData>(num_vars); break;
  //case EIGEN_BASIS:
  //  return std::make_shared<SharedEigenBasisApproxData>(num_vars);   break;
  default:
    PCerr << "Error: SharedBasisApproxData type " << basis_type
	  << " not available." << std::endl;
    return std::shared_ptr<SharedBasisApproxData>(); break;
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
  // Set the rep pointer to the appropriate derived type
  dataRep(get_shared_data(basis_type, approx_order, num_vars,
			  ec_options, bc_options, rc_options))
{
  if ( !dataRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize dataRep to the 
    appropriate derived type. */
std::shared_ptr<SharedBasisApproxData> SharedBasisApproxData::
get_shared_data(short basis_type, const UShortArray& approx_order,
		size_t num_vars, const ExpansionConfigOptions& ec_options,
		const BasisConfigOptions& bc_options,
		const RegressionConfigOptions& rc_options)
{
  switch (basis_type) {
  case GLOBAL_NODAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL:
    return std::make_shared<SharedNodalInterpPolyApproxData>
      (basis_type, num_vars, ec_options, bc_options);
    break;
  case GLOBAL_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
  case PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL:
    return std::make_shared<SharedHierarchInterpPolyApproxData>
      (basis_type, num_vars, ec_options, bc_options);
    break;
  case GLOBAL_REGRESSION_ORTHOGONAL_POLYNOMIAL:
  //case PIECEWISE_REGRESSION_ORTHOGONAL_POLYNOMIAL:
    // L1 or L2 regression
    return std::make_shared<SharedRegressOrthogPolyApproxData>
      (basis_type, approx_order, num_vars, ec_options, bc_options, rc_options);
    break;
  case GLOBAL_PROJECTION_ORTHOGONAL_POLYNOMIAL:
  //case PIECEWISE_PROJECTION_ORTHOGONAL_POLYNOMIAL:
    // projection via numerical integration of inner products
    return std::make_shared<SharedProjectOrthogPolyApproxData>
      (basis_type, approx_order, num_vars, ec_options, bc_options);
    break;
  case GLOBAL_ORTHOGONAL_POLYNOMIAL: //case PIECEWISE_ORTHOGONAL_POLYNOMIAL:
    // coefficient import -- no coefficient computation required
    return std::make_shared<SharedOrthogPolyApproxData>
      (basis_type, approx_order, num_vars, ec_options, bc_options);
    break;
  //case FOURIER_BASIS:
  //  return std::make_shared<SharedFourierBasisApproxData>(num_vars); break;
  //case EIGEN_BASIS:
  //  return std::make_shared<SharedEigenBasisApproxData>(num_vars);   break;
  default:
    PCerr << "Error: SharedBasisApproxData type " << basis_type
	  << " not available." << std::endl;
    return std::shared_ptr<SharedBasisApproxData>(); break;
  }
}


/** Copy constructor manages sharing of dataRep. */
SharedBasisApproxData::
SharedBasisApproxData(const SharedBasisApproxData& shared_data):
  dataRep(shared_data.dataRep)
{ /* empty ctor */ }


SharedBasisApproxData SharedBasisApproxData::
operator=(const SharedBasisApproxData& shared_data)
{
  dataRep = shared_data.dataRep;
  return *this; // calls copy constructor since returned by value
}


SharedBasisApproxData::~SharedBasisApproxData()
{ /* empty dtor */ }


void SharedBasisApproxData::
assign_rep(std::shared_ptr<SharedBasisApproxData> data_rep)
{
  dataRep = data_rep;
}

} // namespace Pecos
