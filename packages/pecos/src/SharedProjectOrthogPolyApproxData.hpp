/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedProjectOrthogPolyApproxData
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        Mike Eldred

#ifndef SHARED_PROJECT_ORTHOG_POLY_APPROX_DATA_HPP
#define SHARED_PROJECT_ORTHOG_POLY_APPROX_DATA_HPP

#include "SharedOrthogPolyApproxData.hpp"

namespace Pecos {

class TensorProductDriver;
class CombinedSparseGridDriver;
class IncrementalSparseGridDriver;
class CubatureDriver;


/// Derived approximation class for multivariate orthogonal polynomial
/// approximation with coefficient estimation via numerical integration.

/** The SharedProjectOrthogPolyApproxData class provides a global approximation
    based on multivariate orthogonal polynomials, where the coefficients are
    computed using numerical integration approaches such as quadrature,
    cubature, sparse grids, and random sampling.  It is used primarily for
    polynomial chaos expansion aproaches to UQ. */

class SharedProjectOrthogPolyApproxData: public SharedOrthogPolyApproxData
{
  //
  //- Heading: Friends
  //

  friend class ProjectOrthogPolyApproximation;

public:

  //
  //- Heading: Constructor and destructor
  //

  /// lightweight constructor
  SharedProjectOrthogPolyApproxData(short basis_type,
				    const UShortArray& approx_order,
				    size_t num_vars);
  /// full constructor
  SharedProjectOrthogPolyApproxData(short basis_type,
				    const UShortArray& approx_order,
				    size_t num_vars,
				    const ExpansionConfigOptions& ec_options,
				    const BasisConfigOptions& bc_options);
  /// destructor
  ~SharedProjectOrthogPolyApproxData();

  //
  //- Heading: Virtual function redefinitions
  //

  void allocate_data();
  void increment_data();
  void decrement_data();

  void pre_push_data();
  void post_push_data();

  void pre_finalize_data();
  void post_finalize_data();

  void pre_combine_data();

  //void construct_basis(const MultivariateDistribution& u_dist);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void increment_component_sobol();

  //
  //- Heading: Member functions
  //

  /// return driverRep
  std::shared_ptr<IntegrationDriver> driver();

private:

  //
  //- Heading: Member functions
  //

  /// initialize multi_index using a sparse grid expansion
  void sparse_grid_multi_index(CombinedSparseGridDriver& csg_driver,
			       UShort2DArray& multi_index);
  /// increment multi_index after a grid refinement
  void increment_sparse_grid_multi_index(
    IncrementalSparseGridDriver& isg_driver, UShort2DArray& multi_index);
  /// decrement multi_index after a grid contraction
  void decrement_sparse_grid_multi_index(
    IncrementalSparseGridDriver& isg_driver, UShort2DArray& multi_index);
  /// return multi_index to a previously incremented state
  void push_sparse_grid_multi_index(IncrementalSparseGridDriver& isg_driver,
				    UShort2DArray& multi_index);

  // initialize tp_multi_index from tpMultiIndexMap
  //void map_tensor_product_multi_index(UShort2DArray& tp_multi_index,
  //				        size_t tp_index);

  /// Perform efficient calculation of tensor-product value via Horner's rule
  Real tensor_product_value(const RealVector& x, const RealVector& tp_coeffs,
			    const UShortArray& approx_order,
			    const UShort2DArray& tp_mi,
			    RealVector& accumulator);

  //
  //- Heading: Data
  //

  /// popped instances of approxOrder that were computed but not selected
  std::map<ActiveKey, UShortArrayDeque> poppedApproxOrder;
};


inline SharedProjectOrthogPolyApproxData::
SharedProjectOrthogPolyApproxData(short basis_type,
				  const UShortArray& approx_order,
				  size_t num_vars):
  SharedOrthogPolyApproxData(basis_type, approx_order, num_vars)
{ }


inline SharedProjectOrthogPolyApproxData::
SharedProjectOrthogPolyApproxData(short basis_type,
				  const UShortArray& approx_order,
				  size_t num_vars,
				  const ExpansionConfigOptions& ec_options,
				  const BasisConfigOptions& bc_options):
  SharedOrthogPolyApproxData(basis_type, approx_order, num_vars,
			     ec_options, bc_options)
{ }


inline SharedProjectOrthogPolyApproxData::~SharedProjectOrthogPolyApproxData()
{ }


/*
// Good for self-contained Pecos, but driverRep grid initialization is
// currently integrated with Dakota::NonDIntegration initialization
inline void SharedProjectOrthogPolyApproxData::
construct_basis(const MultivariateDistribution& u_dist)
{
  SharedOrthogPolyApproxData::construct_basis(u_dist);
  driverRep->initialize_grid(polynomialBasis);
}
*/


inline std::shared_ptr<IntegrationDriver>
SharedProjectOrthogPolyApproxData::driver()
{ return driverRep; }

} // namespace Pecos

#endif
