/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedHierarchInterpPolyApproxData
//- Description:  Class for Nodal Interpolation Polynomial Approximation
//-               
//- Owner:        Mike Eldred

#ifndef SHARED_HIERARCH_INTERP_POLY_APPROX_DATA_HPP
#define SHARED_HIERARCH_INTERP_POLY_APPROX_DATA_HPP

#include "SharedInterpPolyApproxData.hpp"
#include "HierarchSparseGridDriver.hpp"

namespace Pecos {


/// Derived approximation class for hierarchical interpolation polynomials
/// (interpolating values and potentially gradients at collocation points).

/** The SharedHierarchInterpPolyApproxData class provides a polynomial
    approximation based on hierarchical interpolation.  Both local and
    global hierarchical basis functions are available.  It is used
    primarily for stochastic collocation approaches to uncertainty
    quantification. */

class SharedHierarchInterpPolyApproxData: public SharedInterpPolyApproxData
{
  //
  //- Heading: Friends
  //

  friend class HierarchInterpPolyApproximation;

public:

  //
  //- Heading: Constructor and destructor
  //

  /// lightweight constructor
  SharedHierarchInterpPolyApproxData(short basis_type, size_t num_vars);
  /// full constructor
  SharedHierarchInterpPolyApproxData(short basis_type, size_t num_vars,
				     const ExpansionConfigOptions& ec_options,
				     const BasisConfigOptions& bc_options);
  /// destructor
  ~SharedHierarchInterpPolyApproxData();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void allocate_component_sobol();
  void increment_component_sobol();

  void set_new_point(const RealVector& x, const UShortArray& basis_index,
		     short order);
  void set_new_point(const RealVector& x, const UShortArray& basis_index,
		     const SizetList& subset_indices, short order);

  size_t barycentric_exact_index(const UShortArray& basis_index);
  size_t barycentric_exact_index(const UShortArray& basis_index,
				 const SizetList& subset_indices);

  unsigned short tensor_product_num_key(size_t i, unsigned short level_i);
  unsigned short tensor_product_max_key(size_t i, unsigned short level_i);
  void precompute_keys(const UShortArray& basis_index);
  void precompute_keys(const UShortArray& basis_index,
		       const SizetList& subset_indices);
  void precompute_max_keys(const UShortArray& basis_index);
  void precompute_max_keys(const UShortArray& basis_index,
			   const SizetList& subset_indices);

private:

  //
  //- Heading: Convenience functions
  //

  /// return driverRep cast to requested derived type
  HierarchSparseGridDriver* hsg_driver();

  //
  //- Heading: Data
  //

  /// Pecos:PIECEWISE_INTERP_POLYNOMIAL or Pecos:PIECEWISE_CUBIC_INTERP
  short polyType;

  /// used for precomputation of the number of hierarchical keys for a
  /// particular basis_index
  UShortArray tpNumKeys;
  /// used for precomputation of the maximum hierarchical key index
  /// for a particular basis_index
  UShortArray tpMaxKeys;
};


inline SharedHierarchInterpPolyApproxData::
SharedHierarchInterpPolyApproxData(short basis_type, size_t num_vars):
  SharedInterpPolyApproxData(basis_type, num_vars)
{ }


inline SharedHierarchInterpPolyApproxData::
SharedHierarchInterpPolyApproxData(short basis_type, size_t num_vars,
				   const ExpansionConfigOptions& ec_options,
				   const BasisConfigOptions& bc_options):
  SharedInterpPolyApproxData(basis_type, num_vars, ec_options, bc_options)
{ }


inline SharedHierarchInterpPolyApproxData::~SharedHierarchInterpPolyApproxData()
{ }


inline void SharedHierarchInterpPolyApproxData::
set_new_point(const RealVector& x, const UShortArray& basis_index, short order)
{
  unsigned short bi_j; UShortArray delta_key;
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  for (size_t j=0; j<numVars; ++j) {
    bi_j = basis_index[j];
    if (bi_j) { // exclusion of pt must be sync'd w/ factors/scalings
      hsg_driver->level_to_delta_key(j, bi_j, delta_key);
      polynomialBasis[bi_j][j].set_new_point(x[j], order, delta_key);
    }
  }
}


inline void SharedHierarchInterpPolyApproxData::
set_new_point(const RealVector& x, const UShortArray& basis_index,
	      const SizetList& subset_indices, short order)
{
  SizetList::const_iterator cit; size_t j; unsigned short bi_j;
  UShortArray delta_key;
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j = *cit; bi_j = basis_index[j];
    if (bi_j) { // exclusion of pt must be sync'd w/ factors/scalings
      hsg_driver->level_to_delta_key(j, bi_j, delta_key);
      polynomialBasis[bi_j][j].set_new_point(x[j], order, delta_key);
    }
  }
}


inline unsigned short SharedHierarchInterpPolyApproxData::
tensor_product_num_key(size_t i, unsigned short level_i)
{
  // for the case of precomputed keys:
  return tpNumKeys[i];

  //HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  //return hsg_driver->level_to_delta_size(i, level_i);
}


inline unsigned short SharedHierarchInterpPolyApproxData::
tensor_product_max_key(size_t i, unsigned short level_i)
{
  // for the case of precomputed keys:
  return tpMaxKeys[i];

  //HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  //return hsg_driver->level_to_delta_pair(i, level_i).second;
}


inline void SharedHierarchInterpPolyApproxData::
precompute_keys(const UShortArray& basis_index)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  if (tpNumKeys.empty()) tpNumKeys.resize(numVars);
  if (tpMaxKeys.empty()) tpMaxKeys.resize(numVars);
  UShortUShortPair key_pr;
  for (size_t i=0; i<numVars; ++i) {
    key_pr = hsg_driver->level_to_delta_pair(i, basis_index[i]);
    tpNumKeys[i] = key_pr.first; tpMaxKeys[i] = key_pr.second;
  }
}


inline void SharedHierarchInterpPolyApproxData::
precompute_keys(const UShortArray& basis_index, const SizetList& subset_indices)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  if (tpNumKeys.empty()) tpNumKeys.resize(numVars);
  if (tpMaxKeys.empty()) tpMaxKeys.resize(numVars);
  SizetList::const_iterator cit; size_t i; UShortUShortPair key_pr;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    i = *cit;
    key_pr = hsg_driver->level_to_delta_pair(i, basis_index[i]);
    tpNumKeys[i] = key_pr.first; tpMaxKeys[i] = key_pr.second;
  }
}


inline void SharedHierarchInterpPolyApproxData::
precompute_max_keys(const UShortArray& basis_index)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  if (tpMaxKeys.empty()) tpMaxKeys.resize(numVars);
  for (size_t i=0; i<numVars; ++i)
    tpMaxKeys[i] = hsg_driver->level_to_delta_pair(i, basis_index[i]).second;
}


inline void SharedHierarchInterpPolyApproxData::
precompute_max_keys(const UShortArray& basis_index,
		    const SizetList& subset_indices)
{
  HierarchSparseGridDriver* hsg_driver = (HierarchSparseGridDriver*)driverRep;
  if (tpMaxKeys.empty()) tpMaxKeys.resize(numVars);
  SizetList::const_iterator cit; size_t i;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    i = *cit;
    tpMaxKeys[i] = hsg_driver->level_to_delta_pair(i, basis_index[i]).second;
  }
}


inline HierarchSparseGridDriver* SharedHierarchInterpPolyApproxData::
hsg_driver()
{ return (HierarchSparseGridDriver*)driverRep; }

} // namespace Pecos

#endif
