/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedNodalInterpPolyApproxData
//- Description:  Class for Nodal Interpolation Polynomial Approximation
//-               
//- Owner:        Mike Eldred

#ifndef SHARED_NODAL_INTERP_POLY_APPROX_DATA_HPP
#define SHARED_NODAL_INTERP_POLY_APPROX_DATA_HPP

#include "SharedInterpPolyApproxData.hpp"
#include "InterpolationPolynomial.hpp"

namespace Pecos {

class TensorProductDriver;
class CombinedSparseGridDriver;


/// Derived approximation class for nodal interpolation polynomials
/// (global approximation interpolating function values and
/// potentially gradients at collocation points).

/** The SharedNodalInterpPolyApproxData class provides a global polynomial
    approximation based on either Lagrange or Hermite interpolation
    polynomials using a nodal basis approach.  It is used primarily
    for stochastic collocation approaches to uncertainty quantification. */

class SharedNodalInterpPolyApproxData: public SharedInterpPolyApproxData
{
  //
  //- Heading: Friends
  //

  friend class NodalInterpPolyApproximation;

public:

  //
  //- Heading: Constructor and destructor
  //

  /// lightweight constructor
  SharedNodalInterpPolyApproxData(short basis_type, size_t num_vars);
  /// full constructor
  SharedNodalInterpPolyApproxData(short basis_type, size_t num_vars,
				  const ExpansionConfigOptions& ec_options,
				  const BasisConfigOptions& bc_options);
  /// destructor
  ~SharedNodalInterpPolyApproxData();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void allocate_data();
  void allocate_component_sobol();
  void increment_component_sobol();

  void set_new_point(const RealVector& x, const UShortArray& basis_index,
		     short order);
  void set_new_point(const RealVector& x, const UShortArray& basis_index,
		     const SizetList& subset_indices, short order);

  size_t barycentric_exact_index(const UShortArray& basis_index);
  size_t barycentric_exact_index(const UShortArray& basis_index,
				 const SizetList& subset_indices);

  //
  //- Heading: Member functions
  //

  /// return driverRep cast to requested derived type
  TensorProductDriver*      tpq_driver();
  /// return driverRep cast to requested derived type
  CombinedSparseGridDriver* csg_driver();

private:

  //
  //- Heading: Convenience functions
  //

  /// update accumulators for barycentric type1 contributions to moment value
  void accumulate_barycentric(RealVector& t1_accumulator,
			      const UShortArray& lev_index,
			      const UShortArray& key_p);
  /// update accumulators for type1 contributions to moment value
  void accumulate_horners(RealVector& t1_accumulator,
			  const UShortArray& lev_index,
			  const UShortArray& key_p, const RealVector& x);
  /// update accumulators for type1 + type2 contributions to moment value
  void accumulate_horners(RealVector& t1_accumulator,
			  RealMatrix& t2_accumulator,
			  const UShortArray& lev_index,
			  const UShortArray& key_p, const RealVector& x);

  /// update accumulators for barycentric type1 contributions to moment gradient
  void accumulate_barycentric_gradient(RealMatrix& t1_accumulator,
				       const UShortArray& lev_index,
				       const UShortArray& key_p,
				       const SizetArray& dvv);
  /// update accumulators for type1 contributions to moment gradient
  void accumulate_horners_gradient(RealMatrix& t1_accumulator,
				   const UShortArray& lev_index,
				   const UShortArray& key_p,
				   const SizetArray& dvv, const RealVector& x);
  /// update accumulators for type1 + type2 contributions to moment gradient
  void accumulate_horners_gradient(RealMatrix& t1_accumulator,
				   RealMatrixArray& t2_accumulators,
				   const UShortArray& lev_index,
				   const UShortArray& key_p,
				   const SizetArray& dvv, const RealVector& x);

  /// update precomputation of nonzero multidimensional integrals of
  /// products of interpolation polynomials
  void update_nonzero_basis_products(const UShort2DArray& sm_multi_index);

  /// evaluate 1D integral of product of interpolation polynomials
  bool basis_product_1d(InterpolationPolynomial* poly_rep_1,
			InterpolationPolynomial* poly_rep_2,
			unsigned short key_1, unsigned short key_2,
			const RealArray& pts, const RealArray& wts, Real& prod);
  /// lookup multidimensional integral of products of interpolation polynomials
  bool basis_product(const UShortArray& lev_index_1, const UShortArray& key_1,
		     const UShortArray& lev_index_2, const UShortArray& key_2,
		     Real& prod);

  /// computes higher-order grid for tensor reinterpolation of the
  /// covariance fn for non-integrated dimensions in all_variables mode
  void reinterpolated_level(const UShortArray& lev_index);

  /// create a unique map key for value() and gradient() calculation reuse
  void update_member_key(const UShortArray& data,
			 const SizetList& member_indices,
			 UShortArray& member_map_key, size_t cntr);

  //
  //- Heading: Data
  //

  /// type of interpolation for all-variables covariance and variance gradient
  short momentInterpType;

  /// special tensor/sparse integration driver for (exactly) computing
  /// expansion moments using sufficiently high-order Gaussian quadrature
  /// rules on the interpolant
  IntegrationDriver expMomentIntDriver;

  /// map from random index to unique nonZerosMapArray
  SizetArray nonZerosMapIndices;
  /// tracks level maxima already populated within nonZerosMap
  UShortArray nonZerosMapMaxLevels;
  /// expectations of products of interpolation polynomials,
  /// precomputed in update_nonzero_basis_products() for efficiency
  std::vector<UShort2DMultiSetRealMap> nonZerosMapArray;
};


inline SharedNodalInterpPolyApproxData::
SharedNodalInterpPolyApproxData(short basis_type, size_t num_vars):
  SharedInterpPolyApproxData(basis_type, num_vars)//,
  //momentInterpType(INTERPOLATION_OF_PRODUCTS)
  //momentInterpType(REINTERPOLATION_OF_PRODUCTS)
  //momentInterpType(PRODUCT_OF_INTERPOLANTS_FULL)
  //momentInterpType(PRODUCT_OF_INTERPOLANTS_FAST)
{ }


inline SharedNodalInterpPolyApproxData::
SharedNodalInterpPolyApproxData(short basis_type, size_t num_vars,
				const ExpansionConfigOptions& ec_options,
				const BasisConfigOptions& bc_options):
  SharedInterpPolyApproxData(basis_type, num_vars, ec_options, bc_options)//,
  // These 4 compile-time options are relevant for all-variables covariance
  // involving expectations over variable subsets.  Covariance for hierarchical
  // interpolants, nodal covariance in the standard view mode, uses of
  // PolynomialApproximation::integrate_moments(), and Sobol' index
  // calculations all employ an INTERPOLATION_OF_PRODUCTS approach, so that
  // setting is the most self-consistent.  Gradient enhancement is also not
  // currently supported for PRODUCT_OF_INTERPOLANTS approaches.
  //momentInterpType(INTERPOLATION_OF_PRODUCTS)
  //momentInterpType(REINTERPOLATION_OF_PRODUCTS)
  //momentInterpType(PRODUCT_OF_INTERPOLANTS_FULL)
  //momentInterpType(PRODUCT_OF_INTERPOLANTS_FAST)
{ }


inline SharedNodalInterpPolyApproxData::~SharedNodalInterpPolyApproxData()
{ }


inline bool SharedNodalInterpPolyApproxData::
basis_product_1d(InterpolationPolynomial* poly_rep_1,
		 InterpolationPolynomial* poly_rep_2,
		 unsigned short key_1, unsigned short key_2,
		 const RealArray& pts, const RealArray& wts, Real& prod)
{
  Real tol = 1.e-12; // consistent with OrthogonalPolynomial triple product tol
  prod = 0.;
  size_t i, num_pts = pts.size();
  for (i=0; i<num_pts; ++i)
    prod += wts[i] * poly_rep_1->type1_value(pts[i], key_1)
                   * poly_rep_2->type1_value(pts[i], key_2);
  return (std::abs(prod) > tol) ? true : false;
}


inline void SharedNodalInterpPolyApproxData::
update_member_key(const UShortArray& data,
		  const SizetList&   member_indices,
		  UShortArray& member_map_key, size_t cntr)
{
  for (SizetList::const_iterator cit=member_indices.begin();
       cit!=member_indices.end(); ++cit, ++cntr)
    member_map_key[cntr] = data[*cit];
}


inline void SharedNodalInterpPolyApproxData::
set_new_point(const RealVector& x, const UShortArray& basis_index, short order)
{
  unsigned short bi_j;
  for (size_t j=0; j<numVars; ++j) {
    bi_j = basis_index[j];
    if (bi_j) // exclusion of pt must be sync'd w/ factors/scalings
      polynomialBasis[bi_j][j].set_new_point(x[j], order);
  }
}


inline void SharedNodalInterpPolyApproxData::
set_new_point(const RealVector& x, const UShortArray& basis_index,
	      const SizetList& subset_indices, short order)
{
  SizetList::const_iterator cit; size_t j; unsigned short bi_j;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j = *cit; bi_j = basis_index[j];
    if (bi_j) // exclusion of pt must be sync'd w/ factors/scalings
      polynomialBasis[bi_j][j].set_new_point(x[j], order);
  }
}


inline size_t SharedNodalInterpPolyApproxData::
barycentric_exact_index(const UShortArray& basis_index)
{
  size_t j, pt_index = 0, prod = 1; unsigned short bi_j;
  for (j=0; j<numVars; ++j) {
    bi_j = basis_index[j];
    // Note: if bi_j == 0, then constant interp with 1 point: we can replace
    // this constant interpolation with the value at the 1 colloc index (ei=0)
    if (bi_j) {
      BasisPolynomial& poly_i = polynomialBasis[bi_j][j];
      pt_index += poly_i.exact_index() * prod;
      prod     *= poly_i.interpolation_size();
    }
  }
  return pt_index;
}


inline size_t SharedNodalInterpPolyApproxData::
barycentric_exact_index(const UShortArray& basis_index,
			const SizetList& subset_indices)
{
  size_t j, pt_index = 0, prod = 1; unsigned short bi_j;
  SizetList::const_iterator cit;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j = *cit; bi_j = basis_index[j];
    // Note: if bi_j == 0, then constant interp with 1 point: we can replace
    // this constant interpolation with the value at the 1 colloc index (ei=0)
    if (bi_j) {
      BasisPolynomial& poly_j = polynomialBasis[bi_j][j];
      pt_index += poly_j.exact_index() * prod;
      prod     *= poly_j.interpolation_size();
    }
  }
  return pt_index;
}


inline TensorProductDriver* SharedNodalInterpPolyApproxData::tpq_driver()
{ return (TensorProductDriver*)driverRep; }


inline CombinedSparseGridDriver* SharedNodalInterpPolyApproxData::csg_driver()
{ return (CombinedSparseGridDriver*)driverRep; }

} // namespace Pecos

#endif
