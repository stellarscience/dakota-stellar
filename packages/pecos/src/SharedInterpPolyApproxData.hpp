/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedInterpPolyApproxData
//- Description:  Class for Lagrange Interpolation Polynomial ApproxData
//-               
//- Owner:        Mike Eldred

#ifndef SHARED_INTERP_POLY_APPROX_HPP
#define SHARED_INTERP_POLY_APPROX_HPP

#include "SharedPolyApproxData.hpp"
#include "BasisPolynomial.hpp"

namespace Pecos {


/// Derived approximation class for interpolation polynomials (global
/// approximation).

/** The SharedInterpPolyApproxData class provides a global approximation
    based on interpolation polynomials.  It is used primarily for
    stochastic collocation approaches to uncertainty quantification. */

class SharedInterpPolyApproxData: public SharedPolyApproxData
{
  //
  //- Heading: Friends
  //

  friend class InterpPolyApproximation;

public:

  //
  //- Heading: Constructor and destructor
  //

  /// standard constructor
  SharedInterpPolyApproxData(short basis_type, size_t num_vars);
  /// alternate constructor
  SharedInterpPolyApproxData(short basis_type, size_t num_vars,
			     const ExpansionConfigOptions& ec_options,
			     const BasisConfigOptions&     bc_options);
  /// destructor
  ~SharedInterpPolyApproxData();

  //
  //- Heading: member functions
  //

  /// define n-D basis types and collocation rules based on u_types and
  /// basis configuration options for use in configuring integration drivers
  /** These basis types may include orthogonal polynomials for purposes of
      computing their Gauss points and weights within integration drivers;
      thus they differ in general from the interpolation polynomial basis
      used for approximation. */
  static bool initialize_driver_types_rules(const ShortArray& u_types,
    const BasisConfigOptions& bc_options, ShortArray& basis_types,
    ShortArray& colloc_rules);
  /// initialize basis types and collocation rules, construct a vector of basis
  /// polynomials for driver usage (not the vector<vector<BasisPolynomial> >
  /// used herein), and initialize distribution parameters within this basis
  static void construct_basis(const ShortArray& u_types,
			      const AleatoryDistParams& adp,
			      const BasisConfigOptions& bc_options,
			      std::vector<BasisPolynomial>& poly_basis);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void allocate_data();
  void increment_data();
  void decrement_data();
  void post_push_data();
  void post_finalize_data();
  size_t pre_combine_data(short combine_type);
  void post_combine_data(short combine_type);
  void store_data(size_t index = _NPOS);
  void restore_data(size_t index = _NPOS);
  void remove_stored_data(size_t index = _NPOS);

  //
  //- Heading: New virtual functions
  //

  /// allocate sobolIndexMap based on collocation keys
  virtual void allocate_component_sobol() = 0;
  /// update sobolIndexMap based on a refinement increment
  virtual void increment_component_sobol() = 0;

  /// set point values within 1D basis polynomials for purposes of
  /// barycentric precomputation
  virtual void set_new_point(const RealVector& x,
			     const UShortArray& basis_index, short order) = 0;
  /// set point values within subset of 1D basis polynomials for purposes of
  /// barycentric precomputation
  virtual void set_new_point(const RealVector& x,
			     const UShortArray& basis_index,
			     const SizetList& subset_indices, short order) = 0;

  /// for an exact point match in all dimensions, return the
  /// tensor-product index of the matching point
  virtual size_t barycentric_exact_index(const UShortArray& basis_index) = 0;
  /// for an exact point match in all dimensions within subset_indices,
  /// return the tensor-product index of the matching point
  virtual size_t barycentric_exact_index(const UShortArray& basis_index,
					 const SizetList& subset_indices) = 0;

  /// if needed for efficiency, precompute the count and max values
  /// for the collocation keys
  virtual void precompute_keys(const UShortArray& basis_index);
  /// if needed for efficiency, precompute the count and max values
  /// for the collocation keys within the subset
  virtual void precompute_keys(const UShortArray& basis_index,
			       const SizetList& subset_indices);
  /// if needed for efficiency, precompute the max values for the
  /// collocation keys
  virtual void precompute_max_keys(const UShortArray& basis_index);
  /// if needed for efficiency, precompute the max values for the
  /// collocation keys within the subset
  virtual void precompute_max_keys(const UShortArray& basis_index,
				   const SizetList& subset_indices);
  /// return the number of collocation keys for i^{th} variable and
  /// level_i index set
  virtual unsigned short tensor_product_num_key(size_t i,
						unsigned short level_i);
  /// return the maximum collocation key value for i^{th} variable and
  /// level_i index set
  virtual unsigned short tensor_product_max_key(size_t i,
						unsigned short level_i);

  //
  //- Heading: Convenience functions
  //

  /// return value of type 1 interpolation polynomial using all dimensions
  Real type1_interpolant_value(const RealVector& x, const UShortArray& key,
			       const UShortArray& basis_index);
  /// return value of type 1 interpolation polynomial using interpolated
  /// (non-integrated) dimension subset
  Real type1_interpolant_value(const RealVector& x, const UShortArray& key,
			       const UShortArray& basis_index,
			       const SizetList& subset_indices);
  /// return gradient of type 1 interpolation polynomial using all dimensions
  Real type1_interpolant_gradient(const RealVector& x, size_t deriv_index,
				  const UShortArray& key,
				  const UShortArray& basis_index);
  /// return gradient of type 1 interpolation polynomial using interpolated
  /// (non-integrated) dimension subset
  Real type1_interpolant_gradient(const RealVector& x, size_t deriv_index,
				  const UShortArray& key,
				  const UShortArray& basis_index,
				  const SizetList& subset_indices);

  /// return value of type 2 interpolation polynomial using all dimensions
  Real type2_interpolant_value(const RealVector& x, size_t interp_index,
			       const UShortArray& key,
			       const UShortArray& basis_index);
  /// return value of type 2 interpolation polynomial using interpolated
  /// (non-integrated) dimension subset
  Real type2_interpolant_value(const RealVector& x, size_t interp_index,
			       const UShortArray& key,
			       const UShortArray& basis_index,
			       const SizetList& subset_indices);
  /// return gradient of type 2 interpolation polynomial using all dimensions
  Real type2_interpolant_gradient(const RealVector& x, size_t deriv_index,
				  size_t interp_index, const UShortArray& key,
				  const UShortArray& basis_index);
  /// return gradient of type 2 interpolation polynomial using interpolated
  /// (non-integrated) dimension subset
  Real type2_interpolant_gradient(const RealVector& x, size_t deriv_index,
				  size_t interp_index, const UShortArray& key,
				  const UShortArray& basis_index,
				  const SizetList& subset_indices);

  /// compute the number of variables that are active for barycentric
  /// interpolation
  size_t barycentric_active_variables(const UShortArray& basis_index);
  /// compute the number of variables within subset_indices that are
  /// active for barycentric interpolation
  size_t barycentric_active_variables(const UShortArray& basis_index,
				      const SizetList& subset_indices);

  /// compute the product of 1D barycentric value factors
  Real barycentric_value_factor(const UShortArray& key,
				const UShortArray& basis_index);
  /// compute the product of a subset of 1D barycentric value factors
  Real barycentric_value_factor(const UShortArray& key,
				const UShortArray& basis_index,
				const SizetList& subset_indices);
  /// compute the product of 1D barycentric weight factor sums for use in
  /// the denominator of the barycentric interpolation formula (second form)
  Real barycentric_value_scaling(const UShortArray& basis_index);
  /// compute the product of a subset of 1D barycentric weight factor
  /// sums for use in the denominator of the barycentric interpolation
  /// formula (second form)
  Real barycentric_value_scaling(const UShortArray& basis_index,
				 const SizetList& subset_indices);

  /// compute the product of 1D barycentric gradient factors
  Real barycentric_gradient_factor(size_t deriv_index, const UShortArray& key,
				   const UShortArray& basis_index);
  /// compute the product of a subset of 1D barycentric gradient factors
  Real barycentric_gradient_factor(size_t deriv_index, const UShortArray& key,
				   const UShortArray& basis_index,
				   const SizetList& subset_indices);
  /// compute the product of 1D barycentric gradient scalings
  Real barycentric_gradient_scaling(const UShortArray& basis_index);
  /// compute the product of a subset of 1D barycentric gradient scalings
  Real barycentric_gradient_scaling(const UShortArray& basis_index,
				    const SizetList& subset_indices);

  /// return type 1 product weight from integration of type 1 interpolation
  /// polynomials using integrated dimension subset
  Real type1_weight(const UShortArray& key, const UShortArray& basis_index, 
		    const SizetList& subset_indices);
  /// return type 1 product weights from integration of type 1 interpolation
  /// polynomials for both member and nonmember sets
  void type1_weight(const UShortArray& key, const UShortArray& basis_index, 
		    const BitArray& member_bits, Real& member_t1_wt_prod,
		    Real& nonmember_t1_wt_prod);

  /// return type 2 product weight from integration of type 1/2 interpolation
  /// polynomials using integrated dimension subset
  Real type2_weight(size_t interp_index, const UShortArray& key,
		    const UShortArray& basis_index,
		    const SizetList& subset_indices);
  /// return type 2 product weight from integration of type 1/2 interpolation
  /// polynomials for both member and nonmember sets
  void type2_weight(size_t interp_index, const UShortArray& key,
		    const UShortArray& basis_index, const BitArray& member_bits,
		    Real& member_t2_wt_prod, Real& nonmember_t2_wt_prod);

  /// compute the value of a tensor interpolant on a tensor grid;
  /// contributes to value(x)
  Real tensor_product_value(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShortArray& basis_index,  const UShort2DArray& key,
    const SizetArray&  colloc_index);
  /// compute the value of a tensor interpolant on a tensor grid over
  /// a subset of the variables; contributes to value(x)
  Real tensor_product_value(const RealVector& x,
    const RealVector& subset_t1_coeffs, const RealMatrix& subset_t2_coeffs,
    const UShortArray& basis_index,  const UShort2DArray& subset_key,
    const SizetArray& subset_colloc_index, const SizetList& subset_indices);

  /// compute the gradient of a tensor interpolant on a tensor grid
  /// with respect to variables that are included in the polynomial
  /// basis; contributes to gradient_basis_variables(x)
  const RealVector& tensor_product_gradient_basis_variables(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShortArray& basis_index,  const UShort2DArray& key,
    const SizetArray&  colloc_index);
  /// compute the gradient of a tensor interpolant on a tensor grid
  /// with respect to variables that are included in the polynomial
  /// basis; contributes to gradient_basis_variables(x)
  const RealVector& tensor_product_gradient_basis_variables(const RealVector& x,
    const RealVector& subset_t1_coeffs, const RealMatrix& subset_t2_coeffs,
    const UShortArray& basis_index,  const UShort2DArray& subset_key,
    const SizetArray& subset_colloc_index, const SizetList& subset_indices);
  /// compute the gradient of a tensor interpolant on a tensor grid
  /// with respect to variables that are included in the polynomial
  /// basis for given DVV; contributes to gradient_basis_variables(x, dvv)
  const RealVector& tensor_product_gradient_basis_variables(const RealVector& x,
    const RealVector& exp_t1_coeffs, const RealMatrix& exp_t2_coeffs,
    const UShortArray& basis_index,  const UShort2DArray& key,
    const SizetArray& colloc_index,  const SizetArray& dvv);
  /// compute the gradient of a tensor interpolant on a tensor grid
  /// with respect to variables that are not included in the
  /// polynomial basis; contributes to gradient_nonbasis_variables(x)
  const RealVector& tensor_product_gradient_nonbasis_variables(
    const RealVector& x,             const RealMatrix& exp_t1_coeff_grads,
    const UShortArray& basis_index,  const UShort2DArray& key,
    const SizetArray& colloc_index);

  /// resize polynomialBasis to accomodate an update in max interpolation level
  void resize_polynomial_basis(unsigned short max_level);
  /// resize polynomialBasis to accomodate an update in interpolation levels
  void resize_polynomial_basis(const UShortArray& lev_index);
  /// update polynomialBasis for a subset of variables after a change
  /// in quadrature order
  void update_tensor_interpolation_basis(const UShortArray& lev_index,
					 const SizetList& subset_indices);
  /// for a particular level, test for equality between basis v2 and basis v1
  bool same_basis(unsigned short level, size_t v1, size_t v2);

  //
  //- Heading: Data
  //

  /// 2D array of one-dimensional basis polynomial objects used in
  /// constructing the multivariate orthogonal/interpolation polynomials.
  /** Each variable (inner array size = numVars) has multiple
      integration orders associated with it (outer array size). */
  std::vector<std::vector<BasisPolynomial> > polynomialBasis;

  /// flag indicating use of barycentric interpolation for global
  /// value-based Lagrange interpolation
  bool barycentricFlag;

private:

  //
  //- Heading: Convenience functions
  //

  /// define the 1D basis type and collocation rule
  void initialize_polynomial_basis_type(short& poly_type_1d, short& rule);

  /// update polynomialBasis after a change in quadrature order
  void update_tensor_interpolation_basis(const UShortArray& lev_index);
  /// update polynomialBasis after a change in sparse grid level
  void update_sparse_interpolation_basis(unsigned short max_level);
  /// update polynomialBasis for a variable index after an update in level
  void update_interpolation_basis(unsigned short lev_index, size_t var_index);
  /// for a particular level, find index of basis v2 that matches basis v1
  bool find_basis(unsigned short level, size_t v1, size_t& v2);

  /// compute the product of 1D barycentric value factors
  bool barycentric_value_factor(BasisPolynomial& pb_lv, unsigned short key_lv,
				Real& prod);

  /// shared utility for barycentric interpolation over a partial variable
  /// subset: define pt_factors, act_v_set, num_act_pts, and return pt_index
  void barycentric_partial_indexing(const UShortArray& basis_index,
				    //UShortList& num_keys,
				    SizetList& pt_factors, SizetList& act_v_set,
				    size_t& num_act_pts, size_t& pt_index);
  /// shared utility for barycentric interpolation over a partial variable
  /// subset: define pt_factors, act_v_set, num_act_pts, and return pt_index
  void barycentric_partial_indexing(const UShortArray& basis_index,
				    const SizetList& subset_indices,
				    //UShortList& num_keys,
				    SizetList& pt_factors, SizetList& act_v_set,
				    size_t& num_act_pts, size_t& pt_index);

  /// shared code for barycentric interpolation over an active variable subset
  void accumulate_barycentric_partial(const RealVector& t1_coeffs,
				      const UShortArray& basis_index,
				      const UShort2DArray& key,
				      const SizetArray& colloc_index,
				      //const UShortList& num_keys,
				      const SizetList& pt_factors,
				      const SizetList& act_v_set,
				      size_t num_act_pts, size_t pt_index,
				      RealVector& accumulator);
  /// shared code for barycentric gradient evaluation for active variable 0
  void accumulate_barycentric_gradient(unsigned short bi0,
				       unsigned short key_i0,
				       size_t ei_0, Real* accum_0,
				       Real t1_coeff, const RealVector& bc_vf_0,
				       const RealVector& bc_gf_0);
  /// shared code for barycentric gradient evaluation for active variables 1:n
  void accumulate_barycentric_gradient(size_t j, unsigned short bij,
				       unsigned short key_ij,
				       BasisPolynomial& poly_j,
				       RealMatrix& accumulator);

  //
  //- Heading: Data
  //

  /// the gradient of a tensor-product interpolant; a contributor to
  /// approxGradient
  RealVector tpGradient;
};


inline SharedInterpPolyApproxData::
SharedInterpPolyApproxData(short basis_type, size_t num_vars):
  SharedPolyApproxData(basis_type, num_vars)
{ }


inline SharedInterpPolyApproxData::
SharedInterpPolyApproxData(short basis_type, size_t num_vars,
			   const ExpansionConfigOptions& ec_options,
			   const BasisConfigOptions&     bc_options):
  SharedPolyApproxData(basis_type, num_vars, ec_options, bc_options)
{ }


inline SharedInterpPolyApproxData::~SharedInterpPolyApproxData()
{ }


inline void SharedInterpPolyApproxData::
construct_basis(const ShortArray& u_types, const AleatoryDistParams& adp,
		const BasisConfigOptions& bc_options,
		std::vector<BasisPolynomial>& poly_basis)
{
  ShortArray basis_types, colloc_rules;
  bool dist_params = initialize_driver_types_rules(u_types, bc_options,
						   basis_types, colloc_rules);
  initialize_polynomial_basis(basis_types, colloc_rules, poly_basis);
  if (dist_params)
    update_basis_distribution_parameters(u_types, adp, poly_basis);
}


inline void SharedInterpPolyApproxData::
resize_polynomial_basis(unsigned short max_level)
{
  size_t i, basis_size = polynomialBasis.size();
  if (max_level >= basis_size) {
    polynomialBasis.resize(max_level+1);
    for (i=basis_size; i<=max_level; ++i)
      polynomialBasis[i].resize(numVars);
  }
}


inline void SharedInterpPolyApproxData::
resize_polynomial_basis(const UShortArray& lev_index)
{
  unsigned short max_level = lev_index[0];
  for (size_t i=1; i<numVars; ++i)
    if (lev_index[i] > max_level)
      max_level = lev_index[i];
  // For tensor quadrature, order range is 1:m; level range is 0:m-1
  resize_polynomial_basis(max_level);
}


inline void SharedInterpPolyApproxData::store_data(size_t index)
{ driverRep->store_grid(index); }


inline void SharedInterpPolyApproxData::restore_data(size_t index)
{ driverRep->restore_grid(index); }


inline void SharedInterpPolyApproxData::remove_stored_data(size_t index)
{ driverRep->remove_stored_grid(index); }


inline Real SharedInterpPolyApproxData::
type1_interpolant_value(const RealVector& x, const UShortArray& key,
			const UShortArray& basis_index)
{
  Real L1 = 1.;
  for (size_t j=0; j<numVars; ++j)
    L1 *= polynomialBasis[basis_index[j]][j].type1_value(x[j], key[j]);
  return L1;
}


/** Combined expansion version. */
inline Real SharedInterpPolyApproxData::
type1_interpolant_value(const RealVector& x, const UShortArray& key,
			const UShortArray& basis_index,
			const SizetList& subset_indices)
{
  Real L1 = 1.; SizetList::const_iterator cit; size_t j;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j   = *cit;
    L1 *= polynomialBasis[basis_index[j]][j].type1_value(x[j], key[j]);
  }
  return L1;
}


inline Real SharedInterpPolyApproxData::
type1_interpolant_gradient(const RealVector& x, size_t deriv_index,
			   const UShortArray& key,
			   const UShortArray& basis_index)
{
  Real L1_grad = 1.;
  for (size_t k=0; k<numVars; ++k)
    L1_grad *= (k == deriv_index) ?
      polynomialBasis[basis_index[k]][k].type1_gradient(x[k], key[k]) :
      polynomialBasis[basis_index[k]][k].type1_value(x[k],    key[k]);
  return L1_grad;
}


/** Combined expansion version. */
inline Real SharedInterpPolyApproxData::
type1_interpolant_gradient(const RealVector& x, size_t deriv_index,
			   const UShortArray& key,
			   const UShortArray& basis_index,
			   const SizetList& subset_indices)
{
  // deriv_index must be contained within subset_indices, else the grad is zero
  bool deriv = false;

  Real L1_grad = 1.; SizetList::const_iterator cit; size_t k;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    k = *cit;
    if (k == deriv_index) {
      L1_grad *= polynomialBasis[basis_index[k]][k].type1_gradient(x[k],key[k]);
      deriv = true;
    }
    else
      L1_grad *= polynomialBasis[basis_index[k]][k].type1_value(x[k], key[k]);
  }
  return (deriv) ? L1_grad : 0.;
}


inline Real SharedInterpPolyApproxData::
type2_interpolant_value(const RealVector& x,    size_t interp_index,
			const UShortArray& key, const UShortArray& basis_index)
{
  Real L2 = 1.;
  for (size_t k=0; k<numVars; ++k)
    L2 *= (interp_index == k) ?
      polynomialBasis[basis_index[k]][k].type2_value(x[k], key[k]) :
      polynomialBasis[basis_index[k]][k].type1_value(x[k], key[k]);
  return L2;
}


/** Combined expansion version. */
inline Real SharedInterpPolyApproxData::
type2_interpolant_value(const RealVector& x,    size_t interp_index,
			const UShortArray& key, const UShortArray& basis_index,
			const SizetList& subset_indices)
{
  Real L2 = 1.; SizetList::const_iterator cit; size_t k;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    k   = *cit;
    L2 *= (interp_index == k) ?
      polynomialBasis[basis_index[k]][k].type2_value(x[k], key[k]) :
      polynomialBasis[basis_index[k]][k].type1_value(x[k], key[k]);
  }
  return L2;
}


inline Real SharedInterpPolyApproxData::
type2_interpolant_gradient(const RealVector& x, size_t deriv_index,
			   size_t interp_index, const UShortArray& key,
			   const UShortArray& basis_index)
{
  //  deriv_index = desired gradient component
  // interp_index = index of gradient component used in type2 interpolation
  Real L2_grad = 1.;
  for (size_t l=0; l<numVars; ++l)
    if (l == deriv_index)
      L2_grad *= (l == interp_index) ?
	polynomialBasis[basis_index[l]][l].type2_gradient(x[l], key[l]) :
	polynomialBasis[basis_index[l]][l].type1_gradient(x[l], key[l]);
    else
      L2_grad *= (l == interp_index) ?
	polynomialBasis[basis_index[l]][l].type2_value(x[l], key[l]) :
	polynomialBasis[basis_index[l]][l].type1_value(x[l], key[l]);
  return L2_grad;
}


/** Combined expansion version. */
inline Real SharedInterpPolyApproxData::
type2_interpolant_gradient(const RealVector& x, size_t deriv_index,
			   size_t interp_index, const UShortArray& key,
			   const UShortArray& basis_index,
			   const SizetList& subset_indices)
{
  //  deriv_index = desired gradient component
  // interp_index = index of gradient component used in type2 interpolation

  // deriv_index must be contained within subset_indices, else the grad is zero
  bool deriv = false;

  Real L2_grad = 1.; SizetList::const_iterator cit; size_t l;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    l = *cit;
    if (l == deriv_index) {
      L2_grad *= (l == interp_index) ?
	polynomialBasis[basis_index[l]][l].type2_gradient(x[l], key[l]) :
	polynomialBasis[basis_index[l]][l].type1_gradient(x[l], key[l]);
      deriv = true;
    }
    else
      L2_grad *= (l == interp_index) ?
	polynomialBasis[basis_index[l]][l].type2_value(x[l], key[l]) :
	polynomialBasis[basis_index[l]][l].type1_value(x[l], key[l]);
  }
  return (deriv) ? L2_grad : 0.;
}


inline size_t SharedInterpPolyApproxData::
barycentric_active_variables(const UShortArray& basis_index)
{
  size_t j, num_act_v = 0; unsigned short bi_j;
  for (j=0; j<numVars; ++j) {
    bi_j = basis_index[j]; // if 0, then constant interp with 1 pt
    if (bi_j && polynomialBasis[bi_j][j].exact_index() == _NPOS) // active
      ++num_act_v;
  }
  return num_act_v;
}


inline size_t SharedInterpPolyApproxData::
barycentric_active_variables(const UShortArray& basis_index,
			     const SizetList& subset_indices)
{
  size_t j, num_act_v = 0; unsigned short bi_j; SizetList::const_iterator cit;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j = *cit; bi_j = basis_index[j];
    if (bi_j && polynomialBasis[bi_j][j].exact_index() == _NPOS) // active
      ++num_act_v;
  }
  return num_act_v;
}


/** Default implementation; overridden by HierarchInterpPolyApproxData. */
inline void SharedInterpPolyApproxData::
precompute_keys(const UShortArray& basis_index)
{ } // no-op


/** Default implementation; overridden by HierarchInterpPolyApproxData. */
inline void SharedInterpPolyApproxData::
precompute_keys(const UShortArray& basis_index, const SizetList& subset_indices)
{ } // no-op


/** Default implementation; overridden by HierarchInterpPolyApproxData. */
inline void SharedInterpPolyApproxData::
precompute_max_keys(const UShortArray& basis_index)
{ } // no-op


/** Default implementation; overridden by HierarchInterpPolyApproxData. */
inline void SharedInterpPolyApproxData::
precompute_max_keys(const UShortArray& basis_index,
		    const SizetList& subset_indices)
{ } // no-op


/** Default implementation; overridden by HierarchInterpPolyApproxData. */
inline unsigned short SharedInterpPolyApproxData::
tensor_product_num_key(size_t i, unsigned short level_i)
{ return polynomialBasis[level_i][i].interpolation_size(); }


/** Default implementation; overridden by HierarchInterpPolyApproxData. */
inline unsigned short SharedInterpPolyApproxData::
tensor_product_max_key(size_t i, unsigned short level_i)
{ return polynomialBasis[level_i][i].interpolation_size() - 1; }


inline void SharedInterpPolyApproxData::
barycentric_partial_indexing(const UShortArray& basis_index,
			     /*UShortList& num_keys,*/ SizetList& pt_factors,
			     SizetList& act_v_set, size_t& num_act_pts,
			     size_t& pt_index)
{
  // define interpolation set and initial pt_index offset
  size_t j, num_pts = 1, ei_j, edi_j;
  unsigned short pts_j, bi_j;
  num_act_pts = 1; pt_index = 0;
  precompute_keys(basis_index); // precompute count/max if needed for efficiency
  for (j=0; j<numVars; ++j) {
    bi_j = basis_index[j];
    if (bi_j) { // else pts_j = 1 and ei_j can be taken to be 0
      BasisPolynomial& poly_j = polynomialBasis[bi_j][j];
      ei_j = poly_j.exact_index(); pts_j = tensor_product_num_key(j, bi_j);
      if (ei_j == _NPOS) { // active for interpolation
	pt_factors.push_back(num_pts);   act_v_set.push_back(j);
	/* num_keys.push_back(pts_j); */ num_act_pts *= pts_j;
      }
      else {             // inactive for interpolation
	edi_j = poly_j.exact_delta_index();
	if (edi_j == _NPOS)         // exactIndex match but not exactDeltaIndex;
	  { pt_index = _NPOS; break; }   // at least 1 dim has value factor = 0.
	else
	  pt_index += num_pts * edi_j;
      }
      num_pts *= pts_j;
    }
  }
}


inline void SharedInterpPolyApproxData::
barycentric_partial_indexing(const UShortArray& basis_index,
			     const SizetList& subset_indices,
			     /*UShortList& num_keys,*/ SizetList& pt_factors,
			     SizetList& act_v_set, size_t& num_act_pts,
			     size_t& pt_index)
{
  // define interpolation set and initial pt_index offset
  size_t j, num_pts = 1, ei_j, edi_j;
  unsigned short pts_j, bi_j; SizetList::const_iterator cit;
  num_act_pts = 1; pt_index = 0;
  precompute_keys(basis_index, subset_indices); // if needed for efficiency
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j = *cit; bi_j = basis_index[j];
    if (bi_j) { // else pts_j = 1 and ei_j can be taken to be 0
      BasisPolynomial& poly_j = polynomialBasis[bi_j][j];
      ei_j = poly_j.exact_index(); pts_j = tensor_product_num_key(j, bi_j);
      if (ei_j == _NPOS) { // active for interpolation
	pt_factors.push_back(num_pts);   act_v_set.push_back(j);
	/* num_keys.push_back(pts_j); */ num_act_pts *= pts_j;
      }
      else {             // inactive for interpolation
	edi_j = poly_j.exact_delta_index();
	if (edi_j == _NPOS)         // exactIndex match but not exactDeltaIndex;
	  { pt_index = _NPOS; break; }   // at least 1 dim has value factor = 0.
	else
	  pt_index += num_pts * edi_j;
      }
      num_pts *= pts_j;
    }
  }
}


inline Real SharedInterpPolyApproxData::
barycentric_value_factor(const UShortArray& key, const UShortArray& basis_index)
{
  Real b1 = 1., vf; unsigned short bi_j;
  for (size_t j=0; j<numVars; ++j) {
    bi_j = basis_index[j];
    if (bi_j) { // exclusion of bi==0 must be sync'd with bc value scaling
      vf = polynomialBasis[bi_j][j].barycentric_value_factor(key[j]);
      if (vf == 0.)    { b1  = 0.; break; } // key[j] != exactIndex
      else if (vf != 1.) b1 *= vf;          // exactIndex == _NPOS
    }
  }
  return b1;
}


inline Real SharedInterpPolyApproxData::
barycentric_value_factor(const UShortArray& key, const UShortArray& basis_index,
			 const SizetList& subset_indices)
{
  Real b1 = 1., vf;
  SizetList::const_iterator cit; size_t j; unsigned short bi_j;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j = *cit; bi_j = basis_index[j];
    if (bi_j) { // exclusion of bi==0 must be sync'd with bc value scaling
      vf = polynomialBasis[bi_j][j].barycentric_value_factor(key[j]);
      if (vf == 0.)    { b1  = 0.; break; } // key[j] != exactIndex
      else if (vf != 1.) b1 *= vf;          // exactIndex == _NPOS
    }
  }
  return b1;
}


inline Real SharedInterpPolyApproxData::
barycentric_value_scaling(const UShortArray& basis_index)
{
  Real scale = 1.; unsigned short bi_j;
  for (size_t j=0; j<numVars; ++j) {
    bi_j = basis_index[j];
    if (bi_j) { // exclusion of bi==0 must be sync'd with bc value factors
      BasisPolynomial& pb_lv = polynomialBasis[bi_j][j];
      if (pb_lv.exact_index() == _NPOS) // this dim contributes to bc denom
	scale *= pb_lv.barycentric_value_factor_sum();
      //else (new pt is exact match) dim does not contribute to bc denom
    }
  }
  return scale;
}


inline Real SharedInterpPolyApproxData::
barycentric_value_scaling(const UShortArray& basis_index,
			  const SizetList& subset_indices)
{
  Real scale = 1.; SizetList::const_iterator cit; size_t j; unsigned short bi_j;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j = *cit; bi_j = basis_index[j];
    if (bi_j) { // exclusion of bi==0 must be sync'd with bc value factors
      BasisPolynomial& pb_lv = polynomialBasis[bi_j][j];
      if (pb_lv.exact_index() == _NPOS) // this dim contributes to bc denom
	scale *= pb_lv.barycentric_value_factor_sum();
      //else (new pt is exact match) dim does not contribute to bc denom
    }
  }
  return scale;
}


inline Real SharedInterpPolyApproxData::
barycentric_gradient_factor(size_t deriv_index, const UShortArray& key,
			    const UShortArray& basis_index)
{
  Real factor = 1., vf; unsigned short bi_k;
  for (size_t k=0; k<numVars; ++k) {
    bi_k = basis_index[k];
    if (bi_k) { // else omit value factor
      BasisPolynomial& pb_lv = polynomialBasis[bi_k][k];
      if (k == deriv_index) // gradient factor, either exact index or not
	factor *= pb_lv.barycentric_gradient_factor(key[k]);
      else {                // value factor, zero if exact_index
	vf = pb_lv.barycentric_value_factor(key[k]);
	if (vf == 0.)    { factor  = 0.; break; } // key[k] != exactIndex
	else if (vf != 1.) factor *= vf;          // exactIndex == _NPOS
      }
    }
    else if (k == deriv_index) // deriv of constant interp is 0.
      { factor = 0.; break; }
  }
  return factor;
}


/** Combined expansion version. */
inline Real SharedInterpPolyApproxData::
barycentric_gradient_factor(size_t deriv_index, const UShortArray& key,
			    const UShortArray& basis_index,
			    const SizetList& subset_indices)
{
  // deriv_index must be contained within subset_indices, else the grad is zero
  bool deriv = false;

  Real factor = 1., vf;
  SizetList::const_iterator cit; size_t k; unsigned short bi_k;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    k = *cit; bi_k = basis_index[k];
    if (bi_k) { // else omit value factor
      BasisPolynomial& pb_lv = polynomialBasis[bi_k][k];
      if (k == deriv_index) // gradient factor, either exact index or not
	{ factor *= pb_lv.barycentric_gradient_factor(key[k]); deriv = true; }
      else {                // value factor, zero if exact_index
	vf = pb_lv.barycentric_value_factor(key[k]);
	if (vf == 0.)    { factor  = 0.; break; } // key[k] != exactIndex
	else if (vf != 1.) factor *= vf;          // exactIndex == _NPOS
      }
    }
    else if (k == deriv_index) // deriv of constant interp is 0.
      { factor = 0.; break; }
  }
  return (deriv) ? factor : 0.;
}


inline Real SharedInterpPolyApproxData::
barycentric_gradient_scaling(/* size_t deriv_index, */
			     const UShortArray& basis_index)
{
  Real scale = 1.; unsigned short bi_j;
  for (size_t j=0; j<numVars; ++j) {
    bi_j = basis_index[j];
    if (bi_j) { // else omit factor
      BasisPolynomial& pb_lv = polynomialBasis[bi_j][j];
      // if new pt is not an exact match, then dimension contributes a
      // difference product, irregardless of derivative index (whether
      // a gradient or value factor was applied).
      if (pb_lv.exact_index() == _NPOS)
	scale *= pb_lv.barycentric_difference_product();
      // if new pt is exact match, then dimension doesn't contribute to scale,
      // irregardless of whether a gradient or value factor was applied.
    }
  }
  return scale;
}


inline Real SharedInterpPolyApproxData::
barycentric_gradient_scaling(/* size_t deriv_index, */
			     const UShortArray& basis_index,
			     const SizetList& subset_indices)
{
  Real scale = 1.; SizetList::const_iterator cit; size_t j; unsigned short bi_j;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j = *cit; bi_j = basis_index[j];
    if (bi_j) { // else omit factor
      BasisPolynomial& pb_lv = polynomialBasis[bi_j][j];
      // if new pt is not an exact match, then dimension contributes a
      // difference product, irregardless of derivative index (whether
      // a gradient or value factor was applied).
      if (pb_lv.exact_index() == _NPOS)
	scale *= pb_lv.barycentric_difference_product();
      // if new pt is exact match, then dimension doesn't contribute to scale,
      // irregardless of whether a gradient or value factor was applied.
    }
  }
  return scale;
}


/** Combined expansion partial weight. */
inline Real SharedInterpPolyApproxData::
type1_weight(const UShortArray& key, const UShortArray& basis_index, 
	     const SizetList& subset_indices)
{
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();
  Real t1_wt_prod = 1.; SizetList::const_iterator cit; size_t j;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit)
    { j = *cit; t1_wt_prod *= t1_wts_1d[basis_index[j]][j][key[j]]; }
  return t1_wt_prod;
}


/** VBD member and non-member partial weights. */
inline void SharedInterpPolyApproxData::
type1_weight(const UShortArray& key, const UShortArray& basis_index, 
	     const BitArray& member_bits, Real& member_t1_wt_prod,
	     Real& nonmember_t1_wt_prod)
{
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();
  member_t1_wt_prod = nonmember_t1_wt_prod = 1.;
  size_t j, num_bits = member_bits.size();
  for (j=0; j<num_bits; ++j)
    if (member_bits[j])
      member_t1_wt_prod    *= t1_wts_1d[basis_index[j]][j][key[j]];
    else
      nonmember_t1_wt_prod *= t1_wts_1d[basis_index[j]][j][key[j]];
}


/** Combined expansion partial weight. */
inline Real SharedInterpPolyApproxData::
type2_weight(size_t interp_index, const UShortArray& key,
	     const UShortArray& basis_index, const SizetList& subset_indices)
{
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();
  const Real3DArray& t2_wts_1d = driverRep->type2_collocation_weights_1d();
  Real t2_wt_prod = 1.; SizetList::const_iterator cit; size_t j;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    j           = *cit;
    t2_wt_prod *= (interp_index == j) ? t2_wts_1d[basis_index[j]][j][key[j]]
                                      : t1_wts_1d[basis_index[j]][j][key[j]];
  }
  return t2_wt_prod;
}


/** Combined expansion partial weight. */
inline void SharedInterpPolyApproxData::
type2_weight(size_t interp_index, const UShortArray& key,
	     const UShortArray& basis_index, const BitArray& member_bits,
	     Real& member_t2_wt_prod, Real& nonmember_t2_wt_prod)
{
  const Real3DArray& t1_wts_1d = driverRep->type1_collocation_weights_1d();
  const Real3DArray& t2_wts_1d = driverRep->type2_collocation_weights_1d();
  member_t2_wt_prod = nonmember_t2_wt_prod = 1.;
  size_t j, num_bits = member_bits.size();
  for (j=0; j<num_bits; ++j)
    if (member_bits[j])
      member_t2_wt_prod *= (interp_index == j) ?
	t2_wts_1d[basis_index[j]][j][key[j]] :
	t1_wts_1d[basis_index[j]][j][key[j]];
    else
      nonmember_t2_wt_prod *= (interp_index == j) ?
	t2_wts_1d[basis_index[j]][j][key[j]] :
	t1_wts_1d[basis_index[j]][j][key[j]];
}

} // namespace Pecos

#endif
