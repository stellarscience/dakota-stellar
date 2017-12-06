/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedOrthogPolyApproxData
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        Mike Eldred

#ifndef SHARED_ORTHOG_POLY_APPROX_DATA_HPP
#define SHARED_ORTHOG_POLY_APPROX_DATA_HPP

#include "SharedPolyApproxData.hpp"
#include "NumericGenOrthogPolynomial.hpp"
#include "CombinedSparseGridDriver.hpp"

namespace Pecos {


/// Derived approximation class for orthogonal polynomials (global
/// approximation).

/** The SharedOrthogPolyApproxData class provides a global approximation
    based on orthogonal polynomials.  It is used primarily for polynomial
    chaos expansions (for stochastic finite element approaches to
    uncertainty quantification). */

class SharedOrthogPolyApproxData: public SharedPolyApproxData
{
  //
  //- Heading: Friends
  //

  friend class OrthogPolyApproximation;

public:

  //
  //- Heading: Constructor and destructor
  //

  /// standard constructor
  SharedOrthogPolyApproxData(short basis_type, const UShortArray& approx_order,
			     size_t num_vars);
  /// alternate constructor
  SharedOrthogPolyApproxData(short basis_type, const UShortArray& approx_order,
			     size_t num_vars,
			     const ExpansionConfigOptions& ec_options,
			     const BasisConfigOptions&     bc_options);
  /// destructor
  ~SharedOrthogPolyApproxData();

  //
  //- Heading: Member functions
  //

  /// alternative form for setting multiIndex (expansion import) and
  /// allocating related arrays
  void allocate_data(const UShort2DArray& mi);
  /// get multiIndex
  const UShort2DArray& multi_index() const;

  /// retrieve size of multiIndex
  size_t expansion_terms() const;

  /// get approxOrder
  const UShortArray& expansion_order() const;
  /// set approxOrder
  void expansion_order(const UShortArray& order);
  /// uniformly increment approxOrder
  void increment_order();

  /// invoke initialize_orthogonal_basis_types_rules(),
  /// initialize_polynomial_basis(), and, if needed,
  /// update_basis_distribution_parameters() using class member data
  void construct_basis(const ShortArray& u_types,
		       const AleatoryDistParams& adp);
  
  /// invoke initialize_orthogonal_basis_types_rules(),
  /// initialize_polynomial_basis(), and, if needed,
  /// update_basis_distribution_parameters() using passed data
  static void construct_basis(const ShortArray& u_types,
			      const AleatoryDistParams& adp,
			      const BasisConfigOptions& bc_options,
			      std::vector<BasisPolynomial>& poly_basis);

  /// invoke initialize_orthogonal_basis_types_rules(),
  /// initialize_polynomial_basis(), and, if needed,
  /// update_basis_distribution_parameters() using passed data
  static void construct_basis(const ShortArray& u_types,
			      const AleatoryDistParams& adp,
			      const BasisConfigOptions& bc_options,
			      std::vector<BasisPolynomial>& poly_basis,
			      ShortArray &basis_types,ShortArray &colloc_rules);

  /// set orthogPolyTypes
  void orthogonal_basis_types(const ShortArray& opb_types);
  /// get orthogPolyTypes
  const ShortArray& orthogonal_basis_types() const;

  /// get polynomialBasis
  const std::vector<BasisPolynomial>& polynomial_basis() const;
  std::vector<BasisPolynomial>& polynomial_basis();
  /// set polynomialBasis
  void polynomial_basis(const std::vector<BasisPolynomial>& poly_basis);

  /// set NumericGenOrthogPolynomial::coeffsNormsFlag
  void coefficients_norms_flag(bool flag);

  /// set NumericGenOrthogPolynomial::coeffsNormsFlag
  static void coefficients_norms_flag(bool flag,
				      ShortArray &poly_types,
				      std::vector<BasisPolynomial> &poly_basis);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void allocate_data();

  void store_data(size_t index = _NPOS);
  void restore_data(size_t index = _NPOS);
  void remove_stored_data(size_t index = _NPOS);
  
  size_t pre_combine_data(short combine_type);
  void post_combine_data(short combine_type);

  //
  //- Heading: Member functions
  //

  /// detect whether current expansion settings are the most refined
  size_t maximal_expansion();
  /// swap current data and the stored data set identified by index
  void swap_data(size_t index);

  /// convert a sparse grid index set and a growth setting to an integrand_order
  void sparse_grid_level_to_expansion_order(
    CombinedSparseGridDriver* csg_driver, const UShortArray& levels,
    UShortArray& exp_order);//, short growth_rate = UNRESTRICTED_GROWTH);
  /// convert quadrature orders to integrand orders using rigorous mappings
  void quadrature_order_to_integrand_order(IntegrationDriver* int_driver,
					   const UShortArray& quad_order,
					   UShortArray& int_order);
  /// convert integrand orders to expansion orders using rigorous mappings
  void integrand_order_to_expansion_order(const UShortArray& int_order,
					  UShortArray& exp_order);

  /// helper function for incrementing that is modular on sparse grid driver
  /// and multi-index
  void increment_trial_set(CombinedSparseGridDriver* csg_driver,
			   UShort2DArray& aggregated_mi);
  /// helper function for decrementing that is modular on trial set
  /// and multi-index
  void decrement_trial_set(const UShortArray& trial_set,
			   UShort2DArray& aggregated_mi, bool save_map = true);
  /// helper function for restoring that is modular on trial set and multi-index
  void pre_push_trial_set(const UShortArray& trial_set,
			     UShort2DArray& aggregated_mi,
			     bool monotonic = true);
  /// helper function for restoring that is modular on trial set and multi-index
  void post_push_trial_set(const UShortArray& trial_set,
			      UShort2DArray& aggregated_mi,
			      bool save_map = true);
  /// helper function for restoring that is modular on trial set and multi-index
  void push_trial_set(const UShortArray& trial_set,
			 UShort2DArray& aggregated_mi,
			 bool monotonic = true, bool save_map = true);

  /// Precompute a maximal order of quadrature rules (based on a
  /// multiIndex) when too expensive to compute on demand
  void precompute_maximal_rules(const UShort2DArray& multi_index);
  /// Precompute a maximal order of quadrature rules (based on an
  /// approxOrder) when too expensive to compute on demand
  void precompute_maximal_rules(const UShortArray& approx_order);

  /// allocate sobolIndexMap from multi_index
  void allocate_component_sobol(const UShort2DArray& multi_index);
  /// update sobolIndexMap using new multi_index terms (from multifidelity
  /// overlay or a new QoI in orthogonal least interpolation)
  void update_component_sobol(const UShort2DArray& multi_index);

  /// append multi-indices from append_mi that do not already appear
  /// in combined_mi
  void append_multi_index(const UShort2DArray& append_mi,
			  UShort2DArray& combined_mi);
  /// append multi-indices from append_mi that do not already appear
  /// in combined_mi
  void append_multi_index(const UShortArraySet& append_mi,
			  UShort2DArray& combined_mi);
  /// append multi-indices from append_mi that do not already appear
  /// in combined_mi; define append_mi_map and append_mi_map_ref
  void append_multi_index(const UShort2DArray& append_mi,
			  UShort2DArray& combined_mi, SizetArray& append_mi_map,
			  size_t& append_mi_map_ref);
  /// append multi-indices from append_mi that do not already appear
  /// in combined_mi; define append_mi_map and append_mi_map_ref
  void append_multi_index(const UShort2DArray& append_mi,
			  UShort2DArray& combined_mi, SizetSet& append_mi_map,
			  size_t& append_mi_map_ref);
  /// append multi-indices from append_mi that do not already appear in
  /// combined_mi (consistent ordering assumed); define append_mi_map
  /// and append_mi_map_ref
  void append_leading_multi_index(const UShort2DArray& append_mi,
				  UShort2DArray& combined_mi,
				  SizetSet& append_mi_map,
				  size_t& append_mi_map_ref);
  /// append multi-indices from append_mi that do not already appear
  /// in combined_mi, using previously defined append_mi_map and
  /// append_mi_map_ref for mapping
  void append_multi_index(const UShort2DArray& append_mi,
			  SizetArray& append_mi_map, size_t& append_mi_map_ref,
			  UShort2DArray& combined_mi);
  /// append multi-indices from append_mi that do not already appear in
  /// combined_mi, updating sparse_indices, exp_coeffs, and exp_coeff_grads
  void append_multi_index(SizetSet& sparse_indices,
			  const UShort2DArray& append_mi,
			  UShort2DArray& combined_mi, RealVector& exp_coeffs,
			  RealMatrix& exp_coeff_grads);

  /// returns the norm-squared of a particular multivariate polynomial,
  /// treating all variables as probabilistic
  Real norm_squared(const UShortArray& indices);
  /// returns the norm-squared of a particular multivariate polynomial,
  /// treating a subset of the variables as probabilistic
  Real norm_squared(const UShortArray& indices, const SizetList& rand_indices);

  /// calculate a particular multivariate orthogonal polynomial value
  /// evaluated at a particular parameter set
  Real multivariate_polynomial(const RealVector& x,const UShortArray& indices);
  /// calculate a particular multivariate orthogonal polynomial value over
  /// the nonrandom variable subset evaluated at a particular parameter set
  Real multivariate_polynomial(const RealVector& x,const UShortArray& indices,
			       const SizetList& non_rand_indices);
  /// calculate a particular multivariate orthogonal polynomial value
  /// evaluated at a particular parameter set
  static Real multivariate_polynomial(const RealVector& x,
				      const UShortArray& indices,
	       std::vector<BasisPolynomial> &polynomial_basis);
  /// calculate a particular multivariate orthogonal polynomial value over
  /// the nonrandom variable subset evaluated at a particular parameter set
  static Real multivariate_polynomial(const RealVector& x,
				      const UShortArray& indices,
				      const SizetList& non_rand_indices,
	       std::vector<BasisPolynomial> &polynomial_basis);

  /// compute multivariate orthogonal polynomial gradient evaluated at x 
  /// for term corresponding to indices and derivative variable deriv_index
  Real multivariate_polynomial_gradient(const RealVector& x, size_t deriv_index,
    const UShortArray& indices);
  /// compute multivariate orthogonal polynomial gradient evaluated at x for
  /// term corresponding to indices, for derivative variable deriv_index, and
  /// for a nonrandom variable subset
  Real multivariate_polynomial_gradient(const RealVector& x, size_t deriv_index,
    const UShortArray& indices, const SizetList& non_rand_indices);

  /// compute multivariate orthogonal polynomial Hessian for term
  /// corresponding to deriv_index, evaluated at x
  Real multivariate_polynomial_hessian(const RealVector& x,
    size_t deriv_index_i, size_t deriv_index_j, const UShortArray& indices);

  /// calculate multivariate orthogonal polynomial gradient vector
  /// evaluated at a particular parameter set
  const RealVector& multivariate_polynomial_gradient_vector(const RealVector& x,
    const UShortArray& indices);
  /// calculate multivariate orthogonal polynomial gradient vector with
  /// respect to specified dvv and evaluated at a particular parameter set
  const RealVector& multivariate_polynomial_gradient_vector(const RealVector& x,
    const UShortArray& indices, const SizetArray& dvv);

  /// calculate multivariate orthogonal polynomial gradient vector
  /// evaluated at a particular parameter set
  const RealSymMatrix& multivariate_polynomial_hessian_matrix(
    const RealVector& x, const UShortArray& indices);

  /// define/update a combined Pareto set with a new multi_index by
  /// omitting terms that are weakly Pareto dominated (more omissions
  /// = smaller resulting set)
  void update_pareto_set(const UShort2DArray& multi_index,
			 UShort2DArray& combined_pareto);
  /// define/update a combined Pareto set for a new multi_index term
  void update_pareto_set(const UShortArray& mi_i,
			 UShort2DArray& combined_pareto);
  /// define/update a leading multi_index frontier omitting points that are
  /// strongly Pareto dominated (fewer omissions = larger resulting set)
  void update_frontier(const UShortArraySet& multi_index,
		       UShortArraySet& mi_frontier);
  /// define/update a leading multi_index frontier omitting points that are
  /// strongly Pareto dominated (fewer omissions = larger resulting set)
  void update_frontier(const UShort2DArray& multi_index,
		       UShortArraySet& mi_frontier);
  /// update/update a leading multi_index frontier for a new multi_index term
  void update_frontier(const UShortArray& mi_i, UShortArraySet& mi_frontier);

  // assess whether new_pareto is dominated by total_pareto
  //bool assess_dominance(const UShort2DArray& pareto,
  //			  const UShort2DArray& combined_pareto);
  /// assess bi-directional weak dominance for a "challenger" polynomial
  /// index set against an "incumbent" polynomial index set
  void assess_dominance(const UShortArray& new_order,
			const UShortArray& existing_order,
			bool& new_dominated, bool& existing_dominated);
  /// assess bi-directional strong dominance between two polynomial index sets
  void assess_strong_dominance(const UShortArray& order_a,
			       const UShortArray& order_b,
			       bool& a_dominated, bool& b_dominated);

  /// test for nonzero indices in random variable subset
  bool zero_random(const UShortArray& indices) const;

  /// Generate the coefficient tag for variable j of given expansion term order
  void get_tag(char* tag, size_t j, unsigned short order) const;

  /// tests 1D gradient computations (active in DEBUG compile mode)
  void gradient_check();

  //
  //- Heading: Data
  //

  /// array of basis types for each one-dimensional orthogonal polynomial:
  /// HERMITE_ORTHOG, LEGENDRE_ORTHOG, LAGUERRE_ORTHOG, JACOBI_ORTHOG,
  /// GEN_LAGUERRE_ORTHOG, CHEBYSHEV_ORTHOG, or NUM_GEN_ORTHOG
  ShortArray orthogPolyTypes;

  /// array of one-dimensional basis polynomial objects which are used in
  /// constructing the multivariate orthogonal/interpolation polynomials
  std::vector<BasisPolynomial> polynomialBasis;

  /// order of orthogonal polynomial expansion
  UShortArray approxOrder;
  /// previous value of approxOrder; used for detecting when a multiIndex
  /// update is needed
  UShortArray approxOrderPrev;

  /// number of exp terms-by-number of vars array for identifying the orders
  /// of the one-dimensional orthogonal polynomials contributing to each
  /// of the multivariate orthogonal polynomials
  UShort2DArray multiIndex;
  /// multi-index that is the result of expansion combination
  UShort2DArray combinedMultiIndex;

  /// array of stored approxOrder's cached in store_coefficients() for use in
  /// combine_coefficients()
  UShort2DArray storedApproxOrder;
  /// array of stored multiIndex's cached in store_coefficients() for use in
  /// combine_coefficients()
  UShort3DArray storedMultiIndex;
  /// mapping of terms when aggregating storedMultiIndex with multiIndex in
  /// pre_combine_data()
  Sizet2DArray storedMultiIndexMap;

  /// numSmolyakIndices-by-numTensorProductPts-by-numVars array for
  /// identifying the orders of the one-dimensional orthogonal polynomials
  /// contributing to each of the multivariate orthogonal polynomials.
  /** For nested rules (GP, CC, or GK), the integration driver's collocKey
      is insufficient and we must track expansion orders separately. */
  UShort3DArray tpMultiIndex;
  /// sparse grid bookkeeping: mapping from num tensor-products by 
  /// tensor-product multi-indices into aggregated multiIndex
  Sizet2DArray tpMultiIndexMap;
  /// sparse grid bookkeeping: reference points for tpMultiIndexMap
  SizetArray tpMultiIndexMapRef;

  /// popped instances of tpMultiIndex that were computed but not selected
  std::deque<UShort2DArray> poppedTPMultiIndex;
  /// popped instances of tpMultiIndexMap that were computed but not selected
  std::deque<SizetArray> poppedTPMultiIndexMap;
  /// popped instances of tpMultiIndexMapRef that were computed but not selected
  std::deque<size_t> poppedTPMultiIndexMapRef;

  /// index into popped sets of data to be restored (stored in this
  /// class for used by each ProjectOrthogPolyApproximation)
  size_t pushIndex;

  /// Data vector for storing the gradients of individual expansion term
  /// polynomials (see multivariate_polynomial_gradient_vector())
  RealVector mvpGradient;
  /// Data matrix for storing the Hessians of individual expansion term
  /// polynomials (see multivariate_polynomial_hessian_matrix())
  RealSymMatrix mvpHessian;

private:

  //
  //- Heading: Member functions
  //

};


inline SharedOrthogPolyApproxData::
SharedOrthogPolyApproxData(short basis_type, const UShortArray& approx_order,
			   size_t num_vars):
  SharedPolyApproxData(basis_type, num_vars), approxOrder(approx_order)
{ }


inline SharedOrthogPolyApproxData::
SharedOrthogPolyApproxData(short basis_type, const UShortArray& approx_order,
			   size_t num_vars,
			   const ExpansionConfigOptions& ec_options,
			   const BasisConfigOptions&     bc_options):
  SharedPolyApproxData(basis_type, num_vars, ec_options, bc_options),
  approxOrder(approx_order)
{ }


inline SharedOrthogPolyApproxData::~SharedOrthogPolyApproxData()
{ }


inline const UShort2DArray& SharedOrthogPolyApproxData::multi_index() const
{ return multiIndex; }


inline size_t SharedOrthogPolyApproxData::expansion_terms() const
{ return multiIndex.size(); }


inline const UShortArray& SharedOrthogPolyApproxData::expansion_order() const
{ return approxOrder; }


inline void SharedOrthogPolyApproxData::
expansion_order(const UShortArray& order)
{ approxOrder = order; } // multiIndex updated in allocate_arrays()


inline void SharedOrthogPolyApproxData::increment_order()
{
  // increment approxOrder (multiIndex updated in allocate_arrays())
  for (size_t i=0; i<numVars; ++i)
    ++approxOrder[i];
}


/** This function is invoked to create orthogPolyTypes and polynomialBasis
    for cases where they have not already been created by an
    IntegrationDriver (i.e., expansion_samples or regression). */
inline void SharedOrthogPolyApproxData::
construct_basis(const ShortArray& u_types, const AleatoryDistParams& adp)
{
  ShortArray colloc_rules;
  BasisConfigOptions bc_options;
  construct_basis(u_types, adp, basisConfigOptions, polynomialBasis,
		  orthogPolyTypes, colloc_rules);		  
}


/** This function is invoked to create orthogPolyTypes and polynomialBasis
    for cases where they have not already been created by an
    IntegrationDriver (i.e., expansion_samples or regression). */
inline void SharedOrthogPolyApproxData::
construct_basis(const ShortArray& u_types, const AleatoryDistParams& adp,
		const BasisConfigOptions& bc_options,
		std::vector<BasisPolynomial>& poly_basis)
{
  ShortArray basis_types, colloc_rules;
  construct_basis(u_types, adp, bc_options, poly_basis,
		  basis_types, colloc_rules);
}

/** This function is invoked to create orthogPolyTypes and polynomialBasis
    for cases where they have not already been created by an
    IntegrationDriver (i.e., expansion_samples or regression). */
inline void SharedOrthogPolyApproxData::
construct_basis(const ShortArray& u_types, const AleatoryDistParams& adp,
		const BasisConfigOptions& bc_options,
		std::vector<BasisPolynomial>& poly_basis,
		ShortArray &basis_types, ShortArray &colloc_rules)
{
  bool dist_params
    = initialize_orthogonal_basis_types_rules(u_types, bc_options,
					      basis_types, colloc_rules);
  initialize_polynomial_basis(basis_types, colloc_rules, poly_basis);
  if (dist_params)
    update_basis_distribution_parameters(u_types, adp, poly_basis);
}


inline void SharedOrthogPolyApproxData::
orthogonal_basis_types(const ShortArray& opb_types)
{ orthogPolyTypes = opb_types; }


inline const ShortArray& SharedOrthogPolyApproxData::
orthogonal_basis_types() const
{ return orthogPolyTypes; }


inline const std::vector<BasisPolynomial>& SharedOrthogPolyApproxData::
polynomial_basis() const
{ return polynomialBasis; }


inline std::vector<BasisPolynomial>& SharedOrthogPolyApproxData::
polynomial_basis()
{ return polynomialBasis; }


inline void SharedOrthogPolyApproxData::
polynomial_basis(const std::vector<BasisPolynomial>& poly_basis)
{
  polynomialBasis = poly_basis;
  size_t i, num_vars = poly_basis.size();
  orthogPolyTypes.resize(num_vars);
  for (i=0; i<num_vars; ++i)
    orthogPolyTypes[i] = poly_basis[i].basis_type();
}


inline void SharedOrthogPolyApproxData::coefficients_norms_flag(bool flag)
{
  coefficients_norms_flag(flag,orthogPolyTypes,polynomialBasis);
}

inline void SharedOrthogPolyApproxData::coefficients_norms_flag(bool flag,
		  ShortArray &poly_types,
		  std::vector<BasisPolynomial> &poly_basis)
{
  //size_t i, num_basis = orthogPolyTypes.size();
  size_t i, num_basis = poly_basis.size();
  for (i=0; i<num_basis; ++i)
    if (poly_types[i] == NUM_GEN_ORTHOG)
      ((NumericGenOrthogPolynomial*)poly_basis[i].polynomial_rep())
	->coefficients_norms_flag(flag);
}


inline void SharedOrthogPolyApproxData::
push_trial_set(const UShortArray& trial_set, UShort2DArray& aggregated_mi,
		  bool monotonic, bool save_map)
{
  pre_push_trial_set(trial_set, aggregated_mi, monotonic);
  post_push_trial_set(trial_set, aggregated_mi, save_map);
}


inline Real SharedOrthogPolyApproxData::
multivariate_polynomial(const RealVector& x, const UShortArray& indices)
{ return multivariate_polynomial(x, indices, polynomialBasis); }


inline Real SharedOrthogPolyApproxData::
multivariate_polynomial(const RealVector& x, const UShortArray& indices, 
			std::vector<BasisPolynomial> &polynomial_basis )
{
  Real mvp = 1.; unsigned short order_1d; int num_vars = x.length();
  for (size_t i=0; i<num_vars; ++i) {
    order_1d = indices[i];
    if (order_1d)
      mvp *= polynomial_basis[i].type1_value(x[i], order_1d);
  }
  return mvp;
}


/** All variables version. */
inline Real SharedOrthogPolyApproxData::
multivariate_polynomial(const RealVector& x, const UShortArray& indices,
			const SizetList& non_rand_indices)
{
  return multivariate_polynomial(x, indices, non_rand_indices, polynomialBasis);
}


/** All variables version. */
inline Real SharedOrthogPolyApproxData::
multivariate_polynomial(const RealVector& x, const UShortArray& indices,
			const SizetList& non_rand_indices,
			std::vector<BasisPolynomial> &polynomial_basis)
{
  Real mvp = 1.; SizetList::const_iterator cit;
  unsigned short order_1d; size_t i;
  for (cit=non_rand_indices.begin(); cit!=non_rand_indices.end(); ++cit) {
    i = *cit; order_1d = indices[i];
    if (order_1d)
      mvp *= polynomial_basis[i].type1_value(x[i], order_1d);
  }
  return mvp;
}


inline Real SharedOrthogPolyApproxData::
multivariate_polynomial_gradient(const RealVector& x, size_t deriv_index,
				 const UShortArray& indices)
{
  Real mvp_grad = 1.; 
  for (size_t k=0; k<numVars; ++k)
    mvp_grad *= (k == deriv_index) ?
      polynomialBasis[k].type1_gradient(x[k], indices[k]) :
      polynomialBasis[k].type1_value(x[k],    indices[k]);
  return mvp_grad;
}


/** All variables version. */
inline Real SharedOrthogPolyApproxData::
multivariate_polynomial_gradient(const RealVector& x, size_t deriv_index,
				 const UShortArray& indices,
				 const SizetList& non_rand_indices)
{
  Real mvp_grad = 1.; SizetList::const_iterator cit; size_t k;
  for (cit=non_rand_indices.begin(); cit!=non_rand_indices.end(); ++cit) {
    k = *cit;
    mvp_grad *= (k == deriv_index) ?
      polynomialBasis[k].type1_gradient(x[k], indices[k]) :
      polynomialBasis[k].type1_value(x[k],    indices[k]);
  }
  return mvp_grad;
}


inline const RealVector& SharedOrthogPolyApproxData::
multivariate_polynomial_gradient_vector(const RealVector& x,
					const UShortArray& indices)
{
  if (mvpGradient.length() != numVars)
    mvpGradient.sizeUninitialized(numVars);
  for (size_t i=0; i<numVars; ++i)
    mvpGradient[i] = multivariate_polynomial_gradient(x, i, indices);
  return mvpGradient;
}


inline const RealVector& SharedOrthogPolyApproxData::
multivariate_polynomial_gradient_vector(const RealVector& x,
					const UShortArray& indices,
					const SizetArray& dvv)
{
  size_t i, j, deriv_index, num_deriv_vars = dvv.size();
  if (mvpGradient.length() != num_deriv_vars)
    mvpGradient.sizeUninitialized(num_deriv_vars);
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // *** requires an "All" view
    mvpGradient[i] = multivariate_polynomial_gradient(x, deriv_index, indices);
  }
  return mvpGradient;
}


inline Real SharedOrthogPolyApproxData::
multivariate_polynomial_hessian(const RealVector& x, size_t deriv_index_i,
				size_t deriv_index_j,
				const UShortArray& indices)
{
  Real mvp_hess = 1.;
  for (size_t k=0; k<numVars; ++k)
    if (k == deriv_index_i)
      mvp_hess *= (k == deriv_index_j) ?
	polynomialBasis[k].type1_hessian(x[k],  indices[k]) :
	polynomialBasis[k].type1_gradient(x[k], indices[k]);
    else
      mvp_hess *= (k == deriv_index_j) ?
	polynomialBasis[k].type1_gradient(x[k], indices[k]) :
	polynomialBasis[k].type1_value(x[k],    indices[k]);
  return mvp_hess;
}


inline const RealSymMatrix& SharedOrthogPolyApproxData::
multivariate_polynomial_hessian_matrix(const RealVector& x,
				       const UShortArray& indices)
{
  if (mvpHessian.numRows() != numVars)
    mvpHessian.shapeUninitialized(numVars);
  size_t r, c;
  for (r=0; r<numVars; ++r)
    for (c=0; c<=r; ++c) // lower triangle
      mvpHessian(r,c) = multivariate_polynomial_hessian(x, r, c, indices);
  // no operator[] provided for SymMatrix:
  //for (c=0; c<numVars; ++c) {
  //  Real* mvp_hess_col = mvpHessian[c];
  //  for (r=c; r<numVars; ++r) // lower triangle
  //    mvp_hess_col[r] = multivariate_polynomial_hessian(x, r, c, indices);
  //}
  return mvpHessian;
}


inline void SharedOrthogPolyApproxData::
update_pareto_set(const UShort2DArray& multi_index,
		  UShort2DArray& combined_pareto)
{
  // This function can be used for update or for initial definition:
  // > supports the case where the incoming array is a multi-index with
  //   dominated terms --> performs conversion to a non-dominated frontier.

  size_t i, num_p = multi_index.size();
  for (i=0; i<num_p; ++i)
    update_pareto_set(multi_index[i], combined_pareto);
}


inline void SharedOrthogPolyApproxData::
update_frontier(const UShort2DArray& multi_index, UShortArraySet& mi_frontier)
{
  // This function can be used for update or for initial definition:
  // > supports the case where the incoming array is a multi-index with
  //   dominated terms --> performs conversion to a non-dominated frontier.

  size_t i, num_p = multi_index.size();
  for (i=0; i<num_p; ++i)
    update_frontier(multi_index[i], mi_frontier);
}


inline void SharedOrthogPolyApproxData::
update_frontier(const UShortArraySet& multi_index, UShortArraySet& mi_frontier)
{
  // This function can be used for update or for initial definition:
  // > supports the case where the incoming array is a multi-index with
  //   dominated terms --> performs conversion to a non-dominated frontier.

  UShortArraySet::const_iterator cit;
  for (cit=multi_index.begin(); cit!=multi_index.end(); ++cit)
    update_frontier(*cit, mi_frontier);
}


inline bool SharedOrthogPolyApproxData::
zero_random(const UShortArray& indices) const
{
  SizetList::const_iterator cit;
  for (cit=randomIndices.begin(); cit!=randomIndices.end(); ++cit)
    if (indices[*cit])
      return false;
  return true;
}


inline Real SharedOrthogPolyApproxData::norm_squared(const UShortArray& indices)
{
  // the norm squared of a particular multivariate polynomial is the product of
  // the norms squared of the numVars univariate polynomials that comprise it.
  Real norm_sq = 1.; unsigned short order_1d;
  for (size_t i=0; i<numVars; ++i) {
    order_1d = indices[i];
    if (order_1d)
      norm_sq *= polynomialBasis[i].norm_squared(order_1d);
  }
  return norm_sq;
}


/** All variables version. */
inline Real SharedOrthogPolyApproxData::
norm_squared(const UShortArray& indices, const SizetList& rand_indices)
{
  // the norm squared of a particular multivariate polynomial is the product of
  // the norms squared of the numVars univariate polynomials that comprise it.
  Real norm_sq = 1.; SizetList::const_iterator cit;
  unsigned short order_1d; size_t i;
  for (cit=rand_indices.begin(); cit!=rand_indices.end(); ++cit) {
    i = *cit; order_1d = indices[i];
    if (order_1d)
      norm_sq *= polynomialBasis[i].norm_squared(order_1d);
  }
  return norm_sq;
}

} // namespace Pecos

#endif
