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
#include "IncrementalSparseGridDriver.hpp"

namespace Pecos {

class MultivariateDistribution;


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
  //- Heading: Virtual function redefinitions
  //

  bool push_available();
  void construct_basis(const MultivariateDistribution& u_dist);
  void update_basis_distribution_parameters(
    const MultivariateDistribution& u_dist);

  /// get polynomialBasis (const)
  const std::vector<BasisPolynomial>& polynomial_basis() const;
  /// get polynomialBasis
  std::vector<BasisPolynomial>& polynomial_basis();
  /// set polynomialBasis
  void polynomial_basis(const std::vector<BasisPolynomial>& poly_basis);

  //
  //- Heading: Member functions
  //

  /// alternative form for setting multiIndex (expansion import) and
  /// allocating related arrays
  void allocate_data(const UShort2DArray& mi);

  /// get active multiIndex
  const UShort2DArray& multi_index() const;
  /// get active multiIndex
  UShort2DArray& multi_index();
  /// get multiIndex[key]
  const UShort2DArray& multi_index(const ActiveKey& key) const;

  /// return multiIndex (all keys)
  const std::map<ActiveKey, UShort2DArray>& multi_index_map();

  /// retrieve size of multiIndex
  size_t expansion_terms() const;

  /// get active approxOrder
  const UShortArray& expansion_order() const;
  /// get approxOrder[key]; fn name avoids ambiguity with set fn below
  const UShortArray& keyed_expansion_order(const ActiveKey& key) const;
  /// set active approxOrder
  void expansion_order(const UShortArray& order);
  /// set active approxOrder
  void expansion_order(unsigned short new_order, bool one_sided = false);

  /// uniformly increment active approxOrder
  void increment_order();
  /// uniformly decrement active approxOrder
  void decrement_order();

  /// invoke initialize_orthogonal_basis_types_rules() and
  /// initialize_polynomial_basis()
  static void construct_basis(const MultivariateDistribution& u_dist,
			      const BasisConfigOptions& bc_options,
			      std::vector<BasisPolynomial>& poly_basis,
			      ShortArray& basis_types,ShortArray& colloc_rules);

  /// set orthogPolyTypes
  void orthogonal_basis_types(const ShortArray& opb_types);
  /// get orthogPolyTypes
  const ShortArray& orthogonal_basis_types() const;

  /// set NumericGenOrthogPolynomial::coeffsNormsFlag
  void coefficients_norms_flag(bool flag);

  /// set NumericGenOrthogPolynomial::coeffsNormsFlag
  static void coefficients_norms_flag(bool flag,
				      std::vector<BasisPolynomial>& poly_basis);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void active_key(const ActiveKey& key);
  void clear_keys();

  void allocate_data();

  void pre_combine_data();
  //void post_combine_data();
  void combined_to_active(bool clear_combined = true);

  void clear_inactive_data();

  //
  //- Heading: Member functions
  //

  /// update {multiIndex,approxOrd}Iter from activeKey
  void update_active_iterators();

  /// detect whether current expansion settings are the most refined
  const ActiveKey& maximal_expansion();
  // swap current shared data with a stored data set, as identified by index
  //void swap_shared_data(size_t index);

  /// convert a sparse grid index set and a growth setting to an integrand_order
  void sparse_grid_level_to_expansion_order(
    CombinedSparseGridDriver& csg_driver, const UShortArray& levels,
    UShortArray& exp_order);//, short growth_rate = UNRESTRICTED_GROWTH);
  /// convert quadrature orders to integrand orders using rigorous mappings
  void quadrature_order_to_integrand_order(IntegrationDriver& int_driver,
					   const UShortArray& quad_order,
					   UShortArray& int_order);
  /// convert integrand orders to expansion orders using rigorous mappings
  void integrand_order_to_expansion_order(const UShortArray& int_order,
					  UShortArray& exp_order);

  /// helper function for incrementing that is modular on sparse grid driver
  /// and multi-index
  void increment_trial_set(CombinedSparseGridDriver& csg_driver,
			   UShort2DArray& aggregated_mi);
  /// helper function for decrementing that is modular on trial set
  /// and multi-index
  void decrement_trial_set(const UShortArray& trial_set,
			   UShort2DArray& aggregated_mi, bool save_map = true);
  /// helper function for restoring that is modular on trial set and multi-index
  void pre_push_trial_set(const UShortArray& trial_set,
			  UShort2DArray& aggregated_mi, bool monotonic = true);
  /// helper function for restoring that is modular on trial set and multi-index
  void post_push_trial_set(const UShortArray& trial_set,
			   UShort2DArray& aggregated_mi, bool save_map = true);
  /// helper function for restoring that is modular on trial set and multi-index
  void push_trial_set(const UShortArray& trial_set,
		      UShort2DArray& aggregated_mi, bool monotonic = true,
		      bool save_map = true);

  /// Precompute a maximal order of quadrature rules (based on a
  /// multiIndex) when too expensive to compute on demand
  void precompute_maximal_rules(const UShort2DArray& multi_index);
  /// Precompute a maximal order of quadrature rules (based on an
  /// approxOrder) when too expensive to compute on demand
  void precompute_maximal_rules(const UShortArray& approx_order);

  /// allocate sobolIndexMap from active multi_index
  void allocate_component_sobol();
  /// allocate sobolIndexMap from multi_index
  void allocate_component_sobol(const UShort2DArray& multi_index);
  /// update sobolIndexMap using new multi_index terms (from multifidelity
  /// overlay or a new QoI in orthogonal least interpolation)
  void update_component_sobol(const UShort2DArray& multi_index);

  /// define multi_index_c based on products of terms contained within
  /// multi_index_a and multi_index_b
  void product_multi_index(const UShort2DArray& multi_index_a,
			   const UShort2DArray& multi_index_b,
			   UShort2DArray& multi_index_c);

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
  std::map<ActiveKey, UShortArray> approxOrder;
  /// iterator pointing to active node in approxOrder
  std::map<ActiveKey, UShortArray>::iterator approxOrdIter;

  /// initial value of approxOrder passed through ctor
  UShortArray approxOrderSpec;
  /// previous value of active approxOrder; used for detecting when an
  /// expansion update is needed in allocate_arrays()
  UShortArray prevApproxOrder;
  /// previous value of active multiIndex for restoration in decrement_data()
  UShort2DArray prevMultiIndex;
  /// previous value of activeKey; used for detecting when an
  /// expansion update is needed in allocate_arrays()
  ActiveKey prevActiveKey;

  /// number of exp terms-by-number of vars array for identifying the orders
  /// of the one-dimensional orthogonal polynomials contributing to each
  /// of the multivariate orthogonal polynomials
  std::map<ActiveKey, UShort2DArray> multiIndex;
  /// iterator pointing to active node in multiIndex
  std::map<ActiveKey, UShort2DArray>::iterator multiIndexIter;

  /// multi-index that is the final result of a sequence of expansion
  /// combinations
  UShort2DArray combinedMultiIndex;
  /// mapping of terms when aggregating multiIndex into combinedMultiIndex
  /// in pre_combine_data() (case ADD_COMBINE)
  Sizet2DArray combinedMultiIndexMap;
  /// sequence of multi-index products defined in pre_combine_data() for case
  /// MULT_COMBINE.  For combinations of more than two levels, this provides
  /// a bridge between the first multi-index for "a" (multiIndex.begin()) and
  /// the final multi-index for "c" (combinedMultiIndex) in c = a * b.
  UShort3DArray combinedMultiIndexSeq;

  /// numSmolyakIndices-by-numTensorProductPts-by-numVars array for
  /// identifying the orders of the one-dimensional orthogonal polynomials
  /// contributing to each of the multivariate orthogonal polynomials.
  /** For nested rules (GP, CC, or GK), the integration driver's collocKey
      is insufficient and we must track expansion orders separately. */
  std::map<ActiveKey, UShort3DArray> tpMultiIndex;
  /// sparse grid bookkeeping: mapping from num tensor-products by 
  /// tensor-product multi-indices into aggregated multiIndex
  std::map<ActiveKey, Sizet2DArray> tpMultiIndexMap;
  /// sparse grid bookkeeping: reference points for tpMultiIndexMap
  std::map<ActiveKey, SizetArray> tpMultiIndexMapRef;

  /// popped instances of either multiIndex or tpMultiIndex (depending
  /// on expansion solution approach) that were computed but not selected
  std::map<ActiveKey, UShort2DArrayDeque> poppedMultiIndex;
  /// popped instances of tpMultiIndexMap that were computed but not selected
  std::map<ActiveKey, SizetArrayDeque> poppedMultiIndexMap;
  /// popped instances of tpMultiIndexMapRef that were computed but not selected
  std::map<ActiveKey, SizetDeque> poppedMultiIndexMapRef;

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
  SharedPolyApproxData(basis_type, num_vars), approxOrdIter(approxOrder.end()),
  approxOrderSpec(approx_order)
{
  update_active_iterators();
  approxOrdIter->second = approx_order;
}


inline SharedOrthogPolyApproxData::
SharedOrthogPolyApproxData(short basis_type, const UShortArray& approx_order,
			   size_t num_vars,
			   const ExpansionConfigOptions& ec_options,
			   const BasisConfigOptions&     bc_options):
  SharedPolyApproxData(basis_type, num_vars, ec_options, bc_options),
  approxOrdIter(approxOrder.end()), approxOrderSpec(approx_order)
{
  update_active_iterators();
  approxOrdIter->second = approx_order;
}


inline SharedOrthogPolyApproxData::~SharedOrthogPolyApproxData()
{ }


inline void SharedOrthogPolyApproxData::update_active_iterators()
{
  if (approxOrdIter != approxOrder.end() && approxOrdIter->first == activeKey)
    return;

  approxOrdIter  = approxOrder.find(activeKey);
  multiIndexIter =  multiIndex.find(activeKey);

  /* So long as we only create new keys and avoid modifying existing ones,
     this deep copy is not needed.
  ActiveKey active_copy; // share 1 deep copy of current active key
  if (approxOrdIter == approxOrder.end() || multiIndexIter == multiIndex.end())
    active_copy = activeKey.copy();
  */

  if (approxOrdIter == approxOrder.end()) {
    std::pair<ActiveKey, UShortArray>
      ua_pair(activeKey/*active_copy*/, approxOrderSpec);//, UShortArray());
    approxOrdIter = approxOrder.insert(ua_pair).first;
  }
  if (multiIndexIter == multiIndex.end()) {
    std::pair<ActiveKey, UShort2DArray>
      u2a_pair(activeKey/*active_copy*/, UShort2DArray());
    multiIndexIter = multiIndex.insert(u2a_pair).first;
    //updateExpForm = true; // multiIndex to be updated in allocate_arrays()
  }
}


inline const UShort2DArray& SharedOrthogPolyApproxData::multi_index() const
{ return multiIndexIter->second; }


inline UShort2DArray& SharedOrthogPolyApproxData::multi_index()
{ return multiIndexIter->second; }


inline const UShort2DArray& SharedOrthogPolyApproxData::
multi_index(const ActiveKey& key) const
{
  std::map<ActiveKey, UShort2DArray>::const_iterator cit
    = multiIndex.find(key);
  if (cit == multiIndex.end()) {
    PCerr << "Error: key not found in SharedOrthogPolyApproxData::"
	  << "multi_index()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const std::map<ActiveKey, UShort2DArray>& SharedOrthogPolyApproxData::
multi_index_map()
{ return multiIndex; }


inline size_t SharedOrthogPolyApproxData::expansion_terms() const
{ return multiIndexIter->second.size(); }


inline const UShortArray& SharedOrthogPolyApproxData::expansion_order() const
{ return approxOrdIter->second; }


/** Use alternate naming since UShortArray overload already used. */
inline const UShortArray& SharedOrthogPolyApproxData::
keyed_expansion_order(const ActiveKey& key) const
{
  std::map<ActiveKey, UShortArray>::const_iterator cit
    = approxOrder.find(key);
  if (cit == approxOrder.end()) {
    PCerr << "Error: key not found in SharedOrthogPolyApproxData::"
	  << "keyed_expansion_order()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline void SharedOrthogPolyApproxData::
expansion_order(const UShortArray& order)
{
  UShortArray& approx_order = approxOrdIter->second;
  if (approx_order != order) {
    approx_order = order;
    //updateExpForm = true; // multiIndex to be updated in allocate_arrays()
  }
}


inline void SharedOrthogPolyApproxData::
expansion_order(unsigned short new_order, bool one_sided)
{
  UShortArray& approx_order = approxOrdIter->second;
  if (approx_order.empty() || !one_sided)
    approx_order.assign(numVars, new_order);
  else {
    size_t i, num_ao = approx_order.size();
    for (i=0; i<num_ao; ++i)
      if (new_order > approx_order[i])
	approx_order[i] = new_order;
  }
}


inline void SharedOrthogPolyApproxData::allocate_component_sobol()
{ allocate_component_sobol(multiIndexIter->second); }


inline void SharedOrthogPolyApproxData::increment_order()
{
  // increment approxOrder (multiIndex updated in allocate_arrays())
  UShortArray& approx_order = approxOrdIter->second;
  for (size_t i=0; i<numVars; ++i)
    ++approx_order[i];
  //updateExpForm = true; // multiIndex to be updated in allocate_arrays()
}


inline void SharedOrthogPolyApproxData::decrement_order()
{
  // decrement approxOrder (multiIndex updated in allocate_arrays())
  UShortArray& approx_order = approxOrdIter->second;
  for (size_t i=0; i<numVars; ++i)
    --approx_order[i];
  //updateExpForm = true; // multiIndex to be updated in allocate_arrays()
}


/** This function is invoked to create orthogPolyTypes and polynomialBasis
    for cases where they have not already been created by an
    IntegrationDriver (i.e., expansion_samples or regression). */
inline void SharedOrthogPolyApproxData::
construct_basis(const MultivariateDistribution& u_dist,
		const BasisConfigOptions& bc_options,
		std::vector<BasisPolynomial>& poly_basis,
		ShortArray& basis_types, ShortArray& colloc_rules)
{
  // Construct time initializations
  initialize_orthogonal_basis_types_rules(u_dist, bc_options,
					  basis_types, colloc_rules);
  initialize_polynomial_basis(basis_types, colloc_rules, poly_basis);

  // The following update now occurs at run time:
  //update_basis_distribution_parameters(u_dist, poly_basis);
}


inline void SharedOrthogPolyApproxData::
construct_basis(const MultivariateDistribution& u_dist)
{
  ShortArray colloc_rules;
  construct_basis(u_dist, basisConfigOptions, polynomialBasis,
		  orthogPolyTypes, colloc_rules);
}


inline void SharedOrthogPolyApproxData::
update_basis_distribution_parameters(const MultivariateDistribution& u_dist)
{
  SharedPolyApproxData::
    update_basis_distribution_parameters(u_dist, polynomialBasis);
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
{ coefficients_norms_flag(flag, polynomialBasis); }


inline void SharedOrthogPolyApproxData::
coefficients_norms_flag(bool flag, std::vector<BasisPolynomial>& poly_basis)
{
  //size_t i, num_basis = orthogPolyTypes.size();
  size_t i, num_basis = poly_basis.size();
  for (i=0; i<num_basis; ++i) {
    BasisPolynomial& poly_basis_i = poly_basis[i];
    if (poly_basis_i.basis_type() == NUM_GEN_ORTHOG)
      std::static_pointer_cast<NumericGenOrthogPolynomial>
	(poly_basis_i.polynomial_rep())->coefficients_norms_flag(flag);
  }
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
