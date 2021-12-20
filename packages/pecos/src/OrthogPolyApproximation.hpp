/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        OrthogPolyApproximation
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        Mike Eldred

#ifndef ORTHOG_POLY_APPROXIMATION_HPP
#define ORTHOG_POLY_APPROXIMATION_HPP

#include "PolynomialApproximation.hpp"
#include "NumericGenOrthogPolynomial.hpp"
#include "SharedOrthogPolyApproxData.hpp"

namespace Pecos {


/// Derived approximation class for orthogonal polynomials (global
/// approximation).

/** The OrthogPolyApproximation class provides a global approximation
    based on orthogonal polynomials.  It is used primarily for polynomial
    chaos expansions (for stochastic finite element approaches to
    uncertainty quantification). */

class OrthogPolyApproximation: public PolynomialApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  OrthogPolyApproximation(const SharedBasisApproxData& shared_data);
  /// destructor
  ~OrthogPolyApproximation();

  //
  //- Heading: Virtual functions
  //

  /// estimate chaos expansion coefficient decay rates for each random
  /// variable dimension using linear least squares in semilog space
  virtual const RealVector& dimension_decay_rates();

  /// unscale the expansion coefficients following computation using scaled
  /// response data
  virtual void unscale_coefficients(RealVector& exp_coeffs,
				    RealMatrix& exp_coeff_grads);

  /// evaluate all pce basis functions at a single point
  static void basis_value(const RealVector& x,
			  std::vector<BasisPolynomial> &polynomial_basis,
			  const UShort2DArray &multi_index,
			  RealVector &basis_values);
  /// evaluate all pce basis functions at a set of points
  static void basis_matrix(const RealMatrix& x,
			   std::vector<BasisPolynomial> &polynomial_basis,
			   const UShort2DArray &multi_index,
			   RealMatrix &basis_values);

  void basis_matrix(const RealMatrix& x, RealMatrix &basis_values);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  size_t expansion_terms() const;

  int min_coefficients() const;

  /*
  void store_coefficients(size_t index = _NPOS);
  void restore_coefficients(size_t index = _NPOS);
  void swap_coefficients(size_t maximal_index);
  void remove_stored_coefficients(size_t index = _NPOS);
  */
  void clear_inactive();

  void combine_coefficients();
  void combined_to_active(bool clear_combined = true);

  void print_coefficients(std::ostream& s, bool normalized);
  void print_coefficients(std::ostream& s, const UShort2DArray& mi,
			  const RealVector& exp_coeffs, bool normalized);

  /// retrieve or form a set of dense coefficients that correspond to
  /// SharedOrthogPolyApproxData::multiIndex
  RealVector approximation_coefficients(bool normalized) const;
  /// set an array of dense coefficients corresponding to
  /// SharedOrthogPolyApproxData::multiIndex
  void approximation_coefficients(const RealVector& approx_coeffs,
				  bool normalized);

  void coefficient_labels(std::vector<std::string>& all_coeff_tags) const;

  /// initialize multiIndex, expansionCoeffs, et al.
  void allocate_arrays();

  /// Performs global sensitivity analysis via variance-based decomposition;
  /// computes component (main and interaction) Sobol' indices
  void compute_component_sobol();
  /// Performs global sensitivity analysis via variance-based decomposition;
  /// computes total Sobol' indices
  void compute_total_sobol();

  Real value(const RealVector& x);
  const RealVector& gradient_basis_variables(const RealVector& x);
  const RealVector& gradient_basis_variables(const RealVector& x,
    const SizetArray& dvv);
  const RealVector& gradient_nonbasis_variables(const RealVector& x);
  const RealSymMatrix& hessian_basis_variables(const RealVector& x);

  Real stored_value(const RealVector& x, const ActiveKey& key);
  const RealVector& stored_gradient_basis_variables(const RealVector& x,
    const ActiveKey& key);
  const RealVector& stored_gradient_basis_variables(const RealVector& x,
    const SizetArray& dvv, const ActiveKey& key);
  const RealVector& stored_gradient_nonbasis_variables(const RealVector& x,
    const ActiveKey& key);
  const RealSymMatrix& stored_hessian_basis_variables(const RealVector& x,
    const ActiveKey& key);

  const RealVector& expansion_moments() const;
  const RealVector& numerical_integration_moments() const;

  Real mean();
  Real mean(const RealVector& x);
  Real covariance(PolynomialApproximation* poly_approx_2);
  Real covariance(const RealVector& x, PolynomialApproximation* poly_approx_2);

  const RealVector& mean_gradient();
  const RealVector& mean_gradient(const RealVector& x, const SizetArray& dvv);
  const RealVector& variance_gradient();
  const RealVector& variance_gradient(const RealVector& x,
				      const SizetArray& dvv);

  Real combined_mean();
  Real combined_mean(const RealVector& x);
  Real combined_covariance(PolynomialApproximation* poly_approx_2);
  Real combined_covariance(const RealVector& x,
			   PolynomialApproximation* poly_approx_2);

  //
  //- Heading: Member functions
  //

  /// update expCoeff{s,Grads}Iter for new activeKey from sharedDataRep
  bool update_active_iterators(const ActiveKey& key);

  /// size expansion{Coeffs,CoeffGrads} based on the shared multiIndex
  void size_expansion();
  /// size exp_coeff{s,_grads} for an updated number of expansion terms
  void size_expansion(size_t num_exp_terms, RealVector& exp_coeffs,
		      RealMatrix& exp_coeff_grads);
  /// resize expansion{Coeffs,CoeffGrads} for an updated multiIndex
  void resize_expansion();
  /// resize exp_coeff{s,_grads} for an updated number of expansion terms
  void resize_expansion(size_t num_exp_terms, RealVector& exp_coeffs,
			RealMatrix& exp_coeff_grads);

  /// compute the expansion value
  Real value(const RealVector& x, const UShort2DArray& mi,
    const RealVector& exp_coeffs);
  /// compute the expansion gradient with respect to the basis variables
  const RealVector& gradient_basis_variables(const RealVector& x,
    const UShort2DArray& mi, const RealVector& exp_coeffs);
  /// compute the expansion gradient with respect to the basis variables for
  /// gioven DVV components
  const RealVector& gradient_basis_variables(const RealVector& x,
    const SizetArray& dvv, const UShort2DArray& mi,
    const RealVector& exp_coeffs);
  /// compute the expansion gradient with respect to auxilliary
  /// non-basis variables
  const RealVector& gradient_nonbasis_variables(const RealVector& x,
    const UShort2DArray& mi, const RealMatrix& exp_coeff_grads);
  /// compute the expansion Hessian with respect to the basis variables
  const RealSymMatrix& hessian_basis_variables(const RealVector& x,
    const UShort2DArray& mi, const RealVector& exp_coeffs);

  /// overlay the passed expansion with the aggregate
  /// expansion{Coeffs,CoeffGrads} as managed by the multi_index_map
  void overlay_expansion(const SizetArray& multi_index_map,
			 const RealVector& exp_coeffs,
			 const RealMatrix& exp_grads, int coeff,
			 RealVector& exp_coeffs_sum, RealMatrix& exp_grads_sum);
  /// multiply current expansion ("a") with incoming expansion ("b")
  /// and store in product expansion ("c")
  void multiply_expansion(const UShort2DArray& multi_index_a,
			  const RealVector&    exp_coeffs_a,
			  const RealMatrix&    exp_grads_a,
			  const UShort2DArray& multi_index_b,
			  const RealVector&    exp_coeffs_b,
			  const RealMatrix&    exp_grads_b,
			  const UShort2DArray& multi_index_c,
			  RealVector& exp_coeffs_c, RealMatrix& exp_grads_c);

  /// update add_val and add_gradient based on failure map from SurrogateData
  void fail_booleans(SizetShortMap::const_iterator& fit, size_t j,
		     bool& add_val, bool& add_grad);

  /// utility function for solving the least squares estimation of decay rates
  void solve_decay_rates(RealVectorArray& A_vectors, RealVectorArray& b_vectors,
			 UShortArray& max_orders);

  //
  //- Heading: Data
  //

  /// the coefficients of the expansion
  std::map<ActiveKey, RealVector> expansionCoeffs;
  /// iterator pointing to active node in expansionCoeffs
  std::map<ActiveKey, RealVector>::iterator expCoeffsIter;

  /// the gradients of the expansion coefficients
  /** may be interpreted as either the gradients of the expansion
      coefficients or the coefficients of expansions for the response
      gradients.  This array is used when sensitivities of moments are
      needed with respect to variables that do not appear in the
      expansion (e.g., with respect to design or epistemic variables
      for an expansion only over probabilistic variables). */
  std::map<ActiveKey, RealMatrix> expansionCoeffGrads;
  /// iterator pointing to active node in expansionCoeffGrads
  std::map<ActiveKey, RealMatrix>::iterator expCoeffGradsIter;

  /*
  /// copies of expansionCoeffs stored in store_coefficients() for use
  /// in combine_coefficients()
  RealVectorArray storedExpCoeffs;
  /// copies of expansionCoeffGrads stored in store_coefficients() for
  /// use in combine_coefficients()
  RealMatrixArray storedExpCoeffGrads;
  */

  /// roll up of expansion coefficients across all keys
  RealVector combinedExpCoeffs;
  /// roll up of expansion coefficient gradients across all keys
  RealMatrix combinedExpCoeffGrads;

  /// spectral coefficient decay rates estimated by LLS on log of
  /// univariate expansion coefficients
  RealVector decayRates;

private:

  //
  //- Heading: Member functions
  //

  /// compute covariance between two sets of expansion coefficients
  /// (standard variables mode)
  Real covariance(const UShort2DArray& mi, const RealVector& exp_coeffs,
		  const RealVector& exp_coeffs_2);
  /// compute covariance between two sets of expansion coefficients
  /// (combined variables mode)
  Real covariance(const RealVector& x, const UShort2DArray& mi,
		  const RealVector& exp_coeffs, const RealVector& exp_coeffs_2);

  // apply normalization to std_coeffs to create normalized_coeffs
  void normalize(const RealVector& std_coeffs,
		 RealVector& normalized_coeffs) const;
  // remove normalization from normalized_coeffs to create std_coeffs
  void denormalize(const RealVector& normalized_coeffs,
		   RealVector& std_coeffs) const;

  //
  //- Heading: Data
  //

};


inline OrthogPolyApproximation::
OrthogPolyApproximation(const SharedBasisApproxData& shared_data):
  PolynomialApproximation(shared_data), expCoeffsIter(expansionCoeffs.end())
{ }


inline OrthogPolyApproximation::~OrthogPolyApproximation()
{ }


inline bool OrthogPolyApproximation::
update_active_iterators(const ActiveKey& key)
{
  // Test for change
  if (expCoeffsIter != expansionCoeffs.end() && expCoeffsIter->first == key)
    return false;

  expCoeffsIter     = expansionCoeffs.find(key);
  expCoeffGradsIter = expansionCoeffGrads.find(key);

  // share 1 deep copy of current active key
  ActiveKey key_copy;
  if (expCoeffsIter     == expansionCoeffs.end() ||
      expCoeffGradsIter == expansionCoeffGrads.end())
    key_copy = key.copy();

  if (expCoeffsIter == expansionCoeffs.end()) {
    std::pair<ActiveKey, RealVector> rv_pair(key_copy, RealVector());
    expCoeffsIter = expansionCoeffs.insert(rv_pair).first;
  }
  if (expCoeffGradsIter == expansionCoeffGrads.end()) {
    std::pair<ActiveKey, RealMatrix> rm_pair(key_copy, RealMatrix());
    expCoeffGradsIter = expansionCoeffGrads.insert(rm_pair).first;
  }

  surrData.active_key(key);

  PolynomialApproximation::update_active_iterators(key);
  return true;
}


inline const RealVector& OrthogPolyApproximation::expansion_moments() const
{ return primaryMomIter->second; }


inline const RealVector& OrthogPolyApproximation::
numerical_integration_moments() const
{ return secondaryMoments; }


/** default implementation if no sparsity (overridden in
    RegressOrthogPolyApproximation for CS) */
inline size_t OrthogPolyApproximation::expansion_terms() const
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  return data_rep->expansion_terms();
}


inline Real OrthogPolyApproximation::value(const RealVector& x)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  return value(x, data_rep->multi_index(), expCoeffsIter->second);
}


inline Real OrthogPolyApproximation::
stored_value(const RealVector& x, const ActiveKey& key)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  return value(x, data_rep->multi_index(key), expansionCoeffs[key]);
}


inline const RealVector& OrthogPolyApproximation::
gradient_basis_variables(const RealVector& x)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  return gradient_basis_variables(x, data_rep->multi_index(),
				  expCoeffsIter->second);
}


inline const RealVector& OrthogPolyApproximation::
stored_gradient_basis_variables(const RealVector& x, const ActiveKey& key)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  return gradient_basis_variables(x, data_rep->multi_index(key),
				  expansionCoeffs[key]);
}


inline const RealVector& OrthogPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  return gradient_basis_variables(x, dvv, data_rep->multi_index(),
				  expCoeffsIter->second);
}


inline const RealVector& OrthogPolyApproximation::
stored_gradient_basis_variables(const RealVector& x, const SizetArray& dvv,
				const ActiveKey& key)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  return gradient_basis_variables(x, dvv, data_rep->multi_index(key),
				  expansionCoeffs[key]);
}


inline const RealVector& OrthogPolyApproximation::
gradient_nonbasis_variables(const RealVector& x)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  return gradient_nonbasis_variables(x, data_rep->multi_index(),
				     expCoeffGradsIter->second);
}


inline const RealVector& OrthogPolyApproximation::
stored_gradient_nonbasis_variables(const RealVector& x, const ActiveKey& key)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  return gradient_nonbasis_variables(x, data_rep->multi_index(key),
				     expansionCoeffGrads[key]);
}


inline const RealSymMatrix& OrthogPolyApproximation::
hessian_basis_variables(const RealVector& x)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  return hessian_basis_variables(x, data_rep->multi_index(),
				 expCoeffsIter->second);
}


inline const RealSymMatrix& OrthogPolyApproximation::
stored_hessian_basis_variables(const RealVector& x, const ActiveKey& key)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  return hessian_basis_variables(x, data_rep->multi_index(key),
				 expansionCoeffs[key]);
}


inline void OrthogPolyApproximation::
normalize(const RealVector& std_coeffs, RealVector& normalized_coeffs) const
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  const UShort2DArray& mi = data_rep->multi_index();
  size_t i, num_mi = mi.size();
  if (normalized_coeffs.length() != num_mi)
    normalized_coeffs.sizeUninitialized(num_mi);
  // basis is divided by norm, so coeff is multiplied by norm
  for (i=0; i<num_mi; ++i)
    normalized_coeffs[i] = std_coeffs[i]
                         * std::sqrt(data_rep->norm_squared(mi[i]));
}


inline void OrthogPolyApproximation::
denormalize(const RealVector& normalized_coeffs, RealVector& std_coeffs) const
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  const UShort2DArray& mi = data_rep->multi_index();
  size_t i, num_mi = mi.size();
  if (std_coeffs.length() != num_mi)
    std_coeffs.sizeUninitialized(num_mi);
  for (i=0; i<num_mi; ++i)
    std_coeffs[i] = normalized_coeffs[i]
                  / std::sqrt(data_rep->norm_squared(mi[i]));
}


inline RealVector OrthogPolyApproximation::
approximation_coefficients(bool normalized) const
{
  RealVector& exp_coeffs = expCoeffsIter->second;
  if (normalized) {
    RealVector normalized_coeffs;
    normalize(exp_coeffs, normalized_coeffs);
    return normalized_coeffs;
  }
  else
    return RealVector(Teuchos::View, exp_coeffs.values(), exp_coeffs.length());
}


inline void OrthogPolyApproximation::
approximation_coefficients(const RealVector& approx_coeffs, bool normalized)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  update_active_iterators(data_rep->activeKey);

  if (normalized) denormalize(approx_coeffs, expCoeffsIter->second);
  else            expCoeffsIter->second = approx_coeffs;

  // allocate arrays in support of external coefficient import (mirrors
  // allocate_arrays() except for redundant size_expansion())
  allocate_total_sobol();
  allocate_component_sobol();
  RealVector& exp_mom = primaryMomIter->second;
  if (exp_mom.length() != 2) exp_mom.sizeUninitialized(2);
}


inline void OrthogPolyApproximation::
print_coefficients(std::ostream& s, bool normalized)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  print_coefficients(s, data_rep->multi_index(), expCoeffsIter->second,
		     normalized);
}


inline void OrthogPolyApproximation::size_expansion()
{
  size_expansion(expansion_terms(), expCoeffsIter->second,
		 expCoeffGradsIter->second);
}


inline void OrthogPolyApproximation::
size_expansion(size_t num_exp_terms, RealVector& exp_coeffs,
	       RealMatrix& exp_coeff_grads)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  if (expansionCoeffFlag) {
    if (exp_coeffs.length() != num_exp_terms)
      exp_coeffs.sizeUninitialized(num_exp_terms);
  }
  if (expansionCoeffGradFlag) {
    size_t num_deriv_vars = surrData.num_derivative_variables();
    if (exp_coeff_grads.numRows() != num_deriv_vars ||
	exp_coeff_grads.numCols() != num_exp_terms)
      exp_coeff_grads.shapeUninitialized(num_deriv_vars, num_exp_terms);
  }
}


inline void OrthogPolyApproximation::resize_expansion()
{
  resize_expansion(expansion_terms(), expCoeffsIter->second,
		   expCoeffGradsIter->second);
}


inline void OrthogPolyApproximation::
resize_expansion(size_t num_exp_terms, RealVector& exp_coeffs,
		 RealMatrix& exp_coeff_grads)
{
  std::shared_ptr<SharedOrthogPolyApproxData> data_rep =
    std::static_pointer_cast<SharedOrthogPolyApproxData>(sharedDataRep);
  if (expansionCoeffFlag)
    exp_coeffs.resize(num_exp_terms); // new terms initialized to 0
  if (expansionCoeffGradFlag) {
    size_t num_deriv_vars = exp_coeff_grads.numRows();
    exp_coeff_grads.reshape(num_deriv_vars, num_exp_terms);//new terms init to 0
  }
}


inline void OrthogPolyApproximation::
fail_booleans(SizetShortMap::const_iterator& fit, size_t j,
	      bool& add_val, bool& add_grad)
{
  if (fit != surrData.failed_response_data().end() && fit->first == j) {
    short fail_bits = fit->second;
    if (fail_bits & 1) add_val  = false;
    if (fail_bits & 2) add_grad = false;
    ++fit;
  }
}

} // namespace Pecos

#endif
