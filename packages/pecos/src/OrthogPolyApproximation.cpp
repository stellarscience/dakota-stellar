/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        OrthogPolyApproximation
//- Description:  Implementation code for OrthogPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "OrthogPolyApproximation.hpp"
#include "pecos_global_defs.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG
//#define DECAY_DEBUG


namespace Pecos {


int OrthogPolyApproximation::min_coefficients() const
{ return 0; } // coefficient import case


void OrthogPolyApproximation::allocate_arrays()
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  update_active_iterators(data_rep->activeKey);

  // expansion formulation has been defined in Shared*OrthogPolyApproxData::
  // allocate_data(), and its results are employed below

  allocate_total_sobol();
  allocate_component_sobol();

  // size expansion even if !update_exp_form due to possibility of change
  // to expansion{Coeff,GradFlag} settings
  size_expansion();

  RealVector& exp_mom = expMomentsIter->second;
  if (exp_mom.length() != 2) exp_mom.sizeUninitialized(2);
}


void OrthogPolyApproximation::clear_inactive()
{
  std::map<UShortArray, RealVector>::iterator ec_it = expansionCoeffs.begin();
  std::map<UShortArray, RealMatrix>::iterator eg_it
    = expansionCoeffGrads.begin();
  while (ec_it != expansionCoeffs.end())
    if (ec_it == expCoeffsIter) // preserve active
      { ++ec_it, ++eg_it; }
    else // clear inactive: postfix increments manage iterator invalidations
      { expansionCoeffs.erase(ec_it++); expansionCoeffGrads.erase(eg_it++); }
}


void OrthogPolyApproximation::combine_coefficients()
{
  // Combine the data stored previously by store_coefficients()

  // SharedOrthogPolyApproxData::pre_combine_data() appends multi-indices
  // SharedOrthogPolyApproxData::post_combine_data() finalizes multiIndex

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;

  // Coefficient combination is not dependent on active state
  //update_active_iterators(data_rep->activeKey);

  // Note: computed bits are also cleared when refineStatsType is changed
  if (data_rep->expConfigOptions.refineStatsType == COMBINED_EXPANSION_STATS)
    clear_computed_bits();

  //allocate_component_sobol(); // size sobolIndices from shared sobolIndexMap

  std::map<UShortArray, RealVector>::iterator ec_it;
  std::map<UShortArray, RealMatrix>::iterator eg_it;
  switch (data_rep->expConfigOptions.combineType) {
  case MULT_COMBINE: {
    // perform the multiplication of level expansions
    const UShort3DArray& combined_mi_seq = data_rep->combinedMultiIndexSeq;
    size_t cntr, num_seq = combined_mi_seq.size();
    const std::map<UShortArray, UShort2DArray>& mi = data_rep->multiIndex;
    std::map<UShortArray, UShort2DArray>::const_iterator mi_cit = ++mi.begin();
    ec_it = ++expansionCoeffs.begin();  eg_it = ++expansionCoeffGrads.begin();
    for (cntr=0; cntr<=num_seq; ++cntr, ++ec_it, ++eg_it, ++mi_cit) {
      const UShort2DArray& multi_index_a = (cntr) ?
	combined_mi_seq[cntr-1] : mi.begin()->second;
      const RealVector&    exp_coeffs_a  = (cntr) ?
	combinedExpCoeffs       : expansionCoeffs.begin()->second;
      const RealMatrix&    exp_grads_a   = (cntr) ?
	combinedExpCoeffGrads   : expansionCoeffGrads.begin()->second;
      const UShort2DArray& multi_index_c = (cntr < num_seq) ?
	combined_mi_seq[cntr]   : data_rep->combinedMultiIndex;
      // internal tmp copies allow a and c arrays to be the same memory
      multiply_expansion(multi_index_a, exp_coeffs_a, exp_grads_a,
			 mi_cit->second, ec_it->second, eg_it->second,
			 multi_index_c, combinedExpCoeffs,
			 combinedExpCoeffGrads);
    }
    break;
  }
  case ADD_MULT_COMBINE:
    //overlay_expansion(data_rep->combinedMultiIndexMap[i], storedExpCoeffs[i],
    //                  storedExpCoeffGrads[i], addCoeffs, addCoeffGrads);
    //multiply_expansion(data_rep->storedMultiIndex[i], storedExpCoeffs[i],
    //                   storedExpCoeffGrads[i], multCoeffs, multCoeffGrads);
    //compute_combine_factors(addCoeffs, multCoeffs);
    //apply_combine_factors();
    PCerr << "Error : additive+multiplicative combination not yet implemented "
	  << "in OrthogPolyApproximation::combine_coefficients()" << std::endl;
    abort_handler(-1);
    break;
  default: { //case ADD_COMBINE: (correction spec not required)
    // Note: would like to preserve tensor indexing (at least for QUADRATURE
    // case) so that Horner's rule performance opt could be used within
    // tensor_product_value()).  However, a tensor result in the overlay
    // will not occur unless one expansion order dominates the other (partial
    // domination results in sum of tensor expansions as for sparse grids).
    // Therefore, stick with the general-purpose expansion overlay and exclude
    // tensor_product_value() usage for combined coefficient sets.

#ifdef DEBUG
    PCout << "\ncombinedMultiIndex:\n" << data_rep->combinedMultiIndex;
#endif // DEBUG

    // resize combinedExp{Coeffs,CoeffGrads} based on combinedMultiIndex
    resize_expansion(data_rep->combinedMultiIndex.size(), combinedExpCoeffs,
		     combinedExpCoeffGrads);
    combinedExpCoeffs = 0.;  combinedExpCoeffGrads = 0.;

    // update combinedExp{Coeffs,CoeffGrads}
    const Sizet2DArray& combined_mi_map = data_rep->combinedMultiIndexMap;
    size_t i = 0;
    for (ec_it =expansionCoeffs.begin(), eg_it =expansionCoeffGrads.begin();
	 ec_it!=expansionCoeffs.end() && eg_it!=expansionCoeffGrads.end();
	 ++ec_it, ++eg_it, ++i) {
#ifdef DEBUG
      PCout << "\ni = " << i << " combinedMultiIndexMap:\n"
	    << combined_mi_map[i] << "coeffs array:\n" << ec_it->second;
#endif // DEBUG
      overlay_expansion(combined_mi_map[i], ec_it->second, eg_it->second, 1,
			combinedExpCoeffs, combinedExpCoeffGrads);
    }
    break;
  }
  }

  if (data_rep->expConfigOptions.outputLevel >= DEBUG_OUTPUT) {
    const std::map<UShortArray, UShort2DArray>& mi = data_rep->multiIndex;
    std::map<UShortArray, UShort2DArray>::const_iterator mi_cit = mi.begin();
    for (ec_it=expansionCoeffs.begin(); ec_it!=expansionCoeffs.end();
	 ++ec_it, ++mi_cit) {
      PCout << "\nLevel coefficients (unnormalized):";
      print_coefficients(PCout, mi_cit->second, ec_it->second, false);
    }
    PCout << "\nCombined coefficients (unnormalized):";
    print_coefficients(PCout, data_rep->combinedMultiIndex, combinedExpCoeffs,
		       false);
  }
}


void OrthogPolyApproximation::combined_to_active(bool clear_combined)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  update_active_iterators(data_rep->activeKey);

  // Note: swap() is conditionally available for Real{Vector,Matrix}
  if (expansionCoeffFlag) {
    if (clear_combined) {
      expCoeffsIter->second.swap(combinedExpCoeffs); // shallow ptr swap
      combinedExpCoeffs.resize(0);
    }
    else // redundant copy
      expCoeffsIter->second = combinedExpCoeffs;     // deep copy
  }
  if (expansionCoeffGradFlag) {
    if (clear_combined) {
      expCoeffGradsIter->second.swap(combinedExpCoeffGrads); // shallow ptr swap
      combinedExpCoeffGrads.reshape(0, 0);
    }
    else // redundant copy
      expCoeffGradsIter->second = combinedExpCoeffGrads;     // deep copy
  }

  allocate_component_sobol();  // size sobolIndices from shared sobolIndexMap

  // If outgoing stats type is active (e.g., as in Dakota::NonDExpansion::
  // multifidelity_expansion()), then previous active stats are invalidated.
  // But if outgoing stats type is combined, then can avoid recomputation
  // and carry over current moment stats from combined to active. 
  // Note: this reuse optimization introduces an order dependency --> updating
  //       stats type from COMBINED to ACTIVE must occur after this function
  if (data_rep->expConfigOptions.refineStatsType == ACTIVE_EXPANSION_STATS)
    clear_computed_bits();
}


void OrthogPolyApproximation::
overlay_expansion(const SizetArray& multi_index_map,
		  const RealVector& exp_coeffs, const RealMatrix& exp_grads,
		  int coeff, RealVector& exp_coeffs_sum,
		  RealMatrix& exp_grads_sum)
{
#ifdef DEBUG
  PCout << "\nmulti_index_map = " << multi_index_map << " Coeff = " << coeff
	<< " Sizes: " << exp_coeffs.length() << ' ' << exp_coeffs_sum.length();
#endif // DEBUG
  size_t i, j, index, num_terms = multi_index_map.size(), 
    num_deriv_v = exp_grads_sum.numRows();
  for (i=0; i<num_terms; ++i) {
    index = multi_index_map[i];
    if (expansionCoeffFlag)
      exp_coeffs_sum[index] += coeff * exp_coeffs[i];
    if (expansionCoeffGradFlag) {
      Real* exp_grad_sum_i = exp_grads_sum[index];
      const Real* grad_i   = exp_grads[i];
      for (j=0; j<num_deriv_v; ++j)
	exp_grad_sum_i[j] += coeff * grad_i[j];
    }
  }
}


/** Implement c = a * b for expansions a, b, and c. */
void OrthogPolyApproximation::
multiply_expansion(const UShort2DArray& multi_index_a,
		   const RealVector&    exp_coeffs_a,
		   const RealMatrix&    exp_grads_a,
		   const UShort2DArray& multi_index_b,
		   const RealVector&    exp_coeffs_b,
		   const RealMatrix&    exp_grads_b,
		   const UShort2DArray& multi_index_c,
		   RealVector& exp_coeffs_c, RealMatrix& exp_grads_c)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;

  size_t i, j, k, v, num_v = sharedDataRep->numVars,
    num_a = multi_index_a.size(), num_b = multi_index_b.size(),
    num_c = multi_index_c.size(), num_deriv_v = exp_grads_a.numRows();

  // precompute 1D basis triple products required
  unsigned short max_a, max_b, max_c; UShortMultiSet max_abc;
  OrthogonalPolynomial* poly_rep_v;
  for (v=0; v<num_v; ++v) {
    max_a = max_b = max_c = 0; max_abc.clear();
    // could track max_abc within combine_coefficients() and pass in, but would
    // need max orders for each dimension for both factors plus their product.
    // Since this would be awkward and only marginally more efficient, just
    // compute them here from the available multi-index arrays.
    for (i=0; i<num_a; ++i)
      if (multi_index_a[i][v] > max_a)
	max_a = multi_index_a[i][v];
    for (i=0; i<num_b; ++i)
      if (multi_index_b[i][v] > max_b)
	max_b = multi_index_b[i][v];
    for (i=0; i<num_c; ++i)
      if (multi_index_c[i][v] > max_c)
	max_c = multi_index_c[i][v];
    max_abc.insert(max_a); max_abc.insert(max_b); max_abc.insert(max_c); 
    poly_rep_v
      = (OrthogonalPolynomial*)(data_rep->polynomialBasis[v].polynomial_rep());
    poly_rep_v->precompute_triple_products(max_abc);
  }

  // For c = a * b, compute coefficient of product expansion as:
  // \Sum_k c_k \Psi_k = \Sum_i \Sum_j a_i b_j \Psi_i \Psi_j
  //    c_k <\Psi_k^2> = \Sum_i \Sum_j a_i b_j <\Psi_i \Psi_j \Psi_k>
  RealVector exp_coeffs_tmp_c;  RealMatrix exp_grads_tmp_c;
  if (expansionCoeffFlag)     exp_coeffs_tmp_c.size(num_c);        // init to 0
  if (expansionCoeffGradFlag) exp_grads_tmp_c.shape(num_deriv_v, num_c);// to 0
  Real trip_prod, trip_prod_v, norm_sq_k; bool non_zero;
  for (k=0; k<num_c; ++k) {
    for (i=0; i<num_a; ++i) {
      for (j=0; j<num_b; ++j) {
	trip_prod = 1.;
	for (v=0; v<num_v; ++v) {
	  poly_rep_v = (OrthogonalPolynomial*)
	    (data_rep->polynomialBasis[v].polynomial_rep());
	  non_zero = poly_rep_v->triple_product(multi_index_a[i][v],
	    multi_index_b[j][v], multi_index_c[k][v], trip_prod_v);
	  if (non_zero) trip_prod *= trip_prod_v;
	  else          break;
	}
	if (non_zero) {
	  if (expansionCoeffFlag)
	    exp_coeffs_tmp_c[k]
	      += exp_coeffs_a[i] * exp_coeffs_b[j] * trip_prod;
	  if (expansionCoeffGradFlag)
	    for (v=0; v<num_deriv_v; ++v)
	      exp_grads_tmp_c(v,k) += (exp_coeffs_a[i] * exp_grads_b(v,j)
		+ exp_coeffs_b[j] * exp_grads_a(v,i)) * trip_prod;
	}
      }
    }
    norm_sq_k = data_rep->norm_squared(multi_index_c[k]);
    if (expansionCoeffFlag)
      exp_coeffs_tmp_c[k] /= norm_sq_k;
    if (expansionCoeffGradFlag)
      for (v=0; v<num_deriv_v; ++v)
	exp_grads_tmp_c(v,k) /= norm_sq_k;
  }
  // tmp arrays allow incoming a and c to be same arrays (for running products)
  exp_coeffs_c = exp_coeffs_tmp_c;
  exp_grads_c  = exp_grads_tmp_c;
}


void OrthogPolyApproximation::
basis_value(const RealVector& x, std::vector<BasisPolynomial>& polynomial_basis,
	    const UShort2DArray& multi_index, RealVector& basis_values)
{
  size_t i, num_exp_terms = multi_index.size();
  for (i=0; i<num_exp_terms; ++i)
    basis_values[i] = SharedOrthogPolyApproxData::
      multivariate_polynomial(x, multi_index[i], polynomial_basis);
}


void OrthogPolyApproximation::
basis_matrix(const RealMatrix& x, RealMatrix &basis_values)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  basis_matrix(x, data_rep->polynomialBasis, data_rep->multi_index(),
	       basis_values);
}


void OrthogPolyApproximation::
basis_matrix(const RealMatrix& x,
	     std::vector<BasisPolynomial>& polynomial_basis,
	     const UShort2DArray& multi_index, RealMatrix& basis_values)
{
  size_t i, j, num_exp_terms = multi_index.size(), num_samples = x.numCols(),
    num_vars = x.numRows();
  basis_values.shapeUninitialized(num_samples,num_exp_terms);
  for (j=0; j<num_exp_terms; ++j){
    for (i=0; i<num_samples; ++i){
      RealVector sample( Teuchos::View, const_cast<double*>(x[i]), num_vars);
      basis_values(i,j) = SharedOrthogPolyApproxData::
	multivariate_polynomial(sample, multi_index[j], polynomial_basis);
    }
  }
}


Real OrthogPolyApproximation::
value(const RealVector& x, const UShort2DArray& mi,
      const RealVector& exp_coeffs)
{
  // Implement value caching here:
  // > sweep over multiIndex to find max orders per dim
  // > pre-compute for x along each dimension
  // > loop over num_exp_terms with fast lookups

  size_t i, num_terms = mi.size();
  if (!expansionCoeffFlag || !num_terms || exp_coeffs.length() != num_terms) {
    PCerr << "Error: expansion coefficients not available in "
	  << "OrthogPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  Real approx_val = 0.;
  for (i=0; i<num_terms; ++i)
    approx_val += exp_coeffs[i] * data_rep->multivariate_polynomial(x, mi[i]);
  return approx_val;
}


const RealVector& OrthogPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort2DArray& mi,
			 const RealVector& exp_coeffs)
{
  // could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  size_t i, j, num_terms = mi.size(), num_v = sharedDataRep->numVars;
  if (!expansionCoeffFlag || !num_terms || exp_coeffs.length() != num_terms) {
    PCerr << "Error: expansion coefficients not available in OrthogPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != num_v) approxGradient.size(num_v); // init to 0
  else                                  approxGradient = 0.;

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  // sum expansion to get response gradient prediction
  for (i=0; i<num_terms; ++i) {
    const RealVector& term_i_grad
      = data_rep->multivariate_polynomial_gradient_vector(x, mi[i]);
    Real coeff_i = exp_coeffs[i];
    for (j=0; j<num_v; ++j)
      approxGradient[j] += coeff_i * term_i_grad[j];
  }
  return approxGradient;
}


const RealVector& OrthogPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv,
			 const UShort2DArray& mi, const RealVector& exp_coeffs)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  size_t i, j, num_v = dvv.size(), num_terms = mi.size();
  if (!expansionCoeffFlag || !num_terms || exp_coeffs.length() != num_terms) {
    PCerr << "Error: expansion coefficients not available in OrthogPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != num_v) approxGradient.size(num_v); // init to 0
  else                                  approxGradient = 0.;

  // sum expansion to get response gradient prediction
  for (i=0; i<num_terms; ++i) {
    const RealVector& term_i_grad
      = data_rep->multivariate_polynomial_gradient_vector(x, mi[i], dvv);
    Real coeff_i = exp_coeffs[i];
    for (j=0; j<num_v; ++j)
      approxGradient[j] += coeff_i * term_i_grad[j];
  }
  return approxGradient;
}


const RealVector& OrthogPolyApproximation::
gradient_nonbasis_variables(const RealVector& x, const UShort2DArray& mi,
			    const RealMatrix& exp_coeff_grads)
{
  size_t i, j, num_v = exp_coeff_grads.numRows(), num_terms = mi.size();
  if (!expansionCoeffGradFlag ||  !num_terms ||
      exp_coeff_grads.numCols() != num_terms) {
    PCerr << "Error: expansion coefficient gradients not available in Orthog"
	  << "PolyApproximation::gradient_nonbasis_variables()" << std::endl;
    abort_handler(-1);
  }
  
  if (approxGradient.length() != num_v) approxGradient.size(num_v); // init to 0
  else                                  approxGradient = 0.;

  // sum expansion to get response gradient prediction
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  for (i=0; i<num_terms; ++i) {
    Real term_i = data_rep->multivariate_polynomial(x, mi[i]);
    const Real* exp_grad_i = exp_coeff_grads[i];
    for (j=0; j<num_v; ++j)
      approxGradient[j] += exp_grad_i[j] * term_i;
  }
  return approxGradient;
}


const RealSymMatrix& OrthogPolyApproximation::
hessian_basis_variables(const RealVector& x, const UShort2DArray& mi,
			const RealVector& exp_coeffs)
{
  size_t i, row, col, num_terms = mi.size(), num_v = sharedDataRep->numVars;
  if (!expansionCoeffFlag || !num_terms || exp_coeffs.length() != num_terms) {
    PCerr << "Error: expansion coefficients not defined in OrthogPoly"
	  << "Approximation::hessian_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  if (approxHessian.numRows() != num_v) approxHessian.shape(num_v); // init to 0
  else                                  approxHessian = 0.;

  // sum expansion to get response hessian prediction
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  for (i=0; i<num_terms; ++i) {
    const RealSymMatrix& term_i_hess
      = data_rep->multivariate_polynomial_hessian_matrix(x, mi[i]);
    Real coeff_i = exp_coeffs[i];
    for (row=0; row<num_v; ++row)
      for (col=0; col<=row; ++col)
        approxHessian(row,col) += coeff_i * term_i_hess(row,col);
    // no operator[] provided for SymMatrix:
    //for (col=0; col<num_v; ++col) {
    //  Real*       approx_hess_col = approxHessian[col];
    //  const Real* term_i_hess_col = term_i_hess[col];
    //  for (row=col; row<num_v; ++row) // lower triangle
    //    approx_hess_col[row] += coeff_i * term_i_hess_col[row];
    //}
  }
  return approxHessian;
}


/** In this case, all expansion variables are random variables and the
    mean of the expansion is simply the first chaos coefficient. */
Real OrthogPolyApproximation::mean()
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compMeanIter->second & 1))
    return expMomentsIter->second[0];

  Real mean = expCoeffsIter->second[0];
  if (std_mode)
    { expMomentsIter->second[0] = mean; compMeanIter->second |= 1; }
  return mean;
}


/** In this case, a subset of the expansion variables are random
    variables and the mean of the expansion involves evaluating the
    expectation over this subset. */
Real OrthogPolyApproximation::mean(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compMeanIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevMean[key]))
    return expMomentsIter->second[0];

  const RealVector& exp_coeffs = expCoeffsIter->second;
  Real mean = exp_coeffs[0];
  const UShort2DArray& mi = data_rep->multi_index();
  size_t i, num_exp_terms = mi.size();
  for (i=1; i<num_exp_terms; ++i)
    // expectations are zero for expansion terms with nonzero random indices
    if (data_rep->zero_random(mi[i])) {
      mean += exp_coeffs[i] *
	data_rep->multivariate_polynomial(x, mi[i], data_rep->nonRandomIndices);
#ifdef DEBUG
      PCout << "Mean estimate inclusion: term index = " << i << " Psi = "
	    << data_rep->
	         multivariate_polynomial(x, mi[i], data_rep->nonRandomIndices)
	    << " PCE coeff = " << exp_coeffs[i] << " total = " << mean
	    << std::endl;
#endif // DEBUG
    }

  if (all_mode) {
    expMomentsIter->second[0] = mean;
    compMeanIter->second |= 1;  xPrevMean[key] = x;
  }
  return mean;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  In
    this case, the derivative of the expectation is the expectation of
    the derivative.  The mixed derivative case (some design variables
    are inserted and some are augmented) requires no special treatment. */
const RealVector& OrthogPolyApproximation::mean_gradient()
{
  // d/ds \mu_R = d/ds \alpha_0 = <dR/ds>

  // Error check for required data
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in "
	  << "OrthogPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (std_mode && (compMeanIter->second & 2))
    return momentGradsIter->second[0];

  RealVector& mean_grad = momentGradsIter->second[0];
  mean_grad = Teuchos::getCol(Teuchos::Copy, expCoeffGradsIter->second, 0);
  if (std_mode) compMeanIter->second |=  2; // activate 2-bit
  else   compMeanIter->second &= ~2; // deactivate 2-bit: protect mixed use
  return mean_grad;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  In this case, the mean of the expansion is the
    expectation over the random subset and the derivative of the mean
    is the derivative of the remaining expansion over the non-random
    subset.  This function must handle the mixed case, where some
    design/state variables are augmented (and are part of the
    expansion: derivatives are evaluated as described above) and some
    are inserted (derivatives are obtained from expansionCoeffGrads). */
const RealVector& OrthogPolyApproximation::
mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  // if already computed, return previous result
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if ( all_mode && (compMeanIter->second & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevMeanGrad[key]) )
    // && dvv == dvvPrev)
    return momentGradsIter->second[0];

  const UShort2DArray& mi = data_rep->multi_index();
  const RealVector& exp_coeffs      = expCoeffsIter->second;
  const RealMatrix& exp_coeff_grads = expCoeffGradsIter->second;
  size_t i, j, deriv_index, num_deriv_v = dvv.size(),
    num_exp_terms = mi.size(),
    cntr = 0; // insertions carried in order within expansionCoeffGrads
  RealVector& mean_grad = momentGradsIter->second[0];
  if (mean_grad.length() != num_deriv_v)
    mean_grad.sizeUninitialized(num_deriv_v);
  for (i=0; i<num_deriv_v; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    bool random = data_rep->randomVarsKey[deriv_index];
    Real& grad_i = mean_grad[i];
    if (random) {// deriv w.r.t. des var insertion
      // Error check for required data
      if (!expansionCoeffGradFlag) {
	PCerr << "Error: expansion coefficient gradients not defined in "
	      << "OrthogPolyApproximation::mean_gradient()." << std::endl;
	abort_handler(-1);
      }
      grad_i = exp_coeff_grads[0][cntr];
    }
    else {
      grad_i = 0.;
      if (!expansionCoeffFlag) { // check for reqd data
	PCerr << "Error: expansion coefficients not defined in "
	      << "OrthogPolyApproximation::mean_gradient()" << std::endl;
	abort_handler(-1);
      }
    }
    for (j=1; j<num_exp_terms; ++j) {
      // expectations are zero for expansion terms with nonzero random indices
      if (data_rep->zero_random(mi[j])) {
	// In both cases below, term to differentiate is alpha_j(s) Psi_j(s)
	// since <Psi_j>_xi = 1 for included terms.  The difference occurs
	// based on whether a particular s_i dependence appears in alpha
	// (for inserted) or Psi (for augmented).
	if (random)
	  // -------------------------------------------
	  // derivative w.r.t. design variable insertion
	  // -------------------------------------------
	  grad_i += exp_coeff_grads[j][cntr] *
	    data_rep->multivariate_polynomial(x, mi[j],
	      data_rep->nonRandomIndices);
	else
	  // ----------------------------------------------
	  // derivative w.r.t. design variable augmentation
	  // ----------------------------------------------
	  grad_i += exp_coeffs[j] *
	    data_rep->multivariate_polynomial_gradient(x, deriv_index, mi[j],
	      data_rep->nonRandomIndices);
      }
    }
    if (random) // deriv w.r.t. des var insertion
      ++cntr;
  }
  if (all_mode) { compMeanIter->second |=  2; xPrevMeanGrad[key] = x; }
  else   compMeanIter->second &= ~2; // deactivate 2-bit: protect mixed use
  return mean_grad;
}


Real OrthogPolyApproximation::
covariance(const UShort2DArray& mi, const RealVector& exp_coeffs,
	   const RealVector& exp_coeffs_2)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  size_t i, num_exp_terms = mi.size();
  Real covar = 0.;
  for (i=1; i<num_exp_terms; ++i)
    covar += exp_coeffs[i] * exp_coeffs_2[i] * data_rep->norm_squared(mi[i]);
  return covar;
}


Real OrthogPolyApproximation::
covariance(PolynomialApproximation* poly_approx_2)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  OrthogPolyApproximation* opa_2 = (OrthogPolyApproximation*)poly_approx_2;
  bool same = (opa_2 == this), std_mode = data_rep->nonRandomIndices.empty();

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !opa_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if (same && std_mode && (compVarIter->second & 1))
    return expMomentsIter->second[1];
  else {
    Real covar = covariance(data_rep->multi_index(), expCoeffsIter->second,
			    opa_2->expCoeffsIter->second);
    if (same && std_mode)
      { expMomentsIter->second[1] = covar; compVarIter->second |= 1; }
    return covar;
  }
}


Real OrthogPolyApproximation::
covariance(const RealVector& x, const UShort2DArray& mi,
	   const RealVector& exp_coeffs, const RealVector& exp_coeffs_2)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const SizetList&  rand_ind = data_rep->randomIndices;
  const SizetList& nrand_ind = data_rep->nonRandomIndices;
  size_t i1, i2, num_mi = mi.size();
  Real covar = 0.;
  for (i1=1; i1<num_mi; ++i1) {
    // For r = random_vars and nr = non_random_vars,
    // sigma^2_R(nr) = < (R(r,nr) - \mu_R(nr))^2 >_r
    // -> only include terms from R(r,nr) which don't appear in \mu_R(nr)
    const UShortArray& mi1 = mi[i1];
    if (!data_rep->zero_random(mi1)) {
      Real coeff_norm_poly = exp_coeffs[i1] * 
	data_rep->norm_squared(mi1, rand_ind) *
	data_rep->multivariate_polynomial(x, mi1, nrand_ind);
      for (i2=1; i2<num_mi; ++i2) {
	// random polynomial part must be identical to contribute to variance
	// (else orthogonality drops term).  Note that it is not necessary to
	// collapse terms with the same random basis subset, since cross term
	// in (a+b)(a+b) = a^2+2ab+b^2 gets included.  If terms were collapsed
	// (following eval of non-random portions), the nested loop could be
	// replaced with a single loop to evaluate (a+b)^2.
	const UShortArray& mi2 = mi[i2];
	if (data_rep->match_random_key(mi1, mi2))
	  covar += coeff_norm_poly * exp_coeffs_2[i2] *
	    data_rep->multivariate_polynomial(x, mi2, nrand_ind);
      }
    }
  }
  return covar;
}


Real OrthogPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  OrthogPolyApproximation* opa_2 = (OrthogPolyApproximation*)poly_approx_2;
  bool same = (this == opa_2), all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !opa_2->expansionCoeffFlag )) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if ( same && all_mode && (compVarIter->second & 1) &&
       data_rep->match_nonrandom_vars(x, xPrevVar[key]) )
    return expMomentsIter->second[1];

  Real covar = covariance(x, data_rep->multi_index(), expCoeffsIter->second,
			  opa_2->expCoeffsIter->second);
  if (same && all_mode) {
    expMomentsIter->second[1] = covar;
    compVarIter->second |= 1;  xPrevVar[key] = x;
  }
  return covar;
}


/** In this case, all expansion variables are random variables and the
    mean of the expansion is simply the first chaos coefficient. */
Real OrthogPolyApproximation::combined_mean()
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compMeanIter->second & 1))
    return expMomentsIter->second[0];

  Real mean = combinedExpCoeffs[0];
  if (std_mode)
    { expMomentsIter->second[0] = mean; compMeanIter->second |= 1; }
  return mean;
}


/** In this case, a subset of the expansion variables are random
    variables and the mean of the expansion involves evaluating the
    expectation over this subset. */
Real OrthogPolyApproximation::combined_mean(const RealVector& x)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compMeanIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevMean[key]))
    return expMomentsIter->second[0];

  Real mean = combinedExpCoeffs[0];
  const UShort2DArray& comb_mi = data_rep->combinedMultiIndex;
  size_t i, num_exp_terms = comb_mi.size();
  for (i=1; i<num_exp_terms; ++i)
    // expectations are zero for expansion terms with nonzero random indices
    if (data_rep->zero_random(comb_mi[i]))
      mean += combinedExpCoeffs[i] * data_rep->
	multivariate_polynomial(x, comb_mi[i], data_rep->nonRandomIndices);

  if (all_mode) {
    expMomentsIter->second[0] = mean;
    compMeanIter->second |= 1;  xPrevMean[key] = x;
  }
  return mean;
}


Real OrthogPolyApproximation::
combined_covariance(PolynomialApproximation* poly_approx_2)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  OrthogPolyApproximation* opa_2 = (OrthogPolyApproximation*)poly_approx_2;
  bool same = (opa_2 == this), std_mode = data_rep->nonRandomIndices.empty();

  if (same && std_mode && (compVarIter->second & 1))
    return expMomentsIter->second[1];

  Real covar = covariance(data_rep->combinedMultiIndex, combinedExpCoeffs,
			  opa_2->combinedExpCoeffs);
  if (same && std_mode)
    { expMomentsIter->second[1] = covar; compVarIter->second |= 1; }
  return covar;
}


Real OrthogPolyApproximation::
combined_covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  OrthogPolyApproximation* opa_2 = (OrthogPolyApproximation*)poly_approx_2;
  bool same = (this == opa_2), all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;

  if ( same && all_mode && (compVarIter->second & 1) &&
       data_rep->match_nonrandom_vars(x, xPrevVar[key]) )
    return expMomentsIter->second[1];

  Real covar = covariance(x, data_rep->combinedMultiIndex, combinedExpCoeffs,
			  opa_2->combinedExpCoeffs);
  if (same && all_mode) {
    expMomentsIter->second[1] = covar;
    compVarIter->second |= 1;  xPrevVar[key] = x;
  }
  return covar;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  The
    mixed derivative case (some design variables are inserted and some
    are augmented) requires no special treatment. */
const RealVector& OrthogPolyApproximation::variance_gradient()
{
  // d/ds \sigma^2_R = Sum_{j=1}^P <Psi^2_j> d/ds \alpha^2_j
  //                 = 2 Sum_{j=1}^P \alpha_j <dR/ds, Psi_j>

  // Error check for required data
  if (!expansionCoeffFlag ||
      !expansionCoeffGradFlag) {
    PCerr << "Error: insufficient expansion coefficient data in "
	  << "OrthogPolyApproximation::variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;
  if (std_mode && (compVarIter->second & 2))
    return momentGradsIter->second[1];

  const UShort2DArray& mi = data_rep->multi_index();
  const RealVector& exp_coeffs      = expCoeffsIter->second;
  const RealMatrix& exp_coeff_grads = expCoeffGradsIter->second;
  size_t i, j, num_deriv_v = exp_coeff_grads.numRows(),
    num_exp_terms = mi.size();
  RealVector& var_grad = momentGradsIter->second[1];
  if (var_grad.length() != num_deriv_v)
    var_grad.sizeUninitialized(num_deriv_v);
  var_grad = 0.;
  for (i=1; i<num_exp_terms; ++i) {
    Real term_i = 2. * exp_coeffs[i] * data_rep->norm_squared(mi[i]);
    for (j=0; j<num_deriv_v; ++j)
      var_grad[j] += term_i * exp_coeff_grads[i][j];
  }
  if (std_mode) compVarIter->second |=  2;
  else   compVarIter->second &= ~2; // deactivate 2-bit: protect mixed use
  return var_grad;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  This function must handle the mixed case, where
    some design/state variables are augmented (and are part of the
    expansion) and some are inserted (derivatives are obtained from
    expansionCoeffGrads). */
const RealVector& OrthogPolyApproximation::
variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::variance_gradient()" << std::endl;
    abort_handler(-1);
  }

  // if already computed, return previous result
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const SizetList& nrand_ind = data_rep->nonRandomIndices;
  bool all_mode = !nrand_ind.empty();
  const UShortArray& key = data_rep->activeKey;
  if ( all_mode && (compVarIter->second & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevVarGrad[key]) )
    // && dvv == dvvPrev)
    return momentGradsIter->second[1];

  const UShort2DArray&   mi = data_rep->multi_index();
  const SizetList& rand_ind = data_rep->randomIndices;
  const RealVector& exp_coeffs      = expCoeffsIter->second;
  const RealMatrix& exp_coeff_grads = expCoeffGradsIter->second;
  size_t i, j, k, deriv_index, num_deriv_v = dvv.size(),
    num_exp_terms = mi.size(),
    cntr = 0; // insertions carried in order within expansionCoeffGrads
  RealVector& var_grad = momentGradsIter->second[1];
  if (var_grad.length() != num_deriv_v)
    var_grad.sizeUninitialized(num_deriv_v);
  var_grad = 0.;

  Real norm_sq_j, poly_j, poly_grad_j, norm_poly_j, coeff_j, coeff_grad_j;
  for (i=0; i<num_deriv_v; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    bool random = data_rep->randomVarsKey[deriv_index];
    if (random && !expansionCoeffGradFlag){
      PCerr << "Error: expansion coefficient gradients not defined in "
	    << "OrthogPolyApproximation::variance_gradient()." << std::endl;
      abort_handler(-1);
    }
    for (j=1; j<num_exp_terms; ++j) {
      const UShortArray& mi_j = mi[j];
      if (!data_rep->zero_random(mi_j)) {
	coeff_j   = exp_coeffs[j];
	norm_sq_j = data_rep->norm_squared(mi_j, rand_ind);
	poly_j    = data_rep->multivariate_polynomial(x, mi_j, nrand_ind);
	if (random) {
	  norm_poly_j  = norm_sq_j * poly_j;
	  coeff_grad_j = exp_coeff_grads[j][cntr];
	}
	else
	  poly_grad_j = data_rep->
	    multivariate_polynomial_gradient(x, deriv_index, mi_j, nrand_ind);
	for (k=1; k<num_exp_terms; ++k) {
	  // random part of polynomial must be identical to contribute to
	  // variance (else orthogonality drops term)
	  const UShortArray& mi_k = mi[k];
	  if (data_rep->match_random_key(mi_j, mi_k)) {
	    // In both cases below, the term to differentiate is
	    // alpha_j(s) alpha_k(s) <Psi_j^2>_xi Psi_j(s) Psi_k(s) and the
	    // difference occurs based on whether a particular s_i dependence
	    // appears in alpha (for inserted) or Psi (for augmented).
	    if (random)
	      // -------------------------------------------
	      // derivative w.r.t. design variable insertion
	      // -------------------------------------------
	      var_grad[i] += norm_poly_j *
		(coeff_j * exp_coeff_grads[k][cntr] +
		 exp_coeffs[k] * coeff_grad_j) *
		data_rep->multivariate_polynomial(x, mi_k, nrand_ind);
	    else
	      // ----------------------------------------------
	      // derivative w.r.t. design variable augmentation
	      // ----------------------------------------------
	      var_grad[i] +=
		coeff_j * exp_coeffs[k] * norm_sq_j *
		// Psi_j * dPsi_k_ds_i + dPsi_j_ds_i * Psi_k
		(poly_j * data_rep->multivariate_polynomial_gradient(x,
		   deriv_index, mi_k, nrand_ind) +
		 poly_grad_j * data_rep->multivariate_polynomial(x, mi_k,
		   nrand_ind));
	  }
	}
      }
    }
    if (random) // deriv w.r.t. des var insertion
      ++cntr;
  }
  if (all_mode) { compVarIter->second |=  2; xPrevVarGrad[key] = x; }
  else   compVarIter->second &= ~2; // deactivate 2-bit: protect mixed use
  return var_grad;
}


void OrthogPolyApproximation::compute_component_sobol()
{
  // sobolIndices are indexed via a bit array, one bit per variable.
  // A bit is turned on for an expansion term if there is a variable
  // dependence (i.e., its multi-index component is non-zero).  Since
  // the Sobol' indices involve a consolidation of variance contributions
  // from the expansion terms, we define a bit array from the multIndex
  // and then use a lookup within sobolIndexMap to assign the expansion
  // term contribution to the correct Sobol' index.

  // iterate through multiIndex and store sensitivities.  Note: sobolIndices[0]
  // (corresponding to constant exp term with no variable dependence) is unused.
  sobolIndices = 0.; // initialize

  // compute and sum the variance contributions for each expansion term.  For
  // all_vars mode, this approach picks up the total expansion variance, which
  // is the desired reference pt for type-agnostic global sensitivity analysis.
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray&           mi = data_rep->multi_index();
  const BitArrayULongMap& index_map = data_rep->sobolIndexMap;
  const RealVector&      exp_coeffs = expCoeffsIter->second;
  size_t i, j, num_exp_terms = mi.size(), num_v = sharedDataRep->numVars;
  Real p_var, sum_p_var = 0.;
  BitArray set(num_v);
  for (i=1; i<num_exp_terms; ++i) {
    const UShortArray& mi_i = mi[i];
    p_var = exp_coeffs(i) * exp_coeffs(i) * data_rep->norm_squared(mi_i);
    sum_p_var += p_var;

    // determine the bit set corresponding to this expansion term
    for (j=0; j<num_v; ++j)
      if (mi_i[j]) set.set(j);   //   activate bit j
      else         set.reset(j); // deactivate bit j

    // lookup the bit set within sobolIndexMap --> increment the correct
    // Sobol' index with the variance contribution from this expansion term.
    BAULMCIter cit = index_map.find(set);
    if (cit != index_map.end()) // may not be found if vbdOrderLimit
      sobolIndices[cit->second] += p_var; // divide by sum_p_var below
  }
  if (sum_p_var > SMALL_NUMBER) // don't attribute variance if zero/negligible
    sobolIndices.scale(1./sum_p_var);

#ifdef DEBUG
  PCout << "In OrthogPolyApproximation::compute_component_sobol(), "
	<< "sobolIndices =\n" << sobolIndices;
#endif // DEBUG
}


void OrthogPolyApproximation::compute_total_sobol() 
{
  totalSobolIndices = 0.;

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  size_t k, num_v = sharedDataRep->numVars;
  const UShort2DArray&      mi = data_rep->multi_index();
  const RealVector& exp_coeffs = expCoeffsIter->second;
  if (data_rep->expConfigOptions.vbdOrderLimit) {
    // all component indices may not be available, so compute total indices
    // from scratch by computing and summing variance contributions for each
    // expansion term
    size_t i, j, num_exp_terms = mi.size();
    Real p_var, sum_p_var = 0., ratio_i;
    for (i=1; i<num_exp_terms; ++i) {
      const UShortArray& mi_i = mi[i];
      p_var = exp_coeffs(i) * exp_coeffs(i) * data_rep->norm_squared(mi_i);
      sum_p_var += p_var;
      // for any constituent variable j in exansion term i, the expansion
      // term contributes to the total sensitivity of variable j
      for (j=0; j<num_v; ++j)
	if (mi_i[j])
	  totalSobolIndices[j] += p_var; // divide by sum_p_var outside loop
    }
    // if negligible variance (e.g., a deterministic test fn), then attribution
    // of this variance is suspect.  Defaulting totalSobolIndices to zero is a
    // good choice since it drops out from anisotropic refinement based on the
    // response-average of these indices.
    if (sum_p_var > SMALL_NUMBER) // avoid division by zero
      totalSobolIndices.scale(1./sum_p_var);
  }
  else {
    const BitArrayULongMap& index_map = data_rep->sobolIndexMap;
    // all component effects are present, so simply add them up:
    // totalSobolIndices parses the bit sets of each of the sobolIndices
    // and adds them to each matching variable bin
    for (BAULMCIter cit=index_map.begin(); cit!=index_map.end(); ++cit)
      for (k=0; k<num_v; ++k) 
        if (cit->first[k]) // var k is present in this Sobol' index
          totalSobolIndices[k] += sobolIndices[cit->second];
  }

#ifdef DEBUG
  PCout << "In OrthogPolyApproximation::compute_total_sobol(), "
	<< "totalSobolIndices =\n" << totalSobolIndices;
#endif // DEBUG
}


const RealVector& OrthogPolyApproximation::dimension_decay_rates()
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multi_index();
  const RealVector& exp_coeffs = expCoeffsIter->second;
  size_t i, j, num_exp_terms = mi.size(),
    num_v = sharedDataRep->numVars;
  if (decayRates.empty())
    decayRates.sizeUninitialized(num_v);

  // define max_orders for each var for sizing LLS matrices/vectors
  UShortArray max_orders(num_v, 0);
  for (i=0; i<num_exp_terms; ++i)
    for (j=0; j<num_v; ++j)
      if (mi[i][j] > max_orders[j])
	max_orders[j] = mi[i][j];
  // size A_vectors and b_vectors
  RealVectorArray A_vectors(num_v), b_vectors(num_v);
  for (i=0; i<num_v; ++i) {
    A_vectors[i].sizeUninitialized(max_orders[i]);
    b_vectors[i].sizeUninitialized(max_orders[i]);
  }

  // populate A_vectors and b_vectors
  unsigned short order, non_zero, var_index, order_index;
  bool univariate;
  for (i=1; i<num_exp_terms; ++i) {
    univariate = true; non_zero = 0;
    for (j=0; j<num_v; ++j) {
      if (mi[i][j]) {
	++non_zero;
	if (non_zero > 1) { univariate = false; break; }
	else { order = mi[i][j]; var_index = j; order_index = order-1; }
      }
    }
    if (univariate) {
      // find a for y = ax + b with x = term order, y = log(coeff), and
      // b = known intercept for order x = 0
      Real norm =
	std::sqrt(data_rep->polynomialBasis[var_index].norm_squared(order)),
	abs_coeff = std::abs(exp_coeffs[i]);
#ifdef DECAY_DEBUG
      PCout << "Univariate contribution: order = " << order << " coeff = "
	    << abs_coeff << " norm = " << norm << '\n';
#endif // DECAY_DEBUG
      A_vectors[var_index][order_index] = (Real)order;
      b_vectors[var_index][order_index] = (abs_coeff > 1.e-25) ?
	std::log10(abs_coeff * norm) : std::log10(norm) - 25.;
    }
  }
#ifdef DECAY_DEBUG
  PCout << "raw b_vectors:\n";
  for (i=0; i<num_v; ++i)
    PCout << "Variable " << i+1 << '\n' << b_vectors[i];
#endif // DECAY_DEBUG

  solve_decay_rates(A_vectors, b_vectors, max_orders);
  return decayRates;
}


void OrthogPolyApproximation::
solve_decay_rates(RealVectorArray& A_vectors, RealVectorArray& b_vectors,
		  UShortArray& max_orders)
{
  // first coefficient is used in each of the LLS solves
  Real log_coeff0 = std::log10(std::abs(expCoeffsIter->second[0])), tol = -10.;
  short last_index_above = -1, new_size;
  size_t i, j, num_v = sharedDataRep->numVars;;

  for (i=0; i<num_v; ++i) {
    RealVector& A_i = A_vectors[i]; RealVector& b_i = b_vectors[i];

    // Handle case of flatline at numerical precision by ignoring subsequent
    // values below a tolerance (allow first, prune second)
    // > high decay rate will de-emphasize refinement, but consider zeroing
    //   out refinement for a direction that is converged below tol (?)
    // > for now, truncate max_orders and scale back {A,b}_vectors
    unsigned short order = max_orders[i];
    for (j=0; j<order; ++j)
      if (b_i[j] > tol)
	last_index_above = j;
    new_size = last_index_above+2; // include one value below tolerance
    if (new_size < order) {
      max_orders[i] = order = new_size;
      A_i.resize(order); // needed for psuedo-inv but not LAPACK
      b_i.resize(order); // needed for psuedo-inv but not LAPACK
    }

    // subtract intercept b for y = Ax+b  ==>  Ax = y-b
    for (j=0; j<order; ++j)
      b_i[j] -= log_coeff0;

    // employ simple 1-D pseudo inverse for LLS:
    //   A^T A x = A^T(y-b)  ==>  x = A^T(y-b) / A^T A
    // negate negative slope in log space such that large>0 is fast
    // convergence, small>0 is slow convergence, and <0 is divergence
    decayRates[i] = -A_i.dot(b_i) / A_i.dot(A_i);
  }

#ifdef DECAY_DEBUG
  PCout << "Intercept log(abs(coeff0)) = " << log_coeff0
	<< "\nb_vectors after truncation & intercept subtraction:\n";
  for (i=0; i<num_v; ++i)
    PCout << "Variable " << i+1 << '\n' << b_vectors[i];
  PCout << "Individual approximation decay:\n" << decayRates;
#endif // DECAY_DEBUG
}


void OrthogPolyApproximation::
print_coefficients(std::ostream& s, const UShort2DArray& mi,
		   const RealVector& exp_coeffs, bool normalized)
{
  // terms and term identifiers
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  size_t i, j, num_exp_terms = mi.size(), num_v = sharedDataRep->numVars;
  char tag[10];
  for (i=0; i<num_exp_terms; ++i) {
    const UShortArray& mi_i = mi[i];
    s << "\n  " << std::setw(WRITE_PRECISION+7);
    if (normalized) // basis is divided by norm, so coeff is multiplied by norm
      s << exp_coeffs[i] * std::sqrt(data_rep->norm_squared(mi_i));
    else
      s << exp_coeffs[i];
    for (j=0; j<num_v; ++j) {
      data_rep->get_tag(tag, j, mi_i[j]);
      s << std::setw(5) << tag;
    }
  }
  s << '\n';
}


void OrthogPolyApproximation::
coefficient_labels(std::vector<std::string>& coeff_labels) const
{
  // terms and term identifiers
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multi_index();
  size_t i, j, num_exp_terms = mi.size(), num_v = sharedDataRep->numVars;
  coeff_labels.reserve(num_exp_terms);
  char tag[10];
  for (i=0; i<num_exp_terms; ++i) {
    const UShortArray& mi_i = mi[i];
    std::string tags;
    for (j=0; j<num_v; ++j) {
      if (j) tags += ' ';
      data_rep->get_tag(tag, j, mi_i[j]);
      tags += tag;
    }
    coeff_labels.push_back(tags);
  }
}

} // namespace Pecos
