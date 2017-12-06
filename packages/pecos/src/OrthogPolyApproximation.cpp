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
  // expansion formulation has been defined in Shared*OrthogPolyApproxData::
  // allocate_data(), and its results are employed below

  allocate_total_sobol();
  allocate_component_sobol();

  // size expansion even if !update_exp_form due to possibility of change
  // to expansion{Coeff,GradFlag} settings
  size_expansion();

  if (expansionMoments.empty())
    expansionMoments.sizeUninitialized(2);
}


void OrthogPolyApproximation::store_coefficients(size_t index)
{
  // Store the aggregated expansion data.  This is used for multifidelity
  // combination, separate from poppedTP{MultiIndex,Coeffs,CoeffGrads} used
  // for generalized sparse grids.

  size_t stored_len = storedExpCoeffs.size();
  if (index == _NPOS || index == stored_len) { // append
    if (expansionCoeffFlag) storedExpCoeffs.push_back(expansionCoeffs);
    else storedExpCoeffs.push_back(RealVector()); // keep indexing consistent
    if (expansionCoeffGradFlag)
      storedExpCoeffGrads.push_back(expansionCoeffGrads);
    else // keep indexing consistent
      storedExpCoeffGrads.push_back(RealMatrix());
  }
  else if (index < stored_len) { // replace
    if (expansionCoeffFlag) storedExpCoeffs[index] = expansionCoeffs;
    else storedExpCoeffs[index] = RealVector();
    if (expansionCoeffGradFlag)
      storedExpCoeffGrads[index] = expansionCoeffGrads;
    else storedExpCoeffGrads[index] = RealMatrix();
  }
  else {
    PCerr << "Error: bad index (" << index << ") passed in OrthogPoly"
	  << "Approximation::store_coefficients()" << std::endl;
    abort_handler(-1);
  }
}


void OrthogPolyApproximation::restore_coefficients(size_t index)
{
  size_t stored_len = storedExpCoeffs.size();
  if (index == _NPOS) {
    expansionCoeffs     = storedExpCoeffs.back();
    expansionCoeffGrads = storedExpCoeffGrads.back();
  }
  else if (index < stored_len) {
    expansionCoeffs     = storedExpCoeffs[index];
    expansionCoeffGrads = storedExpCoeffGrads[index];
  }
  else {
    PCerr << "Error: bad index (" << index << ") passed in OrthogPoly"
	  << "Approximation::restore_coefficients()" << std::endl;
    abort_handler(-1);
  }
}


void OrthogPolyApproximation::swap_coefficients(size_t maximal_index)
{
  if (expansionCoeffFlag) {
    RealVector tmp_vec(expansionCoeffs);
    expansionCoeffs = storedExpCoeffs[maximal_index];
    storedExpCoeffs[maximal_index] = tmp_vec;
  }
  if (expansionCoeffGradFlag) {
    RealMatrix tmp_mat(expansionCoeffGrads);
    expansionCoeffGrads = storedExpCoeffGrads[maximal_index];
    storedExpCoeffGrads[maximal_index] = tmp_mat;
  }
}


void OrthogPolyApproximation::remove_stored_coefficients(size_t index)
{
  size_t stored_len = storedExpCoeffs.size();
  if (index == _NPOS || index == stored_len)
    { storedExpCoeffs.pop_back(); storedExpCoeffGrads.pop_back(); }
  else if (index < stored_len) {
    RealVectorArray::iterator vit = storedExpCoeffs.begin();
    std::advance(vit, index); storedExpCoeffs.erase(vit);
    RealMatrixArray::iterator mit = storedExpCoeffGrads.begin();
    std::advance(mit, index); storedExpCoeffGrads.erase(mit);
  }
}


void OrthogPolyApproximation::
combine_coefficients(short combine_type, size_t maximal_index)
{
  // based on incoming combine_type, combine the data stored previously
  // by store_coefficients()

  // SharedOrthogPolyApproxData::pre_combine_data() appends multi-indices
  // SharedOrthogPolyApproxData::post_combine_data() finalizes multiIndex

  if (maximal_index != _NPOS) {
    swap_coefficients(maximal_index);
    allocate_component_sobol(); // size sobolIndices from shared sobolIndexMap
  }

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  size_t i, num_stored = storedExpCoeffs.size();
  switch (combine_type) {
  case ADD_COMBINE: {
    // Note: would like to preserve tensor indexing (at least for QUADRATURE
    // case) so that Horner's rule performance opt could be used within
    // tensor_product_value()).  However, a tensor result in the overlay
    // will not occur unless one expansion order dominates the other (partial
    // domination results in sum of tensor expansions as for sparse grids).
    // Therefore, stick with the general-purpose expansion overlay and exclude
    // tensor_product_value() usage for combined coefficient sets.

    // resize expansion{Coeffs,CoeffGrads} based on updated multiIndex
    resize_expansion();
    // update expansion{Coeffs,CoeffGrads}
    for (i=0; i<num_stored; ++i)
      overlay_expansion(data_rep->storedMultiIndexMap[i], storedExpCoeffs[i],
			storedExpCoeffGrads[i], 1);
    break;
  }
  case MULT_COMBINE: {
    // perform the multiplication of current and stored expansions
    for (i=0; i<num_stored; ++i)
      multiply_expansion(data_rep->storedMultiIndex[i], storedExpCoeffs[i],
			 storedExpCoeffGrads[i], data_rep->combinedMultiIndex);
    break;
  }
  case ADD_MULT_COMBINE:
    //overlay_expansion(data_rep->storedMultiIndexMap[i], storedExpCoeffs[i],
    //                  storedExpCoeffGrads[i], addCoeffs, addCoeffGrads);
    //multiply_expansion(data_rep->storedMultiIndex[i], storedExpCoeffs[i],
    //                   storedExpCoeffGrads[i], multCoeffs, multCoeffGrads);
    //compute_combine_factors(addCoeffs, multCoeffs);
    //apply_combine_factors();
    PCerr << "Error : additive+multiplicative combination not yet implemented "
	  << "in OrthogPolyApproximation::combine_coefficients()" << std::endl;
    abort_handler(-1);
    break;
  }

  /* Code moved to ProjectOrthogPolyApproximation::integrate_response_moments()
  if (expansionCoeffFlag)     storedExpCoeffs.clear();
  if (expansionCoeffGradFlag) storedExpCoeffGrads.clear();
  */

  computedMean = computedVariance = 0;
}


void OrthogPolyApproximation::
overlay_expansion(const SizetArray& multi_index_map,
		  const RealVector& exp_coeffs, const RealMatrix& exp_grads,
		  int coeff)
{
  size_t i, j, index, num_terms = multi_index_map.size(), 
    num_deriv_vars = expansionCoeffGrads.numRows();
  for (i=0; i<num_terms; ++i) {
    index = multi_index_map[i];
    if (expansionCoeffFlag)
      expansionCoeffs[index] += coeff * exp_coeffs[i];
    if (expansionCoeffGradFlag) {
      Real*       exp_grad_ndx = expansionCoeffGrads[index];
      const Real* grad_i       = exp_grads[i];
      for (j=0; j<num_deriv_vars; ++j)
	exp_grad_ndx[j] += coeff * grad_i[j];
    }
  }
}


void OrthogPolyApproximation::
multiply_expansion(const UShort2DArray& multi_index_b,
		   const RealVector&    exp_coeffs_b,
		   const RealMatrix&    exp_grads_b,
		   const UShort2DArray& multi_index_c)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;

  const UShort2DArray& multi_index_a = data_rep->multiIndex;
  RealVector exp_coeffs_a = expansionCoeffs; // copy (both expConfigOptions)
  RealMatrix exp_grads_a;
  if (expansionCoeffGradFlag)
    exp_grads_a = expansionCoeffGrads;       // copy (CoeffGrads only)
  size_t i, j, k, v, num_v = sharedDataRep->numVars,
    num_a = multi_index_a.size(), num_b = multi_index_b.size(),
    num_c = multi_index_c.size(), num_deriv_vars = exp_grads_a.numRows();

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
  if (expansionCoeffFlag)
    expansionCoeffs.size(num_c);                      // init to 0
  if (expansionCoeffGradFlag)
    expansionCoeffGrads.shape(num_deriv_vars, num_c); // init to 0
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
	    expansionCoeffs[k] += exp_coeffs_a[i] * exp_coeffs_b[j] * trip_prod;
	  if (expansionCoeffGradFlag)
	    for (v=0; v<num_deriv_vars; ++v)
	      expansionCoeffGrads(v,k) += (exp_coeffs_a[i] * exp_grads_b(v,j)
		+ exp_coeffs_b[j] * exp_grads_a(v,i)) * trip_prod;
	}
      }
    }
    norm_sq_k = data_rep->norm_squared(multi_index_c[k]);
    if (expansionCoeffFlag)
      expansionCoeffs[k] /= norm_sq_k;
    if (expansionCoeffGradFlag)
      for (v=0; v<num_deriv_vars; ++v)
	expansionCoeffGrads(v,k) /= norm_sq_k;
  }
}


Real OrthogPolyApproximation::value(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  // Implement value caching here:
  // > sweep over multiIndex to find max orders per dim
  // > pre-compute for x along each dimension
  // > loop over num_exp_terms with fast lookups

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multiIndex;
  Real approx_val = 0.;
  size_t i, num_exp_terms = mi.size();
  for (i=0; i<num_exp_terms; ++i)
    approx_val += expansionCoeffs[i] *
      data_rep->multivariate_polynomial(x, mi[i]);
  return approx_val;
}

void OrthogPolyApproximation::basis_value(const RealVector& x,
		       std::vector<BasisPolynomial> &polynomial_basis,
					  const UShort2DArray &multi_index,
					  RealVector &basis_values)
{
  size_t i, num_exp_terms = multi_index.size();
  for (i=0; i<num_exp_terms; ++i)
    basis_values[i] = 
      SharedOrthogPolyApproxData::multivariate_polynomial(x, multi_index[i],
							  polynomial_basis);
}

void OrthogPolyApproximation::basis_matrix(const RealMatrix& x,
					   RealMatrix &basis_values){
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  basis_matrix(x,data_rep->polynomialBasis,data_rep->multiIndex,
	       basis_values);
}

void OrthogPolyApproximation::basis_matrix(const RealMatrix& x,
		       std::vector<BasisPolynomial> &polynomial_basis,
					  const UShort2DArray &multi_index,
					  RealMatrix &basis_values)
{
  size_t i, j, num_exp_terms = multi_index.size(), num_samples = x.numCols(),
    num_vars = x.numRows();
  basis_values.shapeUninitialized(num_samples,num_exp_terms);
  for (j=0; j<num_exp_terms; ++j){
    for (i=0; i<num_samples; ++i){
      RealVector sample( Teuchos::View, const_cast<double*>(x[i]), num_vars);
      basis_values(i,j) = 
	SharedOrthogPolyApproxData::multivariate_polynomial(sample, 
							    multi_index[j],
							    polynomial_basis);
    }
  }
}


const RealVector& OrthogPolyApproximation::
gradient_basis_variables(const RealVector& x)
{
  // could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multiIndex;
  size_t i, j, num_exp_terms = mi.size(), num_v = sharedDataRep->numVars;
  if (approxGradient.length() != num_v)
    approxGradient.size(num_v); // init to 0
  else
    approxGradient = 0.;

  // sum expansion to get response gradient prediction
  for (i=0; i<num_exp_terms; ++i) {
    const RealVector& term_i_grad
      = data_rep->multivariate_polynomial_gradient_vector(x, mi[i]);
    Real& coeff_i = expansionCoeffs[i];
    for (j=0; j<num_v; ++j)
      approxGradient[j] += coeff_i * term_i_grad[j];
  }
  return approxGradient;
}


const RealVector& OrthogPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multiIndex;
  size_t i, j, num_deriv_vars = dvv.size(), num_exp_terms = mi.size();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.size(num_deriv_vars); // init to 0
  else
    approxGradient = 0.;

  // sum expansion to get response gradient prediction
  for (i=0; i<num_exp_terms; ++i) {
    const RealVector& term_i_grad
      = data_rep->multivariate_polynomial_gradient_vector(x, mi[i], dvv);
    Real& coeff_i = expansionCoeffs[i];
    for (j=0; j<num_deriv_vars; ++j)
      approxGradient[j] += coeff_i * term_i_grad[j];
  }
  return approxGradient;
}


const RealVector& OrthogPolyApproximation::
gradient_nonbasis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in OrthogPoly"
	  << "Approximation::gradient_coefficient_variables()" << std::endl;
    abort_handler(-1);
  }

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multiIndex;
  size_t i, j, num_deriv_vars = expansionCoeffGrads.numRows(),
    num_exp_terms = mi.size();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.size(num_deriv_vars); // init to 0
  else
    approxGradient = 0.;

  // sum expansion to get response gradient prediction
  for (i=0; i<num_exp_terms; ++i) {
    Real term_i = data_rep->multivariate_polynomial(x, mi[i]);
    const Real* exp_coeff_grad_i = expansionCoeffGrads[i];
    for (j=0; j<num_deriv_vars; ++j)
      approxGradient[j] += exp_coeff_grad_i[j] * term_i;
  }
  return approxGradient;
}


const RealSymMatrix& OrthogPolyApproximation::
hessian_basis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::hessian_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multiIndex;
  size_t i, row, col, num_exp_terms = mi.size(), num_v = sharedDataRep->numVars;
  if (approxHessian.numRows() != num_v)
    approxHessian.shape(num_v); // init to 0
  else
    approxHessian = 0.;

  // sum expansion to get response hessian prediction
  for (i=0; i<num_exp_terms; ++i) {
    const RealSymMatrix& term_i_hess
      = data_rep->multivariate_polynomial_hessian_matrix(x, mi[i]);
    Real& coeff_i = expansionCoeffs[i];
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


Real OrthogPolyApproximation::stored_value(const RealVector& x, size_t index)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& stored_mi = data_rep->storedMultiIndex[index];
  const RealVector& stored_coeffs = storedExpCoeffs[index];

  // Error check for required data
  size_t i, num_stored_terms = stored_mi.size();
  if (!num_stored_terms || stored_coeffs.length() != num_stored_terms) {
    PCerr << "Error: stored expansion coefficients not available in "
	  << "OrthogPolyApproximation::stored_value()" << std::endl;
    abort_handler(-1);
  }

  Real approx_val = 0.;
  for (size_t i=0; i<num_stored_terms; ++i)
    approx_val += stored_coeffs[i] *
      data_rep->multivariate_polynomial(x, stored_mi[i]);
  return approx_val;
}


const RealVector& OrthogPolyApproximation::
stored_gradient_basis_variables(const RealVector& x, size_t index)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& stored_mi = data_rep->storedMultiIndex[index];
  const RealVector& stored_coeffs = storedExpCoeffs[index];

  // Error check for required data
  size_t i, j, num_stored_terms = stored_mi.size(),
    num_v = sharedDataRep->numVars;
  if (!num_stored_terms || stored_coeffs.length() != num_stored_terms) {
    PCerr << "Error: stored expansion coefficients not available in OrthogPoly"
	  << "Approximation::stored_gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != num_v)
    approxGradient.size(num_v); // init to 0
  else
    approxGradient = 0.;

  // sum expansion to get response gradient prediction
  for (i=0; i<num_stored_terms; ++i) {
    const RealVector& term_i_grad
      = data_rep->multivariate_polynomial_gradient_vector(x, stored_mi[i]);
    Real coeff_i = stored_coeffs[i];
    for (j=0; j<num_v; ++j)
      approxGradient[j] += coeff_i * term_i_grad[j];
  }
  return approxGradient;
}


const RealVector& OrthogPolyApproximation::
stored_gradient_nonbasis_variables(const RealVector& x, size_t index)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& stored_mi = data_rep->storedMultiIndex[index];
  const RealMatrix& stored_grads = storedExpCoeffGrads[index];

  // Error check for required data
  size_t i, j, num_stored_terms = stored_mi.size(),
    num_deriv_vars = stored_grads.numRows();
  if (!num_stored_terms || stored_grads.numCols() != num_stored_terms) {
    PCerr << "Error: stored expansion coeff grads not available in OrthogPoly"
	  << "Approximation::stored_gradient_nonbasis_variables()" << std::endl;
    abort_handler(-1);
  }

  if (approxGradient.length() != num_deriv_vars)
    approxGradient.size(num_deriv_vars); // init to 0
  else
    approxGradient = 0.;

  // sum expansion to get response gradient prediction
  for (i=0; i<num_stored_terms; ++i) {
    Real term_i = data_rep->multivariate_polynomial(x, stored_mi[i]);
    const Real* stored_grad_i = stored_grads[i];
    for (j=0; j<num_deriv_vars; ++j)
      approxGradient[j] += stored_grad_i[j] * term_i;
  }
  return approxGradient;
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
  if (std_mode && (computedMean & 1))
    return expansionMoments[0];

  Real mean = expansionCoeffs[0];
  if (std_mode)
    { expansionMoments[0] = mean; computedMean |= 1; }
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
  if (all_mode && (computedMean & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevMean))
    return expansionMoments[0];

  Real mean = expansionCoeffs[0];
  const UShort2DArray& mi = data_rep->multiIndex;
  size_t i, num_exp_terms = mi.size();
  for (i=1; i<num_exp_terms; ++i)
    // expectations are zero for expansion terms with nonzero random indices
    if (data_rep->zero_random(mi[i])) {
      mean += expansionCoeffs[i] *
	data_rep->multivariate_polynomial(x, mi[i], data_rep->nonRandomIndices);
#ifdef DEBUG
      PCout << "Mean estimate inclusion: term index = " << i << " Psi = "
	    << data_rep->multivariate_polynomial(x, mi[i],
						 data_rep->nonRandomIndices)
	    << " PCE coeff = " << expansionCoeffs[i] << " total = " << mean
	    << std::endl;
#endif // DEBUG
    }

  if (all_mode)
    { expansionMoments[0] = mean; computedMean |= 1; xPrevMean = x; }
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
  if (std_mode && (computedMean & 2))
    return meanGradient;

  meanGradient = Teuchos::getCol(Teuchos::Copy, expansionCoeffGrads, 0);
  if (std_mode) computedMean |=  2; //   activate 2-bit
  else          computedMean &= ~2; // deactivate 2-bit: protect mixed usage
  return meanGradient;
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
  if ( all_mode && (computedMean & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevMeanGrad) ) // && dvv == dvvPrev)
    return meanGradient;

  const UShort2DArray& mi = data_rep->multiIndex;
  size_t i, j, deriv_index, num_deriv_vars = dvv.size(),
    num_exp_terms = mi.size(),
    cntr = 0; // insertions carried in order within expansionCoeffGrads
  if (meanGradient.length() != num_deriv_vars)
    meanGradient.sizeUninitialized(num_deriv_vars);
  for (i=0; i<num_deriv_vars; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    bool random = data_rep->randomVarsKey[deriv_index];
    Real& grad_i = meanGradient[i];
    if (random) {// deriv w.r.t. des var insertion
      // Error check for required data
      if (!expansionCoeffGradFlag) {
	PCerr << "Error: expansion coefficient gradients not defined in "
	      << "OrthogPolyApproximation::mean_gradient()." << std::endl;
	abort_handler(-1);
      }
      grad_i = expansionCoeffGrads[0][cntr];
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
	  grad_i += expansionCoeffGrads[j][cntr] *
	    data_rep->multivariate_polynomial(x, mi[j],
	      data_rep->nonRandomIndices);
	else
	  // ----------------------------------------------
	  // derivative w.r.t. design variable augmentation
	  // ----------------------------------------------
	  grad_i += expansionCoeffs[j] *
	    data_rep->multivariate_polynomial_gradient(x, deriv_index, mi[j],
	      data_rep->nonRandomIndices);
      }
    }
    if (random) // deriv w.r.t. des var insertion
      ++cntr;
  }
  if (all_mode) { computedMean |=  2; xPrevMeanGrad = x; }
  else            computedMean &= ~2; // deactivate 2-bit: protect mixed usage
  return meanGradient;
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

  if (same && std_mode && (computedVariance & 1))
    return expansionMoments[1];
  else {
    const UShort2DArray& mi = data_rep->multiIndex;
    size_t i, num_exp_terms = mi.size();
    const RealVector& exp_coeffs_2 = opa_2->expansionCoeffs;
    Real covar = 0.;
    for (i=1; i<num_exp_terms; ++i)
      covar += expansionCoeffs[i] * exp_coeffs_2[i]
	    *  data_rep->norm_squared(mi[i]);
    if (same && std_mode)
      { expansionMoments[1] = covar; computedVariance |= 1; }
    return covar;
  }
}


Real OrthogPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  OrthogPolyApproximation* opa_2 = (OrthogPolyApproximation*)poly_approx_2;
  bool same = (this == opa_2), all_mode = !data_rep->nonRandomIndices.empty();

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !opa_2->expansionCoeffFlag )) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "OrthogPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if ( same && all_mode && (computedVariance & 1) &&
       data_rep->match_nonrandom_vars(x, xPrevVar) )
    return expansionMoments[1];

  const RealVector& exp_coeffs_2 = opa_2->expansionCoeffs;
  const UShort2DArray&    mi = data_rep->multiIndex;
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
      Real coeff_norm_poly = expansionCoeffs[i1] * 
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
  if (same && all_mode)
    { expansionMoments[1] = covar; computedVariance |= 1; xPrevVar = x; }
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
  if (std_mode && (computedVariance & 2))
    return varianceGradient;

  const UShort2DArray& mi = data_rep->multiIndex;
  size_t i, j, num_deriv_vars = expansionCoeffGrads.numRows(),
    num_exp_terms = mi.size();
  if (varianceGradient.length() != num_deriv_vars)
    varianceGradient.sizeUninitialized(num_deriv_vars);
  varianceGradient = 0.;
  for (i=1; i<num_exp_terms; ++i) {
    Real term_i = 2. * expansionCoeffs[i] * data_rep->norm_squared(mi[i]);
    for (j=0; j<num_deriv_vars; ++j)
      varianceGradient[j] += term_i * expansionCoeffGrads[i][j];
  }
  if (std_mode) computedVariance |=  2;
  else          computedVariance &= ~2; // deactivate 2-bit: protect mixed usage
  return varianceGradient;
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
  if ( all_mode && (computedVariance & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevVarGrad) ) // && dvv == dvvPrev)
    return varianceGradient;

  const UShort2DArray&   mi = data_rep->multiIndex;
  const SizetList& rand_ind = data_rep->randomIndices;
  size_t i, j, k, deriv_index, num_deriv_vars = dvv.size(),
    num_exp_terms = mi.size(),
    cntr = 0; // insertions carried in order within expansionCoeffGrads
  if (varianceGradient.length() != num_deriv_vars)
    varianceGradient.sizeUninitialized(num_deriv_vars);
  varianceGradient = 0.;

  Real norm_sq_j, poly_j, poly_grad_j, norm_poly_j, coeff_j, coeff_grad_j;
  for (i=0; i<num_deriv_vars; ++i) {
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
	coeff_j   = expansionCoeffs[j];
	norm_sq_j = data_rep->norm_squared(mi_j, rand_ind);
	poly_j    = data_rep->multivariate_polynomial(x, mi_j, nrand_ind);
	if (random) {
	  norm_poly_j  = norm_sq_j * poly_j;
	  coeff_grad_j = expansionCoeffGrads[j][cntr];
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
	      varianceGradient[i] += norm_poly_j *
		(coeff_j * expansionCoeffGrads[k][cntr] +
		 expansionCoeffs[k] * coeff_grad_j) *
		data_rep->multivariate_polynomial(x, mi_k, nrand_ind);
	    else
	      // ----------------------------------------------
	      // derivative w.r.t. design variable augmentation
	      // ----------------------------------------------
	      varianceGradient[i] +=
		coeff_j * expansionCoeffs[k] * norm_sq_j *
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
  if (all_mode) { computedVariance |=  2; xPrevVarGrad = x; }
  else            computedVariance &= ~2;//deactivate 2-bit: protect mixed usage
  return varianceGradient;
}


void OrthogPolyApproximation::integration_checks()
{
  if (surrData.anchor()) {
    PCerr << "Error: anchor point not supported for numerical integration in "
	  << "SharedOrthogPolyApproxData::integration()." << std::endl;
    abort_handler(-1);
  }
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;

  if (!driver_rep) {
    PCerr << "Error: pointer to integration driver required in "
	  << "SharedOrthogPolyApproxData::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }
  size_t num_data_pts = surrData.points(),
    num_grid_pts = driver_rep->grid_size();
  if (num_data_pts != num_grid_pts) {
    PCerr << "Error: number of current points (" << num_data_pts << ") is "
	  << "not consistent with\n       number of points/weights ("
	  << num_grid_pts << ") from integration driver in\n       "
	  << "SharedOrthogPolyApproxData::compute_coefficients()." << std::endl;
    abort_handler(-1);
  }
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
  const UShort2DArray& mi = data_rep->multiIndex;
  const BitArrayULongMap& index_map = data_rep->sobolIndexMap;
  size_t i, j, num_exp_terms = mi.size(), num_v = sharedDataRep->numVars;
  Real p_var, sum_p_var = 0.;
  BitArray set(num_v);
  for (i=1; i<num_exp_terms; ++i) {
    const UShortArray& mi_i = mi[i];
    p_var = expansionCoeffs(i) * expansionCoeffs(i)
          * data_rep->norm_squared(mi_i);
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
	<< "sobolIndices =\n"; write_data(PCout, sobolIndices);
#endif // DEBUG
}


void OrthogPolyApproximation::compute_total_sobol() 
{
  totalSobolIndices = 0.;

  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  size_t k, num_v = sharedDataRep->numVars;
  const UShort2DArray& mi = data_rep->multiIndex;
  if (data_rep->expConfigOptions.vbdOrderLimit) {
    // all component indices may not be available, so compute total indices
    // from scratch by computing and summing variance contributions for each
    // expansion term
    size_t i, j, num_exp_terms = mi.size();
    Real p_var, sum_p_var = 0., ratio_i;
    for (i=1; i<num_exp_terms; ++i) {
      const UShortArray& mi_i = mi[i];
      p_var = expansionCoeffs(i) * expansionCoeffs(i)
            * data_rep->norm_squared(mi_i);
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
	<< "totalSobolIndices =\n"; write_data(PCout, totalSobolIndices);
#endif // DEBUG
}


const RealVector& OrthogPolyApproximation::dimension_decay_rates()
{
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multiIndex;
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
	abs_coeff = std::abs(expansionCoeffs[i]);
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
    { PCout << "Variable " << i+1 << '\n'; write_data(PCout, b_vectors[i]); }
#endif // DECAY_DEBUG

  solve_decay_rates(A_vectors, b_vectors, max_orders);
  return decayRates;
}


void OrthogPolyApproximation::
solve_decay_rates(RealVectorArray& A_vectors, RealVectorArray& b_vectors,
		  UShortArray& max_orders)
{
  // first coefficient is used in each of the LLS solves
  Real log_coeff0 = std::log10(std::abs(expansionCoeffs[0])), tol = -10.;
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
  PCout << "Intercept log(abs(coeff0)) = " << log_coeff0 << '\n';
  PCout << "b_vectors after truncation & intercept subtraction:\n";
  for (i=0; i<num_v; ++i)
    { PCout << "Variable " << i+1 << '\n'; write_data(PCout, b_vectors[i]); }
  PCout << "Individual approximation decay:\n"; write_data(PCout, decayRates);
#endif // DECAY_DEBUG
}


void OrthogPolyApproximation::
print_coefficients(std::ostream& s, bool normalized)
{
  // terms and term identifiers
  SharedOrthogPolyApproxData* data_rep
    = (SharedOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multiIndex;
  size_t i, j, num_exp_terms = mi.size(), num_v = sharedDataRep->numVars;
  char tag[10];
  for (i=0; i<num_exp_terms; ++i) {
    const UShortArray& mi_i = mi[i];
    s << "\n  " << std::setw(WRITE_PRECISION+7);
    if (normalized) // basis is divided by norm, so coeff is multiplied by norm
      s << expansionCoeffs[i] * std::sqrt(data_rep->norm_squared(mi_i));
    else
      s << expansionCoeffs[i];
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
  const UShort2DArray& mi = data_rep->multiIndex;
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
