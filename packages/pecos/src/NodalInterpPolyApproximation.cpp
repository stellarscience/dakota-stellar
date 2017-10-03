/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        NodalInterpPolyApproximation
//- Description:  Implementation code for NodalInterpPolyApproximation class
//-               
//- Owner:        Mike Eldred

#include "NodalInterpPolyApproximation.hpp"
#include "SharedNodalInterpPolyApproxData.hpp"
#include "TensorProductDriver.hpp"
#include "CombinedSparseGridDriver.hpp"
#include "InterpolationPolynomial.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

//#define DEBUG
//#define VBD_DEBUG

namespace Pecos {


void NodalInterpPolyApproximation::allocate_arrays()
{
  InterpPolyApproximation::allocate_arrays();

  size_t num_colloc_pts = surrData.points(),
    num_deriv_vars = surrData.num_derivative_variables();
  if (surrData.anchor()) ++num_colloc_pts;
  if (expansionCoeffFlag) {
    if (expansionType1Coeffs.length() != num_colloc_pts)
      expansionType1Coeffs.sizeUninitialized(num_colloc_pts);
    SharedNodalInterpPolyApproxData* data_rep
      = (SharedNodalInterpPolyApproxData*)sharedDataRep;
    if ( data_rep->basisConfigOptions.useDerivs &&
	 ( expansionType2Coeffs.numRows() != num_deriv_vars ||
	   expansionType2Coeffs.numCols() != num_colloc_pts ) )
      expansionType2Coeffs.shapeUninitialized(num_deriv_vars, num_colloc_pts);
  }
  if ( expansionCoeffGradFlag &&
       ( expansionType1CoeffGrads.numRows() != num_deriv_vars ||
	 expansionType1CoeffGrads.numCols() != num_colloc_pts ) )
    expansionType1CoeffGrads.shapeUninitialized(num_deriv_vars, num_colloc_pts);

  // checking num_colloc_pts is insufficient due to anisotropy --> changes in
  // anisotropic weights could move points around without changing the total.
  //bool update_exp_form =
  //  ( (expansionCoeffFlag && expansionType1Coeffs.length() != num_colloc_pts)
  // || (expansionCoeffGradFlag &&
  //     expansionType1CoeffGrads.numCols() != num_colloc_pts ) );
}


void NodalInterpPolyApproximation::compute_expansion_coefficients()
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  size_t index = 0, num_colloc_pts = surrData.points(), offset = 0;
  if (surrData.anchor()) {
    offset = 1; ++num_colloc_pts;
    if (expansionCoeffFlag) {
      expansionType1Coeffs[0] = surrData.anchor_function();
      if (data_rep->basisConfigOptions.useDerivs)
	Teuchos::setCol(surrData.anchor_gradient(), 0, expansionType2Coeffs);
    }
    if (expansionCoeffGradFlag)
      Teuchos::setCol(surrData.anchor_gradient(), 0, expansionType1CoeffGrads);
  }

  for (int i=offset; i<num_colloc_pts; ++i, ++index) {
    if (expansionCoeffFlag) {
      expansionType1Coeffs[i] = surrData.response_function(index);
      // Note: gradients from DAKOTA already scaled in u-space Recast
      if (data_rep->basisConfigOptions.useDerivs)
	Teuchos::setCol(surrData.response_gradient(index), i,
			expansionType2Coeffs);
    }
    if (expansionCoeffGradFlag)
      Teuchos::setCol(surrData.response_gradient(index), i,
		      expansionType1CoeffGrads);
  }

  computedMean = computedVariance = 0;
}


void NodalInterpPolyApproximation::store_coefficients(size_t index)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;

  size_t stored_len = storedExpType1Coeffs.size();
  if (index == _NPOS || index == stored_len) { // append
    if (expansionCoeffFlag) {
      storedExpType1Coeffs.push_back(expansionType1Coeffs);
      storedExpType2Coeffs.push_back(expansionType2Coeffs);
    }
    else { // keep indexing consistent
      storedExpType1Coeffs.push_back(RealVector());
      storedExpType2Coeffs.push_back(RealMatrix());
    }
    if (expansionCoeffGradFlag)
      storedExpType1CoeffGrads.push_back(expansionType1CoeffGrads);
    else // keep indexing consistent
      storedExpType1CoeffGrads.push_back(RealMatrix());
  }
  else if (index < stored_len) { // replace
    if (expansionCoeffFlag) {
      storedExpType1Coeffs[index] = expansionType1Coeffs;
      storedExpType2Coeffs[index] = expansionType2Coeffs;
    }
    else { // keep indexing consistent
      storedExpType1Coeffs[index] = RealVector();
      storedExpType2Coeffs[index] = RealMatrix();
    }
    if (expansionCoeffGradFlag)
      storedExpType1CoeffGrads[index] = expansionType1CoeffGrads;
    else // keep indexing consistent
      storedExpType1CoeffGrads[index] = RealMatrix();
  }
  else {
    PCerr << "Error: bad index (" << index << ") passed in NodalInterpPoly"
	  << "Approximation::store_coefficients()" << std::endl;
    abort_handler(-1);
  }
}


void NodalInterpPolyApproximation::restore_coefficients(size_t index)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;

  size_t stored_len = storedExpType1Coeffs.size();
  if (index == _NPOS) {
    expansionType1Coeffs = storedExpType1Coeffs.back();
    expansionType2Coeffs = storedExpType2Coeffs.back();
    expansionType1CoeffGrads = storedExpType1CoeffGrads.back();
  }
  else if (index < stored_len) {
    expansionType1Coeffs = storedExpType1Coeffs[index];
    expansionType2Coeffs = storedExpType2Coeffs[index];
    expansionType1CoeffGrads = storedExpType1CoeffGrads[index];
  }
  else {
    PCerr << "Error: bad index (" << index << ") passed in NodalInterpPoly"
	  << "Approximation::restore_coefficients()" << std::endl;
    abort_handler(-1);
  }
}


void NodalInterpPolyApproximation::swap_coefficients(size_t index)
{
  if (expansionCoeffFlag) {
    RealVector tmp_vec(expansionType1Coeffs);
    expansionType1Coeffs = storedExpType1Coeffs[index];
    storedExpType1Coeffs[index] = tmp_vec;
    SharedNodalInterpPolyApproxData* data_rep
      = (SharedNodalInterpPolyApproxData*)sharedDataRep;
    if (data_rep->basisConfigOptions.useDerivs) {
      RealMatrix tmp_mat(expansionType2Coeffs);
      expansionType2Coeffs = storedExpType2Coeffs[index];
      storedExpType2Coeffs[index] = tmp_mat;
    }
  }
  if (expansionCoeffGradFlag) {
    RealMatrix tmp_mat(expansionType1CoeffGrads);
    expansionType1CoeffGrads = storedExpType1CoeffGrads[index];
    storedExpType1CoeffGrads[index] = tmp_mat;
  }
}


void NodalInterpPolyApproximation::remove_stored_coefficients(size_t index)
{
  size_t stored_len = storedExpType1Coeffs.size();
  if (index == _NPOS || index == stored_len) {
    storedExpType1Coeffs.pop_back(); storedExpType2Coeffs.pop_back();
    storedExpType1CoeffGrads.pop_back();
  }
  else if (index < stored_len) {
    RealVectorArray::iterator vit = storedExpType1Coeffs.begin();
    std::advance(vit, index); storedExpType1Coeffs.erase(vit);
    RealMatrixArray::iterator mit = storedExpType2Coeffs.begin();
    std::advance(mit, index); storedExpType2Coeffs.erase(mit);
    mit = storedExpType1CoeffGrads.begin();
    std::advance(mit, index); storedExpType1CoeffGrads.erase(mit);
  }
}


void NodalInterpPolyApproximation::
combine_coefficients(short combine_type, size_t maximal_index)
{
#ifdef DEBUG
  PCout << "Original type1 expansion coefficients prior to combination:\n";
  write_data(PCout, expansionType1Coeffs);
#endif // DEBUG

  // SharedNodalInterpPolyApproxData::pre_combine_data() has already swapped
  if (maximal_index != _NPOS) {
    swap_coefficients(maximal_index);
    allocate_component_sobol(); // size sobolIndices from shared sobolIndexMap
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;

  // update expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads} by adding or
  // multiplying stored expansion evaluated at current collocation points
  size_t i, v, s, t, offset = 0, num_pts = surrData.points(),
    num_stored = storedExpType1Coeffs.size();
  bool anchor_pt = surrData.anchor();
  if (anchor_pt) { offset = 1; ++num_pts; }
  Real curr_val;
  RealVector stored_vals(num_stored, false);
  //RealVectorArray stored_grads(num_stored);
  for (i=0; i<num_pts; ++i) {
    const RealVector& c_vars = (anchor_pt && i == 0) ?
      surrData.anchor_continuous_variables() :
      surrData.continuous_variables(i-offset);
    if (combine_type == MULT_COMBINE) { // eval once for both Coeffs/CoeffGrads
      curr_val = expansionType1Coeffs[i]; // copy prior to update
      for (s=0; s<num_stored; ++s)
	stored_vals[s] = stored_value(c_vars, s);
      //if (data_rep->basisConfigOptions.useDerivs) 
      //  for (s=0; s<num_stored; ++s)
      //    stored_grads[s] = stored_gradient_basis_variables(c_vars, s);
    }
    if (expansionCoeffFlag) {
      // split up type1/type2 contribs so increments are performed properly
      if (combine_type == ADD_COMBINE)
	for (s=0; s<num_stored; ++s)
	  expansionType1Coeffs[i] += stored_value(c_vars, s);
      else if (combine_type == MULT_COMBINE)
	for (s=0; s<num_stored; ++s)
	  expansionType1Coeffs[i] *= stored_vals[s];
      if (data_rep->basisConfigOptions.useDerivs) {
	Real* exp_t2_coeffs_i = expansionType2Coeffs[i];
	size_t num_deriv_vars = expansionType2Coeffs.numRows();
	if (combine_type == ADD_COMBINE) {
	  for (s=0; s<num_stored; ++s) {
	    const RealVector& stored_grad
	      = stored_gradient_basis_variables(c_vars, s);
	    for (v=0; v<num_deriv_vars; ++v)
	      exp_t2_coeffs_i[v] += stored_grad[v];
	  }
	}
	else if (combine_type == MULT_COMBINE) {
	  // hf = curr*stored --> dhf/dx = dcurr/dx*stored + curr*dstored/dx
	  Real prod = 1.;
	  // first term: dcurr/dx * st * ... * st
	  for (s=0; s<num_stored; ++s)	   prod *= stored_vals[s];
	  for (v=0; v<num_deriv_vars; ++v) exp_t2_coeffs_i[v] *= prod;
	  // second term: curr * st * ... * st * dstored/dx * st * ... * st
	  for (s=0; s<num_stored; ++s) {
	    const RealVector& stored_grad
	      = stored_gradient_basis_variables(c_vars, s);
	    prod = curr_val;
	    for (t=0; t<num_stored; ++t)
	      if (t != s)
		prod *= stored_vals[t];
	    for (v=0; v<num_deriv_vars; ++v)
	      exp_t2_coeffs_i[v] += stored_grad[v] * prod;
	  }
	}
      }
    }
    if (expansionCoeffGradFlag) {
      Real*   exp_t1_grad_i = expansionType1CoeffGrads[i];
      size_t num_deriv_vars = expansionType1CoeffGrads.numRows();
      if (combine_type == ADD_COMBINE) {
	for (s=0; s<num_stored; ++s) {
	  const RealVector& stored_grad
	    = stored_gradient_nonbasis_variables(c_vars, s);
	  for (v=0; v<num_deriv_vars; ++v)
	    exp_t1_grad_i[v] += stored_grad[v];
	}
      }
      else if (combine_type == MULT_COMBINE) {
	// hf = curr*stored --> dhf/dx = dcurr/dx*stored + curr*dstored/dx
	Real prod = 1.;
	// first term: dcurr/dx * st * ... * st
	for (s=0; s<num_stored; ++s)	 prod *= stored_vals[s];
	for (v=0; v<num_deriv_vars; ++v) exp_t1_grad_i[v] *= prod;
	// second term: curr * st * ... * st * dstored/dx * st * ... * st
	for (s=0; s<num_stored; ++s) {
	  const RealVector& stored_grad
	    = stored_gradient_nonbasis_variables(c_vars, s);
	  prod = curr_val;
	  for (t=0; t<num_stored; ++t)
	    if (t != s)
	      prod *= stored_vals[t];
	  for (v=0; v<num_deriv_vars; ++v)
	    exp_t1_grad_i[v] += stored_grad[v] * prod;
	}
      }
    }
  }
#ifdef DEBUG
  PCout << "Updated type1 expansion coefficients following combination:\n";
  write_data(PCout, expansionType1Coeffs);
#endif // DEBUG

  // clear stored data now that it has been combined
  if (expansionCoeffFlag) {
    storedExpType1Coeffs.clear();
    if (data_rep->basisConfigOptions.useDerivs)
      storedExpType2Coeffs.clear();
  }
  if (expansionCoeffGradFlag)
    storedExpType1CoeffGrads.clear();

  computedMean = computedVariance = 0;
}


void NodalInterpPolyApproximation::push_coefficients()
{
  size_t index, offset = 0, old_colloc_pts, new_colloc_pts = surrData.points();
  if (surrData.anchor())
    { offset = 1; ++new_colloc_pts; }

  if (expansionCoeffFlag) {
    old_colloc_pts = expansionType1Coeffs.length();
    expansionType1Coeffs.resize(new_colloc_pts);
    SharedNodalInterpPolyApproxData* data_rep
      = (SharedNodalInterpPolyApproxData*)sharedDataRep;
    if (data_rep->basisConfigOptions.useDerivs) {
      size_t num_deriv_vars = expansionType2Coeffs.numRows();
      expansionType2Coeffs.reshape(num_deriv_vars, new_colloc_pts);
    }
    index = old_colloc_pts - offset;
    for (int i=old_colloc_pts; i<new_colloc_pts; ++i, ++index) {
      expansionType1Coeffs[i] = surrData.response_function(index);
      if (data_rep->basisConfigOptions.useDerivs)
	Teuchos::setCol(surrData.response_gradient(index), i,
			expansionType2Coeffs);
    }
  }
  if (expansionCoeffGradFlag) {
    old_colloc_pts = expansionType1CoeffGrads.numCols();
    size_t num_deriv_vars = expansionType1CoeffGrads.numRows();
    expansionType1CoeffGrads.reshape(num_deriv_vars, new_colloc_pts);
    index = old_colloc_pts - offset;
    for (int i=old_colloc_pts; i<new_colloc_pts; ++i, ++index)
      Teuchos::setCol(surrData.response_gradient(index), i,
		      expansionType1CoeffGrads);
  }

  computedMean = computedVariance = 0;
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
Real NodalInterpPolyApproximation::
tensor_product_mean(const RealVector&    x,   const UShortArray& lev_index,
		    const UShort2DArray& key, const SizetArray&  colloc_index)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::tensor_product_mean()" << std::endl;
    abort_handler(-1);
  }

  size_t i, j, num_colloc_pts = key.size(), num_v = sharedDataRep->numVars;
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  if (data_rep->barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    data_rep->set_new_point(x, lev_index, data_rep->nonRandomIndices, 1);

    // Since we need to include weights for integrated dimensions, loop over
    // all collocation points even if some interpolated dimensions are inactive.
    RealVector accumulator(num_v); // init to 0.
    const Real3DArray& t1_wts_1d = driver_rep->type1_collocation_weights_1d();
    unsigned short       li_0 = lev_index[0], li_j;
    BasisPolynomial&   poly_0 = data_rep->polynomialBasis[li_0][0];
    const RealArray& t1_wts_0 = t1_wts_1d[li_0][0];
    const RealVector& bc_vf_0 = poly_0.barycentric_value_factors();
    size_t c_index, ei0 = poly_0.exact_index(), eij;
    unsigned short key_i0, key_ij, max0 = poly_0.interpolation_size() - 1;
    bool rand_0 = data_rep->randomVarsKey[0];
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i]; key_i0 = key_i[0];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      if (li_0 == 0)          // single integration/interpolation weight = 1
	accumulator[0]  = expansionType1Coeffs[c_index];
      else if (rand_0)        // integration
	accumulator[0] += expansionType1Coeffs[c_index] * t1_wts_0[key_i0];
      else if (ei0 == _NPOS)  // interpolation without exact match
	accumulator[0] += expansionType1Coeffs[c_index] *  bc_vf_0[key_i0];
      else if (ei0 == key_i0) // interpolation with exact match
	accumulator[0]  = expansionType1Coeffs[c_index];
      if (key_i0 == max0)
	data_rep->accumulate_barycentric(accumulator, lev_index, key_i);
    }
    return accumulator[num_v-1] /
      data_rep->barycentric_value_scaling(lev_index,data_rep->nonRandomIndices);
  }
  else if (data_rep->basisConfigOptions.useDerivs) {
    // Horner's rule approach:
    //PCout << "tensor_product_mean(): Horner's useDerivs." << std::endl;
    // Note: while possible to use a single accumulator vector and loop over
    // it n+1 times, value() and gradient() data would then have to be stored.
    RealVector t1_accumulator(num_v);          // init to 0.
    RealMatrix t2_accumulator(num_v, num_v); // init to 0.
    Real *t2_accum_0 = t2_accumulator[0], t1_val, x0 = x[0];
    unsigned short       li_0 = lev_index[0];
    BasisPolynomial&   poly_0 = data_rep->polynomialBasis[li_0][0];
    const Real3DArray& t1_wts_1d = driver_rep->type1_collocation_weights_1d();
    const Real3DArray& t2_wts_1d = driver_rep->type2_collocation_weights_1d();
    const RealArray& t1_wts_0 = t1_wts_1d[li_0][0];
    const RealArray& t2_wts_0 = t2_wts_1d[li_0][0];
    size_t k, c_index; bool rand_0 = data_rep->randomVarsKey[0];
    unsigned short key_i0, max0 = poly_0.interpolation_size() - 1;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i]; key_i0 = key_i[0];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      Real* t2_coeffs_i = expansionType2Coeffs[c_index];
      if (rand_0) { // integration
	if (li_0 == 0) {
	  t1_accumulator[0]  = expansionType1Coeffs[c_index]; // t1 weight is 1
	  //t2_accum_0[0]    = 0.;                            // t2 weight is 0
	  for (k=1; k<num_v; ++k)
	    t2_accum_0[k]    = t2_coeffs_i[k];                // t1 weight is 1
	}
	else {
	  t1_accumulator[0] += expansionType1Coeffs[c_index] * t1_wts_0[key_i0];
	  t2_accum_0[0]     += t2_coeffs_i[0] * t2_wts_0[key_i0];
	  for (k=1; k<num_v; ++k)
	    t2_accum_0[k]   += t2_coeffs_i[k] * t1_wts_0[key_i0];
	}
      }
      else {        // interpolation
	if (li_0 == 0) {
	  t1_accumulator[0]  = expansionType1Coeffs[c_index];  // t1 value is 1
	  t2_accum_0[0]      = t2_coeffs_i[0] * poly_0.type2_value(x0, key_i0);
	  for (k=1; k<num_v; ++k)
	    t2_accum_0[k]    = t2_coeffs_i[k];                 // t1 value is 1
	}
	else {
	  t1_val = poly_0.type1_value(x0, key_i0);
          t1_accumulator[0] += expansionType1Coeffs[c_index] * t1_val;
	  t2_accum_0[0]     += t2_coeffs_i[0] * poly_0.type2_value(x0, key_i0);
	  for (k=1; k<num_v; ++k)
	    t2_accum_0[k]   += t2_coeffs_i[k] * t1_val;
	}
      }
      if (key_i0 == max0)
	data_rep->accumulate_horners(t1_accumulator, t2_accumulator, lev_index,
				     key_i, x);
    }
    Real  tp_mean  = t1_accumulator[num_v-1];
    Real* t2_accum = t2_accumulator[num_v-1];
    for (j=0; j<num_v; ++j)
      tp_mean += t2_accum[j];
    return tp_mean;

    /*
    // Simpler but more expensive approach:
    //PCout << "tensor_product_mean(): Original useDerivs." << std::endl;
    Real tp_mean = 0.; size_t c_index;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      tp_mean += expansionType1Coeffs[c_index]
	      *  type1_interpolant_value(x, key_i, lev_index,
	                                 data_rep->nonRandomIndices)
	      *  type1_weight(key_i, lev_index, data_rep->randomIndices);
      const Real *t2_coeff_i = expansionType2Coeffs[c_index];
      for (j=0; j<num_v; ++j)
	tp_mean += t2_coeff_i[j]
	  * type2_interpolant_value(x, j, key_i, lev_index,
	                            data_rep->nonRandomIndices)
	  * type2_weight(j, key_i, lev_index, data_rep->randomIndices);
    }
    return tp_mean;
    */
  }
  else {
    // Horner's rule approach:
    //PCout << "tensor_product_mean(): Horner's no derivs." << std::endl;
    RealVector accumulator(num_v); // init to 0.
    const Real3DArray& t1_wts_1d = driver_rep->type1_collocation_weights_1d();
    unsigned short       li_0 = lev_index[0];
    BasisPolynomial&   poly_0 = data_rep->polynomialBasis[li_0][0];
    const RealArray& t1_wts_0 = t1_wts_1d[li_0][0];
    size_t c_index; bool rand_0 = data_rep->randomVarsKey[0]; Real x0 = x[0];
    unsigned short key_i0, max0 = poly_0.interpolation_size() - 1;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i]; key_i0 = key_i[0];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      if (li_0 == 0)   // single integration/interpolation weight = 1
	accumulator[0]  = expansionType1Coeffs[c_index];
      else if (rand_0) // integration
	accumulator[0] += expansionType1Coeffs[c_index] * t1_wts_0[key_i0];
      else             // interpolation
	accumulator[0] += expansionType1Coeffs[c_index]
	               *  poly_0.type1_value(x0, key_i0);
      if (key_i0 == max0)
	data_rep->accumulate_horners(accumulator, lev_index, key_i, x);
    }
    return accumulator[num_v-1];

    /*
    // Simpler but more expensive approach:
    //PCout << "tensor_product_mean(): Original no derivs." << std::endl;
    Real tp_mean = 0.;
    if (colloc_index.empty())
      for (i=0; i<num_colloc_pts; ++i)
	tp_mean += expansionType1Coeffs[i]
	  * type1_interpolant_value(x, key[i], lev_index,
	                            data_rep->nonRandomIndices)
	  * type1_weight(key[i], lev_index, data_rep->randomIndices);
      else
	for (i=0; i<num_colloc_pts; ++i)
	  tp_mean += expansionType1Coeffs[colloc_index[i]]
	    * type1_interpolant_value(x, key[i], lev_index,
	                              data_rep->nonRandomIndices)
	    * type1_weight(key[i], lev_index, data_rep->randomIndices);
    return tp_mean;
    */
  }
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_mean_gradient(const RealVector& x, const UShortArray& lev_index,
			     const UShort2DArray& key,
			     const SizetArray& colloc_index,
			     const SizetArray& dvv)
{
  // ---------------------------------------------------------------------
  // For xi = ran vars, sa = augmented des vars, si = inserted design vars
  // Active variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
  //   mu(sa, si)    = Sum_i r_i(sa, si) wt_prod_i
  //   dmu/ds        = Sum_i dr_i/ds wt_prod_i
  // All variable expansion:
  //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
  //   mu(sa, si)    = Sum_i r_i(si) Lsa_i wt_prod_i
  //   dmu/dsa       = Sum_i r_i(si) dLsa_i/dsa wt_prod_i
  //   dmu/dsi       = Sum_i dr_i/dsi Lsa_i wt_prod_i
  // ---------------------------------------------------------------------
  size_t p, d, v, c_index, deriv_index, insert_cntr = 0, 
    num_deriv_vars = dvv.size(), num_colloc_pts = key.size(),
    num_v = sharedDataRep->numVars;
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  if (tpMeanGrad.length() != num_deriv_vars)
    tpMeanGrad.sizeUninitialized(num_deriv_vars);
  tpMeanGrad = 0.;

  // screen for insertion and augmentation and perform error checks
  bool insert = false, augment = false;
  for (d=0; d<num_deriv_vars; ++d) {
    deriv_index = dvv[d] - 1; // OK since we are in an "All" view
    if (data_rep->randomVarsKey[deriv_index]) insert  = true;
    else                                      augment = true;
  }
  if (insert && !expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	  << "InterpPolyApproximation::tensor_product_mean_gradient()."
	  << std::endl;
    abort_handler(-1);
  }
  if (augment && !expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::tensor_product_mean_gradient()" << std::endl;
    abort_handler(-1);
  }

  // Four cases to manage:
  // 1:    random var,  inserted deriv: exp_t1_coeff_grads * type1_weight
  // 2:    random var, augmented deriv: exp_t1_coeffs      * type1_weight
  // 3: nonrandom var,  inserted deriv: exp_t1_coeff_grads * type1_value
  // 4: nonrandom var, augmented deriv: exp_t1_coeffs      * type1_gradient

  if (data_rep->barycentricFlag) {

    // For barycentric interpolation: track x != newPoint within 1D basis
    short order = (augment) ? 3 : 1;
    data_rep->set_new_point(x, lev_index, data_rep->nonRandomIndices, order);

    unsigned short  key_p0, li_0 = lev_index[0];
    BasisPolynomial&      poly_0 = data_rep->polynomialBasis[li_0][0];
    const RealVector&    bc_vf_0 = poly_0.barycentric_value_factors();
    const RealVector&    bc_gf_0 = poly_0.barycentric_gradient_factors();
    const Real3DArray& t1_wts_1d = driver_rep->type1_collocation_weights_1d();
    const RealArray&   t1_wts_0  = t1_wts_1d[li_0][0];
    size_t ei_0 = poly_0.exact_index();
    unsigned short max0 = poly_0.interpolation_size() - 1;
    RealMatrix accumulator(num_deriv_vars, num_v); // init to 0.
    Real *accum_0 = accumulator[0], *t1_coeff_grad, t1_coeff, t1_wt_00,
      bc_vf_00, prod;
    bool rand_0 = data_rep->randomVarsKey[0];
    for (p=0; p<num_colloc_pts; ++p) {
      const UShortArray& key_p = key[p]; key_p0 = key_p[0];
      if (augment)
	t1_coeff = (colloc_index.empty()) ? expansionType1Coeffs[p] :
	  expansionType1Coeffs[colloc_index[p]];
      if (insert)
	t1_coeff_grad = (colloc_index.empty()) ? expansionType1CoeffGrads[p] :
	  expansionType1CoeffGrads[colloc_index[p]];
      if (rand_0) { // integration: deriv of expect = expect of deriv
	if (li_0 == 0) // t1 weight is 1
	  for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	    deriv_index = dvv[d] - 1;
	    accum_0[d] = (data_rep->randomVarsKey[deriv_index]) ?
	      t1_coeff_grad[insert_cntr++] : t1_coeff;            // case 1:2
	  }
	else {
	  t1_wt_00 = t1_wts_0[key_p0];
	  if (augment) prod = t1_coeff * t1_wt_00;
	  for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	    deriv_index = dvv[d] - 1;
	    accum_0[d] += (data_rep->randomVarsKey[deriv_index]) ?
	      t1_coeff_grad[insert_cntr++] * t1_wt_00 : prod;    // case 1:2
	  }
	}
      }
      else { // interpolation: deriv of interpolation of mean
	if (li_0 == 0)  { // grad factor is 0., value factor is omitted
	  for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	    deriv_index = dvv[d] - 1;
	    if (data_rep->randomVarsKey[deriv_index])               // case 3
	      accum_0[d] = t1_coeff_grad[insert_cntr++];
	    else if (deriv_index)                                   // case 4
	      accum_0[d] = t1_coeff;
	  }
	}
	else {
	  bc_vf_00 = bc_vf_0[key_p0];
	  if (augment && ei_0 == _NPOS) prod = t1_coeff * bc_vf_00;
	  for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	    deriv_index = dvv[d] - 1;
	    if (data_rep->randomVarsKey[deriv_index])               // case 3
	      accum_0[d] += t1_coeff_grad[insert_cntr++] * bc_vf_00;
	    else if (deriv_index == 0)                              // case 4
	      accum_0[d] += t1_coeff * bc_gf_0[key_p0];
	    else if (ei_0 == _NPOS)                                 // case 4
	      accum_0[d] += prod;
	    else if (ei_0 == key_p0)                                // case 4
	      accum_0[d] += t1_coeff;
	  }
	}
      }
      if (key_p0 == max0)
	data_rep->accumulate_barycentric_gradient(accumulator, lev_index,
						  key_p, dvv);
    }
    Real scale = data_rep->
      barycentric_gradient_scaling(lev_index, data_rep->nonRandomIndices);
    Real* accum = accumulator[num_v-1];
    for (d=0; d<num_deriv_vars; ++d)
      tpMeanGrad[d] = accum[d] * scale;
  }
  else if (data_rep->basisConfigOptions.useDerivs) {
    if (insert) {
      PCerr << "Error: combination of coefficient gradients and "
            << "use_derivatives in NodalInterpPolyApproximation::"
	    << "tensor_product_mean_gradient()" << std::endl;
      abort_handler(-1);
    }

    //PCout << "tensor_product_mean_gradient(): Horner's useDerivs."<<std::endl;
    unsigned short key_p0, li_0 = lev_index[0];
    BasisPolynomial&      poly_0 = data_rep->polynomialBasis[li_0][0];
    const Real3DArray& t1_wts_1d = driver_rep->type1_collocation_weights_1d();
    const Real3DArray& t2_wts_1d = driver_rep->type2_collocation_weights_1d();
    const RealArray&   t1_wts_0  = t1_wts_1d[li_0][0];
    const RealArray&   t2_wts_0  = t2_wts_1d[li_0][0];
    unsigned short max0 = poly_0.interpolation_size() - 1;
    RealMatrix t1_accumulator(num_deriv_vars, num_v); // init to 0.
    RealMatrixArray t2_accumulators(num_deriv_vars);
    for (v=0; v<num_deriv_vars; ++v)
      t2_accumulators[v].shape(num_v, num_v); // init to 0.
    Real *t1_accum_0 = t1_accumulator[0], *t2_accum_0, t1_coeff, *t2_coeffs_p,
      t1_wt_00, t1_val, t2_val, t1_grad, prod1, prod2, x0 = x[0];
    bool rand_0 = data_rep->randomVarsKey[0];
    for (p=0; p<num_colloc_pts; ++p) {
      const UShortArray& key_p = key[p]; key_p0 = key_p[0];
      c_index     = (colloc_index.empty()) ? p : colloc_index[p];
      t1_coeff    = expansionType1Coeffs[c_index];
      t2_coeffs_p = expansionType2Coeffs[c_index];
      if (rand_0) { // integration: deriv of expect = expect of deriv
	if (li_0) {                                                    // case 2
	  t1_wt_00 = t1_wts_0[key_p0];
	  prod1 = t1_coeff       * t1_wt_00;
	  prod2 = t2_coeffs_p[0] * t2_wts_0[key_p0];
	  for (d=0; d<num_deriv_vars; ++d) {
	    t1_accum_0[d]            += prod1;
	    t2_accumulators[d][0][0] += prod2;
	  }
	  for (v=1; v<num_v; ++v) {
	    prod2 = t2_coeffs_p[v] * t1_wt_00;
	    for (d=0; d<num_deriv_vars; ++d)
	      t2_accumulators[d][0][v] += prod2;
	  }
	}
	else { // t1 weight is 1, t2 weight is 0
	  for (d=0; d<num_deriv_vars; ++d) {
	    t1_accum_0[d] = t1_coeff;                                  // case 2
	    t2_accum_0    = t2_accumulators[d][0];
	    //t2_accum[0] = 0.;
	    for (v=1; v<num_v; ++v)
	      t2_accum_0[v] = t2_coeffs_p[v];
	  }
	}
      }
      else { // interpolation: deriv of interpolation of mean
	if (li_0) {
	  t1_val = poly_0.type1_value(x0, key_p0); prod1 = t1_coeff * t1_val;
	  t2_val = poly_0.type2_value(x0, key_p0);
	  for (d=0; d<num_deriv_vars; ++d) {                           // case 4
	    t2_accum_0 = t2_accumulators[d][0];
	    if (dvv[d] == 1) {
	      t1_grad = poly_0.type1_gradient(x0, key_p0);
	      t1_accum_0[d]   += t1_coeff * t1_grad;
	      t2_accum_0[0]   += t2_coeffs_p[0] *
		poly_0.type2_gradient(x0, key_p0);
	      for (v=1; v<num_v; ++v)
		t2_accum_0[v] += t2_coeffs_p[v] * t1_grad;
	    }
	    else {
	      t1_accum_0[d]   += prod1;
	      t2_accum_0[0]   += t2_coeffs_p[0] * t2_val;
	      for (v=1; v<num_v; ++v)
		t2_accum_0[v] += t2_coeffs_p[v] * t1_val;
	    }
	  }
	}
	else { // t1 value is 1, t1 grad is 0., t2 grad is 1
	  t2_val = poly_0.type2_value(x0, key_p0);
	  for (d=0; d<num_deriv_vars; ++d) {                           // case 4
	    t2_accum_0 = t2_accumulators[d][0];
	    if (dvv[d] == 1) {
	      //t1_accum_0[d] = 0.;
	      t2_accum_0[0] = t2_coeffs_p[0];
	      //for (v=1; v<num_v; ++v)
	      //  t2_accum_0[v] = 0.;
	    }
	    else {
	      t1_accum_0[d]   = t1_coeff;
	      t2_accum_0[0]   = t2_coeffs_p[0] * t2_val;
	      for (v=1; v<num_v; ++v)
		t2_accum_0[v] = t2_coeffs_p[v];
	    }
	  }
	}
      }
      if (key_p0 == max0)
	data_rep->accumulate_horners_gradient(t1_accumulator, t2_accumulators,
					      lev_index, key_p, dvv, x);
    }
    Real *t1_accum = t1_accumulator[num_v-1], *t2_accum;
    for (d=0; d<num_deriv_vars; ++d) {
      tpMeanGrad[d] = t1_accum[d];
      t2_accum = t2_accumulators[d][num_v-1];
      for (v=0; v<num_v; ++v)
	tpMeanGrad[d] += t2_accum[v];
    }

    /*
    // Simpler but less efficient approach:
    //PCout << "tensor_product_mean_gradient(): Original useDerivs."<<std::endl;
    Real t1_coeff, *t2_coeff_p, t1_wt;
    for (p=0; p<num_colloc_pts; ++p) {
      const UShortArray& key_p = key[p];
      c_index = (colloc_index.empty()) ? p : colloc_index[p];
      t1_coeff   = expansionType1Coeffs[c_index];
      t2_coeff_p = expansionType2Coeffs[c_index];
      t1_wt = type1_weight(key_p, lev_index, data_rep->randomIndices);
      for (d=0; d<num_deriv_vars; ++d) {
	deriv_index = dvv[d] - 1; // OK since we are in an "All" view
	// -------------------------------------------------------------------
	// deriv of All var expansion w.r.t. nonrand var (design augmentation)
	// -------------------------------------------------------------------
	Real& tp_mean_grad_d = tpMeanGrad[d];
	tp_mean_grad_d += t1_coeff * t1_wt *
	  type1_interpolant_gradient(x, deriv_index, key_p, lev_index,
				     data_rep->nonRandomIndices);
	for (v=0; v<num_v; ++v)
	  tp_mean_grad_d += t2_coeff_p[v]
	    * type2_weight(v, key_p, lev_index, data_rep->randomIndices)
	    * type2_interpolant_gradient(x, deriv_index, v, key_p, lev_index,
	                                 data_rep->nonRandomIndices);
      }
    }
    */
  }
  else {
    //PCout << "tensor_product_mean_gradient(): Horner's no derivs."<<std::endl;
    // Horner's rule approach:
    unsigned short key_p0, li_0 = lev_index[0];
    BasisPolynomial&      poly_0 = data_rep->polynomialBasis[li_0][0];
    const Real3DArray& t1_wts_1d = driver_rep->type1_collocation_weights_1d();
    const RealArray&   t1_wts_0  = t1_wts_1d[li_0][0];
    unsigned short max0 = poly_0.interpolation_size() - 1;
    RealMatrix accumulator(num_deriv_vars, num_v); // init to 0.
    Real *accum_0 = accumulator[0], *t1_coeff_grad, t1_coeff, t1_wt_00,
      prod, t1_val, x0 = x[0];
    bool rand_0 = data_rep->randomVarsKey[0];
    for (p=0; p<num_colloc_pts; ++p) {
      const UShortArray& key_p = key[p]; key_p0 = key_p[0];
      if (augment)
	t1_coeff = (colloc_index.empty()) ? expansionType1Coeffs[p] :
	  expansionType1Coeffs[colloc_index[p]];
      if (insert)
	t1_coeff_grad = (colloc_index.empty()) ? expansionType1CoeffGrads[p] :
	  expansionType1CoeffGrads[colloc_index[p]];
      if (rand_0) { // integration: deriv of expect = expect of deriv
	if (li_0 == 0) // t1 weight is 1
	  for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	    deriv_index = dvv[d] - 1;
	    if (data_rep->randomVarsKey[deriv_index])
	      accum_0[d] = t1_coeff_grad[insert_cntr++];            // case 1
	    else
	      accum_0[d] = t1_coeff;                                // case 2
	  }
	else {
	  t1_wt_00 = t1_wts_0[key_p0];
	  if (augment) prod = t1_coeff * t1_wt_00;
	  for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	    deriv_index = dvv[d] - 1;
	    if (data_rep->randomVarsKey[deriv_index])
	      accum_0[d] += t1_coeff_grad[insert_cntr++] * t1_wt_00;// case 1
	    else
	      accum_0[d] += prod;                                   // case 2
	  }
	}
      }
      else { // interpolation: deriv of interpolation of mean
	if (li_0 == 0)  { // t1 grad is 0., t1 value is 1
	  for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	    deriv_index = dvv[d] - 1;
	    if (data_rep->randomVarsKey[deriv_index])               // case 3
	      accum_0[d] = t1_coeff_grad[insert_cntr++];
	    else if (deriv_index)                                   // case 4
	      accum_0[d] = t1_coeff;
	  }
	}
	else {
	  t1_val = poly_0.type1_value(x0, key_p0);
	  if (augment) prod = t1_coeff * t1_val;
	  for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	    deriv_index = dvv[d] - 1;
	    if (data_rep->randomVarsKey[deriv_index])               // case 3
	      accum_0[d] += t1_coeff_grad[insert_cntr++] * t1_val;
	    else if (deriv_index == 0)                              // case 4
	      accum_0[d] += t1_coeff * poly_0.type1_gradient(x0, key_p0);
	    else                                                    // case 4
	      accum_0[d] += prod;
	  }
	}
      }
      if (key_p0 == max0)
	data_rep->accumulate_horners_gradient(accumulator, lev_index, key_p,
					      dvv, x);
    }
    copy_data(accumulator[num_v-1], (int)num_deriv_vars, tpMeanGrad);

    /*
    // Simpler but less efficient approach:
    //PCout << "tensor_product_mean_gradient(): Original no derivs."<<std::endl;
    Real t1_coeff, *t1_coeff_grad, t1_wt, t1_val;
    for (p=0; p<num_colloc_pts; ++p) {
      const UShortArray& key_p = key[p];
      c_index = (colloc_index.empty()) ? p : colloc_index[p];
      t1_wt   = type1_weight(key_p, lev_index, data_rep->randomIndices);
      if (insert) {
	t1_coeff_grad = expansionType1CoeffGrads[c_index];
	t1_val = type1_interpolant_value(x, key_p, lev_index,
	                                 data_rep->nonRandomIndices);
      }
      else
	t1_coeff = expansionType1Coeffs[c_index];	
      for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	deriv_index = dvv[d] - 1; // OK since we are in an "All" view
	tpMeanGrad[d] += (data_rep->randomVarsKey[deriv_index]) ?
	  // ------------------------------------------------------------------
	  // derivative of All var expansion w.r.t. rand var (design insertion)
	  // ------------------------------------------------------------------
	  t1_coeff_grad[insert_cntr++] * t1_wt * t1_val :
	  // -------------------------------------------------------------------
	  // deriv of All var expansion w.r.t. nonrand var (design augmentation)
	  // -------------------------------------------------------------------
	  t1_coeff * t1_wt * type1_interpolant_gradient(x, deriv_index, key_p,
	    lev_index, data_rep->nonRandomIndices);
      }
    }
    */
  }

  return tpMeanGrad;
}


/** Covariance of response functions for a matched tensor product grid.
    Supports all_variables mode and either interpolation of products or
    product of interpolants formulations.  For the latter, recursive
    usage only provides the "diagonal" contributions from matched
    tensor products (PRODUCT_OF_INTERPOLANTS_FAST); an exact estimation
    (PRODUCT_OF_INTERPOLANTS_FULL) requires augmentation with mixed tensor
    product contributions using the overloaded form of this function. */
Real NodalInterpPolyApproximation::
tensor_product_covariance(const RealVector& x, const UShortArray& lev_index,
			  const UShort2DArray& key,
			  const SizetArray& colloc_index,
			  NodalInterpPolyApproximation* nip_approx_2)
{
  // Error check for required data
  if (!expansionCoeffFlag || !nip_approx_2->expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::tensor_product_covariance()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  bool same = (this == nip_approx_2);
  Real tp_covar = 0., t1_coeff_1_mm1, t1_coeff_2_mm2, mean_1, mean_2;
  if (data_rep->momentInterpType == PRODUCT_OF_INTERPOLANTS_FAST) {
    mean_1 = tensor_product_mean(x, lev_index, key, colloc_index);
    mean_2 = (same) ? mean_1 :
      nip_approx_2->tensor_product_mean(x, lev_index, key, colloc_index);
  }
  else                 // PRODUCT_OF_INTERPOLANTS_FULL
    { mean_1 = mean(x); mean_2 = (same) ? mean_1 : nip_approx_2->mean(x); }
  size_t i, c_index_i, num_colloc_pts = key.size(),
    num_v = sharedDataRep->numVars;

  switch (data_rep->momentInterpType) {
  case INTERPOLATION_OF_PRODUCTS: {
    const RealVector& t1_coeffs_2 = nip_approx_2->expansionType1Coeffs;
    const Real3DArray&  t1_wts_1d = driver_rep->type1_collocation_weights_1d();
    if (data_rep->barycentricFlag) {

      // For barycentric interpolation: track x != newPoint within 1D basis
      data_rep->set_new_point(x, lev_index, data_rep->nonRandomIndices, 1);

      // Since we need to include weights for integrated dimensions, loop over
      // all collocation pts even if some interpolated dimensions are inactive.
      RealVector accumulator(num_v); // init to 0.
      unsigned short       li_0 = lev_index[0];
      BasisPolynomial&   poly_0 = data_rep->polynomialBasis[li_0][0];
      const RealArray& t1_wts_0 = t1_wts_1d[li_0][0];
      const RealVector& bc_vf_0 = poly_0.barycentric_value_factors();
      size_t ei0 = poly_0.exact_index();
      bool rand_0 = data_rep->randomVarsKey[0];
      unsigned short key_i0, max0 = poly_0.interpolation_size() - 1;
      for (i=0; i<num_colloc_pts; ++i) {
	c_index_i = (colloc_index.empty()) ? i : colloc_index[i];
	t1_coeff_1_mm1 = expansionType1Coeffs[c_index_i] - mean_1;
	t1_coeff_2_mm2 = (same) ? t1_coeff_1_mm1 :
	  t1_coeffs_2[c_index_i] - mean_2;
	const UShortArray& key_i = key[i]; key_i0 = key_i[0];
	if (li_0 == 0)   // single integration/interpolation weight = 1
	  accumulator[0]  = t1_coeff_1_mm1 * t1_coeff_2_mm2;
	else if (rand_0) // integration
	  accumulator[0] += t1_coeff_1_mm1 * t1_coeff_2_mm2 * t1_wts_0[key_i0];
	else if (ei0 == key_i0) // interpolation w/ single pt or exact match
	  accumulator[0]  = t1_coeff_1_mm1 * t1_coeff_2_mm2;
	else if (ei0 == _NPOS)  // interpolation w/o exact match
	  accumulator[0] += t1_coeff_1_mm1 * t1_coeff_2_mm2 * bc_vf_0[key_i0];
	if (key_i0 == max0)
	  data_rep->accumulate_barycentric(accumulator, lev_index, key_i);
      }
      return accumulator[num_v-1] / data_rep->
	barycentric_value_scaling(lev_index, data_rep->nonRandomIndices);
    }
    else if (data_rep->basisConfigOptions.useDerivs) {
      // Horner's rule approach:
      // Note: while possible to use a single accumulator vector and loop over
      // it n+1 times, value() and gradient() data would then have to be stored.
      RealVector t1_accumulator(num_v);          // init to 0.
      RealMatrix t2_accumulator(num_v, num_v); // init to 0.
      Real *t2_accum_0 = t2_accumulator[0], t1_val, x0 = x[0];
      const RealMatrix& t2_coeffs_2 = nip_approx_2->expansionType2Coeffs;
      const Real3DArray& t2_wts_1d = driver_rep->type2_collocation_weights_1d();
      unsigned short       li_0 = lev_index[0];
      BasisPolynomial&   poly_0 = data_rep->polynomialBasis[li_0][0];
      const RealArray& t1_wts_0 = t1_wts_1d[li_0][0]; // li==li for integrated
      const RealArray& t2_wts_0 = t2_wts_1d[li_0][0]; // li==li for integrated
      size_t k; bool rand_0 = data_rep->randomVarsKey[0];
      unsigned short key_i0, max0 = poly_0.interpolation_size() - 1;
      for (i=0; i<num_colloc_pts; ++i) {
	c_index_i = (colloc_index.empty()) ? i : colloc_index[i];
	t1_coeff_1_mm1 = expansionType1Coeffs[c_index_i] - mean_1;
	t1_coeff_2_mm2 = (same) ? t1_coeff_1_mm1 :
	  t1_coeffs_2[c_index_i] - mean_2;
	const Real *t2_coeff_1 = expansionType2Coeffs[c_index_i],
	           *t2_coeff_2 = (same) ? t2_coeff_1 : t2_coeffs_2[c_index_i];
	const UShortArray& key_i = key[i]; key_i0 = key_i[0];
	if (rand_0) { // integration
	  if (li_0 == 0) {
	    t1_accumulator[0]  = t1_coeff_1_mm1 * t1_coeff_2_mm2; // t1 wt is 1
	    t2_accum_0[0]      = 0.;                              // t2 wt is 0
	    for (k=1; k<num_v; ++k)
	      t2_accum_0[k]    = t1_coeff_1_mm1 * t2_coeff_2[k]
		               + t1_coeff_2_mm2 * t2_coeff_1[k];  // t1 wt is 1
	  }
	  else {
	    t1_accumulator[0] += t1_coeff_1_mm1 * t1_coeff_2_mm2
	                      *  t1_wts_0[key_i0];
	    t2_accum_0[0]     += (t1_coeff_1_mm1 * t2_coeff_2[0] +
				  t1_coeff_2_mm2 * t2_coeff_1[0])
	                      *   t2_wts_0[key_i0];
	    for (k=1; k<num_v; ++k)
	      t2_accum_0[k]   += (t1_coeff_1_mm1 * t2_coeff_2[k] +
				  t1_coeff_2_mm2 * t2_coeff_1[k])
		              *   t1_wts_0[key_i0];
	  }
	}
	else {        // interpolation
	  if (li_0 == 0) {
	    t1_accumulator[0]  = t1_coeff_1_mm1 * t1_coeff_2_mm2; // t1 val is 1
	    t2_accum_0[0]      = (t1_coeff_1_mm1 * t2_coeff_2[0] +
				  t1_coeff_2_mm2 * t2_coeff_1[0])
	                       * poly_0.type2_value(x0, key_i0);
	    for (k=1; k<num_v; ++k)
	      t2_accum_0[k]    = t1_coeff_1_mm1 * t2_coeff_2[k]
		               + t1_coeff_2_mm2 * t2_coeff_1[k];  // t1 val is 1
	  }
	  else {
	    t1_val = poly_0.type1_value(x0, key_i0);
	    t1_accumulator[0] += t1_coeff_1_mm1 * t1_coeff_2_mm2 * t1_val;
	    t2_accum_0[0]     += (t1_coeff_1_mm1 * t2_coeff_2[0] +
				  t1_coeff_2_mm2 * t2_coeff_1[0])
	                      *  poly_0.type2_value(x0, key_i0);
	    for (k=1; k<num_v; ++k)
	      t2_accum_0[k]   += (t1_coeff_1_mm1 * t2_coeff_2[k] +
				  t1_coeff_2_mm2 * t2_coeff_1[k]) * t1_val;
	  }
	}
	if (key_i0 == max0)
	  data_rep->accumulate_horners(t1_accumulator, t2_accumulator,
				       lev_index, key_i, x);
      }
      Real  tp_covar = t1_accumulator[num_v-1];
      Real* t2_accum = t2_accumulator[num_v-1];
      for (k=0; k<num_v; ++k)
	tp_covar += t2_accum[k];
      return tp_covar;

      /*
      // Simpler but more expensive approach:
      Real tp_covar = 0.; size_t j;
      for (i=0; i<num_colloc_pts; ++i) {
	const UShortArray& key_i = key[i];
	t1_coeff_1_mm1 = expansionType1Coeffs[c_index_i] - mean_1;
	t1_coeff_2_mm2 = (same) ? t1_coeff_1_mm1 :
	  t1_coeffs_2[c_index_i] - mean_2;
	tp_covar += t1_coeff_1_mm1 * t1_coeff_2_mm2
	  * data_rep->type1_interpolant_value(x, key_i, lev_index,
				              data_rep->nonRandomIndices)
	  * data_rep->type1_weight(key_i, lev_index, data_rep->randomIndices);
	const Real *t2_coeff_1 = expansionType2Coeffs[c_index_i],
	           *t2_coeff_2 = (same) ? t2_coeff_1 : t2_coeffs_2[c_index_i];
	for (j=0; j<num_v; ++j)
	  tp_covar +=
	    (t1_coeff_1_mm1 * t2_coeff_2[j] + t1_coeff_2_mm2 * t2_coeff_1[j])
	    * data_rep->type2_interpolant_value(x, j, key_i, lev_index,
				                data_rep->nonRandomIndices)
	    * data_rep->type2_weight(j, key_i, lev_index,
	                             data_rep->randomIndices);
      }
      return tp_covar;
      */
    }
    else {
      // Horner's rule approach:
      RealVector accumulator(num_v); // init to 0.
      unsigned short       li_0 = lev_index[0];
      BasisPolynomial&   poly_0 = data_rep->polynomialBasis[li_0][0];
      const RealArray& t1_wts_0 = t1_wts_1d[li_0][0]; // li==li for integrated
      unsigned short key_i0, max0 = poly_0.interpolation_size() - 1;
      bool rand_0 = data_rep->randomVarsKey[0]; Real x0 = x[0];
      for (i=0; i<num_colloc_pts; ++i) {
	c_index_i = (colloc_index.empty()) ? i : colloc_index[i];
	t1_coeff_1_mm1 = expansionType1Coeffs[c_index_i] - mean_1;
	t1_coeff_2_mm2 = (same) ? t1_coeff_1_mm1 :
	  t1_coeffs_2[c_index_i] - mean_2;
	const UShortArray& key_i = key[i]; key_i0 = key_i[0];
	if (li_0 == 0)   // single integration/interpolation weight = 1
	  accumulator[0]  = t1_coeff_1_mm1 * t1_coeff_2_mm2;
	else if (rand_0) // integration
	  accumulator[0] += t1_coeff_1_mm1 * t1_coeff_2_mm2 * t1_wts_0[key_i0];
	else             // interpolation
	  accumulator[0] += t1_coeff_1_mm1 * t1_coeff_2_mm2 *
	    poly_0.type1_value(x0, key_i0);
	if (key_i0 == max0)
	  data_rep->accumulate_horners(accumulator, lev_index, key_i, x);
      }
      return accumulator[num_v-1];

      /*
      // Simpler but more expensive approach:
      Real tp_covar = 0.;
      for (i=0; i<num_colloc_pts; ++i) {
	const UShortArray& key_i = key[i];
	c_index_i = (colloc_index.empty()) ? i : colloc_index[i];
	t1_coeff_1_mm1 = expansionType1Coeffs[c_index_i] - mean_1;
	t1_coeff_2_mm2 = (same) ? t1_coeff_1_mm1 :
	  t1_coeffs_2[c_index_i] - mean_2;
	tp_covar += t1_coeff_1_mm1 * t1_coeff_2_mm2
	  * data_rep->type1_interpolant_value(x, key_i, lev_index,
				              data_rep->nonRandomIndices)
	  * data_rep->type1_weight(key_i, lev_index, data_rep->randomIndices);
      }
      return tp_covar;
      */
    }
    break;
  }
  case REINTERPOLATION_OF_PRODUCTS: {
    // For this case, integration can occur on the original grid, but the
    // interpolation portion needs to use a higher order grid to resolve
    // (R-\mu)^2 to comparable accuracy.  This is "bootstrapped" by using
    // the interpolant R-hat from the original grid to evaluate (R-\mu)^2 
    // on the higher order grid.
    const UShortArray&   reinterp_lev_index
      = driver_rep->reinterpolated_level_index();
    const RealMatrix&    reinterp_var_sets
      = driver_rep->reinterpolated_variable_sets();
    const UShort2DArray& reinterp_colloc_key
      = driver_rep->reinterpolated_collocation_key();

    // integrate/interpolate over the new grid
    size_t reinterp_colloc_pts = reinterp_colloc_key.size();
    if (data_rep->barycentricFlag) {

      // For barycentric interpolation: track x != newPoint within 1D basis
      data_rep->set_new_point(x, reinterp_lev_index,
			      data_rep->nonRandomIndices, 1);

      // Since we need to include weights for integrated dimensions, loop over
      // all collocation pts even if some interpolated dimensions are inactive.
      RealVector accumulator(num_v); // init to 0.
      const Real3DArray& t1_wts_1d = driver_rep->type1_collocation_weights_1d();
      unsigned short      rli_0 = reinterp_lev_index[0];
      BasisPolynomial&   poly_0 = data_rep->polynomialBasis[rli_0][0];
      const RealArray& t1_wts_0 = t1_wts_1d[rli_0][0]; // rli==li for integrated
      const RealVector& bc_vf_0 = poly_0.barycentric_value_factors();
      size_t ei0 = poly_0.exact_index();
      bool rand_0 = data_rep->randomVarsKey[0];
      unsigned short key_i0, max0 = poly_0.interpolation_size() - 1;
      for (i=0; i<reinterp_colloc_pts; ++i) {
	RealVector c_vars(Teuchos::View,
			  const_cast<Real*>(reinterp_var_sets[i]), num_v);
	t1_coeff_1_mm1 = value(c_vars) - mean_1;
	t1_coeff_2_mm2 = (same) ? t1_coeff_1_mm1 :
	  nip_approx_2->value(c_vars) - mean_2;
	const UShortArray& key_i = reinterp_colloc_key[i]; key_i0 = key_i[0];
	if (rli_0 == 0)   // single integration/interpolation weight = 1
	  accumulator[0]  = t1_coeff_1_mm1 * t1_coeff_2_mm2;
	else if (rand_0) // integration
	  accumulator[0] += t1_coeff_1_mm1 * t1_coeff_2_mm2 * t1_wts_0[key_i0];
	else if (ei0 == key_i0) // interpolation w/ single pt or exact match
	  accumulator[0]  = t1_coeff_1_mm1 * t1_coeff_2_mm2;
	else if (ei0 == _NPOS)  // interpolation w/o exact match
	  accumulator[0] += t1_coeff_1_mm1 * t1_coeff_2_mm2 * bc_vf_0[key_i0];
	if (key_i0 == max0)
	  data_rep->accumulate_barycentric(accumulator, reinterp_lev_index,
					   key_i);
      }
      return accumulator[num_v-1] / data_rep->
	barycentric_value_scaling(lev_index, data_rep->nonRandomIndices);
    }
    else if (data_rep->basisConfigOptions.useDerivs) {
      // Horner's rule approach:
      // Note: while possible to use a single accumulator vector and loop over
      // it n+1 times, value() and gradient() data would then have to be stored.
      RealVector t1_accumulator(num_v);          // init to 0.
      RealMatrix t2_accumulator(num_v, num_v); // init to 0.
      Real *t2_accum_0 = t2_accumulator[0], t1_val, x0 = x[0];
      const Real3DArray& t1_wts_1d = driver_rep->type1_collocation_weights_1d();
      const Real3DArray& t2_wts_1d = driver_rep->type2_collocation_weights_1d();
      unsigned short      rli_0 = reinterp_lev_index[0];
      BasisPolynomial&   poly_0 = data_rep->polynomialBasis[rli_0][0];
      const RealArray& t1_wts_0 = t1_wts_1d[rli_0][0]; // rli==li for integrated
      const RealArray& t2_wts_0 = t2_wts_1d[rli_0][0]; // rli==li for integrated
      size_t k; bool rand_0 = data_rep->randomVarsKey[0];
      unsigned short key_i0, max0 = poly_0.interpolation_size() - 1;
      for (i=0; i<reinterp_colloc_pts; ++i) {
	RealVector c_vars(Teuchos::View,
			  const_cast<Real*>(reinterp_var_sets[i]), num_v);
	t1_coeff_1_mm1 = value(c_vars) - mean_1;
	t1_coeff_2_mm2 = (same) ? t1_coeff_1_mm1 :
	  nip_approx_2->value(c_vars) - mean_2;
	const RealVector& t2_coeff_1 = gradient_basis_variables(c_vars);
	const RealVector& t2_coeff_2 = (same) ? t2_coeff_1 :
	  nip_approx_2->gradient_basis_variables(c_vars);
	const UShortArray& key_i = reinterp_colloc_key[i]; key_i0 = key_i[0];
	if (rand_0) { // integration
	  if (rli_0 == 0) {
	    t1_accumulator[0]  = t1_coeff_1_mm1 * t1_coeff_2_mm2; // t1 wt is 1
	    t2_accum_0[0]      = 0.;                              // t2 wt is 0
	    for (k=1; k<num_v; ++k)
	      t2_accum_0[k]    = t1_coeff_1_mm1 * t2_coeff_2[k]
		               + t1_coeff_2_mm2 * t2_coeff_1[k];  // t1 wt is 1
	  }
	  else {
	    t1_accumulator[0] += t1_coeff_1_mm1 * t1_coeff_2_mm2
	                      *  t1_wts_0[key_i0];
	    t2_accum_0[0]     += (t1_coeff_1_mm1 * t2_coeff_2[0] +
				  t1_coeff_2_mm2 * t2_coeff_1[0])
	                      *   t2_wts_0[key_i0];
	    for (k=1; k<num_v; ++k)
	      t2_accum_0[k]   += (t1_coeff_1_mm1 * t2_coeff_2[k] +
				  t1_coeff_2_mm2 * t2_coeff_1[k])
		              *   t1_wts_0[key_i0];
	  }
	}
	else {        // interpolation
	  if (rli_0 == 0) {
	    t1_accumulator[0]  = t1_coeff_1_mm1 * t1_coeff_2_mm2; // t1 val is 1
	    t2_accum_0[0]      = (t1_coeff_1_mm1 * t2_coeff_2[0] +
				  t1_coeff_2_mm2 * t2_coeff_1[0])
	                       * poly_0.type2_value(x0, key_i0);
	    for (k=1; k<num_v; ++k)
	      t2_accum_0[k]    = t1_coeff_1_mm1 * t2_coeff_2[k]
		               + t1_coeff_2_mm2 * t2_coeff_1[k];  // t1 val is 1
	  }
	  else {
	    t1_val = poly_0.type1_value(x0, key_i0);
	    t1_accumulator[0] += t1_coeff_1_mm1 * t1_coeff_2_mm2 * t1_val;
	    t2_accum_0[0]     += (t1_coeff_1_mm1 * t2_coeff_2[0] +
				  t1_coeff_2_mm2 * t2_coeff_1[0])
	                      *  poly_0.type2_value(x0, key_i0);
	    for (k=1; k<num_v; ++k)
	      t2_accum_0[k]   += (t1_coeff_1_mm1 * t2_coeff_2[k] +
				  t1_coeff_2_mm2 * t2_coeff_1[k]) * t1_val;
	  }
	}
	if (key_i0 == max0)
	  data_rep->accumulate_horners(t1_accumulator, t2_accumulator,
				       reinterp_lev_index, key_i, x);
      }
      Real  tp_covar = t1_accumulator[num_v-1];
      Real* t2_accum = t2_accumulator[num_v-1];
      for (k=0; k<num_v; ++k)
	tp_covar += t2_accum[k];
      return tp_covar;

      /*
      // Simpler but more expensive approach:
      Real tp_covar = 0.; size_t j;
      for (i=0; i<reinterp_colloc_pts; ++i) {
	const UShortArray& key_i = reinterp_colloc_key[i];
	RealVector c_vars(Teuchos::View,
			  const_cast<Real*>(reinterp_var_sets[i]), num_v);
	t1_coeff_1_mm1 = value(c_vars) - mean_1; // tensor_product_value() ?
	t1_coeff_2_mm2 = (same) ? t1_coeff_1_mm1 :
	  nip_approx_2->value(c_vars) - mean_2;
	tp_covar += t1_coeff_1_mm1 * t1_coeff_2_mm2
	  * data_rep->type1_interpolant_value(x, key_i, reinterp_lev_index,
				              data_rep->nonRandomIndices)
	  * data_rep->type1_weight(key_i, reinterp_lev_index,
	                           data_rep->randomIndices);
	const RealVector& t2_coeff_1 = gradient_basis_variables(c_vars);
	const RealVector& t2_coeff_2 = (same) ? t2_coeff_1 :
	  nip_approx_2->gradient_basis_variables(c_vars);
	for (j=0; j<num_v; ++j)
	  tp_covar +=
	    (t1_coeff_1_mm1 * t2_coeff_2[j] + t1_coeff_2_mm2 * t2_coeff_1[j])
	    * data_rep->type2_interpolant_value(x, j, key_i, reinterp_lev_index,
				                data_rep->nonRandomIndices)
	    * data_rep->type2_weight(j, key_i, reinterp_lev_index,
	                             data_rep->randomIndices);
      }
      return tp_covar;
      */
    }
    else {
      // Horner's rule approach:
      RealVector accumulator(num_v); // init to 0.
      const Real3DArray& t1_wts_1d = driver_rep->type1_collocation_weights_1d();
      unsigned short      rli_0 = reinterp_lev_index[0];
      BasisPolynomial&   poly_0 = data_rep->polynomialBasis[rli_0][0];
      const RealArray& t1_wts_0 = t1_wts_1d[rli_0][0]; // rli==li for integrated
      unsigned short key_i0, max0 = poly_0.interpolation_size() - 1;
      bool rand_0 = data_rep->randomVarsKey[0]; Real x0 = x[0];
      for (i=0; i<reinterp_colloc_pts; ++i) {
	RealVector c_vars(Teuchos::View,
			  const_cast<Real*>(reinterp_var_sets[i]), num_v);
	t1_coeff_1_mm1 = value(c_vars) - mean_1;
	t1_coeff_2_mm2 = (same) ? t1_coeff_1_mm1 :
	  nip_approx_2->value(c_vars) - mean_2;
	const UShortArray& key_i = reinterp_colloc_key[i]; key_i0 = key_i[0];
	if (rli_0 == 0)   // single integration/interpolation weight = 1
	  accumulator[0]  = t1_coeff_1_mm1 * t1_coeff_2_mm2;
	else if (rand_0) // integration
	  accumulator[0] += t1_coeff_1_mm1 * t1_coeff_2_mm2 * t1_wts_0[key_i0];
	else             // interpolation
	  accumulator[0] += t1_coeff_1_mm1 * t1_coeff_2_mm2 *
	    poly_0.type1_value(x0, key_i0);
	if (key_i0 == max0)
	  data_rep->
	    accumulate_horners(accumulator, reinterp_lev_index, key_i, x);
      }
      return accumulator[num_v-1];

      /*
      // Simpler but more expensive approach:
      Real tp_covar = 0.;
      for (i=0; i<reinterp_colloc_pts; ++i) {
	const UShortArray& key_i = reinterp_colloc_key[i];
	RealVector c_vars(Teuchos::View,
			  const_cast<Real*>(reinterp_var_sets[i]), num_v);
	t1_coeff_1_mm1 = value(c_vars) - mean_1; // tensor_product_value() ?
	t1_coeff_2_mm2 = (same) ? t1_coeff_1_mm1 :
	  nip_approx_2->value(c_vars) - mean_2;
	tp_covar += t1_coeff_1_mm1 * t1_coeff_2_mm2
	  * data_rep->type1_interpolant_value(x, key_i, reinterp_lev_index,
				              data_rep->nonRandomIndices)
	  * data_rep->type1_weight(key_i, reinterp_lev_index,
	                           data_rep->randomIndices);
      }
      return tp_covar;
      */
    }
    break;
  }
  case PRODUCT_OF_INTERPOLANTS_FAST: case PRODUCT_OF_INTERPOLANTS_FULL: {
    size_t j, c_index_j; Real t1_wt_Ls_prod_i;
    const RealVector& t1_coeffs_2 = nip_approx_2->expansionType1Coeffs;
    const RealMatrix& t2_coeffs_2 = nip_approx_2->expansionType2Coeffs;
    for (i=0; i<num_colloc_pts; ++i) {
      const UShortArray& key_i = key[i];
      c_index_i = (colloc_index.empty()) ? i : colloc_index[i];
      t1_coeff_1_mm1  = expansionType1Coeffs[c_index_i] - mean_1;
      t1_wt_Ls_prod_i
	= data_rep->type1_weight(key_i, lev_index, data_rep->randomIndices)
	* data_rep->type1_interpolant_value(x, key_i, lev_index,
					    data_rep->nonRandomIndices);
      for (j=0; j<num_colloc_pts; ++j) {
	const UShortArray& key_j = key[j];
	// to include the ij-th term,  basis i must be the same as basis j for
	// the random var subset.  In this case, wt_prod_i may be reused.  Note
	// that it is not necessary to collapse terms with the same random basis
	// subset, since cross term in (a+b)(a+b) = a^2+2ab+b^2 gets included.
	// If terms were collapsed (following eval of non-random portions), the
	// nested loop could be replaced with a single loop to evaluate (a+b)^2.
	if (data_rep->match_random_key(key_i, key_j)) {
	  c_index_j = (colloc_index.empty()) ? j : colloc_index[j];
	  t1_coeff_2_mm2 = t1_coeffs_2[c_index_j] - mean_2;
	  tp_covar += t1_coeff_1_mm1 * t1_coeff_2_mm2 * t1_wt_Ls_prod_i *
	    data_rep->type1_interpolant_value(x, key_j, lev_index,
					      data_rep->nonRandomIndices);
	/* TO DO
	if (data_rep->basisConfigOptions.useDerivs) {
	  const Real *t2_coeff_1i = expansionType2Coeffs[c_index_i],
	             *t2_coeff_2i = t2_coeffs_2[c_index_i];
	  for (j=0; j<num_v; ++j)
	    tp_covar += (t1_coeff_1i_mm1 * t2_coeff_2i[j] +
	                 t1_coeff_2i_mm2 * t2_coeff_1i[j])
	      * data_rep->type2_interpolant_value(x, j, key_i, lev_index,
	                                          data_rep->nonRandomIndices)
	      * data_rep->type2_weight(j, key_i, lev_index,
	                               data_rep->randomIndices);
	}
	*/
	}
      }
    }
    break;
  }
  }

  return tp_covar;
}


/** Covariance of response functions for differing tensor products
    using a PRODUCT_OF_INTERPOLANTS_FULL approach.  Needed for sparse
    interpolants in all_variables mode. */
Real NodalInterpPolyApproximation::
tensor_product_covariance(const RealVector& x, const UShortArray& lev_index_1,
			  const UShort2DArray& key_1,
			  const SizetArray& colloc_index_1,
			  const UShortArray& lev_index_2,
			  const UShort2DArray& key_2,
			  const SizetArray& colloc_index_2,
			  NodalInterpPolyApproximation* nip_approx_2)
{
#ifdef DEBUG
  PCout << "tensor_product_covariance for lev_index_1 =\n" << lev_index_1
	<< "lev_index_2 =\n" << lev_index_2 << "key_1 =\n" << key_1
	<< "key_2 =\n" << key_2 << std::endl;
#endif // DEBUG

  // Error check for required data
  if (!expansionCoeffFlag || !nip_approx_2->expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::tensor_product_covariance()" << std::endl;
    abort_handler(-1);
  }
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  if (data_rep->momentInterpType != PRODUCT_OF_INTERPOLANTS_FULL) {
    PCerr << "Error: mixed tensor product covariance only required for full "
	  << "products of interpolants. " << std::endl;
    abort_handler(-1);
  }

  bool same = (this == nip_approx_2);
  size_t i, j, k, c_index_1i, c_index_2j, num_pts_1 = key_1.size(),
    num_pts_2 = key_2.size();
  const RealVector& t1_coeffs_2 = nip_approx_2->expansionType1Coeffs;
  const RealMatrix& t2_coeffs_2 = nip_approx_2->expansionType2Coeffs;
  Real tp_covar = 0., t1_coeff_1i_mm1, t1_coeff_2j_mm2, t1_Ls_1i,
    basis_prod, basis_prod_k, mean_1 = mean(x),
    mean_2 = (same) ? mean_1 : nip_approx_2->mean(x);

  bool non_zero; unsigned short l1k, l2k; SizetList::iterator it;
  for (i=0; i<num_pts_1; ++i) {
    const UShortArray& key_1i = key_1[i];
    c_index_1i = (colloc_index_1.empty()) ? i : colloc_index_1[i];
    t1_coeff_1i_mm1 = expansionType1Coeffs[c_index_1i] - mean_1;
    t1_Ls_1i = data_rep->type1_interpolant_value(x, key_1i, lev_index_1,
						 data_rep->nonRandomIndices);
    for (j=0; j<num_pts_2; ++j) {
      const UShortArray& key_2j = key_2[j];
      // to include ij-th term, expectation of the product of the interpolation
      // polynomials must be nonzero
      if (data_rep->basis_product(lev_index_1, key_1i,
				  lev_index_2, key_2j, basis_prod)) {
	c_index_2j = (colloc_index_2.empty()) ? j : colloc_index_2[j];
	t1_coeff_2j_mm2 = t1_coeffs_2[c_index_2j] - mean_2;
	tp_covar += basis_prod * t1_coeff_1i_mm1 * t1_coeff_2j_mm2 * t1_Ls_1i *
	  data_rep->type1_interpolant_value(x, key_2j, lev_index_2,
					    data_rep->nonRandomIndices);
	/* TO DO
	if (data_rep->basisConfigOptions.useDerivs) {
	  const Real *t2_coeff_1i = expansionType2Coeffs[c_index_i],
	             *t2_coeff_2i = t2_coeffs_2[c_index_i];
	  for (j=0; j<num_v; ++j)
	    tp_covar += (t1_coeff_1i_mm1 * t2_coeff_2i[j] +
	                 t1_coeff_2i_mm2 * t2_coeff_1i[j])
	      * data_rep->type2_interpolant_value(x, j, key_i, lev_index,
	                                          data_rep->nonRandomIndices)
	      * data_rep->type2_weight(j, key_i, lev_index,
	                               data_rep->randomIndices);
	}
	*/
      }
    }
  }

  return tp_covar;
}


/** Overloaded all_variables version supporting interpolation of
    products or product of interpolants approaches. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_variance_gradient(const RealVector& x,
				 const UShortArray& lev_index,
				 const UShort2DArray& key,
				 const SizetArray& colloc_index,
				 const SizetArray& dvv)
{
  size_t d, insert_cntr, deriv_index, num_deriv_vars = dvv.size(),
    num_v = sharedDataRep->numVars;
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  if (tpVarianceGrad.length() != num_deriv_vars)
    tpVarianceGrad.sizeUninitialized(num_deriv_vars);
  tpVarianceGrad = 0.;
  Real mean_1 = (data_rep->momentInterpType == PRODUCT_OF_INTERPOLANTS_FAST) ?
    tensor_product_mean(x, lev_index, key, colloc_index) : mean(x);
  const RealVector& mean_grad =
    (data_rep->momentInterpType == PRODUCT_OF_INTERPOLANTS_FAST) ?
    tensor_product_mean_gradient(x, lev_index, key, colloc_index, dvv) :
    mean_gradient(x, dvv);

  // screen for insertion and augmentation and perform error checks
  bool insert = false, augment = false;
  for (d=0; d<num_deriv_vars; ++d) {
    deriv_index = dvv[d] - 1; // OK since we are in an "All" view
    if (data_rep->randomVarsKey[deriv_index]) insert  = true;
    else                                      augment = true;
  }
  if (insert && !expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	  << "InterpPolyApproximation::tensor_product_variance_gradient()."
	  << std::endl;
    abort_handler(-1);
  }
  if (augment && !expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::tensor_product_variance_gradient()" << std::endl;
    abort_handler(-1);
  }

  // compute variance gradient, managing mixture of insert/augment modes
  switch (data_rep->momentInterpType) {
  case INTERPOLATION_OF_PRODUCTS: {
    // -------------------------------------------------------------------
    // Mixed variable key:
    //   xi = ran vars, sa = augmented des vars, si = inserted design vars
    //   This assumes that the basis is fixed (ignores parameterized bases).
    // Active variable expansion (interpolation of products):
    //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
    //   covar(sa, si) = Sum_i r1_i(sa, si) r2_i(sa, si) wt_prod_xi_i - mu^2
    //   dcovar/ds     = Sum_i (r1_i dr2_i/ds + dr1_i/ds r2_i) wt_prod_xi_i
    //                 - 2 mu dmu/ds
    // All variable expansion (interpolation of products):
    //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
    //   covar(sa, si) = Sum_i r1_i(si) r2_i(si) Lsa_i wt_prod_xi_i - mu^2
    //   dcovar/dsa    = Sum_i r1_i(si) r2_i(si) dLsa_i/dsa wt_prod_xi_i
    //                 - 2 mu dmu/dsa
    //   dcovar/dsi    = Sum_i (r1_i dr2_j/dsi + dr1_i/dsi r2_j) Lsa_i
    //                   wt_prod_xi_i - 2 mu dmu/dsi
    // -------------------------------------------------------------------

    // Four cases to manage in barycentric/Horner looping:
    // 1:    random exp var,  inserted deriv: variance grad * type1_weight
    // 2:    random exp var, augmented deriv: variance      * type1_weight
    // 3: nonrandom exp var,  inserted deriv: variance grad * type1_value
    // 4: nonrandom exp var, augmented deriv: variance      * type1_gradient

    const Real3DArray& t1_wts_1d = driver_rep->type1_collocation_weights_1d();
    unsigned short  key_p0, li_0 = lev_index[0];
    BasisPolynomial&      poly_0 = data_rep->polynomialBasis[li_0][0];
    const RealArray&    t1_wts_0 = t1_wts_1d[li_0][0];
    size_t    p, c_index_p, ei_0 = poly_0.exact_index(),
                  num_colloc_pts = key.size();
    unsigned short         max_0 = poly_0.interpolation_size() - 1;
    bool                  rand_0 = data_rep->randomVarsKey[0];

    if (data_rep->barycentricFlag) {
      // For barycentric interpolation: track x != newPoint within 1D basis
      short order = (augment) ? 3 : 1;
      data_rep->set_new_point(x, lev_index, data_rep->nonRandomIndices, order);

      const RealVector&   bc_vf_0 = poly_0.barycentric_value_factors();
      const RealVector&   bc_gf_0 = poly_0.barycentric_gradient_factors();
      RealMatrix accumulator(num_deriv_vars, num_v); // init to 0.
      Real *accum_0 = accumulator[0], t1_coeff_p_mm, t1_wt_00, bc_vf_00,
	prod1, prod2, prod3;
      for (p=0; p<num_colloc_pts; ++p) {
	const UShortArray& key_p = key[p]; key_p0 = key_p[0];
	c_index_p = (colloc_index.empty()) ? p : colloc_index[p];
	t1_coeff_p_mm = expansionType1Coeffs[c_index_p] - mean_1;
	const Real* t1_coeff_grad_p = (insert) ?
	  expansionType1CoeffGrads[c_index_p] : NULL;
	if (rand_0) { // integration: deriv of expect = expect of deriv
	  if (li_0) {
	    t1_wt_00 = t1_wts_0[key_p0];
	    if (insert)  prod1 = 2. * t1_coeff_p_mm * t1_wt_00;
	    if (augment) prod2 = t1_coeff_p_mm * t1_coeff_p_mm * t1_wt_00;
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {       // case 1:2
	      deriv_index = dvv[d] - 1;
	      accum_0[d] += (data_rep->randomVarsKey[deriv_index]) ? prod1 *
		(t1_coeff_grad_p[insert_cntr++] - mean_grad[d]) : prod2;
	    }
	  }
	  else {
	    if (insert)  prod1 = 2. * t1_coeff_p_mm;
	    if (augment) prod2 = t1_coeff_p_mm * t1_coeff_p_mm;
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {       // case 1:2
	      deriv_index = dvv[d] - 1;
	      accum_0[d]  = (data_rep->randomVarsKey[deriv_index]) ? prod1 *
		(t1_coeff_grad_p[insert_cntr++] - mean_grad[d]) : prod2;
	    }
	  }
	}
	else { // interpolation: deriv of interpolation of mean
	  if (li_0)  { // grad factor is 0., value factor is omitted
	    bc_vf_00 = bc_vf_0[key_p0];
	    if (insert) prod1  = 2. * t1_coeff_p_mm * bc_vf_00;
	    if (augment) {
	      prod2 = t1_coeff_p_mm * t1_coeff_p_mm;
	      if (ei_0 == _NPOS) prod3 = prod2 * bc_vf_00;
	    }
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	      deriv_index = dvv[d] - 1;
	      if (data_rep->randomVarsKey[deriv_index])               // case 3
		accum_0[d] += prod1 *
		  (t1_coeff_grad_p[insert_cntr++] - mean_grad[d]);
	      else if (deriv_index == 0)                              // case 4
		accum_0[d] += prod2 * bc_gf_0[key_p0];
	      else if (ei_0 == _NPOS)                                 // case 4
		accum_0[d] += prod3;
	      else if (ei_0 == key_p0)                                // case 4
		accum_0[d] += prod2;
	    }
	  }
	  else {
	    if (insert)  prod1 = 2. * t1_coeff_p_mm;
	    if (augment) prod2 = t1_coeff_p_mm * t1_coeff_p_mm;
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	      deriv_index = dvv[d] - 1;
	      if (data_rep->randomVarsKey[deriv_index])               // case 3
		accum_0[d] = prod1 *
		  (t1_coeff_grad_p[insert_cntr++] - mean_grad[d]);
	      else if (deriv_index)                                   // case 4
		accum_0[d] = prod2;
	    }
	  }
	}
	if (key_p0 == max_0) // accumulate sums over vars with max key value
	  data_rep->accumulate_barycentric_gradient(accumulator, lev_index,
						    key_p, dvv);
      }
      Real scale = data_rep->
	barycentric_gradient_scaling(lev_index, data_rep->nonRandomIndices);
      Real* accum = accumulator[num_v-1];
      for (d=0; d<num_deriv_vars; ++d)
	tpVarianceGrad[d] = accum[d] * scale;
    }
    else if (data_rep->basisConfigOptions.useDerivs) {
      if (insert) {
	PCerr << "Error: combination of coefficient gradients and "
	      << "use_derivatives in NodalInterpPolyApproximation::"
	      << "tensor_product_mean_gradient()" << std::endl;
	abort_handler(-1);
      }

      //PCout << "tensor_product_variance_gradient(): Horner's useDerivs."
      //      << std::endl;
      const Real3DArray& t2_wts_1d = driver_rep->type2_collocation_weights_1d();
      const RealArray&   t2_wts_0  = t2_wts_1d[li_0][0];
      RealMatrix t1_accumulator(num_deriv_vars, num_v); // init to 0.
      RealMatrixArray t2_accumulators(num_deriv_vars); size_t v;
      for (v=0; v<num_deriv_vars; ++v)
	t2_accumulators[v].shape(num_v, num_v); // init to 0.
      Real *t1_accum_0 = t1_accumulator[0], *t2_accum_0, t1_coeff_p_mm,
	t1_wt_00, t1_val, t2_val, t1_grad, prod1, prod2, prod3, prod4,
	prod5, prod6, x0 = x[0];
      for (p=0; p<num_colloc_pts; ++p) {
	const UShortArray& key_p = key[p]; key_p0 = key_p[0];
	c_index_p = (colloc_index.empty()) ? p : colloc_index[p];
	t1_coeff_p_mm = expansionType1Coeffs[c_index_p] - mean_1;
	// no DVV: need all grad components for R to expand type2 for (R-mu)^2
	const Real* t2_coeff_p = expansionType2Coeffs[c_index_p];
	prod1 = t1_coeff_p_mm * t1_coeff_p_mm;
	prod2 = 2. * t1_coeff_p_mm;
	if (rand_0) { // integration: deriv of expect = expect of deriv
	  if (li_0) {                                                 // case 2
	    t1_wt_00 = t1_wts_0[key_p0];
	    prod3 = prod1 * t1_wt_00;
	    prod4 = prod2 * t2_coeff_p[0] * t2_wts_0[key_p0];
	    for (d=0; d<num_deriv_vars; ++d) {
	      t1_accum_0[d]            += prod3;
	      t2_accumulators[d][0][0] += prod4;
	    }
	    prod3 = prod2 * t1_wt_00;
	    for (v=1; v<num_v; ++v) {
	      prod4 = prod3 * t2_coeff_p[v];
	      for (d=0; d<num_deriv_vars; ++d)
		t2_accumulators[d][0][v] += prod4;
	    }
	  }
	  else { // t1 weight is 1, t2 weight is 0
	    for (d=0; d<num_deriv_vars; ++d)
	      t1_accum_0[d] = prod1;
	    for (v=1; v<num_v; ++v) {
	      prod4 = prod2 * t2_coeff_p[v];
	      for (d=0; d<num_deriv_vars; ++d)
		t2_accumulators[d][0][v] += prod4;
	    }
	  }
	}
	else if (li_0) { // interpolation: deriv of interpolation of mean
	  t1_val = poly_0.type1_value(x0, key_p0);
	  t2_val = poly_0.type2_value(x0, key_p0);
	  prod3 = prod1 * t1_val; prod4  = prod2 * t1_val;
	  prod5 = prod2 * t2_val;
	  for (d=0; d<num_deriv_vars; ++d) {                         // case 4
	    t2_accum_0 = t2_accumulators[d][0];
	    if (dvv[d] == 1) {
	      t1_grad = poly_0.type1_gradient(x0, key_p0);
	      t1_accum_0[d]   += prod1 * t1_grad;
	      t2_accum_0[0]   += prod2 * t2_coeff_p[0] *
		poly_0.type2_gradient(x0, key_p0);
	      prod6 = prod2 * t1_grad;
	      for (v=1; v<num_v; ++v)
		t2_accum_0[v] += prod6 * t2_coeff_p[v];
	    }
	    else {
	      t1_accum_0[d]   += prod3;
	      t2_accum_0[0]   += prod5 * t2_coeff_p[0];
	      for (v=1; v<num_v; ++v)
		t2_accum_0[v] += prod4 * t2_coeff_p[v];
	    }
	  }
	}
	else { // t1 value is 1, t1 grad is 0., t2 grad is 1
	  prod5 = prod2 * poly_0.type2_value(x0, key_p0);
	  for (d=0; d<num_deriv_vars; ++d) {                         // case 4
	    t2_accum_0 = t2_accumulators[d][0];
	    if (dvv[d] == 1) {
	      //t1_accum_0[d] = 0.;
	      t2_accum_0[0] = prod2 * t2_coeff_p[0];
	      //for (v=1; v<num_v; ++v)
	      //  t2_accum_0[v] = 0.;
	    }
	    else {
	      t1_accum_0[d]   = prod1;
	      t2_accum_0[0]   = prod5 * t2_coeff_p[0];
	      for (v=1; v<num_v; ++v)
		t2_accum_0[v] = prod2 * t2_coeff_p[v];
	    }
	  }
	}
	if (key_p0 == max_0)
	  data_rep->accumulate_horners_gradient(t1_accumulator, t2_accumulators,
						lev_index, key_p, dvv, x);
      }
      Real *t1_accum = t1_accumulator[num_v-1], *t2_accum;
      for (d=0; d<num_deriv_vars; ++d) {
	tpVarianceGrad[d] = t1_accum[d];
	t2_accum = t2_accumulators[d][num_v-1];
	for (v=0; v<num_v; ++v)
	  tpVarianceGrad[d] += t2_accum[v];
      }

      /*
      Real t1_coeff_p_mm;
      for (p=0; p<num_colloc_pts; ++p) {
	const UShortArray& key_p = key[p];
	c_index_p = (colloc_index.empty()) ? p : colloc_index[p];
	t1_coeff_p_mm = expansionType1Coeffs[c_index_p] - mean_1;
	// no DVV: need all grad components for R to expand type2 for (R-mu)^2
	const Real* t2_coeff_p = expansionType2Coeffs[c_index_p];

	for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	  deriv_index = dvv[d] - 1; // OK since we are in an "All" view
	  // ---------------------------------------------------------------
	  // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	  // ---------------------------------------------------------------
	  Real& grad_d = tpVarianceGrad[d];
	  grad_d += t1_coeff_p_mm * t1_coeff_p_mm
	    * type1_interpolant_gradient(x, deriv_index, key_p,
					 lev_index, data_rep->nonRandomIndices)
	    * type1_weight(key_p, lev_index, data_rep->randomIndices);
	  for (v=0; v<num_v; ++v)
	    grad_d += 2. * t1_coeff_p_mm * t2_coeff_p[v]
	      * type2_interpolant_gradient(x, deriv_index, v, key_p, lev_index,
	                                   data_rep->nonRandomIndices)
	      * type2_weight(v, key_p, lev_index, data_rep->randomIndices);
	}
      }
      */
    }
    else {
      //PCout << "tensor_product_variance_gradient(): Horner's no derivs."
      //      << std::endl;
      // Horner's rule approach:
      RealMatrix accumulator(num_deriv_vars, num_v); // init to 0.
      Real *accum_0 = accumulator[0], t1_coeff_p_mm, t1_wt_00, t1_val,
	prod1, prod2, prod3, x0 = x[0];
      for (p=0; p<num_colloc_pts; ++p) {
	const UShortArray& key_p = key[p]; key_p0 = key_p[0];
	c_index_p = (colloc_index.empty()) ? p : colloc_index[p];
	t1_coeff_p_mm = expansionType1Coeffs[c_index_p] - mean_1;
	const Real* t1_coeff_grad_p = (insert) ?
	  expansionType1CoeffGrads[c_index_p] : NULL;
	if (rand_0) { // integration: deriv of expect = expect of deriv
	  if (li_0) { // t1 weight is 1
	    t1_wt_00 = t1_wts_0[key_p0];
	    if (insert)  prod1 = 2. * t1_coeff_p_mm * t1_wt_00;
	    if (augment) prod2 = t1_coeff_p_mm * t1_coeff_p_mm * t1_wt_00;
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {       // case 1:2
	      deriv_index = dvv[d] - 1;
	      accum_0[d] += (data_rep->randomVarsKey[deriv_index]) ? prod1 *
		(t1_coeff_grad_p[insert_cntr++] - mean_grad[d]) : prod2;
	    }
	  }
	  else {
	    if (insert)  prod1 = 2. * t1_coeff_p_mm;
	    if (augment) prod2 = t1_coeff_p_mm * t1_coeff_p_mm;
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	      deriv_index = dvv[d] - 1;
	      accum_0[d]  = (data_rep->randomVarsKey[deriv_index]) ? prod1 *
		(t1_coeff_grad_p[insert_cntr++] - mean_grad[d]) : prod2;
	    }
	  }
	}
	else { // interpolation: deriv of interpolation of mean
	  if (li_0) { // t1 grad is 0., t1 value is 1
	    t1_val = poly_0.type1_value(x0, key_p0);
	    if (insert) prod1  = 2. * t1_coeff_p_mm * t1_val;
	    if (augment)
	      { prod2 = t1_coeff_p_mm * t1_coeff_p_mm; prod3 = prod2 * t1_val; }
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	      deriv_index = dvv[d] - 1;
	      if (data_rep->randomVarsKey[deriv_index])               // case 3
		accum_0[d] += prod1 *
		  (t1_coeff_grad_p[insert_cntr++] - mean_grad[d]);
	      else if (deriv_index == 0)                              // case 4
		accum_0[d] += prod2 * poly_0.type1_gradient(x0, key_p0);
	      else                                                    // case 4
		accum_0[d] += prod3;
	    }
	  }
	  else {
	    if (insert)  prod1 = 2. * t1_coeff_p_mm;
	    if (augment) prod2 = t1_coeff_p_mm * t1_coeff_p_mm;
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	      deriv_index = dvv[d] - 1;
	      if (data_rep->randomVarsKey[deriv_index])               // case 3
		accum_0[d] = prod1 *
		  (t1_coeff_grad_p[insert_cntr++] - mean_grad[d]);
	      else if (deriv_index)                                   // case 4
		accum_0[d] = prod2;
	    }
	  }
	}
	if (key_p0 == max_0)
	  data_rep->accumulate_horners_gradient(accumulator, lev_index, key_p,
						dvv, x);
      }
      copy_data(accumulator[num_v-1], (int)num_deriv_vars, tpVarianceGrad);

      /*
      Real t1_coeff_p_mm;
      size_t p, v, num_colloc_pts = key.size();
      for (p=0; p<num_colloc_pts; ++p) {
	const UShortArray& key_p = key[p];
	c_index_p = (colloc_index.empty()) ? p : colloc_index[p];
	t1_coeff_p_mm = expansionType1Coeffs[c_index_p] - mean_1;
	const Real* t1_coeff_grad_p = (insert) ?
	  expansionType1CoeffGrads[c_index_p] : NULL;
	for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	  deriv_index = dvv[d] - 1; // OK since we are in an "All" view
	  tpVarianceGrad[d] += (data_rep->randomVarsKey[deriv_index]) ?
	    // ---------------------------------------------------------------
	    // deriv of All var expansion w.r.t. random var (design insertion)
	    // ---------------------------------------------------------------
	    // d/dx[(R-mu)^2] = 2(R-mu)(dR/dx - dmu/dx)
	    2. * t1_coeff_p_mm * (t1_coeff_grad_p[insert_cntr++] - mean_grad[d])
	      * type1_interpolant_value(x, key_p, lev_index,
					data_rep->nonRandomIndices)
	      * type1_weight(key_p, lev_index, data_rep->randomIndices) :
	    // ---------------------------------------------------------------
	    // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	    // ---------------------------------------------------------------
	    t1_coeff_p_mm * t1_coeff_p_mm
	      * type1_interpolant_gradient(x, deriv_index, key_p, lev_index,
	                                   data_rep->nonRandomIndices)
	      * type1_weight(key_p, lev_index, data_rep->randomIndices);
	}
      }
      */
    }
    break;
  }
  case REINTERPOLATION_OF_PRODUCTS: {
    const UShortArray&   reinterp_lev_index
      = driver_rep->reinterpolated_level_index();
    const RealMatrix&    reinterp_var_sets
      = driver_rep->reinterpolated_variable_sets();
    const UShort2DArray& reinterp_colloc_key
      = driver_rep->reinterpolated_collocation_key();

    unsigned short key_p0, rli_0 = reinterp_lev_index[0];
    BasisPolynomial&      poly_0 = data_rep->polynomialBasis[rli_0][0];
    const Real3DArray& t1_wts_1d = driver_rep->type1_collocation_weights_1d();
    const RealArray&   t1_wts_0  = t1_wts_1d[rli_0][0];
    size_t ei_0 = poly_0.exact_index(), p, v, d,
      reinterp_colloc_pts = reinterp_colloc_key.size();
    bool rand_0 = data_rep->randomVarsKey[0];
    unsigned short max_0 = poly_0.interpolation_size() - 1;

    if (data_rep->barycentricFlag) {
      // For barycentric interpolation: track x != newPoint within 1D basis
      short order = (augment) ? 3 : 1;
      data_rep->set_new_point(x, reinterp_lev_index,
			      data_rep->nonRandomIndices, order);

      const RealVector& bc_vf_0 = poly_0.barycentric_value_factors();
      const RealVector& bc_gf_0 = poly_0.barycentric_gradient_factors();
      RealMatrix accumulator(num_deriv_vars, num_v); // init to 0.
      Real *accum_0 = accumulator[0], t1_coeff_p_mm, t1_wt_00, bc_vf_00,
	prod1, prod2, prod3;
      RealVector empty_rv;
      for (p=0; p<reinterp_colloc_pts; ++p) {
	const UShortArray& key_p = reinterp_colloc_key[p]; key_p0 = key_p[0];
	RealVector c_vars(Teuchos::View,
	  const_cast<Real*>(reinterp_var_sets[p]), num_v);
	t1_coeff_p_mm = value(c_vars) - mean_1;
	const RealVector& t1_coeff_grad_p = (insert) ?
	  gradient_nonbasis_variables(c_vars) : empty_rv;
	if (rand_0) { // integration: deriv of expect = expect of deriv
	  if (rli_0) {
	    t1_wt_00 = t1_wts_0[key_p0];
	    if (insert)  prod1 = 2. * t1_coeff_p_mm * t1_wt_00;
	    if (augment) prod2 = t1_coeff_p_mm * t1_coeff_p_mm * t1_wt_00;
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {       // case 1:2
	      deriv_index = dvv[d] - 1;
	      accum_0[d] += (data_rep->randomVarsKey[deriv_index]) ? prod1 *
		(t1_coeff_grad_p[insert_cntr++] - mean_grad[d]) : prod2;
	    }
	  }
	  else {
	    if (insert)  prod1 = 2. * t1_coeff_p_mm;
	    if (augment) prod2 = t1_coeff_p_mm * t1_coeff_p_mm;
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {       // case 1:2
	      deriv_index = dvv[d] - 1;
	      accum_0[d]  = (data_rep->randomVarsKey[deriv_index]) ? prod1 *
		(t1_coeff_grad_p[insert_cntr++] - mean_grad[d]) : prod2;
	    }
	  }
	}
	else { // interpolation: deriv of interpolation of mean
	  if (rli_0)  { // grad factor is 0., value factor is omitted
	    bc_vf_00 = bc_vf_0[key_p0];
	    if (insert) prod1  = 2. * t1_coeff_p_mm * bc_vf_00;
	    if (augment) {
	      prod2 = t1_coeff_p_mm * t1_coeff_p_mm;
	      if (ei_0 == _NPOS) prod3 = prod2 * bc_vf_00;
	    }
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	      deriv_index = dvv[d] - 1;
	      if (data_rep->randomVarsKey[deriv_index])               // case 3
		accum_0[d] += prod1 *
		  (t1_coeff_grad_p[insert_cntr++] - mean_grad[d]);
	      else if (deriv_index == 0)                              // case 4
		accum_0[d] += prod2 * bc_gf_0[key_p0];
	      else if (ei_0 == _NPOS)                                 // case 4
		accum_0[d] += prod3;
	      else if (ei_0 == key_p0)                                // case 4
		accum_0[d] += prod2;
	    }
	  }
	  else {
	    if (insert)  prod1 = 2. * t1_coeff_p_mm;
	    if (augment) prod2 = t1_coeff_p_mm * t1_coeff_p_mm;
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	      deriv_index = dvv[d] - 1;
	      if (data_rep->randomVarsKey[deriv_index])               // case 3
		accum_0[d] = prod1 *
		  (t1_coeff_grad_p[insert_cntr++] - mean_grad[d]);
	      else if (deriv_index)                                   // case 4
		accum_0[d] = prod2;
	    }
	  }
	}
	if (key_p0 == max_0) // accumulate sums over vars with max key value
	  data_rep->accumulate_barycentric_gradient(accumulator,
						    reinterp_lev_index,
						    key_p, dvv);
      }
      Real scale = data_rep->
	barycentric_gradient_scaling(reinterp_lev_index,
				     data_rep->nonRandomIndices);
      Real* accum = accumulator[num_v-1];
      for (d=0; d<num_deriv_vars; ++d)
	tpVarianceGrad[d] = accum[d] * scale;
    }
    else if (data_rep->basisConfigOptions.useDerivs) {
      if (insert) {
	PCerr << "Error: combination of coefficient gradients and "
	      << "use_derivatives in NodalInterpPolyApproximation::"
	      << "tensor_product_mean_gradient()" << std::endl;
	abort_handler(-1);
      }

      //PCout << "tensor_product_variance_gradient(): Horner's useDerivs."
      //      << std::endl;
      const Real3DArray& t2_wts_1d = driver_rep->type2_collocation_weights_1d();
      const RealArray&   t2_wts_0  = t2_wts_1d[rli_0][0];
      RealMatrix t1_accumulator(num_deriv_vars, num_v); // init to 0.
      RealMatrixArray t2_accumulators(num_deriv_vars);
      for (v=0; v<num_deriv_vars; ++v)
	t2_accumulators[v].shape(num_v, num_v); // init to 0.
      Real *t1_accum_0 = t1_accumulator[0], *t2_accum_0, t1_coeff_p_mm,
	t1_wt_00, t1_val, t2_val, t1_grad, prod1, prod2, prod3, prod4,
	prod5,prod6, x0 = x[0];
      for (p=0; p<reinterp_colloc_pts; ++p) {
	const UShortArray& key_p = reinterp_colloc_key[p]; key_p0 = key_p[0];
	RealVector c_vars(Teuchos::View,
	  const_cast<Real*>(reinterp_var_sets[p]), num_v);
	t1_coeff_p_mm = value(c_vars) - mean_1;
	// no DVV: need all grad components for R to expand type2 for (R-mu)^2
	const RealVector& t2_coeff_p = gradient_basis_variables(c_vars);
	prod1 = t1_coeff_p_mm * t1_coeff_p_mm;
	prod2 = 2. * t1_coeff_p_mm;
	if (rand_0) { // integration: deriv of expect = expect of deriv
	  if (rli_0) {                                                 // case 2
	    t1_wt_00 = t1_wts_0[key_p0];
	    prod3 = prod1 * t1_wt_00;
	    prod4 = prod2 * t2_coeff_p[0] * t2_wts_0[key_p0];
	    for (d=0; d<num_deriv_vars; ++d) {
	      t1_accum_0[d]            += prod3;
	      t2_accumulators[d][0][0] += prod4;
	    }
	    prod3 = prod2 * t1_wt_00;
	    for (v=1; v<num_v; ++v) {
	      prod4 = prod3 * t2_coeff_p[v];
	      for (d=0; d<num_deriv_vars; ++d)
		t2_accumulators[d][0][v] += prod4;
	    }
	  }
	  else { // t1 weight is 1, t2 weight is 0
	    for (d=0; d<num_deriv_vars; ++d)
	      t1_accum_0[d] = prod1;
	    for (v=1; v<num_v; ++v) {
	      prod4 = prod2 * t2_coeff_p[v];
	      for (d=0; d<num_deriv_vars; ++d)
		t2_accumulators[d][0][v] += prod4;
	    }
	  }
	}
	else if (rli_0) { // interpolation: deriv of interpolation of mean
	  t1_val = poly_0.type1_value(x0, key_p0);
	  t2_val = poly_0.type2_value(x0, key_p0);
	  prod3 = prod1 * t1_val; prod4  = prod2 * t1_val;
	  prod5 = prod2 * t2_val;
	  for (d=0; d<num_deriv_vars; ++d) {                         // case 4
	    t2_accum_0 = t2_accumulators[d][0];
	    if (dvv[d] == 1) {
	      t1_grad = poly_0.type1_gradient(x0, key_p0);
	      t1_accum_0[d]   += prod1 * t1_grad;
	      t2_accum_0[0]   += prod2 * t2_coeff_p[0] *
		poly_0.type2_gradient(x0, key_p0);
	      prod6 = prod2 * t1_grad;
	      for (v=1; v<num_v; ++v)
		t2_accum_0[v] += prod6 * t2_coeff_p[v];
	    }
	    else {
	      t1_accum_0[d]   += prod3;
	      t2_accum_0[0]   += prod5 * t2_coeff_p[0];
	      for (v=1; v<num_v; ++v)
		t2_accum_0[v] += prod4 * t2_coeff_p[v];
	    }
	  }
	}
	else { // t1 value is 1, t1 grad is 0., t2 grad is 1
	  prod5 = prod2 * poly_0.type2_value(x0, key_p0);
	  for (d=0; d<num_deriv_vars; ++d) {                         // case 4
	    t2_accum_0 = t2_accumulators[d][0];
	    if (dvv[d] == 1) {
	      //t1_accum_0[d] = 0.;
	      t2_accum_0[0] = prod2 * t2_coeff_p[0];
	      //for (v=1; v<num_v; ++v)
	      //  t2_accum_0[v] = 0.;
	    }
	    else {
	      t1_accum_0[d]   = prod1;
	      t2_accum_0[0]   = prod5 * t2_coeff_p[0];
	      for (v=1; v<num_v; ++v)
		t2_accum_0[v] = prod2 * t2_coeff_p[v];
	    }
	  }
	}
	if (key_p0 == max_0)
	  data_rep->accumulate_horners_gradient(t1_accumulator, t2_accumulators,
						reinterp_lev_index, key_p, dvv,
						x);
      }
      Real *t1_accum = t1_accumulator[num_v-1], *t2_accum;
      for (d=0; d<num_deriv_vars; ++d) {
	tpVarianceGrad[d] = t1_accum[d];
	t2_accum = t2_accumulators[d][num_v-1];
	for (v=0; v<num_v; ++v)
	  tpVarianceGrad[d] += t2_accum[v];
      }

      /*
      Real t1_coeff_p_mm;
      for (p=0; p<reinterp_colloc_pts; ++p) {
	const UShortArray& key_p = reinterp_colloc_key[p];
	RealVector c_vars(Teuchos::View,
	  const_cast<Real*>(reinterp_var_sets[p]), num_v);
	t1_coeff_p_mm = value(c_vars) - mean_1;
	// no DVV: need all grad components for R to expand type2 for (R-mu)^2
	const RealVector& t2_coeff_p = gradient_basis_variables(c_vars);

	for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	  deriv_index = dvv[d] - 1; // OK since we are in an "All" view
	  // ---------------------------------------------------------------
	  // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	  // ---------------------------------------------------------------
	  Real& grad_d = tpVarianceGrad[d];
	  grad_d += t1_coeff_p_mm * t1_coeff_p_mm
	    * data_rep->type1_interpolant_gradient(x, deriv_index, key_p,
	                                           reinterp_lev_index,
					           data_rep->nonRandomIndices)
	    * data_rep->type1_weight(key_p, reinterp_lev_index,
	                             data_rep->randomIndices);
	  for (v=0; v<num_v; ++v)
	    grad_d += 2. * t1_coeff_p_mm * t2_coeff_p[v]
	      * data_rep->type2_interpolant_gradient(x, deriv_index, v, key_p,
					             reinterp_lev_index,
					             data_rep->nonRandomIndices)
	      * data_rep->type2_weight(v, key_p, reinterp_lev_index,
	                               data_rep->randomIndices);
	}
      }
      */
    }
    else {
      //PCout << "tensor_product_variance_gradient(): Horner's no derivs."
      //      << std::endl;
      // Horner's rule approach:
      RealMatrix accumulator(num_deriv_vars, num_v); // init to 0.
      Real *accum_0 = accumulator[0], t1_coeff_p_mm, t1_wt_00, t1_val,
	prod1, prod2, prod3, x0 = x[0];
      RealVector empty_rv;
      for (p=0; p<reinterp_colloc_pts; ++p) {
	const UShortArray& key_p = reinterp_colloc_key[p]; key_p0 = key_p[0];
	RealVector c_vars(Teuchos::View,
	  const_cast<Real*>(reinterp_var_sets[p]), num_v);
	t1_coeff_p_mm = value(c_vars) - mean_1;
	const RealVector& t1_coeff_grad_p = (insert) ?
	  gradient_nonbasis_variables(c_vars) : empty_rv;
	if (rand_0) { // integration: deriv of expect = expect of deriv
	  if (rli_0) { // t1 weight is 1
	    t1_wt_00 = t1_wts_0[key_p0];
	    if (insert)  prod1 = 2. * t1_coeff_p_mm * t1_wt_00;
	    if (augment) prod2 = t1_coeff_p_mm * t1_coeff_p_mm * t1_wt_00;
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {       // case 1:2
	      deriv_index = dvv[d] - 1;
	      accum_0[d] += (data_rep->randomVarsKey[deriv_index]) ? prod1 *
		(t1_coeff_grad_p[insert_cntr++] - mean_grad[d]) : prod2;
	    }
	  }
	  else {
	    if (insert)  prod1 = 2. * t1_coeff_p_mm;
	    if (augment) prod2 = t1_coeff_p_mm * t1_coeff_p_mm;
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	      deriv_index = dvv[d] - 1;
	      accum_0[d]  = (data_rep->randomVarsKey[deriv_index]) ? prod1 *
		(t1_coeff_grad_p[insert_cntr++] - mean_grad[d]) : prod2;
	    }
	  }
	}
	else { // interpolation: deriv of interpolation of mean
	  if (rli_0) { // t1 grad is 0., t1 value is 1
	    t1_val = poly_0.type1_value(x0, key_p0);
	    if (insert) prod1  = 2. * t1_coeff_p_mm * t1_val;
	    if (augment)
	      { prod2 = t1_coeff_p_mm * t1_coeff_p_mm; prod3 = prod2 * t1_val; }
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	      deriv_index = dvv[d] - 1;
	      if (data_rep->randomVarsKey[deriv_index])               // case 3
		accum_0[d] += prod1 *
		  (t1_coeff_grad_p[insert_cntr++] - mean_grad[d]);
	      else if (deriv_index == 0)                              // case 4
		accum_0[d] += prod2 * poly_0.type1_gradient(x0, key_p0);
	      else                                                    // case 4
		accum_0[d] += prod3;
	    }
	  }
	  else {
	    if (insert)  prod1 = 2. * t1_coeff_p_mm;
	    if (augment) prod2 = t1_coeff_p_mm * t1_coeff_p_mm;
	    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	      deriv_index = dvv[d] - 1;
	      if (data_rep->randomVarsKey[deriv_index])               // case 3
		accum_0[d] = prod1 *
		  (t1_coeff_grad_p[insert_cntr++] - mean_grad[d]);
	      else if (deriv_index)                                   // case 4
		accum_0[d] = prod2;
	    }
	  }
	}
	if (key_p0 == max_0)
	  data_rep->accumulate_horners_gradient(accumulator, reinterp_lev_index,
						key_p, dvv, x);
      }
      copy_data(accumulator[num_v-1], (int)num_deriv_vars, tpVarianceGrad);

      /*
      RealVector empty_rv; Real t1_coeff_p_mm;
      for (p=0; p<reinterp_colloc_pts; ++p) {
	const UShortArray& key_p = reinterp_colloc_key[p];
	RealVector c_vars(Teuchos::View,
	  const_cast<Real*>(reinterp_var_sets[p]), num_v);
	t1_coeff_p_mm = value(c_vars) - mean_1;
	const RealVector& t1_coeff_grad_p = (insert) ?
	  gradient_nonbasis_variables(c_vars) : empty_rv;
	for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	  deriv_index = dvv[d] - 1; // OK since we are in an "All" view
	  tpVarianceGrad[d] += (data_rep->randomVarsKey[deriv_index]) ?
	    // ---------------------------------------------------------------
	    // deriv of All var expansion w.r.t. random var (design insertion)
	    // ---------------------------------------------------------------
	    // d/dx[(R-mu)^2] = 2(R-mu)(dR/dx - dmu/dx)
	    2. * t1_coeff_p_mm * (t1_coeff_grad_p[insert_cntr++] - mean_grad[d])
	      * data_rep->type1_interpolant_value(x, key_p, reinterp_lev_index,
					          data_rep->nonRandomIndices)
	      * data_rep->type1_weight(key_p, reinterp_lev_index,
	                               data_rep->randomIndices) :
	    // ---------------------------------------------------------------
	    // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	    // ---------------------------------------------------------------
	    t1_coeff_p_mm * t1_coeff_p_mm
	      * data_rep->type1_interpolant_gradient(x, deriv_index, key_p,
					             reinterp_lev_index,
						     data_rep->nonRandomIndices)
	      * data_rep->type1_weight(key_p, reinterp_lev_index,
	                               data_rep->randomIndices);
	}
      }
      */
    }
    break;
  }
  case PRODUCT_OF_INTERPOLANTS_FAST: {
    // -------------------------------------------------------------------
    // Mixed variable key:
    //   xi = ran vars, sa = augmented des vars, si = inserted design vars
    //   This assumes that the basis is fixed (ignores parameterized bases).
    // Active variable expansion (product of interpolants):
    //   R(xi, sa, si) = Sum_i r_i(sa, si) L_i(xi)
    //   covar(sa, si) = Sum_i Sum_j r1_i(sa, si) r2_j(sa, si) wt_prod_xi_ij
    //                 - mu^2
    //   dcovar/ds     = Sum_i Sum_j (r1_i dr2_j/ds+dr1_i/ds r2_j) wt_prod_xi_ij
    //                 - 2 mu dmu/ds
    // All variable expansion (product of interpolants):
    //   R(xi, sa, si) = Sum_i r_i(si) L_i(xi, sa)
    //   covar(sa, si) = Sum_i Sum_j r1_i(si) r2_j(si) Lsa_i Lsa_j wt_prod_xi_ij
    //                 - mu^2
    //   dcovar/dsa    = Sum_i Sum_j r1_i(si) r2_j(si) (Lsa_i dLsa_j/dsa
    //                 + dLsa_i/dsa Lsa_j) wt_prod_xi_ij - 2 mu dmu/dsa
    //   dcovar/dsi    = Sum_i Sum_j (r1_i dr2_j/dsi + dr1_i/dsi r2_j)
    //                   Lsa_i Lsa_j wt_prod_xi_ij - 2 mu dmu/dsi
    // -------------------------------------------------------------------
    size_t j, k, c_index_j, c_index_k, num_colloc_pts = key.size();
    Real wt_prod_j, Lsa_j, Lsa_k, dLsa_j_dsa_d, dLsa_k_dsa_d,
      t1_coeff_j_mm, t1_coeff_k_mm;
    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
      deriv_index = dvv[d] - 1; // OK since we are in an "All" view
      Real& grad_d = tpVarianceGrad[d];
      // first loop of double sum
      for (j=0; j<num_colloc_pts; ++j) {
	const UShortArray& key_j = key[j];
	c_index_j = (colloc_index.empty()) ? j : colloc_index[j];
	t1_coeff_j_mm = expansionType1Coeffs[c_index_j] - mean_1;
	// compute wt_prod_j and Lsa_j
	wt_prod_j
	  = data_rep->type1_weight(key_j, lev_index, data_rep->randomIndices);
	Lsa_j = data_rep->type1_interpolant_value(x, key_j, lev_index,
						  data_rep->nonRandomIndices);
	dLsa_j_dsa_d
	  = data_rep->type1_interpolant_gradient(x, deriv_index, key_j,
						 lev_index,
						 data_rep->nonRandomIndices);
	// second loop of double sum
	for (k=0; k<num_colloc_pts; ++k) {
	  const UShortArray& key_k = key[k];
	  c_index_k = (colloc_index.empty()) ? k : colloc_index[k];
	  // to include jk-th term, colloc pts xi_j must be the same as xi_k
	  // for random var subset.  In this case, wt_prod_j may be reused.
	  if (data_rep->match_random_key(key_j, key_k)) {
	    t1_coeff_k_mm = expansionType1Coeffs[c_index_k] - mean_1;
	    Lsa_k =
	      data_rep->type1_interpolant_value(x, key_k, lev_index,
						data_rep->nonRandomIndices);
	    if (data_rep->randomVarsKey[deriv_index])
	      // ---------------------------------------------------------
	      // deriv of All var exp w.r.t. random var (design insertion)
	      // ---------------------------------------------------------
	      grad_d += wt_prod_j * Lsa_j * Lsa_k *
		( t1_coeff_j_mm * ( expansionType1CoeffGrads(insert_cntr,
				    c_index_k) - mean_grad[d] )
		+ t1_coeff_k_mm * ( expansionType1CoeffGrads(insert_cntr,
				    c_index_j) - mean_grad[d] ) );
	    else {
	      // ---------------------------------------------------------------
	      // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	      // ---------------------------------------------------------------
	      dLsa_k_dsa_d = data_rep->
		type1_interpolant_gradient(x, deriv_index, key_k, lev_index,
					   data_rep->nonRandomIndices);
	      grad_d += wt_prod_j * t1_coeff_j_mm * t1_coeff_k_mm *
		(Lsa_j * dLsa_k_dsa_d + dLsa_j_dsa_d * Lsa_k);
	    }
	  }
	}
      }
      if (data_rep->randomVarsKey[deriv_index])// deriv w.r.t. des var insertion
	++insert_cntr;
    }
    break;
  }
  case PRODUCT_OF_INTERPOLANTS_FULL:
    // TO DO
    break;
  }

  return tpVarianceGrad;
}


Real NodalInterpPolyApproximation::value(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response prediction
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    SizetArray colloc_index; // empty -> default indexing
    return data_rep->
      tensor_product_value(x, expansionType1Coeffs, expansionType2Coeffs,
			   tpq_driver->level_index(),
			   tpq_driver->collocation_key(), colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    Real approx_val = 0.;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	approx_val += sm_coeffs[i] * data_rep->
	  tensor_product_value(x, expansionType1Coeffs, expansionType2Coeffs,
			       sm_mi[i], colloc_key[i], colloc_index[i]);
    return approx_val;
    break;
  }
  }
}


/** Special case used for sparse grid interpolation on variable sub-sets
    defined from partial integration. */
Real NodalInterpPolyApproximation::
value(const RealVector& x, const RealVectorArray& t1_coeffs,
      const RealMatrixArray& t2_coeffs, const UShort3DArray& colloc_key,
      const SizetList& subset_indices)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  // Smolyak recursion of anisotropic tensor products
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
  const UShort2DArray&      sm_mi      = csg_driver->smolyak_multi_index();
  const IntArray&           sm_coeffs  = csg_driver->smolyak_coefficients();
  size_t i, num_smolyak_indices = sm_coeffs.size(); SizetArray colloc_index;
  Real approx_val = 0.;
  for (i=0; i<num_smolyak_indices; ++i)
    if (sm_coeffs[i])
      approx_val += sm_coeffs[i] * data_rep->
	tensor_product_value(x, t1_coeffs[i], t2_coeffs[i], sm_mi[i],
			     colloc_key[i], colloc_index, subset_indices);
  return approx_val;
}


const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    SizetArray colloc_index; // empty -> default indexing
    return data_rep->
      tensor_product_gradient_basis_variables(x, expansionType1Coeffs,
					      expansionType2Coeffs,
					      tpq_driver->level_index(),
					      tpq_driver->collocation_key(),
					      colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    size_t num_v = sharedDataRep->numVars;
    if (approxGradient.length() != num_v)
      approxGradient.sizeUninitialized(num_v);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = data_rep->
	  tensor_product_gradient_basis_variables(x, expansionType1Coeffs,
						  expansionType2Coeffs,
						  sm_mi[i], colloc_key[i],
						  colloc_index[i]);
	for (j=0; j<num_v; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const RealVectorArray& t1_coeffs,
			 const RealMatrixArray& t2_coeffs,
			 const UShort3DArray& colloc_key,
			 const SizetList& subset_indices)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  size_t num_v = sharedDataRep->numVars;
  if (approxGradient.length() != num_v)
    approxGradient.sizeUninitialized(num_v);
  approxGradient = 0.;
  // Smolyak recursion of anisotropic tensor products
  CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
  const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
  const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
  size_t i, j, num_smolyak_indices = sm_coeffs.size(); SizetArray colloc_index;
  for (i=0; i<num_smolyak_indices; ++i) {
    int coeff_i = sm_coeffs[i];
    if (coeff_i) {
      const RealVector& tp_grad = data_rep->
	tensor_product_gradient_basis_variables(x, t1_coeffs[i], t2_coeffs[i],
						sm_mi[i], colloc_key[i],
						colloc_index, subset_indices);
      for (j=0; j<num_v; ++j)
	approxGradient[j] += coeff_i * tp_grad[j];
    }
  }
  return approxGradient;
}


/** Special case used for sparse grid interpolation on variable sub-sets
    defined from partial integration. */
const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    SizetArray colloc_index; // empty -> default indexing
    return data_rep->
      tensor_product_gradient_basis_variables(x, expansionType1Coeffs,
					      expansionType2Coeffs,
					      tpq_driver->level_index(),
					      tpq_driver->collocation_key(),
					      colloc_index, dvv);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    size_t num_deriv_vars = dvv.size();
    if (approxGradient.length() != num_deriv_vars)
      approxGradient.sizeUninitialized(num_deriv_vars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = data_rep->
	  tensor_product_gradient_basis_variables(x, expansionType1Coeffs,
						  expansionType2Coeffs,
						  sm_mi[i], colloc_key[i],
						  colloc_index[i], dvv);
	for (j=0; j<num_deriv_vars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x)
{
  // Error check for required data
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in NodalInterp"
	  << "PolyApproximation::gradient_nonbasis_variables()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    SizetArray colloc_index; // empty -> default indexing
    return data_rep->
      tensor_product_gradient_nonbasis_variables(x, expansionType1CoeffGrads,
						 tpq_driver->level_index(),
						 tpq_driver->collocation_key(),
						 colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    size_t num_deriv_vars = expansionType1CoeffGrads.numRows();
    if (approxGradient.length() != num_deriv_vars)
      approxGradient.sizeUninitialized(num_deriv_vars);
    approxGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver   = data_rep->csg_driver();
    const UShort2DArray&      sm_mi        = csg_driver->smolyak_multi_index();
    const IntArray&           sm_coeffs    = csg_driver->smolyak_coefficients();
    const UShort3DArray&      colloc_key   = csg_driver->collocation_key();
    const Sizet2DArray&       colloc_index = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = data_rep->
	  tensor_product_gradient_nonbasis_variables(x,
						     expansionType1CoeffGrads,
						     sm_mi[i], colloc_key[i],
						     colloc_index[i]);
	for (j=0; j<num_deriv_vars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


const RealSymMatrix& NodalInterpPolyApproximation::
hessian_basis_variables(const RealVector& x)
{
  PCerr << "Error: NodalInterpPolyApproximation::hessian_basis_variables() "
	<< "not yet implemented." << std::endl;
  abort_handler(-1);

  return approxHessian;
}


Real NodalInterpPolyApproximation::
stored_value(const RealVector& x, size_t index)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not available in "
	  << "NodalInterpPolyApproximation::stored_value()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  // sum expansion to get response prediction
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    SizetArray colloc_index; // empty -> default indexing
    return data_rep->
      tensor_product_value(x, storedExpType1Coeffs[index],
			   storedExpType2Coeffs[index],
			   tpq_driver->stored_level_index(index),
			   tpq_driver->stored_collocation_key(index),
			   colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
    const IntArray&  sm_coeffs = csg_driver->stored_smolyak_coefficients(index);
    const UShort2DArray& sm_mi = csg_driver->stored_smolyak_multi_index(index);
    const UShort3DArray& colloc_key = csg_driver->stored_collocation_key(index);
    const Sizet2DArray& colloc_index
      = csg_driver->stored_collocation_indices(index);
    // Smolyak recursion of anisotropic tensor products
    size_t i, num_smolyak_indices = sm_coeffs.size();
    Real approx_val = 0.;
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i)
	approx_val += coeff_i * data_rep->
	  tensor_product_value(x, storedExpType1Coeffs[index],
			       storedExpType2Coeffs[index], sm_mi[i],
			       colloc_key[i], colloc_index[i]);
    }
    return approx_val;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
stored_gradient_basis_variables(const RealVector& x, size_t index)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not available in NodalInterpPoly"
	  << "Approximation::stored_gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    SizetArray colloc_index; // empty -> default indexing
    return data_rep->
      tensor_product_gradient_basis_variables(x, storedExpType1Coeffs[index],
	storedExpType2Coeffs[index], tpq_driver->stored_level_index(index),
	tpq_driver->stored_collocation_key(index), colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
    const IntArray&  sm_coeffs = csg_driver->stored_smolyak_coefficients(index);
    const UShort2DArray& sm_mi = csg_driver->stored_smolyak_multi_index(index);
    const UShort3DArray& colloc_key = csg_driver->stored_collocation_key(index);
    const Sizet2DArray& colloc_index
      = csg_driver->stored_collocation_indices(index);
    // Smolyak recursion of anisotropic tensor products
    size_t i, j, num_smolyak_indices = sm_coeffs.size(),
      num_v = sharedDataRep->numVars;
    if (approxGradient.length() != num_v)
      approxGradient.sizeUninitialized(num_v);
    approxGradient = 0.;
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = data_rep->
	  tensor_product_gradient_basis_variables(x,
	    storedExpType1Coeffs[index], storedExpType2Coeffs[index],
	    sm_mi[i], colloc_key[i], colloc_index[i]);
	for (j=0; j<num_v; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
stored_gradient_nonbasis_variables(const RealVector& x, size_t index)
{
  // Error check for required data
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not available in Nodal"
	  << "InterpPolyApproximation::stored_gradient_nonbasis_variables()"
	  << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    SizetArray colloc_index; // empty -> default indexing
    return data_rep->tensor_product_gradient_nonbasis_variables(x,
      storedExpType1CoeffGrads[index], tpq_driver->stored_level_index(index),
      tpq_driver->stored_collocation_key(index), colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
    const IntArray&  sm_coeffs = csg_driver->stored_smolyak_coefficients(index);
    const UShort2DArray& sm_mi = csg_driver->stored_smolyak_multi_index(index);
    const UShort3DArray& colloc_key = csg_driver->stored_collocation_key(index);
    const Sizet2DArray& colloc_index
      = csg_driver->stored_collocation_indices(index);
    // Smolyak recursion of anisotropic tensor products
    size_t i, j, num_smolyak_indices = sm_coeffs.size(),
      num_deriv_vars = storedExpType1CoeffGrads[index].numRows();
    if (approxGradient.length() != num_deriv_vars)
      approxGradient.sizeUninitialized(num_deriv_vars);
    approxGradient = 0.;
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff_i = sm_coeffs[i];
      if (coeff_i) {
	const RealVector& tp_grad = data_rep->
	  tensor_product_gradient_nonbasis_variables(x,
	    storedExpType1CoeffGrads[index], sm_mi[i], colloc_key[i],
	    colloc_index[i]);
	for (j=0; j<num_deriv_vars; ++j)
	  approxGradient[j] += coeff_i * tp_grad[j];
      }
    }
    return approxGradient;
    break;
  }
  }
}


/** In this case, all expansion variables are random variables and the
    mean of the expansion is simply the sum over i of r_i w_i. */
Real NodalInterpPolyApproximation::mean()
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  //IntegrationDriver* driver_rep = data_rep->driverRep;

  // TO DO:
  //if (!driver_rep->track_unique_product_weights()) {
  //  PCerr << "Error: unique product weights required in "
  //	  << "NodalInterpPolyApproximation::mean()" << std::endl;
  //  abort_handler(-1);
  //}

  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedMean & 1))
    return numericalMoments[0];

  Real mean = expectation(expansionType1Coeffs, expansionType2Coeffs);
  if (std_mode)
    { numericalMoments[0] = mean; computedMean |= 1; }
  return mean;
}


/** In this case, a subset of the expansion variables are random
    variables and the mean of the expansion involves integration over
    this subset and evaluation over the subset's complement.  For the
    linear sums of tensor interpolants within a sparse interpolant,
    the expectation can be taken inside the sum and we can simply add
    up the tensor mean contributions. */
Real NodalInterpPolyApproximation::mean(const RealVector& x)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if (all_mode && (computedMean & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevMean))
    return numericalMoments[0];

  Real mean;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    SizetArray colloc_index; // empty -> default indexing
    mean = tensor_product_mean(x, tpq_driver->level_index(),
			       tpq_driver->collocation_key(), colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
    const UShort2DArray& sm_mi        = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs    = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key   = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index = csg_driver->collocation_indices();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    mean = 0.;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	mean += sm_coeffs[i] *
	  tensor_product_mean(x, sm_mi[i], colloc_key[i], colloc_index[i]);
    break;
  }
  }
  if (all_mode)
    { numericalMoments[0] = mean; computedMean |= 1; xPrevMean = x; }

  return mean;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  In
    this case, the derivative of the expectation is the expectation of
    the derivative.  The mixed derivative case (some design variables
    are inserted and some are augmented) requires no special treatment. */
const RealVector& NodalInterpPolyApproximation::mean_gradient()
{
  // d/ds <R> = <dR/ds>

  // Error check for required data
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	  << "InterpPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedMean & 2))
    return meanGradient;

  const RealVector& t1_wts = driver_rep->type1_weight_sets();
  size_t i, j, num_colloc_pts = t1_wts.length(),
    num_deriv_vars = expansionType1CoeffGrads.numRows();
  if (meanGradient.length() != num_deriv_vars)
    meanGradient.sizeUninitialized(num_deriv_vars);
  meanGradient = 0.;
  for (i=0; i<num_colloc_pts; ++i) {
    const Real& t1_wt_i = t1_wts[i];
    for (j=0; j<num_deriv_vars; ++j)
      meanGradient[j] += expansionType1CoeffGrads(j,i) * t1_wt_i;
  }
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
    are inserted (derivatives are obtained from expansionType1CoeffGrads).  
    For the linear sums of tensor interpolants within a sparse 
    interpolant, the expectation can be taken inside the sum and we 
    can simply add up the tensor mean gradient contributions. */
const RealVector& NodalInterpPolyApproximation::
mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  // if already computed, return previous result
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if ( all_mode && (computedMean & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevMeanGrad) ) // && dvv == dvvPrev)
    switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
    case QUADRATURE:           return tpMeanGrad;   break;
    case COMBINED_SPARSE_GRID: return meanGradient; break;
    }

  // compute the gradient of the mean
  if (all_mode) { computedMean |=  2; xPrevMeanGrad = x; }
  else            computedMean &= ~2; // deactivate 2-bit: protect mixed usage
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_mean_gradient(x, tpq_driver->level_index(),
					tpq_driver->collocation_key(),
					colloc_index, dvv);
    break;
  }
  case COMBINED_SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    if (meanGradient.length() != num_deriv_vars)
      meanGradient.sizeUninitialized(num_deriv_vars);
    meanGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      if (coeff) {
	const RealVector& tpm_grad
	  = tensor_product_mean_gradient(x, sm_mi[i], colloc_key[i],
					 colloc_index[i], dvv);
	for (j=0; j<num_deriv_vars; ++j)
	  meanGradient[j] += coeff * tpm_grad[j];
      }
    }
    return meanGradient;
    break;
  }
}


/** In this case, all expansion variables are random variables and the
    variance of the expansion uses an interpolation of central
    products (INTERPOLATION_OF_PRODUCTS is the only option for the
    standard expansion mode).  Since we reinterpolate the central
    products such that we retain the linear sums of tensor
    interpolants within a sparse interpolant, the expectation simply
    involves a summation over the sparse integration weights. */
Real NodalInterpPolyApproximation::
covariance(PolynomialApproximation* poly_approx_2)
{
  NodalInterpPolyApproximation* nip_approx_2
    = (NodalInterpPolyApproximation*)poly_approx_2;
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  bool same = (this == nip_approx_2),
    std_mode = data_rep->nonRandomIndices.empty();

  // Error check for required data
  if ( !expansionCoeffFlag || ( !same && !nip_approx_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  // TO DO:
  //if (!driver_rep->track_unique_product_weights()) {
  //  PCerr << "Error: unique product weights required in "
  //	  << "NodalInterpPolyApproximation::covariance()" << std::endl;
  //  abort_handler(-1);
  //}

  if (same && std_mode && (computedVariance & 1))
    return numericalMoments[1];

  // compute mean_1,mean_2 first, then compute covariance as
  // wt_prod*(coeff1-mean_1)*(coeff2-mean_2) in order to avoid precision
  // loss from computing covariance as <R_i R_j> - \mu_i \mu_j
  // Note: compute_statistics() in dakota/src/NonDExpansion.C orders calls
  //       to reduce repetition in moment calculations.
  Real mean_1 = mean(), mean_2  = (same) ? mean_1 : nip_approx_2->mean();
  const RealVector& t1_coeffs_2 = nip_approx_2->expansionType1Coeffs;
  const RealVector& t1_wts      = driver_rep->type1_weight_sets();
  Real covar = 0.;
  size_t i, j, num_colloc_pts = t1_wts.length(), num_v = sharedDataRep->numVars;
  if (data_rep->basisConfigOptions.useDerivs) {
    const RealMatrix& t2_coeffs_2 = nip_approx_2->expansionType2Coeffs;
    const RealMatrix& t2_wts      = driver_rep->type2_weight_sets();
    for (i=0; i<num_colloc_pts; ++i) {
      // type1 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      Real coeff1_i_mm1 = expansionType1Coeffs[i] - mean_1,
	  coeff1_2i_mm2 = t1_coeffs_2[i]          - mean_2;
      covar += coeff1_i_mm1 * coeff1_2i_mm2 * t1_wts[i];
      // type2 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      // --> interpolated gradients are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
      const Real *coeff2_i  = expansionType2Coeffs[i],
	         *coeff2_2i = t2_coeffs_2[i], *t2_wt_i = t2_wts[i];
      for (j=0; j<num_v; ++j)
	covar  += (coeff1_i_mm1 * coeff2_2i[j] + coeff1_2i_mm2 * coeff2_i[j])
	       *  t2_wt_i[j];
    }
  }
  else
    for (i=0; i<num_colloc_pts; ++i)
      covar += (expansionType1Coeffs[i] - mean_1) * (t1_coeffs_2[i] - mean_2)
	*  t1_wts[i];

  if (same && std_mode)
    { numericalMoments[1] = covar; computedVariance |= 1; }
  return covar;
}


/** In this case, a subset of the expansion variables are random variables and
    the variance of the expansion involves integration over this subset.  Three
    cases are at least partially implemented (type 2 covariance	contributions
    are not currently implemented for PRODUCT_OF_INTERPOLANTS approaches, and
    variance_gradient() is not implemented for PRODUCT_OF_INTERPOLANTS_FULL):
    (1) REINTERPOLATION_OF_PRODUCTS: (R-\mu)^2 is (re)interpolated, allowing
        sparse grid covariance to use a linear sum of tensor contributions.
	This is most consistent with covariance for the standard variables mode.
	The penalty is that the reinterpolated function is higher order.
    (1) INTERPOLATION_OF_PRODUCTS: (R-\mu)^2 is interpolated, allowing sparse
        grid covariance to use a linear sum of tensor contributions.  This is
	most consistent with covariance for the standard variables mode.
    (2) PRODUCT_OF_INTERPOLANTS_FAST: the covariance of the existing sparse
        interpolant is computed without reinterpolation, but only the individual
	covariances from matching tensor grids are included.  Covariance from
	mismatched tensor grids is neglected, and the mean reference values are
	individual tensor means, not the "total" mean.  Due to some fortuitous
	cancellations in the Smolyak construction that have not been fully
	analyzed, this shortcut/fast approach is often quite accurate.
    (3) PRODUCT_OF_INTERPOLANTS_FULL: the covariance of the existing sparse
        interpolant is computed without reinterpolation, the sparse covariance
	among all tensor grid combinations is included, and the total mean is
	used as the reference value for all products.  The covariance among
	mismatched tensor grids requires the evaluation of expectations of
	mismatched interpolation polynomials (matched order basis polynomials
	admit a Kronecker delta simplification).  This option has been the most
	accurate in testing but, despite attempts at pre-computation and fast
	lookup of these products, has been too expensive for large grids. */
Real NodalInterpPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  NodalInterpPolyApproximation* nip_approx_2
    = (NodalInterpPolyApproximation*)poly_approx_2;
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  //IntegrationDriver* driver_rep = data_rep->driverRep;
  bool same = (this == nip_approx_2),
    all_mode = !data_rep->nonRandomIndices.empty();

  // Error check for required data
  if ( !expansionCoeffFlag || ( !same && !nip_approx_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if ( same && all_mode && (computedVariance & 1) &&
       data_rep->match_nonrandom_vars(x, xPrevVar) )
    return numericalMoments[1];

  Real covar;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    SizetArray colloc_index; // empty -> default indexing

    switch (data_rep->momentInterpType) {
    case REINTERPOLATION_OF_PRODUCTS: {
      // compute higher-order covariance grid, if not already computed
      const UShortArray& lev_index = tpq_driver->level_index();
      reinterpolated_level(lev_index);
      // compute TP covariance with reinterpolation on the higher-order grid
      covar = tensor_product_covariance(x, lev_index,
	tpq_driver->collocation_key(), colloc_index, nip_approx_2);
      break;
    }
    case INTERPOLATION_OF_PRODUCTS:
    case PRODUCT_OF_INTERPOLANTS_FAST: case PRODUCT_OF_INTERPOLANTS_FULL:
      covar = tensor_product_covariance(x, tpq_driver->level_index(),
	tpq_driver->collocation_key(), colloc_index, nip_approx_2);
      break;
    }
    break;
  }

  // While we can collapse the Smolyak recursion and combine the weights in the
  // distinct variables case, we cannot do this here for the all_variables case
  // since the non-integrated interpolation polynomial portions are not constant
  // and are coupled with the weight combination.
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices     = sm_coeffs.size();
    covar = 0.;
    switch (data_rep->momentInterpType) {
    // Reinterpolation of covariance function in non-integrated dimensions.
    case REINTERPOLATION_OF_PRODUCTS:
      for (i=0; i<num_smolyak_indices; ++i)
	if (sm_coeffs[i]) {
	  reinterpolated_level(sm_mi[i]);
	  covar += sm_coeffs[i] *
	    tensor_product_covariance(x, sm_mi[i], colloc_key[i],
				      colloc_index[i], nip_approx_2);
	}
      break;
    // For an interpolation of products (which captures product cross-terms),
    // a sum of tensor-product covariances is correct and straightforward.
    // Note that tensor_product_covariance() is overloaded to specialize its
    // logic based on data_rep->momentInterpType.
    case INTERPOLATION_OF_PRODUCTS:
    // For a fast product of interpolants, cross-terms are neglected and
    // results are approximate.  In addition, the mean reference point for
    // each covariance contribution is the tensor mean, not the total mean.
    case PRODUCT_OF_INTERPOLANTS_FAST:
      for (i=0; i<num_smolyak_indices; ++i)
	if (sm_coeffs[i])
	  covar += sm_coeffs[i] * tensor_product_covariance(x, sm_mi[i],
	    colloc_key[i], colloc_index[i], nip_approx_2);
      break;
    // For a full product of interpolants, the rigorous approach must manage
    // the expectation of products of mixed order interpolation polynomials
    // from mixed order tensor products.  These are pre-computed and non-zero
    // products are tabulated for fast lookup, but this approach is still the
    // most expensive and has not been fully implemented for all cases.
    case PRODUCT_OF_INTERPOLANTS_FULL:
      data_rep->update_nonzero_basis_products(sm_mi);
      for (i=0; i<num_smolyak_indices; ++i) {
	int sm_coeff_i = sm_coeffs[i];
	if (sm_coeff_i) {   // one of each diagonal term
#ifdef DEBUG
	  Real tp_covar
	    = tensor_product_covariance(x, sm_mi[i], colloc_key[i],
					colloc_index[i], nip_approx_2);
	  covar += sm_coeff_i * sm_coeff_i * tp_covar;
	  PCout << "Diagonal covar: sm_coeffs[" << i << "] = " << sm_coeff_i
		<< " tp_covar = " << tp_covar << '\n';
#else
	  covar += sm_coeff_i * sm_coeff_i *
	    tensor_product_covariance(x, sm_mi[i], colloc_key[i],
				      colloc_index[i], nip_approx_2);
#endif // DEBUG
	  for (j=0; j<i; ++j)
	    if (sm_coeffs[j]) { // two of each off-diagonal term
#ifdef DEBUG
	      Real tp_covar =
		tensor_product_covariance(x,
		  sm_mi[i], colloc_key[i], colloc_index[i],
		  sm_mi[j], colloc_key[j], colloc_index[j], nip_approx_2);
	      covar += 2. * sm_coeff_i * sm_coeffs[j] * tp_covar;
	      PCout << "Off-diagonal covar: sm_coeffs[" << i << "] = "
		    << sm_coeff_i << " sm_coeffs[" << j << "] = "
		    << sm_coeffs[j] << " tp_covar = " << tp_covar << '\n';
#else
	      covar += 2. * sm_coeff_i * sm_coeffs[j] *
		tensor_product_covariance(x,
		  sm_mi[i], colloc_key[i], colloc_index[i],
		  sm_mi[j], colloc_key[j], colloc_index[j], nip_approx_2);
#endif // DEBUG
	    }
	}
      }
      break;
    }
    break;
  }
  }
  if (same && all_mode)
    { numericalMoments[1] = covar; computedVariance |= 1; xPrevVar = x; }
  return covar;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  The
    mixed derivative case (some design/epistemic variables are
    inserted and some are augmented) requires no special treatment.
    Since we reinterpolate the central products
    (INTERPOLATION_OF_PRODUCTS is the only option for the standard
    expansion mode) such that we retain the linear sums of tensor
    interpolants within a sparse interpolant, the expectation gradient
    simply involves a summation over the sparse integration weights. */
const RealVector& NodalInterpPolyApproximation::variance_gradient()
{
  // Error check for required data
  if (!expansionCoeffFlag || !expansionCoeffGradFlag) {
    PCerr << "Error: insufficient expansion coefficient data in NodalInterp"
	  << "PolyApproximation::variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (computedVariance & 2))
    return varianceGradient;

  const RealVector& t1_wts = driver_rep->type1_weight_sets();
  size_t i, j, num_colloc_pts = t1_wts.length(),
    num_deriv_vars = expansionType1CoeffGrads.numRows();
  if (varianceGradient.length() != num_deriv_vars)
    varianceGradient.sizeUninitialized(num_deriv_vars);
  varianceGradient = 0.;

  Real mean_1 = mean();
  // See Eq. 6.23 in Theory Manual: grad of variance incorporates grad of mean
  for (i=0; i<num_colloc_pts; ++i) {
    Real term_i = 2. * (expansionType1Coeffs[i] - mean_1) * t1_wts[i];
    for (j=0; j<num_deriv_vars; ++j)
      varianceGradient[j] += term_i * expansionType1CoeffGrads(j,i);
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
    expansionType1CoeffGrads). */
const RealVector& NodalInterpPolyApproximation::
variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  //IntegrationDriver* driver_rep = data_rep->driverRep;

  // if already computed, return previous result
  bool all_mode = !data_rep->nonRandomIndices.empty();
  if ( all_mode && (computedVariance & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevVarGrad) ) // && dvv == dvvPrev)
    switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
    case QUADRATURE:           return tpVarianceGrad;   break;
    case COMBINED_SPARSE_GRID: return varianceGradient; break;
    }

  if (all_mode) { computedVariance |=  2; xPrevVarGrad = x; }
  else            computedVariance &= ~2;//deactivate 2-bit: protect mixed usage
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    SizetArray colloc_index; // empty -> default indexing
    switch (data_rep->momentInterpType) {
    case REINTERPOLATION_OF_PRODUCTS: {
      const UShortArray& lev_index = tpq_driver->level_index();
      // compute ord/lev for higher-order covariance grid over non-random vars
      reinterpolated_level(lev_index);
      // compute variance grad with reinterpolation on the higher-order grid
      return tensor_product_variance_gradient(x, lev_index,
        tpq_driver->collocation_key(), colloc_index, dvv);
      break;
    }
    case INTERPOLATION_OF_PRODUCTS:
    case PRODUCT_OF_INTERPOLANTS_FAST: case PRODUCT_OF_INTERPOLANTS_FULL:
      return tensor_product_variance_gradient(x, tpq_driver->level_index(),
        tpq_driver->collocation_key(), colloc_index, dvv);
      break;
    }
    break;
  }
  case COMBINED_SPARSE_GRID: {
    size_t num_deriv_vars = dvv.size();
    if (varianceGradient.length() != num_deriv_vars)
      varianceGradient.sizeUninitialized(num_deriv_vars);
    varianceGradient = 0.;
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
    const UShort2DArray& sm_mi           = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs       = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key      = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index    = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size(); int coeff;
    switch (data_rep->momentInterpType) {
    // Reinterpolation of covariance function in non-integrated dimensions.
    case REINTERPOLATION_OF_PRODUCTS:
      for (i=0; i<num_smolyak_indices; ++i)
	if (coeff = sm_coeffs[i]) {
	  reinterpolated_level(sm_mi[i]);
	  const RealVector& tpv_grad =
	    tensor_product_variance_gradient(x, sm_mi[i], colloc_key[i],
					     colloc_index[i], dvv);
	  for (j=0; j<num_deriv_vars; ++j)
	    varianceGradient[j] += coeff * tpv_grad[j];
	}
      break;
    // For an interpolation of products (which captures product cross-terms),
    // a sum of tensor-product covariances is correct and straightforward.
    case INTERPOLATION_OF_PRODUCTS:
    // For a fast product of interpolants, cross-terms are neglected.
    case PRODUCT_OF_INTERPOLANTS_FAST:
      for (i=0; i<num_smolyak_indices; ++i)
	if (coeff = sm_coeffs[i]) {
	  const RealVector& tpv_grad =
	    tensor_product_variance_gradient(x, sm_mi[i], colloc_key[i],
					     colloc_index[i], dvv);
	  for (j=0; j<num_deriv_vars; ++j)
	    varianceGradient[j] += coeff * tpv_grad[j];
	}
      break;
    case PRODUCT_OF_INTERPOLANTS_FULL:
      // TO DO
      break;
    }
    return varianceGradient;
    break;
  }
  }
}


Real NodalInterpPolyApproximation::
expectation(const RealVector& t1_coeffs, const RealVector& t1_wts,
	    const RealMatrix& t2_coeffs, const RealMatrix& t2_wts)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  //IntegrationDriver* driver_rep = data_rep->driverRep;

  Real integral = 0.;
  size_t num_colloc_pts = t1_coeffs.length();
  if (data_rep->basisConfigOptions.useDerivs) {
    size_t i, j, num_v = t2_coeffs.numRows();
    for (i=0; i<num_colloc_pts; ++i) {
      integral += t1_coeffs[i] * t1_wts[i];
      const Real* coeff2_i = t2_coeffs[i];
      const Real*  t2_wt_i = t2_wts[i];
      for (j=0; j<num_v; ++j)
	integral += coeff2_i[j] * t2_wt_i[j];
    }
  }
  else
    for (size_t i=0; i<num_colloc_pts; ++i)
      integral += t1_coeffs[i] * t1_wts[i];

  return integral;
}


/** Computes the specifics of a higher order grid for reinterpolating
    covariance over dimensions that will not be integrated. */
void NodalInterpPolyApproximation::
reinterpolated_level(const UShortArray& lev_index)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  driver_rep->reinterpolated_tensor_grid(lev_index, data_rep->nonRandomIndices);
  data_rep->
    update_tensor_interpolation_basis(driver_rep->reinterpolated_level_index(),
				      data_rep->nonRandomIndices);
}


void NodalInterpPolyApproximation::
integrate_response_moments(size_t num_moments)
{
  // In this case, we are constrained to use the original collocation points
  // (for which response values are available) within the numerical integration.
  // In the case where interpolatory rules (e.g., Clenshaw-Curtis) have been 
  // used, the integrand accuracy may suffer.  For this reason, the expansion
  // moments (which integrate the interpolant expansion instead of the response,
  // and therefore may employ alternate rules) are preferred.

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in InterpPoly"
	  << "Approximation::integrate_response_moments()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  if (numericalMoments.length() != num_moments)
    numericalMoments.sizeUninitialized(num_moments);
  if (data_rep->basisConfigOptions.useDerivs)
    integrate_moments(expansionType1Coeffs, expansionType2Coeffs,
		      driver_rep->type1_weight_sets(),
		      driver_rep->type2_weight_sets(), numericalMoments);
  else
    integrate_moments(expansionType1Coeffs, driver_rep->type1_weight_sets(),
		      numericalMoments);
}


void NodalInterpPolyApproximation::
integrate_expansion_moments(size_t num_moments)
{
  // In this case, we are integrating the interpolant expansion instead of the
  // response and are *not* constrained to use the original collocation points
  // within the numerical integration.  Thererfore, we replace any interpolatory
  // rules (e.g., Clenshaw-Curtis) with integration rules to improve the
  // integrand accuracy; for this reason, the expansion moments are preferred.
  // Since we are limiting expansion integration to second moments, Gaussian
  // rules of comparable order are generally sufficient (skewness and kurtosis
  // would require higher order rules).

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in InterpPoly"
	  << "Approximation::integrate_expansion_moments()"<< std::endl;
    abort_handler(-1);
  }
  if (expansionMoments.length() != num_moments)
    expansionMoments.sizeUninitialized(num_moments);

  // TO DO: evaluate moments 2/3/4 by evaluating the interpolant of R on higher
  // order grids so that the moment function (R-\mu)^k can be reinterpolated
  // exactly.  Generate higher order rules using driverRep and perform
  // integration here for access to value() and gradient().

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver *alt_driver = data_rep->expMomentIntDriver.driver_rep();
  bool alt_grid = (alt_driver != NULL);

  // Alternate quadrature on interpolant is strictly value-based.  A shared
  // level/order definition could occur in SharedNIPAData::allocate_data(),
  // but we initially compute level/order here on the fly on each QoI.
  if (alt_grid) {
    // synchronize the level/order between alternate and original driver
    if (data_rep->expConfigOptions.expCoeffsSolnApproach == QUADRATURE) {
      TensorProductDriver* tp_driver
	= (TensorProductDriver*)data_rep->driverRep;
      TensorProductDriver* tp_alt_driver = (TensorProductDriver*)alt_driver;
      // match #quad pts: new precision >= old precision
      // Note: Dakota uses quad scalar + dim_pref, Pecos uses aniso quad vector
      tp_alt_driver->quadrature_order(tp_driver->quadrature_order());
    }
    else {
      SparseGridDriver*     sg_driver = (SparseGridDriver*)data_rep->driverRep;
      SparseGridDriver* sg_alt_driver = (SparseGridDriver*)alt_driver;
      // level is matched for now; ignores nonlinear/linear growth mismatch
      sg_alt_driver->level(sg_driver->level());
      sg_alt_driver->anisotropic_weights(sg_driver->anisotropic_weights());
    }
    // compute the points and weights
    RealMatrix alt_pts;
    alt_driver->compute_grid(alt_pts);
    // integrate the interpolant through eval at the quadrature points
    size_t i, num_pts = alt_pts.numCols();
    RealVector t1_exp(num_pts);
    for (i=0; i<num_pts; ++i)
      t1_exp[i] = value(Teuchos::getCol(Teuchos::View, alt_pts, (int)i));
    integrate_moments(t1_exp, alt_driver->type1_weight_sets(),expansionMoments);
  }
  // Native quadrature on interpolant can be value-based or gradient-enhanced.
  // These approaches just replace the values of the response with the values
  // of the interpolant (integrates powers of the interpolant using the same
  // collocation rules/orders used to form the interpolant).
  else {
    size_t i, offset = 0, num_pts = surrData.points();
    bool anchor_pt = surrData.anchor();
    if (anchor_pt) { offset = 1; ++num_pts; }
    RealVector t1_exp(num_pts);
    if (data_rep->basisConfigOptions.useDerivs) { // gradient-enhanced native
      size_t num_v = data_rep->numVars;
      RealMatrix t2_exp(num_v, num_pts);
      for (i=0; i<num_pts; ++i) {
	const RealVector& c_vars = (anchor_pt && i == 0) ?
	  surrData.anchor_continuous_variables() :
	  surrData.continuous_variables(i-offset);
	t1_exp[i] = value(c_vars);
	Teuchos::setCol(gradient_basis_variables(c_vars), (int)i, t2_exp);
      }
      IntegrationDriver* driver_rep = data_rep->driverRep;
      integrate_moments(t1_exp, t2_exp, driver_rep->type1_weight_sets(),
			driver_rep->type2_weight_sets(), expansionMoments);
    }
    else { // value-based native quadrature
      for (i=0; i<num_pts; ++i)
	t1_exp[i] = (anchor_pt && i == 0) ? 
	  value(surrData.anchor_continuous_variables()) :
	  value(surrData.continuous_variables(i-offset));
      integrate_moments(t1_exp, data_rep->driverRep->type1_weight_sets(),
			expansionMoments);
    }
  }
#ifdef DEBUG
  PCout << "Expansion moments type 1 coefficients:\n";
  write_data(PCout, t1_exp);
  PCout << "Expansion moments:\n";
  write_data(PCout, expansionMoments);
#endif // DEBUG
}


/** Computes the variance of component functions. Assumes that all
    subsets of set_value have been computed in advance which will be
    true so long as the partial_variance is called following
    appropriate enumeration of set value  */
void NodalInterpPolyApproximation::
compute_partial_variance(const BitArray& set_value)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  // Perform inner integral over complementary set u' to form new weighted
  // coefficients h; then perform outer integral of h^2 over set u
  Real& variance = partialVariance[data_rep->sobolIndexMap[set_value]];
  variance = member_integral(set_value, 0.);// center = 0
#ifdef VBD_DEBUG
  PCout << "Partial variance = " << variance;
#endif // VBD_DEBUG

  // compute proper subsets and subtract their contributions
  InterpPolyApproximation::compute_partial_variance(set_value);
#ifdef VBD_DEBUG
  PCout << " (raw) " << variance << " (minus subsets)\n";
#endif // VBD_DEBUG
}


void NodalInterpPolyApproximation::compute_total_sobol_indices()
{
  // Compute the total expansion mean and variance.  For standard mode, these
  // are likely already available, as managed by computed{Mean,Variance} data
  // reuse trackers.  For all vars mode, they are computed without partial
  // integration restricted to the random indices.  If negligible variance
  // (deterministic test fn) or negative variance (poor sparse grid resolution),
  // then attribution of this variance is suspect.  Note: zero is a good choice
  // since it drops out from anisotropic refinement based on a response-average
  // of total Sobol' indices.
  Real total_variance = variance();
  if (total_variance <= SMALL_NUMBER)
    { totalSobolIndices = 0.; return; }
  Real total_mean = mean();

  size_t j, num_v = sharedDataRep->numVars;
  BitArray complement_set(num_v);
  for (j=0; j<num_v; ++j) {
    // define complement_set that includes all but index of interest
    complement_set.set(); complement_set.flip(j);

    // Perform inner integral over complementary set u' to form new weighted
    // coeffs h; then perform outer integral of (h-total_mean)^2 over set u
    totalSobolIndices[j] = 1. - member_integral(complement_set, total_mean)
                         / total_variance;
  }
}


/** Forms a new interpolant h over member variables after integrating
    out the non-member variables.  Then finds the variance by
    integrating (h-mean)^2 over the member variables, where the mean
    can be zero for non-central/raw moment cases. */
Real NodalInterpPolyApproximation::
member_integral(const BitArray& member_bits, Real mean)
{
  // Follows Tang, Iaccarino, Eldred (conference paper AIAA-2010-2922)

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  //IntegrationDriver* driver_rep = data_rep->driverRep;

  size_t i, j, num_member_coeffs, num_v = sharedDataRep->numVars;
  SizetList member_indices;
  for (j=0; j<num_v; ++j)
    if (member_bits[j])
      member_indices.push_back(j);

  Real integral = 0.;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = data_rep->tpq_driver();
    const UShortArray&    lev_index = tpq_driver->level_index();
    SizetArray colloc_index; // empty -> default indexing

    // Perform inner integral over complementary set u' to form new weighted
    // coefficients h (stored as member_coeffs)
    RealVector member_t1_coeffs, member_t1_wts;
    RealMatrix member_t2_coeffs, member_t2_wts;
    UShort2DArray member_colloc_key;
    SizetArray member_colloc_index, empty_colloc_index;
    member_coefficients_weights(member_bits, tpq_driver->quadrature_order(),
      lev_index, tpq_driver->collocation_key(), colloc_index, member_t1_coeffs,
      member_t1_wts, member_t2_coeffs, member_t2_wts, member_colloc_key,
      member_colloc_index);

    // Perform outer integral over set u
    num_member_coeffs = member_t1_coeffs.length();
    for (i=0; i<num_member_coeffs; ++i) {
      const RealVector& c_vars
	= surrData.continuous_variables(member_colloc_index[i]);
      Real h_t1_coeff_i_mm = data_rep->
	tensor_product_value(c_vars, member_t1_coeffs, member_t2_coeffs,
			     lev_index, member_colloc_key, empty_colloc_index,
			     member_indices) - mean;
      integral += h_t1_coeff_i_mm * h_t1_coeff_i_mm * member_t1_wts[i];
      if (data_rep->basisConfigOptions.useDerivs) {
	// type2 interpolation of h^2 = 2 h h'
	const RealVector& h_t2_coeffs_i = data_rep->
	  tensor_product_gradient_basis_variables(c_vars, member_t1_coeffs,
	    member_t2_coeffs, lev_index, member_colloc_key, empty_colloc_index,
	    member_indices);
	Real *m_t2_wts_i = member_t2_wts[i];
	for (j=0; j<num_v; ++j)
	  integral += 2. * h_t1_coeff_i_mm * h_t2_coeffs_i[j] * m_t2_wts_i[j];
      }
    }
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = data_rep->csg_driver();
    const IntArray&            sm_coeffs = csg_driver->smolyak_coefficients();
    const UShort2DArray&        sm_index = csg_driver->smolyak_multi_index();
    const UShort3DArray&      colloc_key = csg_driver->collocation_key();
    const Sizet2DArray&     colloc_index = csg_driver->collocation_indices();
    size_t num_smolyak_indices = sm_coeffs.size(),
           num_member_indices  = member_indices.size();
    UShortArray quad_order;

    // Perform inner integral over complementary set u' to form new weighted
    // coefficients h (stored as member {t1,t2} coeffs).  Precompute all
    // coefficients for h since they are needed below for value()/gradient().
    RealVectorArray member_t1_coeffs(num_smolyak_indices),
                    member_t1_wts(num_smolyak_indices);
    RealMatrixArray member_t2_coeffs(num_smolyak_indices),
                    member_t2_wts(num_smolyak_indices);
    UShort3DArray   member_colloc_key(num_smolyak_indices);
    Sizet2DArray    member_colloc_index(num_smolyak_indices);
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i]) {
        csg_driver->level_to_order(sm_index[i], quad_order);
	member_coefficients_weights(member_bits, quad_order, sm_index[i],
	  colloc_key[i], colloc_index[i], member_t1_coeffs[i], member_t1_wts[i],
	  member_t2_coeffs[i], member_t2_wts[i], member_colloc_key[i],
	  member_colloc_index[i]);
      }

    // Perform outer integral for (h-mean)^2.  The key is to evaluate the
    // sparse interpolant for h using value() and gradient_basis_variables(),
    // such that the response {value,gradient} for the reduced-dimension
    // interpolant has the correct aggregated coefficients.  Then we linearly
    // sum the member wts over the Smolyak coeffs (the reduced-dimension
    // sparse grid is not "combined").

    // For efficiency, we eliminate redundant value()/grad() calls using a
    // std::map with unique key created from member variable components of
    // sm_index and member_colloc_key
    std::map<UShortArray, RealRealVectorPair> member_map;
    std::map<UShortArray, RealRealVectorPair>::iterator member_map_it;
    UShortArray member_map_key(2*num_member_indices);
    RealRealVectorPair val_grad_pr;

    Real h_t1_coeff_ij_mm;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i]) {
	num_member_coeffs = member_t1_coeffs[i].length();
	update_member_key(sm_index[i], member_indices, member_map_key, 0);
	UShort2DArray& m_c_key_i   = member_colloc_key[i];
	SizetArray&    m_c_index_i = member_colloc_index[i];
	for (j=0; j<num_member_coeffs; ++j) {
	  const RealVector& c_vars
	    = surrData.continuous_variables(m_c_index_i[j]);
	  // create key for this member coeff; see if val/grad already computed
	  update_member_key(m_c_key_i[j], member_indices, member_map_key,
			    num_member_indices);
	  member_map_it = member_map.find(member_map_key);
	  bool found = (member_map_it != member_map.end());
	  // retrieve or compute type1 coefficient for h-mean
	  if (found)
	    h_t1_coeff_ij_mm = member_map_it->second.first;
	  else
	    val_grad_pr.first = h_t1_coeff_ij_mm =
	      value(c_vars, member_t1_coeffs, member_t2_coeffs,
		    member_colloc_key, member_indices) - mean;
	  // increment integral with type1 contribution of (h-mean)^2
	  integral += h_t1_coeff_ij_mm * h_t1_coeff_ij_mm
	           *  member_t1_wts[i][j] * sm_coeffs[i];
	  if (data_rep->basisConfigOptions.useDerivs) {
	    // retrieve or compute type2 coefficient for h
	    if (!found)
	      val_grad_pr.second = 
		gradient_basis_variables(c_vars, member_t1_coeffs,
		member_t2_coeffs, member_colloc_key, member_indices);
	    const RealVector& h_t2_coeffs_ij = (found) ?
	      member_map_it->second.second : val_grad_pr.second;
	    // increment integral with type2 contribution of (h-mean)^2 = 2 h h'
	    Real *m_t2_wts_ij = member_t2_wts[i][j];
	    for (size_t k=0; k<num_v; ++k)
	      integral += 2. * h_t1_coeff_ij_mm * h_t2_coeffs_ij[k]
		       *  m_t2_wts_ij[k] * sm_coeffs[i];
	  }
	  if (!found)
	    member_map[member_map_key] = val_grad_pr;
	}
      }
    break;
  }
  }

  return integral;
}


void NodalInterpPolyApproximation::
member_coefficients_weights(const BitArray& member_bits,
			    const UShortArray& quad_order,
			    const UShortArray& lev_index,
			    const UShort2DArray& colloc_key,
			    const SizetArray& colloc_index,
			    RealVector& member_t1_coeffs,
			    RealVector& member_t1_wts,
			    RealMatrix& member_t2_coeffs,
			    RealMatrix& member_t2_wts,
			    UShort2DArray& member_colloc_key,
			    SizetArray& member_colloc_index)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  // get number of expansion coeffs in member-variable-only expansion and
  // precompute indexing factors, since they only depend on j
  size_t i, j, num_v = sharedDataRep->numVars,
    num_member_coeffs = 1;  // # exp coeffs in member-var-only exp
  SizetArray indexing_factor; // factors for indexing member coeffs,wts
  for (j=0; j<num_v; ++j)
    if (member_bits[j]) {
      indexing_factor.push_back(num_member_coeffs); // for member_index below
      num_member_coeffs *= quad_order[j];
    }

  // Size vectors to store new coefficients
  member_t1_coeffs.size(num_member_coeffs); // init to 0
  member_t1_wts.size(num_member_coeffs);    // init to 0
  if (data_rep->basisConfigOptions.useDerivs) {
    member_t2_coeffs.shape(num_v, num_member_coeffs); // init to 0
    member_t2_wts.shape(num_v, num_member_coeffs);    // init to 0
  }
  member_colloc_key.resize(num_member_coeffs);
  member_colloc_index.resize(num_member_coeffs);

  // Perform inner integral over complementary set u' (non-member vars) to
  // form new weighted expansion coefficients h (stored as member_t1_coeffs)
  size_t num_tp_pts = colloc_key.size(), member_index, c_index, cntr;
  Real member_wt, nonmember_wt;
  for (i=0; i<num_tp_pts; ++i) {
    // convert key_i to a corresponding index on member coeffs and wts.  A
    // more rigorous approach would define a mapping to a member-variable
    // multi-index/collocation key (from collapsing non-member indices/keys),
    // but the current approach is simpler and sufficiently general for this
    // purpose.  Note: the precise sequence of member indices is not important
    // in its current usage -- it just must enumerate the unique members.
    const UShortArray& key_i = colloc_key[i];
    for (j=0, member_index=0, cntr=0; j<num_v; ++j)
      if (member_bits[j])
	member_index += key_i[j] * indexing_factor[cntr++];

    // integrate out the nonmember dimensions and aggregate with type1 coeffs
    data_rep->
      type1_weight(key_i, lev_index, member_bits, member_wt, nonmember_wt);
    c_index = (colloc_index.empty()) ? i : colloc_index[i];
    member_t1_coeffs[member_index]
      += nonmember_wt * expansionType1Coeffs[c_index];
    // reduced dimension data updated more times than necessary, but tracking
    // this redundancy would be more expensive/complicated.  Note: non-member
    // key and c_vars data may change, but member data should be consistent.
    member_t1_wts[member_index]       = member_wt;
    member_colloc_key[member_index]   = key_i;  // links back to interp poly's
    member_colloc_index[member_index] = c_index;// links back to surrData c_vars

    // now do the same for the type2 coeffs and weights
    if (data_rep->basisConfigOptions.useDerivs) {
      Real *m_t2_coeffs_i = member_t2_coeffs[member_index],
	   *m_t2_wts_i    = member_t2_wts[member_index];
      const Real *t2_coeffs_i = expansionType2Coeffs[c_index];
      for (j=0; j<num_v; ++j) {
	data_rep->type2_weight(j, key_i, lev_index, member_bits,
			       member_wt, nonmember_wt);
	m_t2_coeffs_i[j] += nonmember_wt * t2_coeffs_i[j];
	m_t2_wts_i[j]    =  member_wt;
      }
    }
  }
#ifdef VBD_DEBUG
  PCout << "member_bits: " << member_bits << '\n'; // MSB->LSB: order reversed
  PCout << "member_t1_coeffs:\n"; write_data(PCout, member_t1_coeffs);
  PCout << "member_t1_wts:\n";    write_data(PCout, member_t1_wts);
  if (data_rep->basisConfigOptions.useDerivs) {
    PCout << "member_t2_coeffs:\n";
    write_data(PCout, member_t2_coeffs, false, true, true);
    PCout << "member_t2_wts:\n";
    write_data(PCout, member_t2_wts,    false, true, true);
  }
#endif // VBD_DEBUG
}

} // namespace Pecos
