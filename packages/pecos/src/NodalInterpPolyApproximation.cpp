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
//#define INTERPOLATION_TEST

namespace Pecos {


void NodalInterpPolyApproximation::allocate_arrays()
{
  InterpPolyApproximation::allocate_arrays();

  size_t num_colloc_pts = modSurrData.points(),
    num_deriv_vars = modSurrData.num_derivative_variables();
  if (expansionCoeffFlag) { // else coeff arrays not entered into maps
    RealVector& exp_t1_coeffs = expT1CoeffsIter->second;
    if (exp_t1_coeffs.length() != num_colloc_pts)
      exp_t1_coeffs.sizeUninitialized(num_colloc_pts);
    SharedNodalInterpPolyApproxData* data_rep
      = (SharedNodalInterpPolyApproxData*)sharedDataRep;
    if (data_rep->basisConfigOptions.useDerivs) {
      // find existing or create new
      RealMatrix& exp_t2_coeffs = expT2CoeffsIter->second;
      if ( exp_t2_coeffs.numRows() != num_deriv_vars ||
	   exp_t2_coeffs.numCols() != num_colloc_pts )
	exp_t2_coeffs.shapeUninitialized(num_deriv_vars, num_colloc_pts);
    }
  }
  if (expansionCoeffGradFlag) { // else coeff grads not entered into map
    // find existing or create new
    RealMatrix& exp_t1_coeff_grads = expT1CoeffGradsIter->second;
    if ( exp_t1_coeff_grads.numRows() != num_deriv_vars ||
	 exp_t1_coeff_grads.numCols() != num_colloc_pts )
      exp_t1_coeff_grads.shapeUninitialized(num_deriv_vars, num_colloc_pts);
  }

  // checking num_colloc_pts is insufficient due to anisotropy --> changes in
  // anisotropic weights could move points around without changing the total.
  //bool update_exp_form =
  //  ( (expansionCoeffFlag && exp_t1_coeffs.length() != num_colloc_pts) ||
  //    (expansionCoeffGradFlag &&
  //     exp_t1_coeff_grads.numCols() != num_colloc_pts ) );
}


void NodalInterpPolyApproximation::compute_coefficients()
{
  PolynomialApproximation::compute_coefficients();
  if (!expansionCoeffFlag && !expansionCoeffGradFlag)
    return;

  allocate_arrays();

  const SDRArray& sdr_array = modSurrData.response_data();
  size_t num_colloc_pts = modSurrData.points(); int i;
  if (expansionCoeffFlag) {
    RealVector& exp_t1_coeffs = expT1CoeffsIter->second;
    RealMatrix& exp_t2_coeffs = expT2CoeffsIter->second;
    SharedNodalInterpPolyApproxData* data_rep
      = (SharedNodalInterpPolyApproxData*)sharedDataRep;
    bool derivs = data_rep->basisConfigOptions.useDerivs;
    for (i=0; i<num_colloc_pts; ++i) {
      exp_t1_coeffs[i] = sdr_array[i].response_function();
      // Note: gradients from DAKOTA already scaled in u-space Recast
      if (derivs)
	Teuchos::setCol(sdr_array[i].response_gradient(), i,
			exp_t2_coeffs);
    }
  }
  if (expansionCoeffGradFlag) {
    RealMatrix& exp_t1_coeff_grads = expT1CoeffGradsIter->second;
    for (i=0; i<num_colloc_pts; ++i)
      Teuchos::setCol(sdr_array[i].response_gradient(), i,
		      exp_t1_coeff_grads);
  }

#ifdef INTERPOLATION_TEST
  test_interpolation();
#endif

  clear_computed_bits();
}


void NodalInterpPolyApproximation::pop_coefficients(bool save_data)
{
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  update_active_iterators(data_rep->activeKey);

  size_t new_colloc_pts = modSurrData.points();
  if (expansionCoeffFlag) {
    RealVector& exp_t1_coeffs = expT1CoeffsIter->second;
    exp_t1_coeffs.resize(new_colloc_pts);
    if (data_rep->basisConfigOptions.useDerivs) {
      RealMatrix& exp_t2_coeffs = expT2CoeffsIter->second;
      size_t num_deriv_vars = exp_t2_coeffs.numRows();
      exp_t2_coeffs.reshape(num_deriv_vars, new_colloc_pts);//prune trailing pts
    }
  }
  if (expansionCoeffGradFlag) {
    RealMatrix& exp_t1_coeff_grads = expT1CoeffGradsIter->second;
    size_t num_deriv_vars = exp_t1_coeff_grads.numRows();
    exp_t1_coeff_grads.reshape(num_deriv_vars, new_colloc_pts);// prune trailing
  }

  clear_computed_bits();
}


void NodalInterpPolyApproximation::update_expansion_coefficients()
{
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  update_active_iterators(data_rep->activeKey);

  // TO DO: partial sync for new TP data set, e.g. update_surrogate_data() ?
  synchronize_surrogate_data(); // modSurrData updates required

  size_t old_colloc_pts, new_colloc_pts = modSurrData.points();
  const SDRArray& sdr_array = modSurrData.response_data();
  bool append
    = (data_rep->expConfigOptions.expCoeffsSolnApproach != QUADRATURE);
  if (expansionCoeffFlag) {
    RealVector& exp_t1_coeffs = expT1CoeffsIter->second;
    RealMatrix& exp_t2_coeffs = expT2CoeffsIter->second;
    old_colloc_pts = (append) ? exp_t1_coeffs.length() : 0;
    exp_t1_coeffs.resize(new_colloc_pts);
    if (data_rep->basisConfigOptions.useDerivs) {
      size_t num_deriv_vars = exp_t2_coeffs.numRows();
      exp_t2_coeffs.reshape(num_deriv_vars, new_colloc_pts);
    }
    for (int i=old_colloc_pts; i<new_colloc_pts; ++i) {
      const SurrogateDataResp& sdr = sdr_array[i];
      exp_t1_coeffs[i] = sdr.response_function();
      if (data_rep->basisConfigOptions.useDerivs)
	Teuchos::setCol(sdr.response_gradient(), i, exp_t2_coeffs);
    }
  }
  if (expansionCoeffGradFlag) {
    RealMatrix& exp_t1_coeff_grads = expT1CoeffGradsIter->second;
    old_colloc_pts = (append) ? exp_t1_coeff_grads.numCols() : 0;
    size_t num_deriv_vars = exp_t1_coeff_grads.numRows();
    exp_t1_coeff_grads.reshape(num_deriv_vars, new_colloc_pts);
    for (int i=old_colloc_pts; i<new_colloc_pts; ++i)
      Teuchos::setCol(sdr_array[i].response_gradient(), i,
		      exp_t1_coeff_grads);
  }

  clear_computed_bits();
}


void NodalInterpPolyApproximation::clear_inactive()
{
  std::map<UShortArray, RealVector>::iterator e1c_it
    = expansionType1Coeffs.begin();
  std::map<UShortArray, RealMatrix>::iterator e2c_it
    = expansionType2Coeffs.begin();
  std::map<UShortArray, RealMatrix>::iterator e1g_it
    = expansionType1CoeffGrads.begin();
  while (e1c_it != expansionType1Coeffs.end())
    if (e1c_it == expT1CoeffsIter) // preserve active
      { ++e1c_it; ++e2c_it; ++e1g_it; }
    else { // clear inactive: postfix increments manage iterator invalidations
      expansionType1Coeffs.erase(e1c_it++);
      expansionType2Coeffs.erase(e2c_it++);
      expansionType1CoeffGrads.erase(e1g_it++);
    }
}


void NodalInterpPolyApproximation::combine_coefficients()
{
  // update expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads} by adding or
  // multiplying stored expansion evaluated at current collocation points

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;

  // Coefficient combination is not dependent on active state, but active
  // iterators are used below to avoid additional key lookups in stored_*()
  update_active_iterators(data_rep->activeKey);

  // Note: computed bits are also cleared when refineStatsType is changed
  if (data_rep->expConfigOptions.refineStatsType == COMBINED_EXPANSION_STATS)
    clear_computed_bits();
  
  short combine_type = data_rep->expConfigOptions.combineType;
  std::map<UShortArray, RealVector>::iterator ec1_it;
  std::map<UShortArray, RealMatrix>::iterator ec2_it
    = expansionType2Coeffs.begin(), eg1_it = expansionType1CoeffGrads.begin();
  CombinedSparseGridDriver* csg_driver
    = (CombinedSparseGridDriver*)data_rep->driver();
  const RealMatrix& comb_var_sets = csg_driver->combined_variable_sets();
  size_t p, v, num_pts = comb_var_sets.numCols(),
    num_v   = comb_var_sets.numRows(), num_t2v = ec2_it->second.numRows(),
    num_t1v = eg1_it->second.numRows();
  bool use_derivs = data_rep->basisConfigOptions.useDerivs;
  if (expansionCoeffFlag) {
    if (combinedExpT1Coeffs.length() != num_pts)
      combinedExpT1Coeffs.size(num_pts);
    else combinedExpT1Coeffs = 0.;
    if (use_derivs) {
      if (combinedExpT2Coeffs.numRows() != num_pts ||
	  combinedExpT2Coeffs.numCols() != num_t2v)
	combinedExpT2Coeffs.shape(num_t2v, num_pts);
      else combinedExpT2Coeffs = 0.;
    }
  }
  if (expansionCoeffGradFlag) {
    if (combinedExpT1CoeffGrads.numRows() != num_pts ||
	combinedExpT1CoeffGrads.numCols() != num_t1v)
      combinedExpT1CoeffGrads.shape(num_t1v, num_pts);
    else combinedExpT1CoeffGrads = 0.;
  }

  switch (combine_type) {
  case MULT_COMBINE: { // multiplication of current and stored expansions
    size_t l, ll, num_lev = expansionType1Coeffs.size();
    RealVector t1c_vals(num_lev, false);  Real t1c_prod;
    for (p=0; p<num_pts; ++p) {
      RealVector c_vars(Teuchos::View,
			const_cast<Real*>(comb_var_sets[p]), (int)num_v);
      // pre-compute stored vals once for both Coeffs/CoeffGrads
      for (ec1_it =expansionType1Coeffs.begin(), l=0;
	   ec1_it!=expansionType1Coeffs.end(); ++ec1_it, ++l)
	t1c_vals[l] = (ec1_it == expT1CoeffsIter) ? value(c_vars) :
	  stored_value(c_vars, ec1_it->first);

      if (expansionCoeffFlag) {
	// split up type1/type2 contribs so increments are performed properly
	Real& combined_t1c_p = combinedExpT1Coeffs[p];
	combined_t1c_p = t1c_vals[0];
	for (l=1; l<num_lev; ++l)
	  combined_t1c_p *= t1c_vals[l];
	if (use_derivs) {
	  Real* combined_t2c_p = combinedExpT2Coeffs[p];
	  // hf = curr*stored --> dhf/dx = dcurr/dx*stored + curr*dstored/dx
	  for (ec2_it =expansionType2Coeffs.begin(), l=0;
	       ec2_it!=expansionType2Coeffs.end(); ++ec2_it, ++l) {
	    const RealVector& basis_grad = (ec2_it == expT2CoeffsIter) ?
	      gradient_basis_variables(c_vars) :
	      stored_gradient_basis_variables(c_vars, ec2_it->first);
	    t1c_prod = 1.;
	    for (ll=0; ll<num_lev; ++ll)
	      if (ll != l)
		t1c_prod *= t1c_vals[l];
	    for (v=0; v<num_t2v; ++v)
	      combined_t2c_p[v] += basis_grad[v] * t1c_prod;
	  }
	}
      }
      if (expansionCoeffGradFlag) {
	Real* combined_t1g_p = combinedExpT1CoeffGrads[p];
	// hf = curr*stored --> dhf/dx = dcurr/dx*stored + curr*dstored/dx
	for (eg1_it =expansionType1CoeffGrads.begin(), l=0;
	     eg1_it!=expansionType1CoeffGrads.end(); ++eg1_it, ++l) {
	  const RealVector& nonbasis_grad = (eg1_it == expT1CoeffGradsIter) ?
	    gradient_nonbasis_variables(c_vars) :
	    stored_gradient_nonbasis_variables(c_vars, ec2_it->first);
	  t1c_prod = 1.;
	  for (ll=0; ll<num_lev; ++ll)
	    if (ll != l)
	      t1c_prod *= t1c_vals[l];
	  for (v=0; v<num_t1v; ++v)
	    combined_t1g_p[v] += nonbasis_grad[v] * t1c_prod;
	}
      }
    }
    break;
  }
  default: //case ADD_COMBINE: // addition of current and stored expansions
    for (p=0; p<num_pts; ++p) {
      RealVector c_vars(Teuchos::View,
			const_cast<Real*>(comb_var_sets[p]), (int)num_v);
      if (expansionCoeffFlag) {
	for (ec1_it  = expansionType1Coeffs.begin();
	     ec1_it != expansionType1Coeffs.end(); ++ec1_it)
	  combinedExpT1Coeffs[p] += (ec1_it == expT1CoeffsIter) ?
	    value(c_vars) : stored_value(c_vars, ec1_it->first);
	if (use_derivs) {
	  Real* combined_t2c_p = combinedExpT2Coeffs[p];
	  for (ec2_it  = expansionType2Coeffs.begin();
	       ec2_it != expansionType2Coeffs.end(); ++ec2_it) {
	    const RealVector& basis_grad = (ec2_it == expT2CoeffsIter) ?
	      gradient_basis_variables(c_vars) :
	      stored_gradient_basis_variables(c_vars, ec2_it->first);
	    for (v=0; v<num_t2v; ++v)
	      combined_t2c_p[v] += basis_grad[v];
	  }
	}
      }
      if (expansionCoeffGradFlag) {
	Real* combined_t1g_p = combinedExpT1CoeffGrads[p];
	for (eg1_it  = expansionType1CoeffGrads.begin();
	     eg1_it != expansionType1CoeffGrads.end(); ++eg1_it) {
	  const RealVector& nonbasis_grad = (eg1_it == expT1CoeffGradsIter) ?
	    gradient_nonbasis_variables(c_vars) :
	    stored_gradient_nonbasis_variables(c_vars, ec2_it->first);
	  for (v=0; v<num_t1v; ++v)
	    combined_t1g_p[v] += nonbasis_grad[v];
	}
      }
    }
    break;
  }

#ifdef DEBUG
  PCout << "Combined type1 expansion coefficients:\n" << combinedExpT1Coeffs;
  if (use_derivs) {
    PCout << "Combined type2 expansion coefficients:\n";
    write_data(PCout, combinedExpT2Coeffs, false, true, true);
  }
  if (expansionCoeffGradFlag) {
    PCout << "Combined type1 expansion coefficient gradients:\n";
    write_data(PCout, combinedExpT1CoeffGrads, false, true, true);
  }
#endif // DEBUG
}


void NodalInterpPolyApproximation::combined_to_active(bool clear_combined)
{
  // replace active expansions with combined expansion arrays
  // > clear_inactive() takes care of the auxilliary inactive expansions that
  //   are now assimilated within each new active expansion
  // > swap() is conditionally available for Real{Vector,Matrix}

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  if (expansionCoeffFlag) {
    if (clear_combined) {
      expT1CoeffsIter->second.swap(combinedExpT1Coeffs); // shallow ptr swap
      combinedExpT1Coeffs.resize(0);
    }
    else // redundant copy
      expT1CoeffsIter->second = combinedExpT1Coeffs;     // deep copy

    if (data_rep->basisConfigOptions.useDerivs) {
      if (clear_combined) {
	expT2CoeffsIter->second.swap(combinedExpT2Coeffs); // shallow ptr swap
	combinedExpT2Coeffs.reshape(0, 0);
      }
      else // redundant copy
	expT2CoeffsIter->second = combinedExpT2Coeffs;     // deep copy
    }
  }

  if (expansionCoeffGradFlag) {
    if (clear_combined) {
      expT1CoeffGradsIter->second.swap(combinedExpT1CoeffGrads); // ptr swap
      combinedExpT1CoeffGrads.reshape(0, 0);
    }
    else // redundant copy
      expT1CoeffGradsIter->second = combinedExpT1CoeffGrads;     // deep copy
  }

  // resize sobolIndices to sync with resize of sobolIndexMap
  allocate_component_sobol();

  // Create a dummy modSurrData for the combined-now-active coeffs, for
  // accelerating FINAL_RESULTS (integration, VBD processing, etc.)
  synthetic_surrogate_data(modSurrData); // overwrite data for activeKey

  // If outgoing stats type is active (e.g., as in Dakota::NonDExpansion::
  // multifidelity_expansion()), then previous active stats are invalidated.
  // But if outgoing stats type is combined, then can avoid recomputation
  // and carry over current moment stats from combined to active. 
  // Note: due to this carry-over optimization, updating of stats type from
  //       COMBINED to ACTIVE must follow this function
  if (data_rep->expConfigOptions.refineStatsType == ACTIVE_EXPANSION_STATS)
    clear_computed_bits();
}


void NodalInterpPolyApproximation::
synthetic_surrogate_data(SurrogateData& surr_data)
{
  // Update the active key of surr_data with synthetic data based on the
  // active grid from csg_driver

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;

  // CombinedSparseGridDriver::combined_to_active() transfers all data except
  // collocation indices, which are invalidated by the combination.  In support
  // of the synthetic data to be created, new colloc indices are defined at the
  // end of CSGDriver::combined_to_active(), and ordering is preserved below. 
  const RealMatrix&  var_sets = data_rep->driver()->variable_sets();
  const RealVector& t1_coeffs = expT1CoeffsIter->second;
  const RealMatrix& t2_coeffs = expT2CoeffsIter->second;

  // shallow copy vars array data using CombinedSparseGridDriver::variableSets
  surr_data.clear_all_active();
  size_t num_v = var_sets.numRows(), num_pts = var_sets.numCols();
  bool use_derivs = data_rep->basisConfigOptions.useDerivs;
  short bits = (use_derivs) ? 3 : 1;
  surr_data.resize(num_pts, bits, num_v);
  
  // use interpolant to produce data values that are being interpolated
  // Note: SurrogateData reconstruction follows order of unique variable sets,
  //       also used in defining expansion coeffs (inverse of mapping used in
  //       compute_coefficients())
  SDVArray& sdv_array = surr_data.variables_data(); 
  SDRArray& sdr_array = surr_data.response_data();
  for (size_t pt=0; pt<num_pts; ++pt) {
    RealVector c_vars(Teuchos::View, const_cast<Real*>(var_sets[pt]),
		      (int)num_v);
    sdv_array[pt].continuous_variables(c_vars);
    sdr_array[pt].response_function( value(c_vars, t1_coeffs, t2_coeffs) );
    if (use_derivs)
      sdr_array[pt].response_gradient(
	gradient_basis_variables(c_vars, t1_coeffs, t2_coeffs) );
  }
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
Real NodalInterpPolyApproximation::
tensor_product_mean(const RealVector& x, const RealVector& exp_t1_coeffs,
		    const RealMatrix& exp_t2_coeffs,
		    const UShortArray& lev_index,
		    const UShort2DArray& colloc_key,
		    const SizetArray& colloc_index)
{
  size_t i, j, num_colloc_pts = colloc_key.size(),
    num_v = sharedDataRep->numVars;
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
      const UShortArray& key_i = colloc_key[i]; key_i0 = key_i[0];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      if (li_0 == 0)          // single integration/interpolation weight = 1
	accumulator[0]  = exp_t1_coeffs[c_index];
      else if (rand_0)        // integration
	accumulator[0] += exp_t1_coeffs[c_index] * t1_wts_0[key_i0];
      else if (ei0 == _NPOS)  // interpolation without exact match
	accumulator[0] += exp_t1_coeffs[c_index] *  bc_vf_0[key_i0];
      else if (ei0 == key_i0) // interpolation with exact match
	accumulator[0]  = exp_t1_coeffs[c_index];
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
      const UShortArray& key_i = colloc_key[i]; key_i0 = key_i[0];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      const Real* t2_coeffs_i = exp_t2_coeffs[c_index];
      if (rand_0) { // integration
	if (li_0 == 0) {
	  t1_accumulator[0]  = exp_t1_coeffs[c_index]; // t1 weight is 1
	  //t2_accum_0[0]    = 0.;                     // t2 weight is 0
	  for (k=1; k<num_v; ++k)
	    t2_accum_0[k]    = t2_coeffs_i[k];         // t1 weight is 1
	}
	else {
	  t1_accumulator[0] += exp_t1_coeffs[c_index] * t1_wts_0[key_i0];
	  t2_accum_0[0]     += t2_coeffs_i[0] * t2_wts_0[key_i0];
	  for (k=1; k<num_v; ++k)
	    t2_accum_0[k]   += t2_coeffs_i[k] * t1_wts_0[key_i0];
	}
      }
      else {        // interpolation
	if (li_0 == 0) {
	  t1_accumulator[0]  = exp_t1_coeffs[c_index];  // t1 value is 1
	  t2_accum_0[0]      = t2_coeffs_i[0] * poly_0.type2_value(x0, key_i0);
	  for (k=1; k<num_v; ++k)
	    t2_accum_0[k]    = t2_coeffs_i[k];          // t1 value is 1
	}
	else {
	  t1_val = poly_0.type1_value(x0, key_i0);
          t1_accumulator[0] += exp_t1_coeffs[c_index] * t1_val;
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
      const UShortArray& key_i = colloc_key[i];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      tp_mean += exp_t1_coeffs[c_index]
	      *  type1_interpolant_value(x, key_i, lev_index,
	                                 data_rep->nonRandomIndices)
	      *  type1_weight(key_i, lev_index, data_rep->randomIndices);
      const Real *t2_coeff_i = exp_t2_coeffs[c_index];
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
      const UShortArray& key_i = colloc_key[i]; key_i0 = key_i[0];
      c_index = (colloc_index.empty()) ? i : colloc_index[i];
      if (li_0 == 0)   // single integration/interpolation weight = 1
	accumulator[0]  = exp_t1_coeffs[c_index];
      else if (rand_0) // integration
	accumulator[0] += exp_t1_coeffs[c_index] * t1_wts_0[key_i0];
      else             // interpolation
	accumulator[0] += exp_t1_coeffs[c_index]
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
	tp_mean += exp_t1_coeffs[i]
	  * type1_interpolant_value(x, colloc_key[i], lev_index,
	                            data_rep->nonRandomIndices)
	  * type1_weight(colloc_key[i], lev_index, data_rep->randomIndices);
      else
	for (i=0; i<num_colloc_pts; ++i)
	  tp_mean += exp_t1_coeffs[colloc_index[i]]
	    * type1_interpolant_value(x, colloc_key[i], lev_index,
	                              data_rep->nonRandomIndices)
	    * type1_weight(colloc_key[i], lev_index, data_rep->randomIndices);
    return tp_mean;
    */
  }
}


/** Overloaded all_variables version supporting Smolyak sparse grids. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_mean_gradient(const RealVector& x,
			     const RealVector& exp_t1_coeffs,
			     const RealMatrix& exp_t2_coeffs,
			     const RealMatrix& exp_t1_coeff_grads,
			     const UShortArray& lev_index,
			     const UShort2DArray& colloc_key,
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
    num_deriv_vars = dvv.size(), num_colloc_pts = colloc_key.size(),
    num_v = sharedDataRep->numVars;
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  if (tpMomentGrads.size() != 2) tpMomentGrads.resize(2);
  RealVector& tp_mean_grad = tpMomentGrads[0];
  if (tp_mean_grad.length() != num_deriv_vars)
    tp_mean_grad.sizeUninitialized(num_deriv_vars);
  tp_mean_grad = 0.;

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
    Real *accum_0 = accumulator[0], t1_coeff, t1_wt_00, bc_vf_00, prod;
    const Real *t1_coeff_grad;
    bool rand_0 = data_rep->randomVarsKey[0];
    for (p=0; p<num_colloc_pts; ++p) {
      const UShortArray& key_p = colloc_key[p]; key_p0 = key_p[0];
      if (augment)
	t1_coeff = (colloc_index.empty()) ? exp_t1_coeffs[p] :
	  exp_t1_coeffs[colloc_index[p]];
      if (insert)
	t1_coeff_grad = (colloc_index.empty()) ? exp_t1_coeff_grads[p] :
	  exp_t1_coeff_grads[colloc_index[p]];
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
      tp_mean_grad[d] = accum[d] * scale;
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
    Real *t1_accum_0 = t1_accumulator[0], *t2_accum_0, t1_coeff, t1_wt_00,
      t1_val, t2_val, t1_grad, prod1, prod2, x0 = x[0];
    const Real *t2_coeffs_p;
    bool rand_0 = data_rep->randomVarsKey[0];
    for (p=0; p<num_colloc_pts; ++p) {
      const UShortArray& key_p = colloc_key[p]; key_p0 = key_p[0];
      c_index     = (colloc_index.empty()) ? p : colloc_index[p];
      t1_coeff    = exp_t1_coeffs[c_index];
      t2_coeffs_p = exp_t2_coeffs[c_index];
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
      tp_mean_grad[d] = t1_accum[d];
      t2_accum = t2_accumulators[d][num_v-1];
      for (v=0; v<num_v; ++v)
	tp_mean_grad[d] += t2_accum[v];
    }

    /*
    // Simpler but less efficient approach:
    //PCout << "tensor_product_mean_gradient(): Original useDerivs."<<std::endl;
    Real t1_coeff, *t2_coeff_p, t1_wt;
    for (p=0; p<num_colloc_pts; ++p) {
      const UShortArray& key_p = colloc_key[p];
      c_index = (colloc_index.empty()) ? p : colloc_index[p];
      t1_coeff   = exp_t1_coeffs[c_index];
      t2_coeff_p = exp_t2_coeffs[c_index];
      t1_wt = type1_weight(key_p, lev_index, data_rep->randomIndices);
      for (d=0; d<num_deriv_vars; ++d) {
	deriv_index = dvv[d] - 1; // OK since we are in an "All" view
	// -------------------------------------------------------------------
	// deriv of All var expansion w.r.t. nonrand var (design augmentation)
	// -------------------------------------------------------------------
	Real& tp_mean_grad_d = tp_mean_grad[d];
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
    Real *accum_0 = accumulator[0], t1_coeff, t1_wt_00, prod, t1_val, x0 = x[0];
    const Real *t1_coeff_grad;
    bool rand_0 = data_rep->randomVarsKey[0];
    for (p=0; p<num_colloc_pts; ++p) {
      const UShortArray& key_p = colloc_key[p]; key_p0 = key_p[0];
      if (augment)
	t1_coeff = (colloc_index.empty()) ? exp_t1_coeffs[p] :
	  exp_t1_coeffs[colloc_index[p]];
      if (insert)
	t1_coeff_grad = (colloc_index.empty()) ? exp_t1_coeff_grads[p] :
	  exp_t1_coeff_grads[colloc_index[p]];
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
    copy_data(accumulator[num_v-1], (int)num_deriv_vars, tp_mean_grad);

    /*
    // Simpler but less efficient approach:
    //PCout << "tensor_product_mean_gradient(): Original no derivs."<<std::endl;
    Real t1_coeff, *t1_coeff_grad, t1_wt, t1_val;
    for (p=0; p<num_colloc_pts; ++p) {
      const UShortArray& key_p = colloc_key[p];
      c_index = (colloc_index.empty()) ? p : colloc_index[p];
      t1_wt   = type1_weight(key_p, lev_index, data_rep->randomIndices);
      if (insert) {
	t1_coeff_grad = exp_t1_coeff_grads[c_index];
	t1_val = type1_interpolant_value(x, key_p, lev_index,
	                                 data_rep->nonRandomIndices);
      }
      else
	t1_coeff = exp_t1_coeffs[c_index];	
      for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	deriv_index = dvv[d] - 1; // OK since we are in an "All" view
	tp_mean_grad[d] += (data_rep->randomVarsKey[deriv_index]) ?
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

  return tp_mean_grad;
}


/** Covariance of response functions for a matched tensor product grid.
    Supports all_variables mode and either interpolation of products or
    product of interpolants formulations.  For the latter, recursive
    usage only provides the "diagonal" contributions from matched
    tensor products (PRODUCT_OF_INTERPOLANTS_FAST); an exact estimation
    (PRODUCT_OF_INTERPOLANTS_FULL) requires augmentation with mixed tensor
    product contributions using the overloaded form of this function. */
Real NodalInterpPolyApproximation::
tensor_product_covariance(const RealVector& x, Real mean_1, Real mean_2,
			  const RealVector& exp_t1c_1,
			  const RealMatrix& exp_t2c_1,
			  const RealVector& exp_t1c_2,
			  const RealMatrix& exp_t2c_2,
			  const UShortArray& lev_index,
			  const UShort2DArray& colloc_key,
			  const SizetArray& colloc_index)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  //bool same = (this == nip_approx_2);
  Real t1_coeff_1_mm1, t1_coeff_2_mm2;
  size_t i, c_index_i, num_colloc_pts = colloc_key.size(),
    num_v = sharedDataRep->numVars;

  switch (data_rep->momentInterpType) {
  case INTERPOLATION_OF_PRODUCTS: {
    const Real3DArray& t1_wts_1d = driver_rep->type1_collocation_weights_1d();
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
	t1_coeff_1_mm1 = exp_t1c_1[c_index_i] - mean_1;
	t1_coeff_2_mm2 = exp_t1c_2[c_index_i] - mean_2;
	const UShortArray& key_i = colloc_key[i]; key_i0 = key_i[0];
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
      const Real3DArray& t2_wts_1d = driver_rep->type2_collocation_weights_1d();
      unsigned short       li_0 = lev_index[0];
      BasisPolynomial&   poly_0 = data_rep->polynomialBasis[li_0][0];
      const RealArray& t1_wts_0 = t1_wts_1d[li_0][0]; // li==li for integrated
      const RealArray& t2_wts_0 = t2_wts_1d[li_0][0]; // li==li for integrated
      size_t k; bool rand_0 = data_rep->randomVarsKey[0];
      unsigned short key_i0, max0 = poly_0.interpolation_size() - 1;
      for (i=0; i<num_colloc_pts; ++i) {
	c_index_i = (colloc_index.empty()) ? i : colloc_index[i];
	t1_coeff_1_mm1 = exp_t1c_1[c_index_i] - mean_1;
	t1_coeff_2_mm2 = exp_t1c_2[c_index_i] - mean_2;
	const Real *t2_coeff_1 = exp_t2c_1[c_index_i],
	           *t2_coeff_2 = exp_t2c_2[c_index_i];
	const UShortArray& key_i = colloc_key[i]; key_i0 = key_i[0];
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
	const UShortArray& key_i = colloc_key[i];
	t1_coeff_1_mm1 = exp_t1c_1[c_index_i] - mean_1;
	t1_coeff_2_mm2 = exp_t1c_2[c_index_i] - mean_2;
	tp_covar += t1_coeff_1_mm1 * t1_coeff_2_mm2
	  * data_rep->type1_interpolant_value(x, key_i, lev_index,
				              data_rep->nonRandomIndices)
	  * data_rep->type1_weight(key_i, lev_index, data_rep->randomIndices);
	const Real *t2_coeff_1 = exp_t2c_1[c_index_i],
	           *t2_coeff_2 = exp_t2c_2[c_index_i];
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
	t1_coeff_1_mm1 = exp_t1c_1[c_index_i] - mean_1;
	t1_coeff_2_mm2 = exp_t1c_2[c_index_i] - mean_2;
	const UShortArray& key_i = colloc_key[i]; key_i0 = key_i[0];
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
	const UShortArray& key_i = colloc_key[i];
	c_index_i = (colloc_index.empty()) ? i : colloc_index[i];
	t1_coeff_1_mm1 = exp_t1c_1[c_index_i] - mean_1;
	t1_coeff_2_mm2 = exp_t1c_2[c_index_i] - mean_2;
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
	t1_coeff_1_mm1 = value(c_vars, exp_t1c_1, exp_t2c_1) - mean_1;
	t1_coeff_2_mm2 = //(same) ? t1_coeff_1_mm1 :
	  value(c_vars, exp_t1c_2, exp_t2c_2) - mean_2;
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
	t1_coeff_1_mm1 = value(c_vars, exp_t1c_1, exp_t2c_1) - mean_1;
	t1_coeff_2_mm2 = //(same) ? t1_coeff_1_mm1 :
	  value(c_vars, exp_t1c_2, exp_t2c_2) - mean_2;
	const RealVector& t2_coeff_1 =
	  gradient_basis_variables(c_vars, exp_t1c_1, exp_t2c_1);
	const RealVector& t2_coeff_2 = //(same) ? t2_coeff_1 :
	  gradient_basis_variables(c_vars, exp_t1c_2, exp_t2c_2);
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
	t1_coeff_1_mm1 = value(c_vars, exp_t1c_1, exp_t2c_1) - mean_1; // tensor_product_value() ?
	t1_coeff_2_mm2 = //(same) ? t1_coeff_1_mm1 :
	  value(c_vars, exp_t1c_2, exp_t2c_2) - mean_2;
	tp_covar += t1_coeff_1_mm1 * t1_coeff_2_mm2
	  * data_rep->type1_interpolant_value(x, key_i, reinterp_lev_index,
				              data_rep->nonRandomIndices)
	  * data_rep->type1_weight(key_i, reinterp_lev_index,
	                           data_rep->randomIndices);
	const RealVector& t2_coeff_1 =
	  gradient_basis_variables(c_vars, exp_t1c_1, exp_t2c_1);
	const RealVector& t2_coeff_2 = //(same) ? t2_coeff_1 :
	  gradient_basis_variables(c_vars, exp_t1c_2, exp_t2c_2);
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
	t1_coeff_1_mm1 = value(c_vars, exp_t1c_1, exp_t2c_1) - mean_1;
	t1_coeff_2_mm2 = //(same) ? t1_coeff_1_mm1 :
	  value(c_vars, exp_t1c_2, exp_t2c_2) - mean_2;
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
	t1_coeff_1_mm1 = value(c_vars, exp_t1c_1, exp_t2c_1) - mean_1; // tensor_product_value() ?
	t1_coeff_2_mm2 = //(same) ? t1_coeff_1_mm1 :
	  value(c_vars, exp_t1c_2, exp_t2c_2) - mean_2;
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
  case PRODUCT_OF_INTERPOLANTS_FAST: // Note: used for TPQ; SSG short-cuts
  case PRODUCT_OF_INTERPOLANTS_FULL:
    return product_of_interpolants(x, mean_1, mean_2, exp_t1c_1, exp_t2c_1,
				   exp_t1c_2, exp_t2c_2, lev_index, colloc_key,
				   colloc_index);
    break;
  }
}


/** Overloaded all_variables version supporting interpolation of
    products or product of interpolants approaches. */
const RealVector& NodalInterpPolyApproximation::
tensor_product_variance_gradient(const RealVector& x, Real mean,
				 const RealVector& mean_grad,
				 const RealVector& exp_t1_coeffs,
				 const RealMatrix& exp_t2_coeffs,
				 const RealMatrix& exp_t1_coeff_grads,
				 const UShortArray& lev_index,
				 const UShort2DArray& colloc_key,
				 const SizetArray& colloc_index,
				 const SizetArray& dvv)
{
  size_t d, insert_cntr, deriv_index, num_deriv_vars = dvv.size(),
    num_v = sharedDataRep->numVars;
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  if (tpMomentGrads.size() != 2) tpMomentGrads.resize(2);
  RealVector& tp_var_grad = tpMomentGrads[1];
  if (tp_var_grad.length() != num_deriv_vars)
    tp_var_grad.sizeUninitialized(num_deriv_vars);
  tp_var_grad = 0.;

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
                  num_colloc_pts = colloc_key.size();
    unsigned short         max_0 = poly_0.interpolation_size() - 1;
    bool                  rand_0 = data_rep->randomVarsKey[0];

    if (data_rep->barycentricFlag) {

      // For barycentric interpolation: track x != newPoint within 1D basis
      short order = (augment) ? 3 : 1;
      data_rep->set_new_point(x, lev_index, data_rep->nonRandomIndices, order);

      const RealVector& bc_vf_0 = poly_0.barycentric_value_factors();
      const RealVector& bc_gf_0 = poly_0.barycentric_gradient_factors();
      RealMatrix accumulator(num_deriv_vars, num_v); // init to 0.
      Real *accum_0 = accumulator[0], t1_coeff_p_mm, t1_wt_00, bc_vf_00,
	prod1, prod2, prod3;
      for (p=0; p<num_colloc_pts; ++p) {
	const UShortArray& key_p = colloc_key[p]; key_p0 = key_p[0];
	c_index_p = (colloc_index.empty()) ? p : colloc_index[p];
	t1_coeff_p_mm = exp_t1_coeffs[c_index_p] - mean;
	const Real* t1_coeff_grad_p = (insert) ?
	  exp_t1_coeff_grads[c_index_p] : NULL;
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
	tp_var_grad[d] = accum[d] * scale;
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
	const UShortArray& key_p = colloc_key[p]; key_p0 = key_p[0];
	c_index_p = (colloc_index.empty()) ? p : colloc_index[p];
	t1_coeff_p_mm = exp_t1_coeffs[c_index_p] - mean;
	// no DVV: need all grad components for R to expand type2 for (R-mu)^2
	const Real* t2_coeff_p = exp_t2_coeffs[c_index_p];
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
	tp_var_grad[d] = t1_accum[d];
	t2_accum = t2_accumulators[d][num_v-1];
	for (v=0; v<num_v; ++v)
	  tp_var_grad[d] += t2_accum[v];
      }

      /*
      Real t1_coeff_p_mm;
      for (p=0; p<num_colloc_pts; ++p) {
	const UShortArray& key_p = colloc_key[p];
	c_index_p = (colloc_index.empty()) ? p : colloc_index[p];
	t1_coeff_p_mm = exp_t1_coeffs[c_index_p] - mean;
	// no DVV: need all grad components for R to expand type2 for (R-mu)^2
	const Real* t2_coeff_p = exp_t2_coeffs[c_index_p];

	for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	  deriv_index = dvv[d] - 1; // OK since we are in an "All" view
	  // ---------------------------------------------------------------
	  // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	  // ---------------------------------------------------------------
	  Real& grad_d = tp_var_grad[d];
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
	const UShortArray& key_p = colloc_key[p]; key_p0 = key_p[0];
	c_index_p = (colloc_index.empty()) ? p : colloc_index[p];
	t1_coeff_p_mm = exp_t1_coeffs[c_index_p] - mean;
	const Real* t1_coeff_grad_p = (insert) ?
	  exp_t1_coeff_grads[c_index_p] : NULL;
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
      copy_data(accumulator[num_v-1], (int)num_deriv_vars, tp_var_grad);

      /*
      Real t1_coeff_p_mm;
      size_t p, v, num_colloc_pts = colloc_key.size();
      for (p=0; p<num_colloc_pts; ++p) {
	const UShortArray& key_p = colloc_key[p];
	c_index_p = (colloc_index.empty()) ? p : colloc_index[p];
	t1_coeff_p_mm = exp_t1_coeffs[c_index_p] - mean;
	const Real* t1_coeff_grad_p = (insert) ?
	  exp_t1_coeff_grads[c_index_p] : NULL;
	for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	  deriv_index = dvv[d] - 1; // OK since we are in an "All" view
	  tp_var_grad[d] += (data_rep->randomVarsKey[deriv_index]) ?
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
	t1_coeff_p_mm = value(c_vars, exp_t1_coeffs, exp_t2_coeffs) - mean;
	const RealVector& t1_coeff_grad_p = (insert) ?
	  gradient_nonbasis_variables(c_vars, exp_t1_coeff_grads) : empty_rv;
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
	tp_var_grad[d] = accum[d] * scale;
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
	t1_coeff_p_mm = value(c_vars, exp_t1_coeffs, exp_t2_coeffs) - mean;
	// no DVV: need all grad components for R to expand type2 for (R-mu)^2
	const RealVector& t2_coeff_p
	  = gradient_basis_variables(c_vars, exp_t1_coeffs, exp_t2_coeffs);
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
	tp_var_grad[d] = t1_accum[d];
	t2_accum = t2_accumulators[d][num_v-1];
	for (v=0; v<num_v; ++v)
	  tp_var_grad[d] += t2_accum[v];
      }

      /*
      Real t1_coeff_p_mm;
      for (p=0; p<reinterp_colloc_pts; ++p) {
	const UShortArray& key_p = reinterp_colloc_key[p];
	RealVector c_vars(Teuchos::View,
	  const_cast<Real*>(reinterp_var_sets[p]), num_v);
	t1_coeff_p_mm = value(c_vars, exp_t1_coeffs, exp_t2_coeffs) - mean;
	// no DVV: need all grad components for R to expand type2 for (R-mu)^2
	const RealVector& t2_coeff_p
	  = gradient_basis_variables(c_vars, exp_t1_coeffs, exp_t2_coeffs);

	for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	  deriv_index = dvv[d] - 1; // OK since we are in an "All" view
	  // ---------------------------------------------------------------
	  // deriv of All var exp w.r.t. nonrandom var (design augmentation)
	  // ---------------------------------------------------------------
	  Real& grad_d = tp_var_grad[d];
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
	t1_coeff_p_mm = value(c_vars, exp_t1_coeffs, exp_t2_coeffs) - mean;
	const RealVector& t1_coeff_grad_p = (insert) ?
	  gradient_nonbasis_variables(c_vars, exp_t1_coeff_grads) : empty_rv;
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
      copy_data(accumulator[num_v-1], (int)num_deriv_vars, tp_var_grad);

      /*
      RealVector empty_rv; Real t1_coeff_p_mm;
      for (p=0; p<reinterp_colloc_pts; ++p) {
	const UShortArray& key_p = reinterp_colloc_key[p];
	RealVector c_vars(Teuchos::View,
	  const_cast<Real*>(reinterp_var_sets[p]), num_v);
	t1_coeff_p_mm = value(c_vars, exp_t1_coeffs, exp_t2_coeffs) - mean;
	const RealVector& t1_coeff_grad_p = (insert) ?
	  gradient_nonbasis_variables(c_vars, exp_t1_coeff_grads) : empty_rv;
	for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
	  deriv_index = dvv[d] - 1; // OK since we are in an "All" view
	  tp_var_grad[d] += (data_rep->randomVarsKey[deriv_index]) ?
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
  case PRODUCT_OF_INTERPOLANTS_FULL: // for matched products; mismatched = TO DO
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
    size_t j, k, c_index_j, c_index_k, num_colloc_pts = colloc_key.size();
    Real wt_prod_j, Lsa_j, Lsa_k, dLsa_j_dsa_d, dLsa_k_dsa_d,
      t1_coeff_j_mm, t1_coeff_k_mm;
    for (d=0, insert_cntr=0; d<num_deriv_vars; ++d) {
      deriv_index = dvv[d] - 1; // OK since we are in an "All" view
      Real& grad_d = tp_var_grad[d];
      // first loop of double sum
      for (j=0; j<num_colloc_pts; ++j) {
	const UShortArray& key_j = colloc_key[j];
	c_index_j = (colloc_index.empty()) ? j : colloc_index[j];
	t1_coeff_j_mm = exp_t1_coeffs[c_index_j] - mean;
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
	  const UShortArray& key_k = colloc_key[k];
	  c_index_k = (colloc_index.empty()) ? k : colloc_index[k];
	  // to include jk-th term, colloc pts xi_j must be the same as xi_k
	  // for random var subset.  In this case, wt_prod_j may be reused.
	  if (data_rep->match_random_key(key_j, key_k)) {
	    t1_coeff_k_mm = exp_t1_coeffs[c_index_k] - mean;
	    Lsa_k =
	      data_rep->type1_interpolant_value(x, key_k, lev_index,
						data_rep->nonRandomIndices);
	    if (data_rep->randomVarsKey[deriv_index])
	      // ---------------------------------------------------------
	      // deriv of All var exp w.r.t. random var (design insertion)
	      // ---------------------------------------------------------
	      grad_d += wt_prod_j * Lsa_j * Lsa_k *
		( t1_coeff_j_mm * ( exp_t1_coeff_grads(insert_cntr, c_index_k) -
				    mean_grad[d] )
		+ t1_coeff_k_mm * ( exp_t1_coeff_grads(insert_cntr, c_index_j) -
				    mean_grad[d] ) );
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
  }

  return tp_var_grad;
}


Real NodalInterpPolyApproximation::
product_of_interpolants(const RealVector& x, Real mean_1, Real mean_2,
			const RealVector& exp_t1c_1,
			const RealMatrix& exp_t2c_1,
			const RealVector& exp_t1c_2,
			const RealMatrix& exp_t2c_2,
			const UShortArray& lev_index,
			const UShort2DArray& colloc_key,
			const SizetArray& colloc_index)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  size_t i, j, c_index_i, c_index_j, num_colloc_pts = colloc_key.size(),
    num_v = sharedDataRep->numVars;
  Real tp_covar = 0., t1_wt_Ls_prod_i, t1_coeff_1_mm1, t1_coeff_2_mm2;
  for (i=0; i<num_colloc_pts; ++i) {
    const UShortArray& key_i = colloc_key[i];
    c_index_i = (colloc_index.empty()) ? i : colloc_index[i];
    t1_coeff_1_mm1  = exp_t1c_1[c_index_i] - mean_1;
    t1_wt_Ls_prod_i
      = data_rep->type1_weight(key_i, lev_index, data_rep->randomIndices)
      * data_rep->type1_interpolant_value(x, key_i, lev_index,
					    data_rep->nonRandomIndices);
    for (j=0; j<num_colloc_pts; ++j) {
      const UShortArray& key_j = colloc_key[j];
      // to include the ij-th term,  basis i must be the same as basis j for
      // the random var subset.  In this case, wt_prod_i may be reused.  Note
      // that it is not necessary to collapse terms with the same random basis
      // subset, since cross term in (a+b)(a+b) = a^2+2ab+b^2 gets included.
      // If terms were collapsed (following eval of non-random portions), the
      // nested loop could be replaced with a single loop to evaluate (a+b)^2.
      if (data_rep->match_random_key(key_i, key_j)) {
	c_index_j = (colloc_index.empty()) ? j : colloc_index[j];
	t1_coeff_2_mm2 = exp_t1c_2[c_index_j] - mean_2;
	tp_covar += t1_coeff_1_mm1 * t1_coeff_2_mm2 * t1_wt_Ls_prod_i *
	  data_rep->type1_interpolant_value(x, key_j, lev_index,
					    data_rep->nonRandomIndices);
	/* TO DO
	if (data_rep->basisConfigOptions.useDerivs) {
	  const Real *t2_coeff_1i = exp_t2c_1[c_index_i],
	             *t2_coeff_2i = exp_t2c_2[c_index_i];
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


/** Covariance of response functions for differing tensor products
    using a PRODUCT_OF_INTERPOLANTS_FULL approach.  Needed for sparse
    interpolants in all_variables mode. */
Real NodalInterpPolyApproximation::
product_of_interpolants(const RealVector& x, Real mean_1, Real mean_2,
			const RealVector& exp_t1c_1,
			const RealMatrix& exp_t2c_1,
			const RealVector& exp_t1c_2,
			const RealMatrix& exp_t2c_2,
			const UShortArray& lev_index_1,
			const UShort2DArray& colloc_key_1,
			const SizetArray& colloc_index_1,
			const UShortArray& lev_index_2,
			const UShort2DArray& colloc_key_2,
			const SizetArray& colloc_index_2)
{
#ifdef DEBUG
  PCout << "product_of_interpolants for lev_index_1 =\n" << lev_index_1
	<< "lev_index_2 =\n" << lev_index_2 << "colloc_key_1 =\n"
	<< colloc_key_1 << "colloc_key_2 =\n" << colloc_key_2 << std::endl;
#endif // DEBUG

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  if (data_rep->momentInterpType != PRODUCT_OF_INTERPOLANTS_FULL) {
    PCerr << "Error: mixed tensor product covariance only required for full "
	  << "products of interpolants. " << std::endl;
    abort_handler(-1);
  }

  size_t i, j, k, c_index_1i, c_index_2j, num_pts_1 = colloc_key_1.size(),
    num_pts_2 = colloc_key_2.size();
  Real tp_covar = 0., t1_coeff_1i_mm1, t1_coeff_2j_mm2, t1_Ls_1i,
    basis_prod, basis_prod_k;
  bool non_zero; unsigned short l1k, l2k; SizetList::iterator it;
  for (i=0; i<num_pts_1; ++i) {
    const UShortArray& key_1i = colloc_key_1[i];
    c_index_1i = (colloc_index_1.empty()) ? i : colloc_index_1[i];
    t1_coeff_1i_mm1 = exp_t1c_1[c_index_1i] - mean_1;
    t1_Ls_1i = data_rep->type1_interpolant_value(x, key_1i, lev_index_1,
						 data_rep->nonRandomIndices);
    for (j=0; j<num_pts_2; ++j) {
      const UShortArray& key_2j = colloc_key_2[j];
      // to include ij-th term, expectation of the product of the interpolation
      // polynomials must be nonzero
      if (data_rep->basis_product(lev_index_1, key_1i,
				  lev_index_2, key_2j, basis_prod)) {
	c_index_2j = (colloc_index_2.empty()) ? j : colloc_index_2[j];
	t1_coeff_2j_mm2 = exp_t1c_2[c_index_2j] - mean_2;
	tp_covar += basis_prod * t1_coeff_1i_mm1 * t1_coeff_2j_mm2 * t1_Ls_1i *
	  data_rep->type1_interpolant_value(x, key_2j, lev_index_2,
					    data_rep->nonRandomIndices);
	/* TO DO
	if (data_rep->basisConfigOptions.useDerivs) {
	  const Real *t2_coeff_1i = exp_t2c_1[c_index_i],
	             *t2_coeff_2i = exp_t2c_2[c_index_i];
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


Real NodalInterpPolyApproximation::
value(const RealVector& x, const RealVector& exp_t1_coeffs,
      const RealMatrix& exp_t2_coeffs)
{
  // sum expansion to get response prediction
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    return value(x, exp_t1_coeffs, exp_t2_coeffs, tpq_driver->level_index(),
		 tpq_driver->collocation_key());
    break;
  }
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    return value(x, exp_t1_coeffs, exp_t2_coeffs,
		 csg_driver->smolyak_multi_index(),
		 csg_driver->smolyak_coefficients(),
		 csg_driver->collocation_key(),
		 csg_driver->collocation_indices());
    break;
  }
  }
}


Real NodalInterpPolyApproximation::
stored_value(const RealVector& x, const UShortArray& key)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not available in "
	  << "NodalInterpPolyApproximation::stored_value()" << std::endl;
    abort_handler(-1);
  }

  // sum expansion to get response prediction
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    return value(x, expansionType1Coeffs[key], expansionType2Coeffs[key],
		 tpq_driver->level_index(key),
		 tpq_driver->collocation_key(key));
    break;
  }
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    return value(x, expansionType1Coeffs[key], expansionType2Coeffs[key],
		 csg_driver->smolyak_multi_index(key),
		 csg_driver->smolyak_coefficients(key),
		 csg_driver->collocation_key(key),
		 csg_driver->collocation_indices(key));
    break;
  }
  }
}


Real NodalInterpPolyApproximation::
value(const RealVector& x, const RealVector& exp_t1_coeffs,
      const RealMatrix& exp_t2_coeffs, const UShortArray& lev_index,
      const UShort2DArray& colloc_key)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> default indexing
  return data_rep->
    tensor_product_value(x, exp_t1_coeffs, exp_t2_coeffs, lev_index,
			 colloc_key, colloc_index);
}



Real NodalInterpPolyApproximation::
value(const RealVector& x, const RealVector& exp_t1_coeffs,
      const RealMatrix& exp_t2_coeffs, const UShort2DArray& sm_mi,
      const IntArray& sm_coeffs, const UShort3DArray& colloc_key,
      const Sizet2DArray& colloc_index)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  size_t i, num_smolyak_indices = sm_coeffs.size();
  Real approx_val = 0.;
  for (i=0; i<num_smolyak_indices; ++i)
    if (sm_coeffs[i])
      approx_val += sm_coeffs[i] * data_rep->
	tensor_product_value(x, exp_t1_coeffs, exp_t2_coeffs, sm_mi[i],
			     colloc_key[i], colloc_index[i]);
  return approx_val;
}


/** Special case used for sparse grid interpolation on variable sub-sets
    defined from partial integration. */
Real NodalInterpPolyApproximation::
value(const RealVector& x, const RealVectorArray& t1_coeffs,
      const RealMatrixArray& t2_coeffs, const UShort3DArray& colloc_key,
      const SizetList& subset_indices)
{
  // Smolyak recursion of anisotropic tensor products
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  CombinedSparseGridDriver* csg_driver
    = (CombinedSparseGridDriver*)data_rep->driver();
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
gradient_basis_variables(const RealVector& x, const RealVector& exp_t1_coeffs,
			 const RealMatrix& exp_t2_coeffs)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    return gradient_basis_variables(x, exp_t1_coeffs, exp_t2_coeffs,
				    tpq_driver->level_index(),
				    tpq_driver->collocation_key());
    break;
  }
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    return gradient_basis_variables(x, exp_t1_coeffs, exp_t2_coeffs,
				    csg_driver->smolyak_multi_index(),
				    csg_driver->smolyak_coefficients(),
				    csg_driver->collocation_key(),
				    csg_driver->collocation_indices());
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const RealVector& exp_t1_coeffs,
			 const RealMatrix& exp_t2_coeffs, const SizetArray& dvv)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    return gradient_basis_variables(x, exp_t1_coeffs, exp_t2_coeffs,
				    tpq_driver->level_index(),
				    tpq_driver->collocation_key(), dvv);
    break;
  }
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    return gradient_basis_variables(x, exp_t1_coeffs, exp_t2_coeffs,
				    csg_driver->smolyak_multi_index(),
				    csg_driver->smolyak_coefficients(),
				    csg_driver->collocation_key(),
				    csg_driver->collocation_indices(), dvv);
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
stored_gradient_basis_variables(const RealVector& x, const UShortArray& key)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::stored_gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    return gradient_basis_variables(x, expansionType1Coeffs[key],
				    expansionType2Coeffs[key],
				    tpq_driver->level_index(key),
				    tpq_driver->collocation_key(key));
    break;
  }
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    return gradient_basis_variables(x, expansionType1Coeffs[key],
				    expansionType2Coeffs[key],
				    csg_driver->smolyak_multi_index(key),
				    csg_driver->smolyak_coefficients(key),
				    csg_driver->collocation_key(key),
				    csg_driver->collocation_indices(key));
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
stored_gradient_basis_variables(const RealVector& x, const SizetArray& dvv,
				const UShortArray& key)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::stored_gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    return gradient_basis_variables(x, expansionType1Coeffs[key],
				    expansionType2Coeffs[key],
				    tpq_driver->level_index(key),
				    tpq_driver->collocation_key(key), dvv);
    break;
  }
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    return gradient_basis_variables(x, expansionType1Coeffs[key],
				    expansionType2Coeffs[key],
				    csg_driver->smolyak_multi_index(key),
				    csg_driver->smolyak_coefficients(key),
				    csg_driver->collocation_key(key),
				    csg_driver->collocation_indices(key), dvv);
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const RealVector& exp_t1_coeffs,
			 const RealMatrix& exp_t2_coeffs,
			 const UShortArray& lev_index,
			 const UShort2DArray& colloc_key)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> default indexing
  return data_rep->
    tensor_product_gradient_basis_variables(x, exp_t1_coeffs, exp_t2_coeffs,
					    lev_index, colloc_key,colloc_index);
}


const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const RealVector& exp_t1_coeffs,
			 const RealMatrix& exp_t2_coeffs,
			 const UShortArray& lev_index,
			 const UShort2DArray& colloc_key, const SizetArray& dvv)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  SizetArray colloc_index; // empty -> default indexing
  return data_rep->
    tensor_product_gradient_basis_variables(x, exp_t1_coeffs, exp_t2_coeffs,
					    lev_index, colloc_key,
					    colloc_index, dvv);
}


const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const RealVector& exp_t1_coeffs,
			 const RealMatrix& exp_t2_coeffs,
			 const UShort2DArray& sm_mi, const IntArray& sm_coeffs,
			 const UShort3DArray& colloc_key,
			 const Sizet2DArray& colloc_index)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  size_t num_v = data_rep->numVars;
  if (approxGradient.length() != num_v)
    approxGradient.sizeUninitialized(num_v);
  approxGradient = 0.;
  // Smolyak recursion of anisotropic tensor products
  size_t i, j, num_smolyak_indices = sm_coeffs.size();
  for (i=0; i<num_smolyak_indices; ++i) {
    int coeff_i = sm_coeffs[i];
    if (coeff_i) {
      const RealVector& tp_grad = data_rep->
	tensor_product_gradient_basis_variables(x, exp_t1_coeffs, exp_t2_coeffs,
						sm_mi[i], colloc_key[i],
						colloc_index[i]);
      for (j=0; j<num_v; ++j)
	approxGradient[j] += coeff_i * tp_grad[j];
    }
  }
  return approxGradient;
}


const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const RealVector& exp_t1_coeffs,
			 const RealMatrix& exp_t2_coeffs,
			 const UShort2DArray& sm_mi, const IntArray& sm_coeffs,
			 const UShort3DArray& colloc_key,
			 const Sizet2DArray& colloc_index,
			 const SizetArray& dvv)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  size_t num_deriv_vars = dvv.size();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.sizeUninitialized(num_deriv_vars);
  approxGradient = 0.;
  // Smolyak recursion of anisotropic tensor products
  size_t i, j, num_smolyak_indices = sm_coeffs.size();
  for (i=0; i<num_smolyak_indices; ++i) {
    int coeff_i = sm_coeffs[i];
    if (coeff_i) {
      const RealVector& tp_grad = data_rep->
	tensor_product_gradient_basis_variables(x, exp_t1_coeffs, exp_t2_coeffs,
						sm_mi[i], colloc_key[i],
						colloc_index[i], dvv);
      for (j=0; j<num_deriv_vars; ++j)
	approxGradient[j] += coeff_i * tp_grad[j];
    }
  }
  return approxGradient;
}


/** Special case used for sparse grid interpolation on variable sub-sets
    defined from partial integration. */
const RealVector& NodalInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const RealVectorArray& t1_coeffs,
			 const RealMatrixArray& t2_coeffs,
			 const UShort3DArray& colloc_key,
			 const SizetList& subset_indices)
{
  // this could define a default_dvv and call gradient_basis_variables(x, dvv),
  // but we want this fn to be as fast as possible

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  size_t num_v = sharedDataRep->numVars;
  if (approxGradient.length() != num_v)
    approxGradient.sizeUninitialized(num_v);
  approxGradient = 0.;
  // Smolyak recursion of anisotropic tensor products
  CombinedSparseGridDriver* csg_driver
    = (CombinedSparseGridDriver*)data_rep->driver();
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


const RealVector& NodalInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x,
			    const RealMatrix& exp_t1_coeff_grads)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    return gradient_nonbasis_variables(x, exp_t1_coeff_grads,
				       tpq_driver->level_index(),
				       tpq_driver->collocation_key());
    break;
  }
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    return gradient_nonbasis_variables(x, exp_t1_coeff_grads,
				       csg_driver->smolyak_multi_index(),
				       csg_driver->smolyak_coefficients(),
				       csg_driver->collocation_key(),
				       csg_driver->collocation_indices());
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
stored_gradient_nonbasis_variables(const RealVector& x, const UShortArray& key)
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
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    return gradient_nonbasis_variables(x, expansionType1CoeffGrads[key],
				       tpq_driver->level_index(key),
				       tpq_driver->collocation_key(key));
    break;
  }
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    return gradient_nonbasis_variables(x, expansionType1CoeffGrads[key],
				       csg_driver->smolyak_multi_index(key),
				       csg_driver->smolyak_coefficients(key),
				       csg_driver->collocation_key(key),
				       csg_driver->collocation_indices(key));
    break;
  }
  }
}


const RealVector& NodalInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x,
			    const RealMatrix& exp_t1_coeff_grads,
			    const UShortArray& lev_index,
			    const UShort2DArray& colloc_key)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
  SizetArray colloc_index; // empty -> default indexing
  return data_rep->
    tensor_product_gradient_nonbasis_variables(x, exp_t1_coeff_grads, lev_index,
					       colloc_key, colloc_index);
}


const RealVector& NodalInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x,
			    const RealMatrix& exp_t1_coeff_grads,
			    const UShort2DArray& sm_mi,
			    const IntArray& sm_coeffs,
			    const UShort3DArray& colloc_key,
			    const Sizet2DArray& colloc_index)
{
  size_t num_deriv_vars = exp_t1_coeff_grads.numRows();
  if (approxGradient.length() != num_deriv_vars)
    approxGradient.sizeUninitialized(num_deriv_vars);
  approxGradient = 0.;
  // Smolyak recursion of anisotropic tensor products
  size_t i, j, num_smolyak_indices = sm_coeffs.size();
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  for (i=0; i<num_smolyak_indices; ++i) {
    int coeff_i = sm_coeffs[i];
    if (coeff_i) {
      const RealVector& tp_grad = data_rep->
	tensor_product_gradient_nonbasis_variables(x, exp_t1_coeff_grads,
						   sm_mi[i], colloc_key[i],
						   colloc_index[i]);
      for (j=0; j<num_deriv_vars; ++j)
	approxGradient[j] += coeff_i * tp_grad[j];
    }
  }
  return approxGradient;
}


const RealSymMatrix& NodalInterpPolyApproximation::
hessian_basis_variables(const RealVector& x)
{
  PCerr << "Error: NodalInterpPolyApproximation::hessian_basis_variables() "
	<< "not yet implemented." << std::endl;
  abort_handler(-1);

  return approxHessian;
}


const RealSymMatrix& NodalInterpPolyApproximation::
stored_hessian_basis_variables(const RealVector& x, const UShortArray& key)
{
  PCerr << "Error: NodalInterpPolyApproximation::stored_hessian_basis_"
	<< "variables() not yet implemented." << std::endl;
  abort_handler(-1);

  return approxHessian;
}


/** In this case, a subset of the expansion variables are random
    variables and the mean of the expansion involves integration over
    this subset and evaluation over the subset's complement.  For the
    linear sums of tensor interpolants within a sparse interpolant,
    the expectation can be taken inside the sum and we can simply add
    up the tensor mean contributions. */
Real NodalInterpPolyApproximation::
mean(const RealVector& x, const RealVector& exp_t1_coeffs, 
     const RealMatrix& exp_t2_coeffs)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  Real mean;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    SizetArray colloc_index; // empty -> default indexing
    mean = tensor_product_mean(x, exp_t1_coeffs, exp_t2_coeffs,
			       tpq_driver->level_index(),
			       tpq_driver->collocation_key(), colloc_index);
    break;
  }
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    const UShort2DArray& sm_mi        = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs    = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key   = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index = csg_driver->collocation_indices();
    size_t i, num_smolyak_indices = sm_coeffs.size();
    mean = 0.;
    for (i=0; i<num_smolyak_indices; ++i)
      if (sm_coeffs[i])
	mean += sm_coeffs[i] *
	  tensor_product_mean(x, exp_t1_coeffs, exp_t2_coeffs, sm_mi[i],
			      colloc_key[i], colloc_index[i]);
    break;
  }
  }
  return mean;
}


/** In this function, all expansion variables are random variables and
    any design/state variables are omitted from the expansion.  In
    this case, the derivative of the expectation is the expectation of
    the derivative.  The mixed derivative case (some design variables
    are inserted and some are augmented) requires no special treatment. */
const RealVector& NodalInterpPolyApproximation::
mean_gradient(const RealMatrix& exp_t1_coeff_grads,
	      const RealVector& t1_wts)
{
  // d/ds <R> = <dR/ds>

  size_t i, j, num_colloc_pts = t1_wts.length(),
    num_deriv_vars = exp_t1_coeff_grads.numRows();
  RealVector& mean_grad = momentGradsIter->second[0];
  if (mean_grad.length() == num_deriv_vars) mean_grad = 0.;
  else mean_grad.size(num_deriv_vars);
  Real t1_wt_i;
  for (i=0; i<num_colloc_pts; ++i) {
    t1_wt_i = t1_wts[i];
    for (j=0; j<num_deriv_vars; ++j)
      mean_grad[j] += exp_t1_coeff_grads(j,i) * t1_wt_i;
  }
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
    are inserted (derivatives are obtained from expansionType1CoeffGrads).  
    For the linear sums of tensor interpolants within a sparse 
    interpolant, the expectation can be taken inside the sum and we 
    can simply add up the tensor mean gradient contributions. */
const RealVector& NodalInterpPolyApproximation::
mean_gradient(const RealVector& x, const RealVector& exp_t1_coeffs,
	      const RealMatrix& exp_t2_coeffs,
	      const RealMatrix& exp_t1_coeff_grads, const SizetArray& dvv)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    SizetArray colloc_index; // empty -> default indexing
    return tensor_product_mean_gradient(x, exp_t1_coeffs, exp_t2_coeffs,
      exp_t1_coeff_grads, tpq_driver->level_index(),
      tpq_driver->collocation_key(), colloc_index, dvv);
    break;
  }
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID:
    size_t num_deriv_vars = dvv.size();
    RealVector& mean_grad = momentGradsIter->second[0];
    if (mean_grad.length() != num_deriv_vars)
      mean_grad.sizeUninitialized(num_deriv_vars);
    mean_grad = 0.;
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    const UShort2DArray& sm_mi        = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs    = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key   = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size();
    for (i=0; i<num_smolyak_indices; ++i) {
      int coeff = sm_coeffs[i];
      if (coeff) {
	const RealVector& tpm_grad =
	  tensor_product_mean_gradient(x, exp_t1_coeffs, exp_t2_coeffs,
	    exp_t1_coeff_grads, sm_mi[i], colloc_key[i], colloc_index[i], dvv);
	for (j=0; j<num_deriv_vars; ++j)
	  mean_grad[j] += coeff * tpm_grad[j];
      }
    }
    return mean_grad;
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
covariance(Real mean_1, Real mean_2,    const RealVector& exp_t1c_1,
	   const RealMatrix& exp_t2c_1, const RealVector& exp_t1c_2,
	   const RealMatrix& exp_t2c_2, const RealVector& t1_wts,
	   const RealMatrix& t2_wts)
{
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  // compute mean_1,mean_2 first, then compute covariance as
  // wt_prod*(coeff1-mean_1)*(coeff2-mean_2) in order to avoid precision
  // loss from computing covariance as <R_i R_j> - \mu_i \mu_j
  // Note: compute_statistics() in dakota/src/NonDExpansion.C orders calls
  //       to reduce repetition in moment calculations.
  Real covar = 0.;
  size_t i, j, num_colloc_pts = t1_wts.length(), num_v = sharedDataRep->numVars;
  if (data_rep->basisConfigOptions.useDerivs) {
    for (i=0; i<num_colloc_pts; ++i) {
      // type1 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      Real t1c_1i_mm1 = exp_t1c_1[i] - mean_1,
	   t1c_2i_mm2 = exp_t1c_2[i] - mean_2;
      covar += t1c_1i_mm1 * t1c_2i_mm2 * t1_wts[i];
      // type2 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      // --> interpolated gradients are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
      const Real *t2c_1i  = exp_t2c_1[i], *t2c_2i = exp_t2c_2[i],
	         *t2_wt_i = t2_wts[i];
      for (j=0; j<num_v; ++j)
	covar += (t1c_1i_mm1 * t2c_2i[j] + t1c_2i_mm2 * t2c_1i[j]) * t2_wt_i[j];
    }
  }
  else
    for (i=0; i<num_colloc_pts; ++i)
      covar += (exp_t1c_1[i] - mean_1) * (exp_t1c_2[i] - mean_2) * t1_wts[i];

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
covariance(const RealVector& x, Real mean_1, Real mean_2,
	   const RealVector& exp_t1c_1, const RealMatrix& exp_t2c_1,
	   const RealVector& exp_t1c_2, const RealMatrix& exp_t2c_2)
{
  Real covar;
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    SizetArray colloc_index; // empty -> default indexing

    switch (data_rep->momentInterpType) {
    case REINTERPOLATION_OF_PRODUCTS: {
      // compute higher-order covariance grid, if not already computed
      const UShortArray& lev_index = tpq_driver->level_index();
      reinterpolated_level(lev_index);
      // compute TP covariance with reinterpolation on the higher-order grid
      covar = tensor_product_covariance(x, mean_1, mean_2,
	exp_t1c_1, exp_t2c_1, exp_t1c_2, exp_t2c_2, lev_index,
	tpq_driver->collocation_key(), colloc_index);
      break;
    }
    case PRODUCT_OF_INTERPOLANTS_FAST: // don't recompute TPQ mean (unlike SSG)
    case INTERPOLATION_OF_PRODUCTS: case PRODUCT_OF_INTERPOLANTS_FULL:
      covar = tensor_product_covariance(x, mean_1, mean_2, exp_t1c_1,
	exp_t2c_1, exp_t1c_2, exp_t2c_2, tpq_driver->level_index(),
	tpq_driver->collocation_key(), colloc_index);
      break;
    }
    break;
  }

  // While we can collapse the Smolyak recursion and combine the weights in the
  // distinct variables case, we cannot do this here for the all_variables case
  // since the non-integrated interpolation polynomial portions are not constant
  // and are coupled with the weight combination.
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    const UShort2DArray& sm_mi        = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs    = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key   = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices  = sm_coeffs.size();
    covar = 0.;
    switch (data_rep->momentInterpType) {

    // Reinterpolation of covariance function in non-integrated dimensions.
    case REINTERPOLATION_OF_PRODUCTS:
      for (i=0; i<num_smolyak_indices; ++i)
	if (sm_coeffs[i]) {
	  reinterpolated_level(sm_mi[i]);
	  covar += sm_coeffs[i] * tensor_product_covariance(x, mean_1, mean_2,
	    exp_t1c_1, exp_t2c_1, exp_t1c_2, exp_t2c_2, sm_mi[i], colloc_key[i],
	    colloc_index[i]);
	}
      break;

    // For an interpolation of products (which captures product cross-terms),
    // a sum of tensor-product covariances is correct and straightforward.
    // Note that tensor_product_covariance() is overloaded to specialize its
    // logic based on data_rep->momentInterpType.
    case INTERPOLATION_OF_PRODUCTS:
      for (i=0; i<num_smolyak_indices; ++i)
	if (sm_coeffs[i])
	  covar += sm_coeffs[i] * tensor_product_covariance(x, mean_1, mean_2,
	    exp_t1c_1, exp_t2c_1, exp_t1c_2, exp_t2c_2, sm_mi[i], colloc_key[i],
	    colloc_index[i]);
      break;

    // For a fast product of interpolants, cross-terms are neglected and
    // results are approximate.  In addition, the mean reference point for
    // each covariance contribution is the tensor mean, not the total mean.
    case PRODUCT_OF_INTERPOLANTS_FAST:
      for (i=0; i<num_smolyak_indices; ++i)
	if (sm_coeffs[i]) {
	  const UShortArray&  sm_mi_i = sm_mi[i];
	  const UShort2DArray&  key_i = colloc_key[i];
	  const SizetArray& c_index_i = colloc_index[i];
	  // Replace incoming expansion means with tensor means
	  mean_1 = tensor_product_mean(x, exp_t1c_1, exp_t2c_1, sm_mi_i,
				       key_i, c_index_i);
	  mean_2 = //(same) ? mean_1 :
	    tensor_product_mean(x, exp_t1c_2, exp_t2c_2, sm_mi_i,
				key_i, c_index_i);
	  covar += sm_coeffs[i] *
	    // 10/2018: if not coeff^2, then more like INTERP_OF_PRODUCTS !!
	    product_of_interpolants(x, mean_1, mean_2, exp_t1c_1, exp_t2c_1,
	      exp_t1c_2, exp_t2c_2, sm_mi_i, key_i, c_index_i);
	    // short-cut tensor_product_covariance(...)
	}
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
	  const UShortArray&  sm_mi_i = sm_mi[i];
	  const UShort2DArray&  key_i = colloc_key[i];
	  const SizetArray& c_index_i = colloc_index[i];
#ifdef DEBUG
	  Real tp_covar = product_of_interpolants(x, mean_1, mean_2, exp_t1c_1,
	    exp_t2c_1, exp_t1c_2, exp_t2c_2, sm_mi_i, key_i, c_index_i);
	  covar += sm_coeff_i * sm_coeff_i * tp_covar;
	  PCout << "Diagonal covar: sm_coeffs[" << i << "] = " << sm_coeff_i
		<< " tp_covar = " << tp_covar << '\n';
#else
	  covar += sm_coeff_i * sm_coeff_i * // coeff^2 for prod of interp
	    product_of_interpolants(x, mean_1, mean_2, exp_t1c_1, exp_t2c_1,
	      exp_t1c_2, exp_t2c_2, sm_mi_i, key_i, c_index_i);
	    // short-cut tensor_product_covariance(...)
#endif // DEBUG
	  for (j=0; j<i; ++j)
	    if (sm_coeffs[j]) { // two of each off-diagonal term
#ifdef DEBUG
	      Real tp_covar = product_of_interpolants(x, mean_1, mean_2,
		exp_t1c_1, exp_t2c_1, exp_t1c_2, exp_t2c_2, sm_mi_i, key_i,
		c_index_i, sm_mi[j], colloc_key[j], colloc_index[j]);
	      covar += 2. * sm_coeff_i * sm_coeffs[j] * tp_covar;
	      PCout << "Off-diagonal covar: sm_coeffs[" << i << "] = "
		    << sm_coeff_i << " sm_coeffs[" << j << "] = "
		    << sm_coeffs[j] << " tp_covar = " << tp_covar << '\n';
#else
	      covar += 2. * sm_coeff_i * sm_coeffs[j] *
		product_of_interpolants(x, mean_1, mean_2, exp_t1c_1,
		  exp_t2c_1, exp_t1c_2, exp_t2c_2, sm_mi_i, key_i, c_index_i,
		  sm_mi[j], colloc_key[j], colloc_index[j]);
#endif // DEBUG
	    }
	}
      }
      break;
    }
    break;
  }
  }

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
const RealVector& NodalInterpPolyApproximation::
variance_gradient(Real mean, const RealVector& exp_t1_coeffs,
		  const RealMatrix& exp_t1_coeff_grads,
		  const RealVector& t1_wts)
{
  size_t i, j, num_colloc_pts = t1_wts.length(),
    num_deriv_vars = exp_t1_coeff_grads.numRows();
  RealVector& var_grad = momentGradsIter->second[1];
  if  (var_grad.length() == num_deriv_vars) var_grad = 0.;
  else var_grad.size(num_deriv_vars);

  // See Eq. 6.23 in Theory Manual: grad of variance incorporates grad of mean
  for (i=0; i<num_colloc_pts; ++i) {
    Real term_i = 2. * (exp_t1_coeffs[i] - mean) * t1_wts[i];
    for (j=0; j<num_deriv_vars; ++j)
      var_grad[j] += term_i * exp_t1_coeff_grads(j,i);
  }
  return var_grad;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  This function must handle the mixed case, where
    some design/state variables are augmented (and are part of the
    expansion) and some are inserted (derivatives are obtained from
    expansionType1CoeffGrads). */
const RealVector& NodalInterpPolyApproximation::
variance_gradient(const RealVector& x, Real mean, const RealVector& mean_grad,
		  const RealVector& exp_t1_coeffs,
		  const RealMatrix& exp_t2_coeffs,
		  const RealMatrix& exp_t1_coeff_grads, const SizetArray& dvv)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    SizetArray colloc_index; // empty -> default indexing
    switch (data_rep->momentInterpType) {
    case REINTERPOLATION_OF_PRODUCTS: {
      const UShortArray& lev_index = tpq_driver->level_index();
      // compute ord/lev for higher-order covariance grid over non-random vars
      reinterpolated_level(lev_index);
      // compute variance grad with reinterpolation on the higher-order grid
      return tensor_product_variance_gradient(x, mean, mean_grad,
	exp_t1_coeffs, exp_t2_coeffs, exp_t1_coeff_grads, lev_index,
	tpq_driver->collocation_key(), colloc_index, dvv);
      break;
    }
    case INTERPOLATION_OF_PRODUCTS:
    case PRODUCT_OF_INTERPOLANTS_FAST: case PRODUCT_OF_INTERPOLANTS_FULL:
      return tensor_product_variance_gradient(x, mean, mean_grad, exp_t1_coeffs,
	exp_t2_coeffs, exp_t1_coeff_grads, tpq_driver->level_index(),
	tpq_driver->collocation_key(), colloc_index, dvv);
      break;
    }
    break;
  }
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    size_t num_deriv_vars = dvv.size();
    RealVector& var_grad = momentGradsIter->second[1];
    if  (var_grad.length() == num_deriv_vars) var_grad = 0.;
    else var_grad.size(num_deriv_vars);
    // Smolyak recursion of anisotropic tensor products
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    const UShort2DArray& sm_mi        = csg_driver->smolyak_multi_index();
    const IntArray&      sm_coeffs    = csg_driver->smolyak_coefficients();
    const UShort3DArray& colloc_key   = csg_driver->collocation_key();
    const Sizet2DArray&  colloc_index = csg_driver->collocation_indices();
    size_t i, j, num_smolyak_indices = sm_coeffs.size(); int coeff;
    switch (data_rep->momentInterpType) {

    // Reinterpolation of covariance function in non-integrated dimensions.
    case REINTERPOLATION_OF_PRODUCTS:
      for (i=0; i<num_smolyak_indices; ++i)
	if (coeff = sm_coeffs[i]) {
	  reinterpolated_level(sm_mi[i]);
	  const RealVector& tpv_grad =
	    tensor_product_variance_gradient(x, mean, mean_grad, exp_t1_coeffs,
	    exp_t2_coeffs, exp_t1_coeff_grads, sm_mi[i], colloc_key[i],
	    colloc_index[i], dvv);
	  for (j=0; j<num_deriv_vars; ++j)
	    var_grad[j] += coeff * tpv_grad[j];
	}
      break;

    // For an interpolation of products (which captures product cross-terms),
    // a sum of tensor-product covariances is correct and straightforward.
    case INTERPOLATION_OF_PRODUCTS:
      for (i=0; i<num_smolyak_indices; ++i)
	if (coeff = sm_coeffs[i]) {
	  const RealVector& tpv_grad =
	    tensor_product_variance_gradient(x, mean, mean_grad, exp_t1_coeffs,
	    exp_t2_coeffs, exp_t1_coeff_grads, sm_mi[i], colloc_key[i],
	    colloc_index[i], dvv);
	  for (j=0; j<num_deriv_vars; ++j)
	    var_grad[j] += coeff * tpv_grad[j];
	}

    // For a fast product of interpolants, cross-terms are neglected and
    // mean/mean_grad are for each tensor product grid.
    case PRODUCT_OF_INTERPOLANTS_FAST:
      for (i=0; i<num_smolyak_indices; ++i)
	if (coeff = sm_coeffs[i]) {
	  const UShortArray&  sm_mi_i = sm_mi[i];
	  const UShort2DArray&  key_i = colloc_key[i];
	  const SizetArray& c_index_i = colloc_index[i];
	  // Replace incoming expansion mean/mean_grad with tensor mean_grad
	  Real tpm = tensor_product_mean(x, exp_t1_coeffs, exp_t2_coeffs,
	    sm_mi_i, key_i, c_index_i);
	  const RealVector& tpm_grad =
	    tensor_product_mean_gradient(x, exp_t1_coeffs, exp_t2_coeffs,
	    exp_t1_coeff_grads, sm_mi_i, key_i, c_index_i, dvv);
	  const RealVector& tpv_grad =
	    tensor_product_variance_gradient(x, tpm, tpm_grad, exp_t1_coeffs,
	    exp_t2_coeffs, exp_t1_coeff_grads, sm_mi_i, key_i, c_index_i, dvv);
	  for (j=0; j<num_deriv_vars; ++j)
	    var_grad[j] += coeff * tpv_grad[j];
	}
      break;

    case PRODUCT_OF_INTERPOLANTS_FULL:
      // Note:  matched variance gradient contributions as above
      // TO DO: variance gradient contributions from mis-matched tensor grids
      PCerr << "Error: variance gradient not yet implemented for "
	    << "PRODUCT_OF_INTERPOLANTS_FULL." << std::endl;
      abort_handler(-1);
      break;
    }
    return var_grad;
    break;
  }
  }
}


Real NodalInterpPolyApproximation::
expectation(const RealVector& t1_coeffs, const RealMatrix& t2_coeffs,
	    const RealVector& t1_wts,    const RealMatrix& t2_wts)
{
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  //IntegrationDriver* driver_rep = data_rep->driverRep;

  Real integral = 0.;
  size_t i, num_colloc_pts = t1_coeffs.length();
  if (data_rep->basisConfigOptions.useDerivs) {
    size_t j, num_v = t2_coeffs.numRows();
    for (i=0; i<num_colloc_pts; ++i) {
      integral += t1_coeffs[i] * t1_wts[i];
      const Real *t2c_i = t2_coeffs[i], *t2_wt_i = t2_wts[i];
      for (j=0; j<num_v; ++j)
	integral += t2c_i[j] * t2_wt_i[j];
    }
  }
  else
    for (i=0; i<num_colloc_pts; ++i)
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
integrate_response_moments(size_t num_moments, bool combined_stats)
{
  // In this case, we are constrained to use the original collocation points
  // (for which response values are available) within the numerical integration.
  // In the case where interpolatory rules (e.g., Clenshaw-Curtis) have been 
  // used, the integrand accuracy may suffer.  For this reason, the expansion
  // moments (which integrate the interpolant expansion instead of the response,
  // and therefore may employ alternate rules) are preferred.

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::integrate_response_moments()" << std::endl;
    abort_handler(-1);
  }

  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  IntegrationDriver* driver_rep = data_rep->driverRep;
  RealVector& numer_mom = numMomentsIter->second;
  if (numer_mom.length() != num_moments)
    numer_mom.sizeUninitialized(num_moments);

  // Support combined_stats for completeness
  // > use of combined_to_active() prior to full_stats computation makes
  //   this moot / unused for NIPA
  // > SharedNodalInterpPolyApproxData::pre_combine_data() activates the
  //   maximal grid and NodalInterpPolyApproximation::combine_coefficients()
  //   computes combinedExpT{1,2}Coeffs based on it, so
  //   driver_rep->type{1,2}_weight_sets() are synchronized/sufficient
  //   (no need for driver_rep->combined_type{1,2}_weight_sets()
  if (data_rep->basisConfigOptions.useDerivs) {
    if (combined_stats)
      integrate_moments(combinedExpT1Coeffs, combinedExpT2Coeffs,
	driver_rep->combined_type1_weight_sets(),
	driver_rep->combined_type2_weight_sets(), numer_mom);
    else
      integrate_moments(expT1CoeffsIter->second, expT2CoeffsIter->second,
	driver_rep->type1_weight_sets(), driver_rep->type2_weight_sets(),
	numer_mom);
  }
  else if (combined_stats)
    integrate_moments(combinedExpT1Coeffs,
      driver_rep->combined_type1_weight_sets(), numer_mom);
  else
    integrate_moments(expT1CoeffsIter->second, driver_rep->type1_weight_sets(),
      numer_mom);
}


void NodalInterpPolyApproximation::
integrate_expansion_moments(size_t num_moments, bool combined_stats)
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
    PCerr << "Error: expansion coefficients not defined in NodalInterpPoly"
	  << "Approximation::integrate_expansion_moments()"<< std::endl;
    abort_handler(-1);
  }
  if (combined_stats) {
    // requires use of modular value() and gradient_basis_variables() for
    // combinedExpT{1,2}Coeffs, which are available, but don't bother for now
    PCerr << "Error: combined_stats unavailable.  NodalInterpPolyApproximation"
	  << "::integrate_expansion_moments()\n       currently requires "
	  << "promotion of combined to active." << std::endl;
    abort_handler(-1);
  }
  RealVector& exp_mom = expMomentsIter->second;
  if (exp_mom.length() != num_moments) exp_mom.sizeUninitialized(num_moments);

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
    integrate_moments(t1_exp, alt_driver->type1_weight_sets(), exp_mom);
  }
  /*
  // Native quadrature on interpolant can be value-based or gradient-enhanced.
  // These approaches just replace the values of the response with the values
  // of the interpolant (integrates powers of the interpolant using the same
  // collocation rules/orders used to form the interpolant).
  else if (csg_driver->track_collocation_details() &&
	   csg_driver->collocation_indices().empty()) {// invalidated by combine
    IntegrationDriver* driver_rep = data_rep->driverRep;
    const RealMatrix&  var_sets = driver_rep->variable_sets();
    size_t i, num_pts = var_sets.numCols(), num_v = var_sets.numRows();
    RealVector t1_exp(num_pts);
    if (data_rep->basisConfigOptions.useDerivs) { // gradient-enhanced native
      RealMatrix t2_exp(num_v, num_pts);
      for (i=0; i<num_pts; ++i) {
	RealVector c_vars(Teuchos::View,
			  const_cast<Real*>(var_sets[i]), (int)num_v);
	t1_exp[i] = value(c_vars); // *** requires colloc_indices! ***
	Teuchos::setCol(gradient_basis_variables(c_vars), (int)i, t2_exp);
      }
      integrate_moments(t1_exp, t2_exp, driver_rep->type1_weight_sets(),
			driver_rep->type2_weight_sets(), exp_mom);
    }
    else { // value-based native quadrature
      for (i=0; i<num_pts; ++i) {
	RealVector c_vars(Teuchos::View,
			  const_cast<Real*>(var_sets[i]), (int)num_v);
	t1_exp[i] = value(c_vars); // *** requires colloc_indices! ***
      }
      integrate_moments(t1_exp, data_rep->driverRep->type1_weight_sets(),
			exp_mom);
    }
  }
  */
  else { // use modSurrData: original single-level or synthetic combined
    IntegrationDriver* driver_rep = data_rep->driverRep;
    const SDRArray& sdr_array = modSurrData.response_data();
    size_t i, num_pts = sdr_array.size();
    RealVector t1_exp(num_pts);
    if (data_rep->basisConfigOptions.useDerivs) { // gradient-enhanced native
      RealMatrix t2_exp(data_rep->numVars, num_pts);
      for (i=0; i<num_pts; ++i) {
	t1_exp[i] = sdr_array[i].response_function();
	Teuchos::setCol(sdr_array[i].response_gradient(), (int)i, t2_exp);
      }
      integrate_moments(t1_exp, t2_exp, driver_rep->type1_weight_sets(),
			driver_rep->type2_weight_sets(), exp_mom);
    }
    else { // value-based native quadrature
      for (i=0; i<num_pts; ++i)
	t1_exp[i] = sdr_array[i].response_function();
      integrate_moments(t1_exp, data_rep->driverRep->type1_weight_sets(),
			exp_mom);
    }
  }
#ifdef DEBUG
  PCout << "Expansion moments type 1 coefficients:\n" << t1_exp
	<< "Expansion moments:\n" << exp_mom;
#endif // DEBUG
}


/** Computes the variance of component functions. Assumes that all
    subsets of set_value have been computed in advance which will be
    true so long as the partial_variance is called following
    appropriate enumeration of set value  */
void NodalInterpPolyApproximation::
compute_partial_variance(const BitArray& set_value)
{
  // Perform inner integral over complementary set u' to form new weighted
  // coefficients h; then perform outer integral of h^2 over set u
  SharedNodalInterpPolyApproxData* data_rep
    = (SharedNodalInterpPolyApproxData*)sharedDataRep;
  const BitArrayULongMap& sobol_index_map = data_rep->sobolIndexMap;

  unsigned long pv_index;
  BitArrayULongMap::const_iterator cit = sobol_index_map.find(set_value);
  if (cit == sobol_index_map.end()) {
    PCerr << "Error in compute_partial_variance(): key not found in "
	  << "sobolIndexMap." << std::endl;
    abort_handler(-1);
  }
  else pv_index = cit->second;

  Real& variance = partialVariance[pv_index];
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
  size_t i, j, num_member_coeffs, num_v = sharedDataRep->numVars;
  SizetList member_indices;
  for (j=0; j<num_v; ++j)
    if (member_bits[j])
      member_indices.push_back(j);

  Real integral = 0.;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)data_rep->driver();
    const RealMatrix&      var_sets = tpq_driver->variable_sets();
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
      RealVector c_vars(Teuchos::View,
	const_cast<Real*>(var_sets[member_colloc_index[i]]), (int)num_v);
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
  case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver
      = (CombinedSparseGridDriver*)data_rep->driver();
    const RealMatrix&       var_sets = csg_driver->variable_sets();
    const IntArray&        sm_coeffs = csg_driver->smolyak_coefficients();
    const UShort2DArray&    sm_index = csg_driver->smolyak_multi_index();
    const UShort3DArray&  colloc_key = csg_driver->collocation_key();
    const Sizet2DArray& colloc_index = csg_driver->collocation_indices();
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
	  RealVector c_vars(Teuchos::View,
	    const_cast<Real*>(var_sets[m_c_index_i[j]]), (int)num_v);
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

  const RealVector& exp_t1_coeffs = expT1CoeffsIter->second;
  const RealMatrix& exp_t2_coeffs = expT2CoeffsIter->second;

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
      += nonmember_wt * exp_t1_coeffs[c_index];
    // reduced dimension data updated more times than necessary, but tracking
    // this redundancy would be more expensive/complicated.  Note: non-member
    // key and c_vars data may change, but member data should be consistent.
    member_t1_wts[member_index]       = member_wt;
    member_colloc_key[member_index]   = key_i;  // links back to interp poly's
    member_colloc_index[member_index] = c_index;// links back to Driver varSets

    // now do the same for the type2 coeffs and weights
    if (data_rep->basisConfigOptions.useDerivs) {
      Real *m_t2_coeffs_i = member_t2_coeffs[member_index],
	   *m_t2_wts_i    = member_t2_wts[member_index];
      const Real *t2_coeffs_i = exp_t2_coeffs[c_index];
      for (j=0; j<num_v; ++j) {
	data_rep->type2_weight(j, key_i, lev_index, member_bits,
			       member_wt, nonmember_wt);
	m_t2_coeffs_i[j] += nonmember_wt * t2_coeffs_i[j];
	m_t2_wts_i[j]    =  member_wt;
      }
    }
  }
#ifdef VBD_DEBUG
  PCout << "member_bits: "         << member_bits // MSB->LSB: order reversed
	<< "\nmember_t1_coeffs:\n" << member_t1_coeffs
	<< "member_t1_wts:\n"      << member_t1_wts;
  if (data_rep->basisConfigOptions.useDerivs) {
    PCout << "member_t2_coeffs:\n";
    write_data(PCout, member_t2_coeffs, false, true, true);
    PCout << "member_t2_wts:\n";
    write_data(PCout, member_t2_wts,    false, true, true);
  }
#endif // VBD_DEBUG
}

} // namespace Pecos
