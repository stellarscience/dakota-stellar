/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        RegressOrthogPolyApproximation
//- Description:  Implementation code for RegressOrthogPolyApproximation class
//-               
//- Owner:        John Jakeman

#include "RegressOrthogPolyApproximation.hpp"
#include "pecos_global_defs.hpp"
#include "pecos_math_util.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

// headers necessary for cross validation
#include "math_tools.hpp"
#include "CrossValidation.hpp"

//#define DEBUG

namespace Pecos {


int RegressOrthogPolyApproximation::min_coefficients() const
{
  // return the minimum number of data instances required to build the 
  // surface of multiple dimensions
  if (expansionCoeffFlag || expansionCoeffGradFlag)
    // Now that L1-regression has been implemented. There is no longer a need 
    // to enforce a lower bound on the number of data instances.
    return 1;
    // multiIndex is computed by allocate_arrays() in compute_coefficients(),
    // which is too late for use of this fn by ApproximationInterface::
    // minimum_samples() in DataFitSurrModel::build_global(), so number of
    // expansion terms must be calculated.
    //return total_order_terms(approxOrder);
  else
    return 0;
}


/** In this case, regression is used in place of spectral projection.  That
    is, instead of calculating the PCE coefficients using inner products, 
    linear least squares is used to estimate the PCE coefficients which
    best match a set of response samples.  The least squares estimation is
    performed using DGELSS (SVD) or DGGLSE (equality-constrained) from
    LAPACK, based on anchor point and derivative data availability. */
void RegressOrthogPolyApproximation::select_solver(bool cv_active)
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  bool fn_constrained_lls
    = (data_rep->basisConfigOptions.useDerivs && faultInfo.constr_eqns &&
       faultInfo.constr_eqns < data_rep->expansion_terms()); // candidate exp
  bool eq_con
    = (fn_constrained_lls || faultInfo.anchor_fn || faultInfo.anchor_grad);

  // **************************
  // Algorithm selection logic:
  // **************************

  // Assign defaults:
  short ec_options_solver = data_rep->expConfigOptions.expCoeffsSolnApproach;
  if (ec_options_solver == DEFAULT_REGRESSION) {
    if (faultInfo.under_determined)
      CSOpts.solver = LASSO_REGRESSION;
    else
      CSOpts.solver = (eq_con && !cv_active) ? EQ_CON_LEAST_SQ_REGRESSION :
	SVD_LEAST_SQ_REGRESSION;
  }
  else if (ec_options_solver == DEFAULT_LEAST_SQ_REGRESSION)
    CSOpts.solver = (eq_con && !faultInfo.under_determined && !cv_active) ?
      EQ_CON_LEAST_SQ_REGRESSION : SVD_LEAST_SQ_REGRESSION;
  else // Assign user's selection (not a default)
    CSOpts.solver = ec_options_solver;

  // Special case error traps (Note: only an error for user request, as CV
  // solver may select EQ_CON_LEAST_SQ_REGRESSION on the fly...)
  if ( cv_active && ( ec_options_solver == EQ_CON_LEAST_SQ_REGRESSION ||
		      ec_options_solver == ORTHOG_LEAST_INTERPOLATION ) )
    throw( std::runtime_error("Cross validation does not currently support equality constrained least squares or orthogonal least interpolation") );

  // Manage overrides/settings for BP, BPDN, equality-constrained least sq:
  switch (CSOpts.solver) {
  case BASIS_PURSUIT: case BASIS_PURSUIT_DENOISING:
    if ( !faultInfo.under_determined ) {  // no BP/BPDN for over-determined
      PCerr << "Warning: Could not perform BP/BPDN for over-determined system."
	    << "\n         Using SVD least squares regression instead.\n";
      CSOpts.solver = SVD_LEAST_SQ_REGRESSION; // force override to least sq
    }
    break;
  case EQ_CON_LEAST_SQ_REGRESSION:
    if ( faultInfo.under_determined ) {
      PCerr << "Warning: Could not perform equality constrained least-squares."
	    << "\n         Using LASSO regression instead.\n";
      CSOpts.solver = LASSO_REGRESSION;
    }
    else if ( !eq_con ) {
      PCerr << "Warning: Could not perform equality constrained least-squares."
	    << "\n         Using SVD least squares regression instead.\n";
      CSOpts.solver = SVD_LEAST_SQ_REGRESSION;
    }
    break;
  }

  // Set solver parameters
  RealVector noise_tols = data_rep->regressConfigOptions.noiseTols; // copy
  if ( CSOpts.solver == EQ_CON_LEAST_SQ_REGRESSION )
    CSOpts.numFunctionSamples = modSurrData.points();
  if ( CSOpts.solver == LASSO_REGRESSION )
    CSOpts.delta = data_rep->regressConfigOptions.l2Penalty;
  if ( noise_tols.empty() ) {
    noise_tols.size( 1 );
    noise_tols[0] = (CSOpts.solver == BASIS_PURSUIT_DENOISING) ? 1.e-3 :
      CSOpts.epsilon;
  }
  else
    CSOpts.epsilon = noise_tols[0];
  CSOpts.solverTolerance = (CSOpts.solver == SVD_LEAST_SQ_REGRESSION) ? -1. :
    data_rep->expConfigOptions.convergenceTol;
  CSOpts.verbosity = std::max(0, data_rep->expConfigOptions.outputLevel - 1);
  if ( data_rep->expConfigOptions.maxSolverIterations > 0 )
    CSOpts.maxNumIterations = data_rep->expConfigOptions.maxSolverIterations;

  // define global flags
  // Previous default assignments + method overrides simplifies this logic
  switch (CSOpts.solver) {
  case ORTHOG_LEAST_INTERPOLATION: case ORTHOG_MATCH_PURSUIT:
  case LASSO_REGRESSION:           case LEAST_ANGLE_REGRESSION:
  case BASIS_PURSUIT:              case BASIS_PURSUIT_DENOISING:
    sparseSoln =  true; break;
  //case BASIS_PURSUIT: case BASIS_PURSUIT_DENOISING:
  //  sparseSoln = faultInfo.under_determined; break; // handled by override
  default:
    sparseSoln = false; break;
  }
}


void RegressOrthogPolyApproximation::allocate_arrays()
{
  if (sparseSoln) {
    SharedRegressOrthogPolyApproxData* data_rep
      = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
    update_active_iterators(data_rep->activeKey);
    allocate_total_sobol(); // no dependencies

    // defer allocations until sparsity is known
    if (data_rep->expConfigOptions.vbdFlag && 
	data_rep->expConfigOptions.vbdOrderLimit == 1)
      allocate_component_sobol(); // no dependence on multiIndex for order 1
    // else defer until sparse recovery

    //size_expansion(); // defer until sparse recovery

    RealVector& exp_mom = expMomentsIter->second;
    if (exp_mom.length() != 2) exp_mom.sizeUninitialized(2);
  }
  else // pre-allocate total-order expansion
    OrthogPolyApproximation::allocate_arrays();
}


void RegressOrthogPolyApproximation::compute_coefficients()
{
  PolynomialApproximation::compute_coefficients();
  if (!expansionCoeffFlag && !expansionCoeffGradFlag)
    return;

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;

  // check data set for gradients/constraints/faults to determine settings
  modSurrData.data_checks();
#ifdef DEBUG
  data_rep->gradient_check();
#endif // DEBUG
  // multiIndex size (from Shared allocate_data()) used to set
  // faultInfo.under_determined.  Note: OLI's multiIndex is empty,
  // but does not use under_determined flag.
  set_fault_info();
  // initial solver settings (may be updated for CV folds)
  select_solver(data_rep->regressConfigOptions.crossValidation);
  // array allocations dependent on solver type
  allocate_arrays();

  switch (data_rep->expConfigOptions.expBasisType) {
  case DEFAULT_BASIS: // least interpolation case
  case TOTAL_ORDER_BASIS: case TENSOR_PRODUCT_BASIS:
    run_regression(); // solve for PCE coefficients, optionally with cross valid
    break;
  case ADAPTED_BASIS_GENERALIZED: case ADAPTED_BASIS_EXPANDING_FRONT:
    adapt_regression(); // adapt for the best multiIndex for a fixed data set
    break;
  }

  clear_computed_bits();
}


void RegressOrthogPolyApproximation::increment_coefficients()
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  update_active_iterators(data_rep->activeKey); // redundant but needed for
                                                // prevExp prior to compute

  // for use in pop_coefficients()
  prevExpCoeffs     = expCoeffsIter->second;     // copy
  prevExpCoeffGrads = expCoeffGradsIter->second; // copy
  prevSparseIndices = sparseIndIter->second;     // copy

  compute_coefficients(); // from scratch
}


void RegressOrthogPolyApproximation::pop_coefficients(bool save_data)
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const UShortArray& key = data_rep->activeKey;

  // likely overkill, but multilevel roll up after increment modifies and
  // then restores active key
  update_active_iterators(key);

  RealVector& exp_coeffs = expCoeffsIter->second;
  RealMatrix& exp_grads  = expCoeffGradsIter->second;
  SizetSet&   sparse_ind = sparseIndIter->second;

  // store the incremented coeff state for possible push'
  if (save_data) {
    poppedExpCoeffs[key].push_back(exp_coeffs);
    poppedExpCoeffGrads[key].push_back(exp_grads);
    poppedSparseInd[key].push_back(sparse_ind);
  }

  // reset expansion{Coeffs,CoeffGrads}
  exp_coeffs = prevExpCoeffs;  exp_grads = prevExpCoeffGrads;
  sparse_ind = prevSparseIndices;

  clear_computed_bits();
}


void RegressOrthogPolyApproximation::push_coefficients()
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const UShortArray& key = data_rep->activeKey;

  // synchronize expansionCoeff{s,Grads} and approxData
  update_active_iterators(key);

  // SharedPolyApproxData::candidate_index() currently returns 0 for
  // all cases other than generalized sparse grids
  size_t p_index = data_rep->push_index();

  // store current state for use in pop_coefficients()
  prevExpCoeffs     = expCoeffsIter->second;     // copy
  prevExpCoeffGrads = expCoeffGradsIter->second; // copy
  prevSparseIndices = sparseIndIter->second;     // copy

  // retrieve a previously popped state
  std::map<UShortArray, RealVectorDeque>::iterator prv_it
    = poppedExpCoeffs.find(key);
  std::map<UShortArray, RealMatrixDeque>::iterator prm_it
    = poppedExpCoeffGrads.find(key);
  std::map<UShortArray, SizetSetDeque>::iterator pss_it
    = poppedSparseInd.find(key);
  RealVectorDeque::iterator rv_it;  RealMatrixDeque::iterator rm_it;
  SizetSetDeque::iterator   ss_it;
  if (prv_it != poppedExpCoeffs.end()) {
    rv_it = prv_it->second.begin();     std::advance(rv_it, p_index);
    expCoeffsIter->second = *rv_it;     prv_it->second.erase(rv_it);
  }
  if (prm_it != poppedExpCoeffGrads.end()) {
    rm_it = prm_it->second.begin();     std::advance(rm_it, p_index);
    expCoeffGradsIter->second = *rm_it; prm_it->second.erase(rm_it);
  }
  if (pss_it != poppedSparseInd.end()) {
    ss_it = pss_it->second.begin();     std::advance(ss_it, p_index);
    sparseIndIter->second = *ss_it;     pss_it->second.erase(ss_it);
  }

  clear_computed_bits();
}


void RegressOrthogPolyApproximation::combine_coefficients()
{
  // Combine the data stored previously by store_coefficients()
  bool sparse = false;
  std::map<UShortArray, SizetSet>::iterator sp_it;
  if (!sparseIndices.empty())
    for (sp_it=sparseIndices.begin(); sp_it!=sparseIndices.end(); ++sp_it)
      if (!sp_it->second.empty())
	{ sparse = true; break; }

  if (!sparse)
    { OrthogPolyApproximation::combine_coefficients(); return; }

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  // Coefficient combination is not dependent on active state
  //update_active_iterators(data_rep->activeKey);

  // Note: computed bits are also cleared when refineStatsType is changed
  if (data_rep->expConfigOptions.refineStatsType == COMBINED_EXPANSION_STATS)
    clear_computed_bits();

  //allocate_component_sobol(); // size sobolIndices from shared sobolIndexMap

  // support mixed case using simple approach of populating the missing
  // sparse indices arrays (not optimal for performance but a lot less code),
  // prior to expansion aggregation.  Note: sparseSobolIndexMap is updated
  // following aggregation.
  const std::map<UShortArray, UShort2DArray>& mi = data_rep->multiIndex;
  std::map<UShortArray, UShort2DArray>::const_iterator mi_cit = mi.begin();
  for (sp_it =sparseIndices.begin();
       sp_it!=sparseIndices.end() && mi_cit!=mi.end(); ++sp_it, ++mi_cit) {
    SizetSet& sparse_ind = sp_it->second;
    if (sparse_ind.empty())
      inflate(sparse_ind, mi_cit->second.size());
  }

  std::map<UShortArray, RealVector>::iterator ec_it;
  std::map<UShortArray, RealMatrix>::iterator eg_it;
  switch (data_rep->expConfigOptions.combineType) {
  case MULT_COMBINE: {
    // perform the multiplication of level expansions
    const UShort3DArray& combined_mi_seq = data_rep->combinedMultiIndexSeq;
    size_t cntr = 0, num_seq = combined_mi_seq.size();
    /*
    ec_it = expansionCoeffs.begin();  eg_it = expansionCoeffGrads.begin();
    sp_it = sparseIndices.begin();   mi_cit = mi.begin();
    const UShort2DArray& multi_index_a = mi_cit->second; ++mi_cit;
    const RealVector&    exp_coeffs_a  =  ec_it->second; ++ec_it;
    const RealMatrix&    exp_grads_a   =  eg_it->second; ++eg_it;
    const SizetSet&      sparse_ind_a  =  sp_it->second; ++sp_it;
    const UShort2DArray& multi_index_c = (num_seq) ?
      combined_mi_seq[cntr] : data_rep->combinedMultiIndex;
    multiply_expansion(multi_index_a, sparse_ind_a, exp_coeffs_a, exp_grads_a,
		       mi_cit->second, sp_it->second, ec_it->second,
		       eg_it->second, multi_index_c, combinedSparseIndices,
		       combinedExpCoeffs, combinedExpCoeffGrads);
    for (cntr=1; cntr<=num_seq; ++cntr, ++sp_it, ++ec_it, ++eg_it, ++mi_cit) {
      const UShort2DArray& multi_index_c = (cntr < num_seq) ?
	combined_mi_seq[cntr] : data_rep->combinedMultiIndex;
      multiply_expansion(combined_mi_seq[cntr-1], combinedSparseIndices,
			 combinedExpCoeffs, combinedExpCoeffGrads,
			 mi_cit->second, sp_it->second, ec_it->second,
			 eg_it->second, multi_index_c, combinedSparseIndices,
			 combinedExpCoeffs, combinedExpCoeffGrads);
    }
    */
    ec_it = ++expansionCoeffs.begin();  eg_it = ++expansionCoeffGrads.begin();
    sp_it = ++sparseIndices.begin();   mi_cit = ++mi.begin();
    for (cntr=0; cntr<=num_seq; ++cntr, ++sp_it, ++ec_it, ++eg_it, ++mi_cit) {
      const UShort2DArray& multi_index_a = (cntr) ?
	combined_mi_seq[cntr-1] : mi.begin()->second;
      const RealVector&    exp_coeffs_a  = (cntr) ?
	combinedExpCoeffs       : expansionCoeffs.begin()->second;
      const RealMatrix&    exp_grads_a   = (cntr) ?
	combinedExpCoeffGrads   : expansionCoeffGrads.begin()->second;
      const SizetSet&      sparse_ind_a  = (cntr) ?
	combinedSparseIndices   : sparseIndices.begin()->second;
      const UShort2DArray& multi_index_c = (cntr < num_seq) ?
	combined_mi_seq[cntr]   : data_rep->combinedMultiIndex;
      multiply_expansion(multi_index_a, sparse_ind_a, exp_coeffs_a, exp_grads_a,
			 mi_cit->second, sp_it->second, ec_it->second,
			 eg_it->second, multi_index_c, combinedSparseIndices,
			 combinedExpCoeffs, combinedExpCoeffGrads);
    }
    break;
  }
  case ADD_MULT_COMBINE:
    //overlay_expansion(data_rep->storedMultiIndex[i], storedExpCoeffs[i],
    //                  storedExpCoeffGrads[i], addCoeffs, addCoeffGrads);
    //multiply_expansion(data_rep->storedMultiIndex[i], storedExpCoeffs[i],
    //                   storedExpCoeffGrads[i], multCoeffs, multCoeffGrads);
    //compute_combine_factors(addCoeffs, multCoeffs);
    //apply_combine_factors();
    PCerr << "Error : additive+multiplicative combination not yet "
	  << "implemented in OrthogPolyApproximation::combine_coefficients()"
	  << std::endl;
    abort_handler(-1);
    break;
  default: { //case ADD_COMBINE:
    const Sizet2DArray& combined_mi_map = data_rep->combinedMultiIndexMap;
    size_t i, num_combine = combined_mi_map.size();

    // overlay/add the level expansions, reindexing the sparse indices

    /*
    // Incurs a performance penalty by not assuming form of combined_mi_map[0]:
    combinedSparseIndices.clear(); // combined coeffs,grads cleared in overlay
    for (i=0, sp_it = sparseIndices.begin(), ec_it = expansionCoeffs.begin(),
	      eg_it = expansionCoeffGrads.begin();
	 i<num_combine; ++i, ++sp_it, ++ec_it, ++eg_it)
      overlay_expansion(sp_it->second, combined_mi_map[i], ec_it->second,
			eg_it->second, 1,  combinedSparseIndices,
			combinedExpCoeffs, combinedExpCoeffGrads);
    */

    // More efficient for (coarse/large) leading expansion but fragile in that
    // it ignores combined_mi_map[0] (but from SharedOrthogPolyApproxData::
    // pre_combine_data(), we know combined_mi_map[0] will sequence the
    // leading terms of combinedMultiIndex)
    sp_it = sparseIndices.begin();
    ec_it = expansionCoeffs.begin();  eg_it = expansionCoeffGrads.begin();
    // avoid overhead of unnecessary reindexing on first overlay
    combinedSparseIndices = sp_it->second; combinedExpCoeffs = ec_it->second;
    combinedExpCoeffGrads = eg_it->second; ++sp_it; ++ec_it; ++eg_it;
    // overlay with reindexing of sparse indices
    for (i=1; i<num_combine; ++i, ++sp_it, ++ec_it, ++eg_it)
      overlay_expansion(sp_it->second, combined_mi_map[i], ec_it->second,
			eg_it->second, 1,  combinedSparseIndices,
			combinedExpCoeffs, combinedExpCoeffGrads);
    break;
  }
  }
}


void RegressOrthogPolyApproximation::combined_to_active(bool clear_combined)
{
  // swap coefficients, reset keys and computed flags
  OrthogPolyApproximation::combined_to_active(clear_combined);

  // update sparseIndices and sparseSobolIndexMap
  if (!combinedSparseIndices.empty()) {
    sparseIndIter->second = combinedSparseIndices;
    // update sparseSobolIndexMap
    // Note 1: if sobol indices used in future combined roll-ups (e.g., for
    //   anisotropic sparse grid adaptation with total Sobol' indices), then
    //   will need to track combinedSparseSobolIndexMap and promote here.
    //   For now, regenerate sparseSobolIndexMap after roll ups are complete. 
    // Note 2: SharedOrthogPolyApproxData::combined_to_active() has promoted
    //   combinedMultiIndex to active and updated sobolIndexMap.
    SharedRegressOrthogPolyApproxData* data_rep
      = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
    update_sparse_sobol(combinedSparseIndices, data_rep->multi_index(),
			data_rep->sobolIndexMap);
    if (clear_combined)
      combinedSparseIndices.clear();
  }
}


void RegressOrthogPolyApproximation::run_regression()
{
  // Assume all function values are stored in top block of matrix in rows
  // 0 to num_surr_data_pts-1. Gradient information will be stored
  // in the bottom block of the matrix in rows num_surr_data_pts to
  // num_surr_data_pts + num_data_pts_grad * num_v. All the gradient 
  // information of point 0 will be stored consecutively then all the gradient
  // data of point 1, and so on.

  // Currently nothing is done  to modify the regression linear system matrices
  // A and B if modSurrData.anchor() is true, as currently modSurrData.anchor()
  // is always false. If in the future modSurrData.anchor() is enabled then
  // A must be adjusted to include the extra constraint information associated
  // with the anchor data. That is, if using EQ_CON_LEAST_SQUARES C matrix 
  // (top block of A ) must contain the fn and grad data of anchor point.
  // This will violate the first assumption discussed above and effect cross
  // validation. For this reason no modification is made as yet.

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;

  // Perform cross validation loop over degrees here.
  // Current cross validation will not work for equality 
  // constrained least squares
  if (data_rep->regressConfigOptions.crossValidation) // CV: MULTIPLE CS SOLVES
    run_cross_validation_expansion(); // updates all global bookkeeping
                                      // multiple RHS not currently supported
  else {
    RealMatrix A, B, points;
    build_linear_system( A, B, points );
    IntVector index_mapping;
    if ( data_rep->expConfigOptions.expCoeffsSolnApproach == 
	 ORTHOG_LEAST_INTERPOLATION ) { // SINGLE OLI SOLVE
      remove_faulty_data( A, B, points, index_mapping, faultInfo,
			  modSurrData.failed_response_data() );
      //faultInfo.under_determined = false;
      PCout << "Forming least interpolant for " << points.numCols()
	    << " points.\n";
      least_interpolation( points, B ); // updates all global bookkeeping
                                        // multiple RHS not currently supported
    }
    else { // SINGLE CS SOLVE
      RealMatrix points_dummy;
      remove_faulty_data( A, B, points_dummy, index_mapping, faultInfo,
			  modSurrData.failed_response_data() );
      //faultInfo.under_determined = A.numRows() < A.numCols();
      PCout << "Applying regression to compute " << data_rep->expansion_terms()
	    << " chaos coefficients using " << A.numRows() << " equations.\n";
      compressed_sensing(A, B); // updates all global bookkeeping
                                // includes support for multiple RHS
    }
  }
}


void RegressOrthogPolyApproximation::adapt_regression()
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  UShort2DArray& mi = data_rep->multi_index();
  Real rel_delta_star, abs_conv_tol = DBL_EPSILON, // for now
    rel_conv_tol = data_rep->expConfigOptions.convergenceTol;
  short basis_type = data_rep->expConfigOptions.expBasisType;

  // For a level 0 reference multi-index, we could safely omit evaluating the
  // initial CV error and always perform one cycle of select_best_active_mi(),
  // since the CV error would presumably not be low enough to trigger immediate
  // termination.  However, this requires care in selecting the lowest CV error
  // on the first cycle (largest CV change would need to utilize worst - best;
  // DBL_MAX cannot be used as an error ref due to loss of precision).  We also
  // need to support general reference points, since adaptive algorithms can
  // stall without good starting points.  Therefore, go ahead and compute the
  // CV err for the reference candidate basis.
  bestAdaptedMultiIndex = mi;
  SizetSet& sparse_ind = sparseIndIter->second;
  cvErrorRef
    = run_cross_validation_solver(bestAdaptedMultiIndex, expCoeffsIter->second,
				  sparse_ind);
  PCout << "<<<<< Cross validation error reference = " << cvErrorRef << '\n';

  // absolute error instead of delta error for initial tolerance check
  unsigned short soft_conv_limit = data_rep->expConfigOptions.softConvLimit,
    soft_conv_count = (cvErrorRef > abs_conv_tol) ? 0 : 1;

  if (soft_conv_count < soft_conv_limit) {
    adaptedMultiIndex = bestAdaptedMultiIndex; // starting point for adaptation
    adaptedSparseIndices = sparse_ind; // needed when using restriction
    switch (basis_type) {
    case ADAPTED_BASIS_GENERALIZED:
      data_rep->lsgDriver.initialize_sets(); // initialize the active sets
      break;
    //case ADAPTED_BASIS_EXPANDING_FRONT:
    //  break;
    }
  }

  while (soft_conv_count < soft_conv_limit) {
    // invoke run_cross_validation_solver() for each candidate and select best
    rel_delta_star = (basis_type == ADAPTED_BASIS_GENERALIZED) ?
      select_best_active_multi_index() : select_best_basis_expansion();
    // Several possible convergence assessment approaches:
    // > capture when no CV error reduction is obtained for any candidate:
    //   no reduction --> rel_delta_star <= 0 --> increment soft conv count.
    // > relative change (rel_delta_star): a good conv_tol might be O(1);
    //   the Dakota default 1.e-4 may be too tight.
    // > absolute change (delta_star): a good conv_tol might be O(10^{-16});
    //   the Dakota 1.e-4 default is definitely too loose.
    if (rel_delta_star > rel_conv_tol)
      soft_conv_count = 0;
    else
      ++soft_conv_count;
  }

  // different finalize: don't add in any remaining evaluated sets; rather,
  // we need to backtrack and restore the best solution with lowest CV error
  adaptedMultiIndex.clear(); adaptedSparseIndices.clear();
  data_rep->clear_adapted();

  // Once done for this QoI, append adaptedMultiIndex to shared multiIndex,
  // update sparseIndices (which corresponds to bestAdaptedMultiIndex) to
  // point into shared multiIndex, reorder expansionCoeff{s,Grads} as needed,
  // and clear bestAdaptedMultiIndex
  data_rep->
    append_sparse_multi_index(sparse_ind, bestAdaptedMultiIndex, mi,
			      expCoeffsIter->second, expCoeffGradsIter->second);
  bestAdaptedMultiIndex.clear();

  // now update sobolIndexMap and sparseSobolIndexMap
  data_rep->update_component_sobol(mi);
  update_sparse_sobol(sparse_ind, mi, data_rep->sobolIndexMap);
}


Real RegressOrthogPolyApproximation::select_best_active_multi_index()
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  LightweightSparseGridDriver* lsg_driver = &data_rep->lsgDriver;
  const UShortArray& key = data_rep->activeKey;

  // Perform a (coarse-grained) restriction operation that updates
  // oldMultiIndex and activeMultiIndex. This requires the same heavyweight
  // update of bestAdaptedMultiIndex used for expanding front (bottom of fn).
  SizetSet save_tp; // set of tensor products to retain
  bool mi_restricted = data_rep->
    set_restriction(adaptedMultiIndex, adaptedSparseIndices, save_tp);
  if (mi_restricted)
    lsg_driver->prune_sets(save_tp); // update from scratch

  const UShortArraySet& active_mi = lsg_driver->active_multi_index();
  UShortArraySet::const_iterator cit, cit_star;
  Real curr_cv_err, cv_err_star, rel_delta, rel_delta_star = -DBL_MAX;
  RealVector curr_exp_coeffs; SizetSet curr_sparse_ind;
  // Reevaluate the effect of every active set every time
  for (cit=active_mi.begin(); cit!=active_mi.end(); ++cit) {

    // increment grid with current candidate
    const UShortArray& trial_set = *cit;
    PCout << "\n>>>>> Evaluating trial index set:\n" << trial_set;
    //lsg_driver->push_trial_set(trial_set);

    // trial index set -> tpMultiIndex -> append to (local) adaptedMultiIndex
    if (data_rep->push_available())
      // cannot assume monotonicity in bookkeeping due to restriction above
      data_rep->push_trial_set(trial_set, adaptedMultiIndex, false);
    else
      data_rep->increment_trial_set(trial_set, adaptedMultiIndex);

    // Solve CS with cross-validation applied to solver settings (e.g., noise
    // tolerance), but not expansion order (since we are manually adapting it).
    // The CV error is L2 (sum of squares of mismatch) and is non-negative.
    curr_cv_err = run_cross_validation_solver(adaptedMultiIndex,
					      curr_exp_coeffs, curr_sparse_ind);

    // Smallest absolute error is largest decrease from consistent ref.
    // rel_delta is a signed quantity in order to detect when best case error
    // increases relative to ref -> terminate or increment soft conv counter
    // (allow some number of successive increases before abandoning hope).
    rel_delta = (cvErrorRef > 0.) ? 1. - curr_cv_err / cvErrorRef
                                  : cvErrorRef - curr_cv_err;
    if (data_rep->regressConfigOptions.normalizeCV) {
      // Normalize rel_delta based on size of candidate basis expansion
      // (number of unique points added is equivalent to number of candidate
      // expansion terms added for Gauss quadrature, but not other cases)
      //int new_terms = lsg_driver->unique_trial_points();
      size_t new_terms = adaptedMultiIndex.size()
	               - data_rep->tpMultiIndexMapRef[key].back();
      rel_delta /= new_terms;
    }
    PCout << "\n<<<<< Trial set refinement metric = " << rel_delta
	  << " (relative error reduction)\n";

    // track best increment evaluated thus far
    if (rel_delta > rel_delta_star) {
      cit_star = cit; rel_delta_star = rel_delta;     // best for this iteration
      adaptedSparseIndices = curr_sparse_ind; // best for this iteration
      if (rel_delta_star > 0.) {
	cv_err_star           = curr_cv_err;     // best overall
	expCoeffsIter->second = curr_exp_coeffs; // best overall
	sparseIndIter->second = curr_sparse_ind; // best overall
      }
    }

    // restore previous state (destruct order is reversed from construct order)
    data_rep->decrement_trial_set(trial_set, adaptedMultiIndex);
    //lsg_driver->pop_set();
  }
  const UShortArray& best_set = *cit_star;
  PCout << "\n<<<<< Evaluation of active index sets completed.\n"
	<< "\n<<<<< Index set selection:\n" << best_set;

  // apply best increment; can assume monotonicity in bookkeeping so long as
  // Map,MapRef are updated for all candidates by push/increment above
  data_rep->push_trial_set(best_set, adaptedMultiIndex, true);
  lsg_driver->update_sets(best_set); // if restriction: rebuilt on next iter

  // update reference points and current best soln if CV error has been reduced
  // Note: bestAdaptedMultiIndex assignment is heavier weight than currently
  // required, prior to use of restriction at top of this fn.
  if (rel_delta_star > 0.) {
    cvErrorRef = cv_err_star; bestAdaptedMultiIndex = adaptedMultiIndex;
    PCout << "<<<<< New cross validation error reference = " << cvErrorRef
	  << '\n';
  }

  return rel_delta_star;
}


Real RegressOrthogPolyApproximation::select_best_basis_expansion()
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;

  // restrict adaptedMultiIndex and then advance this frontier.  Notes:
  // > candidateBasisExp is an array of sets that store multi-index increments;
  //   thus, the adaptedMultiIndex reference must be preserved to recover best.
  // > input/output from this function is a candidate multi_index along with 
  //   its sparse solution (expansionCoeffs + sparseIndices key); multi-index
  //   frontiers and advancement sets are only used internal to this fn.
  // > for soft convergence to be meaningful (allowing adaptation to step over
  //   dead spots), the restriction must be applied to the selected increment
  //   from the previous iteration, even if not an improvement in CV error.
  //   Thus, we use adaptedSparseIndices below instead of sparseIndices.
  UShortArraySetArray candidate_basis_exp;
  if (data_rep->regressConfigOptions.advanceByFrontier) {
    // define and advance a frontier without gaps
    frontier_restriction(adaptedMultiIndex, adaptedSparseIndices);
    advance_multi_index_front(adaptedMultiIndex, candidate_basis_exp);
  }
  else { // advance complete adaptedMultiIndex, including any gaps
    sparse_restriction(adaptedMultiIndex, adaptedSparseIndices);
    advance_multi_index(adaptedMultiIndex, candidate_basis_exp);
  }

  // Evaluate the effect of each candidate basis expansion
  size_t i, i_star = 0, num_candidates = candidate_basis_exp.size(), 
    size_star = adaptedMultiIndex.size();
  Real curr_cv_err, cv_err_star, rel_delta, rel_delta_star = -DBL_MAX;
  RealVector curr_exp_coeffs; SizetSet curr_sparse_ind;
  for (i=0; i<num_candidates; ++i) {

    // append candidate expansion to adaptedMultiIndex.  advance_multi_index()
    // already provides unique additions --> append_multi_index() is inefficient
    // since duplication checks are redundant.
    const UShortArraySet& candidate_exp = candidate_basis_exp[i];
    //append_multi_index(candidate_exp, adaptedMultiIndex);
    PCout << "\n>>>>> Evaluating trial basis " << i+1 << " expanded from "
	  << adaptedMultiIndex.size() << " to ";
    adaptedMultiIndex.insert(adaptedMultiIndex.end(), candidate_exp.begin(),
			     candidate_exp.end());
    PCout << adaptedMultiIndex.size() << " terms\n";

    // Solve CS with cross-validation applied to solver settings (e.g., noise
    // tolerance), but not expansion order (since we are manually adapting it).
    // The CV error is L2 (sum of squares of mismatch) and is non-negative.
    curr_cv_err = run_cross_validation_solver(adaptedMultiIndex,
					      curr_exp_coeffs, curr_sparse_ind);

    // Smallest absolute error is largest decrease from consistent ref.
    // rel_delta is a signed quantity in order to detect when best case error
    // increases relative to ref -> terminate or increment soft conv counter
    // (allow some number of successive increases before abandoning hope).
    rel_delta = (cvErrorRef > 0.) ? 1. - curr_cv_err / cvErrorRef
                                  : cvErrorRef - curr_cv_err;
    if (data_rep->regressConfigOptions.normalizeCV) // not as well justified
      // Normalize rel_delta based on size of candidate basis expansion
      // (number of unique points added is equivalent to number of candidate
      // expansion terms added for Gauss quadrature, but not other cases)
      rel_delta /= candidate_exp.size();
    PCout << "\n<<<<< Trial set refinement metric = " << rel_delta
	  << " (relative error reduction)\n";

    // track best increment evaluated thus far.  Note: an increment is always
    // selected even if no improvement in CV error --> rely on soft convergence
    // and best solution tracking.
    if (rel_delta > rel_delta_star) {
      i_star = i; rel_delta_star = rel_delta; // best for this iteration
      size_star = adaptedMultiIndex.size();   // best for this iteration
      adaptedSparseIndices = curr_sparse_ind; // best for this iteration
      if (rel_delta_star > 0.) {              // reduction in best CV error
	cv_err_star           = curr_cv_err;        // best overall
	expCoeffsIter->second = curr_exp_coeffs;    // best overall
	sparseIndIter->second = curr_sparse_ind;    // best overall
      }
    }
  }

  //const UShort2DArray& best_set = candidate_basis_exp[i_star];
  size_t id_star = i_star + 1;
  PCout << "\n<<<<< Evaluation of candidate basis expansions completed.\n"
	<< "\n<<<<< Selection of basis expansion set " << id_star << '\n';
  // rewind adaptedMultiIndex to best candidate
  if (id_star != num_candidates)
    adaptedMultiIndex.resize(size_star);
  // update reference points for best solution if CV error has been reduced
  if (rel_delta_star > 0.) {
    cvErrorRef = cv_err_star; bestAdaptedMultiIndex = adaptedMultiIndex;
    PCout << "<<<<< New cross validation error reference = " << cvErrorRef
	  << '\n';
  }

  return rel_delta_star;
}


void RegressOrthogPolyApproximation::
overlay_expansion(const SizetSet& sparse_ind, const SizetArray& multi_index_map,
		  const RealVector& exp_coeffs, const RealMatrix& exp_grads,
		  int coeff, SizetSet& sparse_ind_sum,
		  RealVector& exp_coeffs_sum, RealMatrix& exp_grads_sum)
{
  // update sparse_ind_sum w/ new contributions relative to updated multi-index
  // generated by SharedOrthogPolyApproxData::pre_combine_data():
  // 1. previous sparse_ind_sum into the updated multi-index are still valid
  //    (new mi entries are appended without reordering such that previous mi
  //    ordering is preserved)
  // 2. new sparse_ind into previous multi-index terms may be interwoven with
  //    sparse_ind_sum due to ordered set insertion -> insert combined index
  //    from multi_index_map (reject duplicates).  This case invalidates the
  //    previous sparse mapping and requires re-indexing (see mapping from
  //    sparse_ind_sum to combined_sparse_ind below)
  // 3. new sparse_ind into appended mi are new -> insert combined index from
  //    multi_index_map
  StSCIter cit;
  SizetSet combined_sparse_ind = sparse_ind_sum; // copy
  for (cit=sparse_ind.begin(); cit!=sparse_ind.end(); ++cit)
    combined_sparse_ind.insert(multi_index_map[*cit]); // duplicates rejected

  size_t i, j, combined_index, num_deriv_v,
    num_combined_terms = combined_sparse_ind.size();
  // initialize combined_exp_coeff{s,_grads}.  previous expansionCoeff{s,Grads}
  // cannot simply be resized since they must correspond to the combined sparse
  // indices, which can be reordered (see note 2. above)
  RealVector combined_exp_coeffs; RealMatrix combined_exp_grads;
  if (expansionCoeffFlag)
    combined_exp_coeffs.size(num_combined_terms); // init to 0
  if (expansionCoeffGradFlag) {
    num_deriv_v = exp_grads_sum.numRows();
    combined_exp_grads.shape(num_deriv_v, num_combined_terms); // init to 0
  }

  for (i=0, cit=sparse_ind_sum.begin(); cit!=sparse_ind_sum.end(); ++i, ++cit) {
    combined_index = find_index(combined_sparse_ind, *cit);
    if (expansionCoeffFlag)
      combined_exp_coeffs[combined_index] = exp_coeffs_sum[i];
    if (expansionCoeffGradFlag)
      copy_data(exp_grads_sum[i], num_deriv_v,
		combined_exp_grads[combined_index]);
  }

  for (i=0, cit=sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit) {
    combined_index = find_index(combined_sparse_ind, multi_index_map[*cit]);
    if (expansionCoeffFlag)
      combined_exp_coeffs[combined_index] += coeff * exp_coeffs[i];
    if (expansionCoeffGradFlag) {
      Real*       combined_exp_grad = combined_exp_grads[combined_index];
      const Real* grad_i            = exp_grads[i];
      for (j=0; j<num_deriv_v; ++j)
	combined_exp_grad[j] += coeff * grad_i[j];
    }
  }

  // overlay is complete; can now overwrite previous state
  sparse_ind_sum = combined_sparse_ind;
  if (expansionCoeffFlag)     exp_coeffs_sum = combined_exp_coeffs;
  if (expansionCoeffGradFlag) exp_grads_sum  = combined_exp_grads;
}


void RegressOrthogPolyApproximation::
multiply_expansion(const UShort2DArray& multi_index_a,
		   const SizetSet&      sparse_ind_a,
		   const RealVector&    exp_coeffs_a,
		   const RealMatrix&    exp_grads_a,
		   const UShort2DArray& multi_index_b,
		   const SizetSet&      sparse_ind_b,
		   const RealVector&    exp_coeffs_b,
		   const RealMatrix&    exp_grads_b,
		   const UShort2DArray& multi_index_c, SizetSet& sparse_ind_c,
		   RealVector& exp_coeffs_c, RealMatrix& exp_grads_c)
{
  // sparsity in exp_{coeffs,grads}_c determined *after* all multi_index_c
  // projection terms have been computed

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  size_t i, j, k, v, si, sj, num_v = sharedDataRep->numVars,
    num_deriv_v = exp_grads_a.numRows(), num_c = multi_index_c.size();

  // precompute 1D basis triple products required
  unsigned short max_a, max_b, max_c; UShortMultiSet max_abc;
  OrthogonalPolynomial* poly_rep_v;   StSCIter ait, bit;
  for (v=0; v<num_v; ++v) {
    max_a = max_b = max_c = 0; max_abc.clear();
    // could track max_abc within combine_coefficients() and pass in, but would
    // need max orders for each dimension for both factors plus their product.
    // Since this would be awkward and only marginally more efficient, just
    // compute them here from the available multi-index arrays.
    for (ait=++sparse_ind_a.begin(); ait!=sparse_ind_a.end(); ++ait)
      if (multi_index_a[*ait][v] > max_a)
	max_a = multi_index_a[*ait][v];
    for (bit=++sparse_ind_b.begin(); bit!=sparse_ind_b.end(); ++bit)
      if (multi_index_b[*bit][v] > max_b)
	max_b = multi_index_b[*bit][v];
    for (k=1; k<num_c; ++k) // sparse_ind_c to be defined from product results
      if (multi_index_c[k][v] > max_c)
	max_c = multi_index_c[k][v];
    max_abc.insert(max_a); max_abc.insert(max_b); max_abc.insert(max_c); 
    poly_rep_v
      = (OrthogonalPolynomial*)(data_rep->polynomialBasis[v].polynomial_rep());
    poly_rep_v->precompute_triple_products(max_abc);
  }

  // For c = a * b, compute coefficient of product expansion as:
  // \Sum_k c_k \Psi_k = \Sum_i \Sum_j a_i b_j \Psi_i \Psi_j
  //    c_k <\Psi_k^2> = \Sum_i \Sum_j a_i b_j <\Psi_i \Psi_j \Psi_k>
  RealVector exp_coeffs_tmp_c; RealVectorArray exp_grads_tmp_c;
  if (expansionCoeffFlag)
    exp_coeffs_tmp_c.size(num_c);     // init to 0
  if (expansionCoeffGradFlag) {
    exp_grads_tmp_c.resize(num_deriv_v);
    for (v=0; v<num_deriv_v; ++v)
      exp_grads_tmp_c[v].size(num_c); // init to 0
  }
  Real trip_prod, trip_prod_v, norm_sq_k; bool non_zero;
  for (k=0; k<num_c; ++k) {
    for (i=0, ait=sparse_ind_a.begin(); ait!=sparse_ind_a.end(); ++i, ++ait) {
      si = *ait;
      for (j=0, bit=sparse_ind_b.begin(); bit!=sparse_ind_b.end(); ++j, ++bit) {
	sj = *bit;
	trip_prod = 1.;
	for (v=0; v<num_v; ++v) {
	  poly_rep_v = (OrthogonalPolynomial*)
	    (data_rep->polynomialBasis[v].polynomial_rep());
	  non_zero = poly_rep_v->triple_product(multi_index_a[si][v],
	    multi_index_b[sj][v], multi_index_c[k][v], trip_prod_v);
	  if (non_zero) trip_prod *= trip_prod_v;
	  else          break;
	}
	if (non_zero) {
	  if (expansionCoeffFlag)
	    exp_coeffs_tmp_c[k]
	      += exp_coeffs_a[i] * exp_coeffs_b[j] * trip_prod;
	  if (expansionCoeffGradFlag) {
	    const Real* exp_grads_a_i = exp_grads_a[i];
	    const Real* exp_grads_b_j = exp_grads_b[j];
	    for (v=0; v<num_deriv_v; ++v)
	      exp_grads_tmp_c[v][k] += (exp_coeffs_a[i] * exp_grads_b_j[v]
		+ exp_coeffs_b[j] * exp_grads_a_i[v]) * trip_prod;
	  }
	}
      }
    }
    norm_sq_k = data_rep->norm_squared(multi_index_c[k]);
    if (expansionCoeffFlag)
      exp_coeffs_tmp_c[k] /= norm_sq_k;
    if (expansionCoeffGradFlag)
      for (v=0; v<num_deriv_v; ++v)
	exp_grads_tmp_c[v][k] /= norm_sq_k;
  }

  // update sparse bookkeeping based on nonzero terms in dense tmp arrays
  sparse_ind_c.clear();
  if (expansionCoeffFlag)
    update_sparse_indices(exp_coeffs_tmp_c.values(), num_c, sparse_ind_c);
  if (expansionCoeffGradFlag)
    for (v=0; v<num_deriv_v; ++v)
      update_sparse_indices(exp_grads_tmp_c[v].values(), num_c, sparse_ind_c);
  // update exp_{coeffs,grads}_c from dense tmp arrays & updated sparse indices
  if (expansionCoeffFlag)
    update_sparse_coeffs(exp_coeffs_tmp_c.values(), exp_coeffs_c, sparse_ind_c);
  if (expansionCoeffGradFlag) {
    exp_grads_c.shapeUninitialized(num_deriv_v, sparse_ind_c.size());
    for (v=0; v<num_deriv_v; ++v)
      update_sparse_coeff_grads(exp_grads_tmp_c[v].values(), v, exp_grads_c,
				sparse_ind_c);
  }
}


Real RegressOrthogPolyApproximation::
value(const RealVector& x, const UShort2DArray& mi,
      const RealVector& exp_coeffs, const SizetSet& sparse_ind)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "RegressOrthogPolyApproximation::value()" << std::endl;
    abort_handler(-1);
  }

  Real approx_val = 0.;
  size_t i; StSCIter cit;
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  for (i=0, cit=sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit)
    approx_val += exp_coeffs[i]
      * data_rep->multivariate_polynomial(x, mi[*cit]);
  return approx_val;
}


const RealVector& RegressOrthogPolyApproximation::
gradient_basis_variables(const RealVector& x, const UShort2DArray& mi,
			 const RealVector& exp_coeffs,
			 const SizetSet& sparse_ind)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in RegressOrthogPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t i, j, num_v = sharedDataRep->numVars;  StSCIter cit;
  if (approxGradient.length() != num_v) approxGradient.size(num_v); // init to 0
  else                                  approxGradient = 0.;

  // sum expansion to get response gradient prediction
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  for (i=0, cit=sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit) {
    const RealVector& term_i_grad
      = data_rep->multivariate_polynomial_gradient_vector(x, mi[*cit]);
    Real coeff_i = exp_coeffs[i];
    for (j=0; j<num_v; ++j)
      approxGradient[j] += coeff_i * term_i_grad[j];
  }
  return approxGradient;
}


const RealVector& RegressOrthogPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv,
			 const UShort2DArray& mi, const RealVector& exp_coeffs,
			 const SizetSet& sparse_ind)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in RegressOrthogPoly"
	  << "Approximation::gradient_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t i, j, num_v = dvv.size();  StSCIter cit;
  if (approxGradient.length() != num_v) approxGradient.size(num_v); // init to 0
  else                                  approxGradient = 0.;

  // sum expansion to get response gradient prediction
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  for (i=0, cit=sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit) {
    const RealVector& term_i_grad
      = data_rep->multivariate_polynomial_gradient_vector(x, mi[*cit], dvv);
    Real coeff_i = exp_coeffs[i];
    for (j=0; j<num_v; ++j)
      approxGradient[j] += coeff_i * term_i_grad[j];
  }
  return approxGradient;
}


const RealVector& RegressOrthogPolyApproximation::
gradient_nonbasis_variables(const RealVector& x, const UShort2DArray& mi,
			    const RealMatrix& exp_coeff_grads,
			    const SizetSet& sparse_ind)
{
  // Error check for required data
  if (!expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in RegressOrth"
	  << "ogPolyApproximation::gradient_nonbasis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t i, j, num_v = exp_coeff_grads.numRows();  StSCIter cit;
  if (approxGradient.length() != num_v) approxGradient.size(num_v); // init to 0
  else                                  approxGradient = 0.;

  // sum expansion to get response gradient prediction
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  for (i=0, cit=sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit) {
    Real term_i = data_rep->multivariate_polynomial(x, mi[*cit]);
    const Real* exp_coeff_grad_i = exp_coeff_grads[i];
    for (j=0; j<num_v; ++j)
      approxGradient[j] += exp_coeff_grad_i[j] * term_i;
  }
  return approxGradient;
}


const RealSymMatrix& RegressOrthogPolyApproximation::
hessian_basis_variables(const RealVector& x, const UShort2DArray& mi,
			const RealVector& exp_coeffs,
			const SizetSet& sparse_ind)
{
  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in RegressOrthogPoly"
	  << "Approximation::hessian_basis_variables()" << std::endl;
    abort_handler(-1);
  }

  size_t i, row, col, num_v = sharedDataRep->numVars;  StSCIter cit;
  if (approxHessian.numRows() != num_v) approxHessian.shape(num_v); // init to 0
  else                                  approxHessian = 0.;

  // sum expansion to get response hessian prediction
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  for (i=0, cit=sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit) {
    const RealSymMatrix& term_i_hess
      = data_rep->multivariate_polynomial_hessian_matrix(x, mi[*cit]);
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


/** In this case, a subset of the expansion variables are random
    variables and the mean of the expansion involves evaluating the
    expectation over this subset. */
Real RegressOrthogPolyApproximation::mean(const RealVector& x)
{
  if (sparseIndIter == sparseIndices.end() || sparseIndIter->second.empty())
    return OrthogPolyApproximation::mean(x);

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "RegressOrthogPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const SizetList& nrand_ind = data_rep->nonRandomIndices;
  bool all_mode = !nrand_ind.empty();
  const UShortArray& key = data_rep->activeKey;
  if (all_mode && (compMeanIter->second & 1) &&
      data_rep->match_nonrandom_vars(x, xPrevMean[key]))
    return expMomentsIter->second[0];

  const UShort2DArray& mi = data_rep->multi_index();
  const SizetSet& sparse_ind = sparseIndIter->second;
  const RealVector& exp_coeffs = expCoeffsIter->second;
  Real mean = exp_coeffs[0];
  size_t i; StSCIter cit;
  for (i=1, cit=++sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit) {
    const UShortArray& mi_i = mi[*cit];
    // expectations are zero for expansion terms with nonzero random indices
    if (data_rep->zero_random(mi_i)) {
      mean += exp_coeffs[i] *
	data_rep->multivariate_polynomial(x, mi_i, nrand_ind);
#ifdef DEBUG
      PCout << "Mean estimate inclusion: term index = " << i << " Psi = "
	    << data_rep->multivariate_polynomial(x, mi_i, nrand_ind)
	    << " PCE coeff = " << exp_coeffs[i] << " total = " << mean
	    << std::endl;
#endif // DEBUG
    }
  }

  if (all_mode) {
    expMomentsIter->second[0] = mean;
    compMeanIter->second |= 1;  xPrevMean[key] = x;
  }
  return mean;
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
const RealVector& RegressOrthogPolyApproximation::
mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  if (sparseIndIter == sparseIndices.end() || sparseIndIter->second.empty())
    return OrthogPolyApproximation::mean_gradient(x, dvv);

  // if already computed, return previous result
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const SizetList& nrand_ind = data_rep->nonRandomIndices;
  bool all_mode = !nrand_ind.empty();
  const UShortArray& key = data_rep->activeKey;
  if ( all_mode && (compMeanIter->second & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevMeanGrad[key]) )
    // && dvv == dvvPrev)
    return momentGradsIter->second[0];

  size_t i, j, deriv_index, num_deriv_v = dvv.size(),
    cntr = 0; // insertions carried in order within expansionCoeffGrads
  RealVector& mean_grad = momentGradsIter->second[0];
  if (mean_grad.length() != num_deriv_v)
    mean_grad.sizeUninitialized(num_deriv_v);
  const UShort2DArray& mi = data_rep->multi_index();
  const RealVector& exp_coeffs      = expCoeffsIter->second;
  const RealMatrix& exp_coeff_grads = expCoeffGradsIter->second;
  const SizetSet&   sparse_ind      = sparseIndIter->second;  StSCIter cit;
  for (i=0; i<num_deriv_v; ++i) {
    Real& grad_i = mean_grad[i];
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    bool random = data_rep->randomVarsKey[deriv_index];
    if (random) { // deriv w.r.t. des var insertion
      if (!expansionCoeffGradFlag) { // error check for required data
	PCerr << "Error: expansion coefficient gradients not defined in Regress"
	      << "OrthogPolyApproximation::mean_gradient()." << std::endl;
	abort_handler(-1);
      }
      grad_i = exp_coeff_grads[0][cntr];
    }
    else {
      grad_i = 0.;
      if (!expansionCoeffFlag) { // check for reqd data
	PCerr << "Error: expansion coefficients not defined in RegressOrthog"
	      << "PolyApproximation::mean_gradient()" << std::endl;
	abort_handler(-1);
      }
    }
    for (j=1, cit=++sparse_ind.begin(); cit!=sparse_ind.end(); ++j, ++cit) {
      const UShortArray& mi_j = mi[*cit];
      // expectations are zero for expansion terms with nonzero random indices
      if (data_rep->zero_random(mi_j)) {
	// In both cases below, term to differentiate is alpha_j(s) Psi_j(s)
	// since <Psi_j>_xi = 1 for included terms.  The difference occurs
	// based on whether a particular s_i dependence appears in alpha
	// (for inserted) or Psi (for augmented).
	if (random)
	  // -------------------------------------------
	  // derivative w.r.t. design variable insertion
	  // -------------------------------------------
	  grad_i += exp_coeff_grads[j][cntr] * data_rep->
	    multivariate_polynomial(x, mi_j, nrand_ind);
	else
	  // ----------------------------------------------
	  // derivative w.r.t. design variable augmentation
	  // ----------------------------------------------
	  grad_i += exp_coeffs[j] * data_rep->
	    multivariate_polynomial_gradient(x, deriv_index, mi_j, nrand_ind);
      }
    }
    if (random) // deriv w.r.t. des var insertion
      ++cntr;
  }
  if (all_mode) { compMeanIter->second |=  2; xPrevMeanGrad[key] = x; }
  else   compMeanIter->second &= ~2; // deactivate 2-bit: protect mixed use
  return mean_grad;
}


Real RegressOrthogPolyApproximation::
variance(const UShort2DArray& mi, const RealVector& exp_coeffs,
	 const SizetSet& sparse_ind)
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  size_t i; StSCIter cit;
  Real var = 0.;
  for (i=1, cit=++sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit)
    var += exp_coeffs[i] * exp_coeffs[i] * data_rep->norm_squared(mi[*cit]);
  return var;
}


Real RegressOrthogPolyApproximation::
covariance(const UShort2DArray& mi,    const RealVector& exp_coeffs,
	   const SizetSet& sparse_ind, const RealVector& exp_coeffs_2,
	   const SizetSet& sparse_ind_2)
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  size_t i1, i2, si1, si2; StSCIter cit1, cit2;
  Real covar = 0.;
  if (sparse_ind.empty()) // mixed: one set is sparse
    for (i2=1, cit2=++sparse_ind_2.begin(); cit2!=sparse_ind_2.end();
	 ++i2, ++cit2) {
      si2 = *cit2;
      covar += exp_coeffs[si2] * exp_coeffs_2[i2]
	    *  data_rep->norm_squared(mi[si2]);
    }
  else if (sparse_ind_2.empty()) // mixed: one set is sparse
    for (i1=1, cit1=++sparse_ind.begin(); cit1!=sparse_ind.end();
	 ++i1, ++cit1) {
      si1 = *cit1;
      covar += exp_coeffs[i1] * exp_coeffs_2[si1]
	    *  data_rep->norm_squared(mi[si1]);
    }
  else { // both are (distinctly) sparse
    i1=1; i2=1; cit1=++sparse_ind.begin(); cit2=++sparse_ind_2.begin();
    while (cit1!=sparse_ind.end() && cit2!=sparse_ind_2.end()) {
      si1 = *cit1; si2 = *cit2;
      if (si1 == si2) {
	covar += exp_coeffs[i1] * exp_coeffs_2[i2]
	      *  data_rep->norm_squared(mi[si1]);
	++i1; ++cit1; ++i2; ++cit2;
      }
      else if (si1 < si2) { ++i1; ++cit1; }
      else                { ++i2; ++cit2; }
    }
  }
  return covar;
}


Real RegressOrthogPolyApproximation::
covariance(PolynomialApproximation* poly_approx_2)
{
  RegressOrthogPolyApproximation* ropa_2
    = (RegressOrthogPolyApproximation*)poly_approx_2;
  if ( ( sparseIndIter == sparseIndices.end() ||
	 sparseIndIter->second.empty()) &&
       ( ropa_2->sparseIndIter == ropa_2->sparseIndices.end() ||
	 ropa_2->sparseIndIter->second.empty()) )
    return OrthogPolyApproximation::covariance(poly_approx_2);

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  bool same = (ropa_2 == this);

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !ropa_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "RegressOrthogPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if (same) {
    bool std_mode = data_rep->nonRandomIndices.empty();
    if (std_mode && (compVarIter->second & 1))
      return expMomentsIter->second[1];
    else {
      Real var = variance(data_rep->multi_index(), expCoeffsIter->second,
			  sparseIndIter->second);
      if (std_mode)
	{ expMomentsIter->second[1] = var; compVarIter->second |= 1; }
      return var;
    }
  }
  else
    return covariance(data_rep->multi_index(), expCoeffsIter->second,
		      sparseIndIter->second, ropa_2->expCoeffsIter->second,
		      ropa_2->sparseIndIter->second);
}


Real RegressOrthogPolyApproximation::
combined_covariance(PolynomialApproximation* poly_approx_2)
{
  RegressOrthogPolyApproximation* ropa_2
    = (RegressOrthogPolyApproximation*)poly_approx_2;
  if (combinedSparseIndices.empty() && ropa_2->combinedSparseIndices.empty())
    return OrthogPolyApproximation::covariance(poly_approx_2);

  bool same = (ropa_2 == this);

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !ropa_2->expansionCoeffFlag ) ) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "RegressOrthogPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  if (same) {
    bool std_mode = data_rep->nonRandomIndices.empty();
    if (std_mode && (compVarIter->second & 1))
      return expMomentsIter->second[1];
    else {
      Real var = variance(data_rep->combinedMultiIndex, combinedExpCoeffs,
			  combinedSparseIndices);
      if (std_mode)
	{ expMomentsIter->second[1] = var; compVarIter->second |= 1; }
      return var;
    }
  }
  else
    return covariance(data_rep->combinedMultiIndex, combinedExpCoeffs,
		      combinedSparseIndices, ropa_2->combinedExpCoeffs,
		      ropa_2->combinedSparseIndices);
}


Real RegressOrthogPolyApproximation::
covariance(const RealVector& x, const UShort2DArray& mi,
	   const RealVector& exp_coeffs,   const SizetSet& sparse_ind,
	   const RealVector& exp_coeffs_2, const SizetSet& sparse_ind_2)
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const SizetList&  rand_ind = data_rep->randomIndices;
  const SizetList& nrand_ind = data_rep->nonRandomIndices;

  Real covar = 0.;
  if (sparse_ind.empty()) { // mixed mode
    size_t i1, i2, num_mi = mi.size(); StSCIter cit2;
    for (i1=1; i1<num_mi; ++i1) {
      const UShortArray& mi1 = mi[i1];
      if (!data_rep->zero_random(mi1)) {
	Real coeff_norm_poly = exp_coeffs[i1] * 
	  data_rep->norm_squared(mi1, rand_ind) *
	  data_rep->multivariate_polynomial(x, mi1, nrand_ind);
	for (cit2=++sparse_ind_2.begin(), i2=1; cit2!=sparse_ind_2.end();
	     ++cit2, ++i2) {
	  const UShortArray& mi2 = mi[*cit2];
	  if (data_rep->match_random_key(mi1, mi2))
	    covar += coeff_norm_poly * exp_coeffs_2[i2] * 
	      data_rep->multivariate_polynomial(x, mi2, nrand_ind);
	}
      }
    }
  }
  else if (sparse_ind_2.empty()) { // mixed mode
    size_t i1, i2, num_mi = mi.size(); StSCIter cit1;
    for (cit1=++sparse_ind.begin(), i1=1; cit1!=sparse_ind.end(); ++cit1,++i1) {
      const UShortArray& mi1 = mi[*cit1];
      if (!data_rep->zero_random(mi1)) {
	Real coeff_norm_poly = exp_coeffs[i1] * 
	  data_rep->norm_squared(mi1, rand_ind) *
	  data_rep->multivariate_polynomial(x, mi1, nrand_ind);
	for (i2=1; i2<num_mi; ++i2) {
	  const UShortArray& mi2 = mi[i2];
	  if (data_rep->match_random_key(mi1, mi2))
	    covar += coeff_norm_poly * exp_coeffs_2[i2] * 
	      data_rep->multivariate_polynomial(x, mi2, nrand_ind);
	}
      }
    }
  }
  else {
    size_t i1, i2; StSCIter cit1, cit2;
    for (cit1=++sparse_ind.begin(), i1=1; cit1!=sparse_ind.end(); ++cit1,++i1) {
      // For r = random_vars and nr = non_random_vars,
      // sigma^2_R(nr) = < (R(r,nr) - \mu_R(nr))^2 >_r
      // -> only include terms from R(r,nr) which don't appear in \mu_R(nr)
      const UShortArray& mi1 = mi[*cit1];
      if (!data_rep->zero_random(mi1)) {
	Real coeff_norm_poly = exp_coeffs[i1] * 
	  data_rep->norm_squared(mi1, rand_ind) *
	  data_rep->multivariate_polynomial(x, mi1, nrand_ind);
	for (cit2=++sparse_ind_2.begin(), i2=1; cit2!=sparse_ind_2.end();
	     ++cit2, ++i2) {
	  const UShortArray& mi2 = mi[*cit2];
	  // random polynomial part must be identical to contribute to variance
	  // (else orthogonality drops term).  Note that it is not necessary to
	  // collapse terms with the same random basis subset, since cross term
	  // in (a+b)(a+b) = a^2+2ab+b^2 gets included.  If terms were collapsed
	  // (following eval of non-random portions), the nested loop could be
	  // replaced with a single loop to evaluate (a+b)^2.
	  if (data_rep->match_random_key(mi1, mi2))
	    covar += coeff_norm_poly * exp_coeffs_2[i2] * 
	      data_rep->multivariate_polynomial(x, mi2, nrand_ind);
	}
      }
    }
  }
  return covar;
}


Real RegressOrthogPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  RegressOrthogPolyApproximation* ropa_2
    = (RegressOrthogPolyApproximation*)poly_approx_2;
  if ( ( sparseIndIter == sparseIndices.end() ||
	 sparseIndIter->second.empty()) &&
       ( ropa_2->sparseIndIter == ropa_2->sparseIndices.end() ||
	 ropa_2->sparseIndIter->second.empty()) )
    return OrthogPolyApproximation::covariance(x, poly_approx_2);

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  bool same = (this == ropa_2), all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !ropa_2->expansionCoeffFlag )) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "RegressOrthogPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if ( same && all_mode && (compVarIter->second & 1) &&
       data_rep->match_nonrandom_vars(x, xPrevVar[key]) )
    return expMomentsIter->second[1];

  Real covar = covariance(x, data_rep->multi_index(), expCoeffsIter->second,
			  sparseIndIter->second, ropa_2->expCoeffsIter->second,
			  ropa_2->sparseIndIter->second);
  if (same && all_mode) {
    expMomentsIter->second[1] = covar;
    compVarIter->second |= 1;  xPrevVar[key] = x;
  }
  return covar;
}


Real RegressOrthogPolyApproximation::
combined_covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  RegressOrthogPolyApproximation* ropa_2
    = (RegressOrthogPolyApproximation*)poly_approx_2;
  if (combinedSparseIndices.empty() && ropa_2->combinedSparseIndices.empty())
    return OrthogPolyApproximation::covariance(x, poly_approx_2);

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  bool same = (this == ropa_2), all_mode = !data_rep->nonRandomIndices.empty();
  const UShortArray& key = data_rep->activeKey;

  // Error check for required data
  if ( !expansionCoeffFlag ||
       ( !same && !ropa_2->expansionCoeffFlag )) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "RegressOrthogPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }

  if ( same && all_mode && (compVarIter->second & 1) &&
       data_rep->match_nonrandom_vars(x, xPrevVar[key]) )
    return expMomentsIter->second[1];

  Real covar = covariance(x, data_rep->combinedMultiIndex, combinedExpCoeffs,
			  combinedSparseIndices, ropa_2->combinedExpCoeffs,
			  ropa_2->combinedSparseIndices);
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
const RealVector& RegressOrthogPolyApproximation::variance_gradient()
{
  if (sparseIndIter == sparseIndices.end() || sparseIndIter->second.empty())
    return OrthogPolyApproximation::variance_gradient();

  // d/ds \sigma^2_R = Sum_{j=1}^P <Psi^2_j> d/ds \alpha^2_j
  //                 = 2 Sum_{j=1}^P \alpha_j <dR/ds, Psi_j>

  // Error check for required data
  if (!expansionCoeffFlag || !expansionCoeffGradFlag) {
    PCerr << "Error: insufficient expansion coefficient data in RegressOrthog"
	  << "PolyApproximation::variance_gradient()." << std::endl;
    abort_handler(-1);
  }

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  bool std_mode = data_rep->nonRandomIndices.empty();
  if (std_mode && (compVarIter->second & 2))
    return momentGradsIter->second[1];

  const RealVector& exp_coeffs      = expCoeffsIter->second;
  const RealMatrix& exp_coeff_grads = expCoeffGradsIter->second;
  size_t i, j, num_deriv_v = exp_coeff_grads.numRows();
  RealVector& var_grad = momentGradsIter->second[1];
  if (var_grad.length() != num_deriv_v)
    var_grad.sizeUninitialized(num_deriv_v);
  var_grad = 0.;
  const UShort2DArray& mi = data_rep->multi_index();
  const SizetSet& sparse_ind = sparseIndIter->second;  StSCIter cit;
  for (i=1, cit=++sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit) {
    Real term_i = 2. * exp_coeffs[i] * data_rep->norm_squared(mi[*cit]);
    for (j=0; j<num_deriv_v; ++j)
      var_grad[j] += term_i * exp_coeff_grads[i][j];
  }
  if (std_mode) compVarIter->second |=  2;
  else          compVarIter->second &= ~2;// deactivate 2-bit: protect mixed use
  return var_grad;
}


/** In this function, a subset of the expansion variables are random
    variables and any augmented design/state variables (i.e., not
    inserted as random variable distribution parameters) are included
    in the expansion.  This function must handle the mixed case, where
    some design/state variables are augmented (and are part of the
    expansion) and some are inserted (derivatives are obtained from
    expansionCoeffGrads). */
const RealVector& RegressOrthogPolyApproximation::
variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  if (sparseIndIter == sparseIndices.end() || sparseIndIter->second.empty())
    return OrthogPolyApproximation::variance_gradient(x, dvv);

  // Error check for required data
  if (!expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "RegressOrthogPolyApproximation::variance_gradient()" << std::endl;
    abort_handler(-1);
  }

  // if already computed, return previous result
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const SizetList& nrand_ind = data_rep->nonRandomIndices;
  bool all_mode = !nrand_ind.empty();
  const UShortArray& key = data_rep->activeKey;
  if ( all_mode && (compVarIter->second & 2) &&
       data_rep->match_nonrandom_vars(x, xPrevVarGrad[key]) )
    // && dvv == dvvPrev)
    return momentGradsIter->second[1];

  size_t i, j, k, deriv_index, num_deriv_v = dvv.size(),
    cntr = 0; // insertions carried in order within expansionCoeffGrads
  RealVector& var_grad = momentGradsIter->second[1];
  if (var_grad.length() != num_deriv_v)
    var_grad.sizeUninitialized(num_deriv_v);
  var_grad = 0.;

  const UShort2DArray& mi           = data_rep->multi_index();
  const RealVector& exp_coeffs      = expCoeffsIter->second;
  const RealMatrix& exp_coeff_grads = expCoeffGradsIter->second;
  const SizetSet&   sparse_ind      = sparseIndIter->second;  StSCIter cit;
  const SizetList&    rand_ind      = data_rep->randomIndices;
  Real norm_sq_j, poly_j, poly_grad_j, norm_poly_j, coeff_j, coeff_grad_j;
  StSCIter jit, kit;
  for (i=0; i<num_deriv_v; ++i) {
    deriv_index = dvv[i] - 1; // OK since we are in an "All" view
    bool random = data_rep->randomVarsKey[deriv_index];
    if (random && !expansionCoeffGradFlag){
      PCerr << "Error: expansion coefficient gradients not defined in Regress"
	    << "OrthogPolyApproximation::variance_gradient()." << std::endl;
      abort_handler(-1);
    }
    for (j=1, jit=++sparse_ind.begin(); jit!=sparse_ind.end(); ++j, ++jit) {
      const UShortArray& mi_j = mi[*jit];
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
	for (k=1, kit=++sparse_ind.begin(); kit!=sparse_ind.end(); ++k, ++kit) {
	  // random part of polynomial must be identical to contribute to
	  // variance (else orthogonality drops term)
	  const UShortArray& mi_k = mi[*kit];
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
		( coeff_j * exp_coeff_grads[k][cntr] +
		  exp_coeffs[k] * coeff_grad_j ) *
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
  else compVarIter->second &= ~2;//deactivate 2-bit: protect mixed use
  return var_grad;
}


void RegressOrthogPolyApproximation::set_fault_info()
{
  size_t constr_eqns, num_data_pts_fn,
    num_data_pts_grad, total_eqns, num_surr_data_pts;
  bool under_determined = false, reuse_solver_data,
    anchor_fn = false, anchor_grad = false;

  // compute order of data contained within modSurrData
  short data_order = (expansionCoeffFlag) ? 1 : 0;
  if (modSurrData.num_gradient_variables())  data_order |= 2;
  //if (modSurrData.num_hessian_variables()) data_order |= 4;

  // verify support for basisConfigOptions.useDerivs, which indicates usage of
  // derivative data with respect to expansion variables (aleatory or combined)
  // within the expansion coefficient solution process, which must be
  // distinguished from usage of derivative data with respect to non-expansion
  // variables (the expansionCoeffGradFlag case).
  bool config_err = false;
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  if (data_rep->basisConfigOptions.useDerivs) {
    if (!(data_order & 2)) {
      PCerr << "Error: useDerivs configuration option lacks data support in "
	    << "RegressOrthogPolyApproximation::regression()" << std::endl;
      config_err = true;
    }
    if (expansionCoeffGradFlag) {
      PCerr << "Error: useDerivs configuration option conflicts with gradient "
	    << "expansion request in RegressOrthogPolyApproximation::"
	    << "regression()" << std::endl;
      config_err = true;
    }
    //if (data_order & 4)
    //  PCerr << "Warning: useDerivs configuration option does not yet support "
    //	      << "Hessian data in RegressOrthogPolyApproximation::regression()"
    //	      << std::endl;
  }
  if (config_err)
    abort_handler(-1);

  // compute data counts
  const SizetShortMap& failed_resp_data = modSurrData.failed_response_data();
  size_t num_failed_surr_fn = 0, num_failed_surr_grad = 0,
    num_v = sharedDataRep->numVars;
  SizetShortMap::const_iterator fit; bool faults_differ = false;
  for (fit=failed_resp_data.begin(); fit!=failed_resp_data.end(); ++fit) {
    short fail_bits = fit->second;
    if (fail_bits & 1) ++num_failed_surr_fn;
    if (fail_bits & 2) ++num_failed_surr_grad;
    // if failure omissions are not consistent, manage differing Psi matrices
    if ( (fail_bits & data_order) != data_order ) faults_differ = true;
  }
  num_surr_data_pts = modSurrData.points();
  num_data_pts_fn   = num_surr_data_pts - num_failed_surr_fn;
  num_data_pts_grad = num_surr_data_pts - num_failed_surr_grad;
  if (modSurrData.anchor()) {
    short failed_anchor_data = modSurrData.failed_anchor_data();
    if ((data_order & 1) && !(failed_anchor_data & 1)) anchor_fn   = true;
    if ((data_order & 2) && !(failed_anchor_data & 2)) anchor_grad = true;
  }

  // detect underdetermined system of equations (following fault omissions)
  // for either expansion coeffs or coeff grads (switch logic treats together)
  reuse_solver_data
    = (expansionCoeffFlag && expansionCoeffGradFlag && !faults_differ);
  constr_eqns = 0;
  if (expansionCoeffFlag) {
    constr_eqns = num_data_pts_fn;
    total_eqns = (data_rep->basisConfigOptions.useDerivs) ?
      constr_eqns + num_data_pts_grad * num_v : constr_eqns;
    if (total_eqns < data_rep->expansion_terms()) // candidate expansion size
      under_determined = true;
  }
  if (expansionCoeffGradFlag) {
    total_eqns = num_data_pts_grad;
    if (total_eqns < data_rep->expansion_terms()) // candidate expansion size
      under_determined = true;
  }

  faultInfo.set_info( constr_eqns, anchor_fn, anchor_grad,
		      under_determined, num_data_pts_fn, num_data_pts_grad,
		      reuse_solver_data, total_eqns, num_surr_data_pts,
		      num_v, data_rep->basisConfigOptions.useDerivs,
		      modSurrData.num_derivative_variables() );
}


void RegressOrthogPolyApproximation::
build_linear_system( RealMatrix &A, const UShort2DArray& multi_index)
{
  size_t i, j, a_cntr = 0, num_surr_data_pts = modSurrData.points(),
    num_v = sharedDataRep->numVars,  a_grad_cntr = 0;
  int num_rows_A, num_cols_A = multi_index.size(), // candidate expansion size
    num_data_pts_fn   = num_surr_data_pts, // failed data is removed downstream
    num_data_pts_grad = num_surr_data_pts; // failed data is removed downstream
  bool add_val, add_grad;
  const SDVArray& sdv_array = modSurrData.variables_data();

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  if (expansionCoeffFlag) {
    // matrix/vector sizing
    num_rows_A = (data_rep->basisConfigOptions.useDerivs) ?
      num_data_pts_fn + num_data_pts_grad * num_v : num_data_pts_fn;

    A.shapeUninitialized(num_rows_A, num_cols_A);
    Real *A_matrix = A.values();
    // The "A" matrix is a contiguous block of memory packed in column-major
    // ordering as required by F77 for the GELSS subroutine from LAPACK.  For
    // example, the 6 elements of A(2,3) are stored in the order A(1,1),
    // A(2,1), A(1,2), A(2,2), A(1,3), A(2,3).
    for (i=0; i<num_cols_A; ++i) {
      a_cntr = num_rows_A*i;
      a_grad_cntr = a_cntr + num_data_pts_fn;
      const UShortArray& mi_i = multi_index[i];
      for (j=0; j<num_surr_data_pts; ++j) {
	add_val = true; add_grad = data_rep->basisConfigOptions.useDerivs;
	data_rep->pack_polynomial_data(sdv_array[j].continuous_variables(),
				       mi_i, add_val, A_matrix, a_cntr,
				       add_grad, A_matrix, a_grad_cntr);
      }
    }
  }
  else if (expansionCoeffGradFlag) {
    num_rows_A = num_data_pts_grad;
    A.shapeUninitialized(num_rows_A, num_cols_A);
    Real *A_matrix = A.values();

    // repack "A" matrix with different Psi omissions
    a_cntr = 0;
    for (i=0; i<num_cols_A; ++i) {
      const UShortArray& mi_i = multi_index[i];
      for (j=0; j<num_surr_data_pts; ++j) {
	//add_val = false; add_grad = true;
	//if (add_grad) {
	  A_matrix[a_cntr] = data_rep->
	    multivariate_polynomial(sdv_array[j].continuous_variables(), mi_i);
	  ++a_cntr;
	//}
      }
    }
  }
}


void RegressOrthogPolyApproximation::
build_linear_system( RealMatrix &A, RealMatrix &B,
		     const UShort2DArray& multi_index)
{
  size_t i, j, b_cntr = 0, num_surr_data_pts = modSurrData.points(),
    num_deriv_v = modSurrData.num_derivative_variables(),
    num_v = sharedDataRep->numVars, b_grad_cntr = 0;
  int num_rows_B, num_coeff_rhs, num_grad_rhs = num_deriv_v, num_rhs,
    num_data_pts_fn   = num_surr_data_pts, // failed data is removed downstream
    num_data_pts_grad = num_surr_data_pts; // failed data is removed downstream
  bool add_val, add_grad;
  const SDRArray& sdr_array = modSurrData.response_data();

  // populate A
  build_linear_system(A, multi_index);

  if (expansionCoeffFlag) {
    
    // matrix/vector sizing
    SharedRegressOrthogPolyApproxData* data_rep
      = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
    num_rows_B = (data_rep->basisConfigOptions.useDerivs) ?
      num_data_pts_fn + num_data_pts_grad * num_v : num_data_pts_fn;
    num_coeff_rhs = 1;
    num_rhs = (expansionCoeffGradFlag) ?
      num_coeff_rhs + num_grad_rhs : num_coeff_rhs;

    B.shapeUninitialized(num_rows_B, num_rhs);
    Real *b_vectors = B.values();
    
    // response data (values/gradients) define the multiple RHS which are
    // matched in the LS soln.  b_vectors is num_data_pts (rows) x num_rhs
    // (cols), arranged in column-major order.
    b_cntr = 0; b_grad_cntr = num_data_pts_fn;
    for (i=0; i<num_surr_data_pts; ++i) {
      add_val = true; add_grad = data_rep->basisConfigOptions.useDerivs;
      data_rep->pack_response_data(sdr_array[i], add_val, b_vectors, b_cntr,
				   add_grad, b_vectors, b_grad_cntr);
    }
  }

  if (expansionCoeffGradFlag) {

    if (!expansionCoeffFlag) {
      num_rows_B = num_data_pts_grad; num_rhs = num_grad_rhs; num_coeff_rhs = 0;
      B.shapeUninitialized(num_rows_B, num_rhs);
    }
    
    // response data (values/gradients) define the multiple RHS which are
    // matched in the LS soln.  b_vectors is num_data_pts (rows) x num_rhs
    // (cols), arranged in column-major order.
    Real *b_vectors = B.values();
    b_cntr = 0;
    for (i=0; i<num_surr_data_pts; ++i) {
      //add_val = false; add_grad = true;
      //if (add_grad) {
	const RealVector& resp_grad = sdr_array[i].response_gradient();
	for (j=0; j<num_grad_rhs; ++j) // i-th point, j-th grad component
	  b_vectors[(j+num_coeff_rhs)*num_data_pts_grad+b_cntr] = resp_grad[j];
	++b_cntr;
      //}
    }
  }
}


void RegressOrthogPolyApproximation::
build_linear_system( RealMatrix &A, RealMatrix &B, RealMatrix &points,
		     const UShort2DArray& multi_index)
{
  // populate A and B
  build_linear_system(A, B, multi_index);

  // populate points
  size_t i, j, num_surr_data_pts = modSurrData.points(),
    num_v = sharedDataRep->numVars;
  const SDVArray& sdv_array = modSurrData.variables_data();
  points.shapeUninitialized( num_v, num_surr_data_pts );
  for (i=0; i<num_surr_data_pts; ++i) {
    const RealVector& c_vars = sdv_array[i].continuous_variables();
    for (j=0; j<num_v; ++j)
      points(j,i) = c_vars[j];
  }
  //points.print(std::cout);
}


void RegressOrthogPolyApproximation::
augment_linear_system( const RealVectorArray& samples, RealMatrix &A,
		       const UShort2DArray& multi_index)
{
  size_t i, j, a_cntr = 0, num_v = sharedDataRep->numVars, a_grad_cntr = 0;
  int num_samp = samples.size(), orig_rows_A = A.numRows(), num_rows_A,
    num_cols_A = multi_index.size(); // candidate expansion size
  bool add_val, add_grad;

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  if (expansionCoeffFlag) {
    // matrix/vector sizing
    num_rows_A += (data_rep->basisConfigOptions.useDerivs) ?
      orig_rows_A + num_samp*(1+num_v) : orig_rows_A + num_samp;

    A.reshape(num_rows_A, num_cols_A);
    Real *A_matrix = A.values();
    // The "A" matrix is a contiguous block of memory packed in column-major
    // ordering as required by F77 for the GELSS subroutine from LAPACK.  For
    // example, the 6 elements of A(2,3) are stored in the order A(1,1),
    // A(2,1), A(1,2), A(2,2), A(1,3), A(2,3).
    for (i=0; i<num_cols_A; ++i) {
      a_cntr = orig_rows_A + num_rows_A*i; // offset by original rows
      a_grad_cntr = a_cntr + num_samp;
      const UShortArray& mi_i = multi_index[i];
      for (j=0; j<num_samp; ++j) {
	add_val = true; add_grad = data_rep->basisConfigOptions.useDerivs;
	data_rep->pack_polynomial_data(samples[j], mi_i, add_val, A_matrix,
				       a_cntr, add_grad, A_matrix, a_grad_cntr);
      }
    }
  }
  else if (expansionCoeffGradFlag) {
    num_rows_A = orig_rows_A + num_samp;
    A.reshape(num_rows_A, num_cols_A);
    Real *A_matrix = A.values();

    // repack "A" matrix with different Psi omissions
    a_cntr = 0;
    for (i=0; i<num_cols_A; ++i) {
      const UShortArray& mi_i = multi_index[i];
      a_cntr += orig_rows_A; // TO DO: verify offset
      for (j=0; j<num_samp; ++j) {
	//add_val = false; add_grad = true;
	//if (add_grad) {
	  A_matrix[a_cntr] = data_rep->multivariate_polynomial(samples[j],mi_i);
	  ++a_cntr;
        //}
      }
    }
  }
}


Real RegressOrthogPolyApproximation::
run_cross_validation_solver(const UShort2DArray& multi_index,
			    RealVector& best_exp_coeffs,
			    SizetSet& best_sparse_indices)
{
  // TO DO: employ this fn as component within primary CV context

  RealMatrix A, B;
  build_linear_system( A, B, multi_index );

  int num_data_pts_fn = modSurrData.points();

  RealVector b( Teuchos::Copy, B.values(), B.numRows() );
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const UShortArray& approx_order = data_rep->expansion_order();
  int num_rhs = B.numCols(), num_dims( approx_order.size() );

  Real best_score = std::numeric_limits<Real>::max();
  int best_basis_parameters_index = 0;
  int num_build_points = num_data_pts_fn;
  bool use_gradients = false;
  if ( A.numRows() > num_data_pts_fn ) use_gradients = true;
  MultipleSolutionLinearModelCrossValidationIterator cv_iterator;
  int cv_seed = data_rep->regressConfigOptions.randomSeed; 
  cv_iterator.set_seed( cv_seed );
  cv_iterator.set_fault_data( faultInfo,
			      modSurrData.failed_response_data() );
  
  data_rep->CSTool.set_linear_solver( CSOpts );
  LinearSolver_ptr linear_solver = data_rep->CSTool.get_linear_solver();
  cv_iterator.set_solver( linear_solver );

  // default to 10 folds and revert to leave-one-out for fewer than 10
  // data points
  int num_folds = std::min(10, num_data_pts_fn);
  //int num_folds = num_data_pts_fn; //HACK
  int max_num_pts_per_fold = num_data_pts_fn / num_folds;
  if ( num_data_pts_fn % num_folds ) ++max_num_pts_per_fold;
  if ( CSOpts.solver != ORTHOG_MATCH_PURSUIT   &&
       CSOpts.solver != LASSO_REGRESSION       &&
       CSOpts.solver != LEAST_ANGLE_REGRESSION &&
       ( num_data_pts_fn - max_num_pts_per_fold < A.numCols() ) &&
       ( num_data_pts_fn >= A.numCols() ) ) {
    // use one at a time cross validation
    PCout << "The linear system will switch from over-determined to " << 
      " under-determined when K = 10 cross validation is employed. " <<
      " Switching to one at a time cross validation to avoid this\n";
    num_folds = num_data_pts_fn;

    if ( num_data_pts_fn == A.numCols() )
      PCout << "Warning: The linear system will switch from exactly " <<
	" determined to under-determined when leave one out cross " << 
	" validation is employed\n.";
  }

  cv_iterator.set_max_num_unique_tolerances( 100 );
  cv_iterator.set_num_folds( num_folds );
  cv_iterator.set_num_points( num_build_points );
  if ( use_gradients )
    cv_iterator.set_num_equations_per_point( sharedDataRep->numVars + 1 );
  else 
    cv_iterator.set_num_equations_per_point( 1 );

  Real score = cv_iterator.run_cross_validation( A, b );
  Real best_tolerance = cv_iterator.get_best_residual_tolerance();

  if ( data_rep->expConfigOptions.outputLevel >= NORMAL_OUTPUT )
    PCout << "Cross validation score: " << score << "\nBest tolerance chosen "
	  << "by cross validation: " << best_tolerance << "\n";

  // CV is complete, now compute final solution with all data points:
  IntVector index_mapping;
  RealMatrix points_dummy;
  remove_faulty_data( A, b, points_dummy, index_mapping,
		      faultInfo, modSurrData.failed_response_data() );
  int num_rows_V = A.numRows(), num_cols_V = A.numCols();
  faultInfo.under_determined = num_rows_V < num_cols_V;
  select_solver(false); // CV no longer active for final soln
  PCout << "Applying regression to compute " << num_cols_V
	<< " chaos coefficients using " << num_rows_V << " equations.\n";
  RealMatrix solutions, metrics;
  linear_solver->set_residual_tolerance( best_tolerance );
  linear_solver->solve( A, b, solutions, metrics );

  // In current usage, global (shared multiIndex/sobolIndexMap) and local
  // (sparseSobolIndexMap) bookkeeping are *not* updated since higher level
  // logic (select_best_*()) determines acceptance of candidate solutions.
  best_sparse_indices.clear();
  Real* dense_coeffs = solutions[solutions.numCols()-1];
  update_sparse_indices(dense_coeffs, A.numCols(), best_sparse_indices);
  update_sparse_coeffs(dense_coeffs, best_exp_coeffs, best_sparse_indices);  

  // JDJ: We always want to only keep sparse indices
  //else { // least sq case does not require sparse indices
  //  best_sparse_indices.clear();
  //  copy_data(solutions[solutions.numCols()-1], A.numCols(), best_exp_coeffs);
  //}

  return score;
}


Real RegressOrthogPolyApproximation::run_cross_validation_expansion()
{
  RealMatrix A, B;
  build_linear_system( A, B );

  int num_data_pts_fn = modSurrData.points();

  RealVector b( Teuchos::Copy, B.values(), B.numRows() );
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const UShortArray& approx_order = data_rep->expansion_order();
  int num_rhs = B.numCols(), num_dims( approx_order.size() );
  // Do cross validation for varing polynomial orders up to 
  // a maximum order defined by approxOrder[0]
  unsigned short ao0 = approx_order[0],
    min_order = (data_rep->regressConfigOptions.crossValidNoiseOnly)
    ? ao0 : std::min((unsigned short)1, ao0);
  //if ( min_order > ao0 ) min_order = ao0;

  Real best_score = std::numeric_limits<Real>::max(), best_tolerance = 0.;
  int best_basis_parameters_index = 0;
  int num_build_points = num_data_pts_fn;
  bool use_gradients = false;
  if ( A.numRows() > num_data_pts_fn ) use_gradients = true;
  MultipleSolutionLinearModelCrossValidationIterator cv_iterator;
  int cv_seed = data_rep->regressConfigOptions.randomSeed; 
  cv_iterator.set_seed( cv_seed );
  cv_iterator.set_fault_data( faultInfo,
			      modSurrData.failed_response_data() );
  
  RealVector basis_scores;
  basis_scores.sizeUninitialized( ao0 - min_order + 1 );
  basis_scores = std::numeric_limits<double>::max();

  data_rep->CSTool.set_linear_solver( CSOpts );
  LinearSolver_ptr linear_solver = data_rep->CSTool.get_linear_solver();
  cv_iterator.set_solver( linear_solver );

  int i = 0;
  bestApproxOrder.size( 1 ); 
  for ( int order = min_order; order <= ao0; order++ )
    {
      if (data_rep->expConfigOptions.outputLevel >= QUIET_OUTPUT)
	PCout << "Testing PCE order " << order << std::endl;
      int num_basis_terms = util::nchoosek( num_dims + order, order );
      RealMatrix vandermonde_submatrix( Teuchos::View, 
					A,
					A.numRows(),
					num_basis_terms, 0, 0 );


      int num_folds = 10;
      int max_num_pts_per_fold = num_data_pts_fn / num_folds;
      if ( num_data_pts_fn % num_folds != 0 ); max_num_pts_per_fold++;
      if ( CSOpts.solver != ORTHOG_MATCH_PURSUIT   &&
	   CSOpts.solver != LASSO_REGRESSION       &&
	   CSOpts.solver != LEAST_ANGLE_REGRESSION &&
	   ( num_data_pts_fn - max_num_pts_per_fold < vandermonde_submatrix.numCols() ) &&
	   ( num_data_pts_fn >= vandermonde_submatrix.numCols() ) ) {
	  // use one at a time cross validation
	PCout << "The linear system will switch from over-determined to " << 
	  " under-determined when K = 10 cross validation is employed. " <<
	  " Switching to one at a time cross validation to avoid this\n";
	num_folds = num_data_pts_fn;

	if ( num_data_pts_fn == vandermonde_submatrix.numCols() )
	  PCout << "Warning: The linear system will switch from exactly " <<
	    " determined to under-determined when leave one out cross " << 
	    " validation is employed\n.";
	}

      cv_iterator.set_max_num_unique_tolerances( 100 );
      cv_iterator.set_num_folds( num_folds );
      cv_iterator.set_num_points( num_build_points );
      if ( use_gradients )
	cv_iterator.set_num_equations_per_point( sharedDataRep->numVars + 1 );
      else 
	cv_iterator.set_num_equations_per_point( 1 );


      Real score = cv_iterator.run_cross_validation( vandermonde_submatrix, b );

      if ( score < best_score )
	{
	  best_score = score;
	  best_tolerance = cv_iterator.get_best_residual_tolerance();
	  best_basis_parameters_index = i;
	  bestApproxOrder[0] = order;
	}
      basis_scores[i] = score;

      //if ( score >= best_score && i - best_basis_parameters_index >= 2 )
      //break;

      if ( data_rep->expConfigOptions.outputLevel >= QUIET_OUTPUT )
	{
	  PCout << "Cross validation error for degree ";
	  PCout << order << ": " << score << "\n";
	}
      i++;
    }

  int num_basis_terms = util::nchoosek( num_dims + bestApproxOrder[0],
					      bestApproxOrder[0] );
  if (data_rep->expConfigOptions.outputLevel >= QUIET_OUTPUT)
    PCout << "Best approximation order: " << bestApproxOrder[0]
	  << "\nBest cross validation error: " << best_score << "\n";
  // set CSOpts so that best PCE can be built. We are assuming num_rhs=1
  RealMatrix vandermonde_submatrix( Teuchos::View, A, A.numRows(),
				    num_basis_terms, 0, 0 );
  IntVector index_mapping;
  RealMatrix points_dummy;
  remove_faulty_data( vandermonde_submatrix, b, points_dummy, index_mapping,
		      faultInfo, modSurrData.failed_response_data() );
  int num_rows_V = vandermonde_submatrix.numRows(),
      num_cols_V = vandermonde_submatrix.numCols();
  faultInfo.under_determined = num_rows_V < num_cols_V;
  PCout << "Applying regression to compute " << num_cols_V
	<< " chaos coefficients using " << num_rows_V << " equations.\n";
  select_solver(false); // CV no longer active for final soln
  RealMatrix solutions, metrics;
  linear_solver->set_residual_tolerance( best_tolerance );
  linear_solver->solve( vandermonde_submatrix, b, solutions, metrics );

  int last_index = solutions.numCols() - 1;
  // exploit CS sparsity for OMP,LASSO,LARS regardless of data size and for
  // BP,BPDN if under-determined (these revert to least sq for over determined)
  if (sparseSoln)
    update_sparse(solutions[last_index], num_basis_terms);
  else {
    copy_data(solutions[last_index], num_basis_terms, expCoeffsIter->second);
    // if best expansion order is less than maximum candidate, define
    // sparseIndices to define active subset within data_rep->multi_index().
    // Note that this requires care in cross-expansion evaluations such as
    // off-diagonal covariance.
    SizetSet& sparse_ind = sparseIndIter->second;
    sparse_ind.clear();
    if (num_basis_terms < data_rep->expansion_terms()) { // candidate exp size
      inflate(sparse_ind, num_basis_terms); // define from leading subset
      update_sparse_sobol(sparse_ind, data_rep->multi_index(),
			  data_rep->sobolIndexMap);
    }
  }

  return best_score;
}


void RegressOrthogPolyApproximation::
compressed_sensing( RealMatrix &A, RealMatrix &B )
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;

  CSOpts.standardizeInputs = false;// false essential when using derivatives

  CompressedSensingOptionsList opts_list;
  RealMatrixArray solutions;
  data_rep->CSTool.solve( A, B, solutions, CSOpts, opts_list );

  // update bookkeeping for sparse solutions
  bool multiple_rhs = (expansionCoeffFlag && expansionCoeffGradFlag);
  int num_expansion_terms = data_rep->expansion_terms();
  if ( expansionCoeffFlag && !multiple_rhs ) {
    if (sparseSoln) // exploit CS sparsity
      update_sparse(solutions[0][0], num_expansion_terms);
    else {                          // retain full solution
      copy_data(solutions[0][0], num_expansion_terms, expCoeffsIter->second);
      if (sparseIndIter != sparseIndices.end()) {
	//sparseIndices.erase(sparseIndIter);
	//sparseIndIter = sparseIndices.end();
	sparseIndIter->second.clear();
      }
    }
  }
  else {
    int i, j, num_grad_rhs = modSurrData.num_derivative_variables(),
      num_coeff_rhs = ( !multiple_rhs && expansionCoeffGradFlag ) ? 0 : 1;
    if (sparseSoln) { // exploit CS sparsity
      // overlay sparse solutions into an aggregated set of sparse indices
      SizetSet& sparse_ind = sparseIndIter->second;
      sparse_ind.clear();
      if (multiple_rhs)
	update_sparse_indices(solutions[0][0], num_expansion_terms,
			      sparse_ind);
      for (i=0; i<num_grad_rhs; ++i)
	update_sparse_indices(solutions[i+num_coeff_rhs][0],
			      num_expansion_terms, sparse_ind);
      // update expansion{Coeffs,CoeffGrads} from sparse indices
      if (multiple_rhs)
	update_sparse_coeffs(solutions[0][0], expCoeffsIter->second,
			     sparse_ind);
      for (i=0; i<num_grad_rhs; ++i)
	update_sparse_coeff_grads(solutions[i+num_coeff_rhs][0], i,
				  expCoeffGradsIter->second, sparse_ind);
      // update sobol index bookkeeping
      update_sparse_sobol(sparse_ind, data_rep->multi_index(),
			  data_rep->sobolIndexMap);
    }
    else { // retain original multiIndex layout
      if (multiple_rhs)
	copy_data(solutions[0][0], num_expansion_terms, expCoeffsIter->second);
      RealMatrix& exp_coeff_grads = expCoeffGradsIter->second;
      for (i=0; i<num_grad_rhs; ++i) {
	Real* dense_coeffs = solutions[i+num_coeff_rhs][0];
	for (j=0; j<num_expansion_terms; ++j)
	  exp_coeff_grads(i,j) = dense_coeffs[j];
      }
      if (sparseIndIter != sparseIndices.end()) {
	//sparseIndices.erase(sparseIndIter);
	//sparseIndIter = sparseIndices.end();
	sparseIndIter->second.clear();
      }
    }
  }
}


void RegressOrthogPolyApproximation::
least_interpolation( RealMatrix &pts, RealMatrix &vals )
{
#ifdef DEBUG
  if ( pts.numCols() != vals.numRows() ) {
    std::string msg
      = "least_interpolation() dimensions of pts and vals are inconsistent";
    throw( std::runtime_error( msg ) );
  }
#endif

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  SizetSet& sparse_ind = sparseIndIter->second;
  // if no sim faults on subsequent QoI and previous QoI interp size matches,
  // then reuse previous factorization.  Detecting non-null fault sets that are
  // consistent is more complicated and would require changes to the Dakota::
  // ApproximationInterface design to propagate fault deltas among the QoI.
  // Note 1: size logic is more general than currently necessary since OLI does
  //         not currently support deriv enhancement.
  // Note 2: multiIndex size check captures first QoI pass as well as reentrancy
  //         (changes in points sets for OUU and mixed UQ) due to clear() in
  //         SharedRegressOrthogPolyApproxData::allocate_data().
  bool faults = !modSurrData.failed_response_data().empty(),
    inconsistent_prev = ( data_rep->multi_index().empty() ||
      modSurrData.active_response_size() != data_rep->pivotHistory.numRows() );
  if (faults || inconsistent_prev) {
    // compute the least factorization that interpolates this data set
    UShort2DArray local_multi_index; IntVector k;
    least_factorization( pts, local_multi_index, data_rep->lowerFactor,
			 data_rep->upperFactor, data_rep->pivotHistory,
			 data_rep->pivotVect, k );
    // update approxOrder (for use in, e.g., combine_coefficients())
    int last_index = k.length() - 1, new_order = k[last_index];
    data_rep->expansion_order(new_order, true);
    // update sparseIndices and shared multiIndex from local_multi_index
    // Note: sparseIndices definition does not involve exp coeffs in this case
    size_t local_mi_ref;
    data_rep->append_leading_multi_index(local_multi_index,
					 data_rep->multi_index(),
					 sparse_ind, local_mi_ref);
    // update shared sobolIndexMap from local_multi_index
    data_rep->update_component_sobol(local_multi_index);
  }
  else // define sparseIndices for this QoI (sparseSobolIndexMap updated below)
    inflate(sparse_ind, data_rep->expansion_terms());

  // define sparseSobolIndexMap from sparseIndices, shared multiIndex,
  // and shared sobolIndexMap
  update_sparse_sobol(sparse_ind, data_rep->multi_index(),
		      data_rep->sobolIndexMap);

  RealMatrix coefficients;
  transform_least_interpolant( data_rep->lowerFactor,  data_rep->upperFactor,
			       data_rep->pivotHistory, data_rep->pivotVect,
			       vals );
}


void RegressOrthogPolyApproximation::
transform_least_interpolant( RealMatrix &L, RealMatrix &U, RealMatrix &H,
			     IntVector &p,  RealMatrix &vals )
{
  int num_pts = vals.numRows(), num_qoi = vals.numCols();
  
  RealMatrix LU_inv;
  IntVector dummy;
  util::lu_inverse( L, U, dummy, LU_inv );

  RealMatrix V_inv( H.numCols(), num_pts );
  V_inv.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, H, LU_inv, 0.0 );

  IntVector p_col;
  util::argsort( p, p_col );
  util::permute_matrix_columns( V_inv, p_col );

  RealMatrix coefficients;
  coefficients.shapeUninitialized( V_inv.numRows(), num_qoi );
  coefficients.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, V_inv, 
			 vals, 0.0 );

  // multiIndex should be consistent across QoI vector
  copy_data(coefficients.values(), (int)sparseIndIter->second.size(),
	    expCoeffsIter->second);
}


void RegressOrthogPolyApproximation::
least_factorization( RealMatrix &pts, UShort2DArray &basis_indices,
		     RealMatrix &l, RealMatrix &u, RealMatrix &H,
		     IntVector &p, IntVector &k )
{
  int num_vars = pts.numRows(), num_pts = pts.numCols();

  util::eye( num_pts, l );
  util::eye( num_pts, u );

  util::range( p, 0, num_pts, 1 );

  //This is just a guess: this vector could be much larger, or much smaller
  RealVector v( 1000 );
  int v_index = 0;

  // Current polynomial degree
  int k_counter = 0;
  // k[q] gives the degree used to eliminate the q'th point
  k.size( num_pts );

  // The current LU row to factor out:
  int lu_row = 0;

  UShort2DArray internal_basis_indices;

  // Current degree is k_counter, and we iterate on this
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  while ( lu_row < num_pts )
    {
      UShort2DArray new_indices;
      if ( basis_indices.size() == 0 )
	{
	  data_rep->total_order_multi_index( (unsigned short)k_counter,
					     (size_t)num_vars, new_indices );
	  int num_indices = internal_basis_indices.size();
	}
      else
	{
	  for ( int i = 0; i < (int)basis_indices.size(); i++ )
	    {
	      if ( l1_norm( basis_indices[i] ) == k_counter )
		{
		  new_indices.push_back( basis_indices[i] );
		}
	    }
	  // If the basis_indices set is very sparse then not all degrees may 
	  // be represented in basis_indices. Thus increase degree counter
	  if ( new_indices.size() == 0 )
	    k_counter++;
	  if ( ( basis_indices.size() == internal_basis_indices.size() ) && 
	       ( new_indices.size() == 0 ) )
	    {
	      std::string msg = "least_factorization() the basis indices ";
	      msg += "specified were insufficient to interpolate the data";
	      throw( std::runtime_error( msg ) ); 
	    }
	}

      if ( (  new_indices.size() > 0 ) )
	{	
      
	  internal_basis_indices.insert( internal_basis_indices.end(), 
					 new_indices.begin(), 
					 new_indices.end() );
      
	  int current_dim = new_indices.size();

	  // Evaluate the basis
	  //RealMatrix W;
	  //basis->value_set( pts, new_indices, W );
	  RealMatrix W( num_pts, current_dim, false);
	  for ( int j = 0; j < num_pts; j++ )
	    {
	      RealVector x( Teuchos::View, pts[j], num_vars );
	      for ( int i = 0; i < (int)new_indices.size(); i++ )
		{ 
		  W(j,i) = data_rep->
		    multivariate_polynomial( x, new_indices[i] );
		}
	    }

	  util::permute_matrix_rows( W, p );
      
	  // Row-reduce W according to previous elimination steps
	  int m = W.numRows(), n = W.numCols();
	  for ( int q = 0; q < lu_row; q++ )
	    {
	      for ( int i = 0; i < n; i++ )
		W(q,i) /= l(q,q);

	      RealMatrix tmp( Teuchos::View, W, m-q-1, n, q+1, 0 ),
		l_col( Teuchos::View, l, m-q-1, 1, q+1, q ),
		W_row( Teuchos::View, W, 1, n, q, 0 ); 
	      tmp.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1.0, l_col, 
			    W_row, 1.0 );
	    }

	  //RealMatrix wm( Teuchos::View, W, m-lu_row, n, lu_row, 0 );
	  RealMatrix wm( n, m-lu_row );
	  for ( int i = 0; i < m-lu_row; i++ )
	    {
	      for ( int j = 0; j < n; j++ )
		wm(j,i) = W(lu_row+i,j);
	    }
	  RealMatrix Q, R;
	  IntVector evec;
	  util::pivoted_qr_factorization( wm, Q, R, evec );
      
	  // Compute rank
	  int rnk = 0;
	  for ( int i = 0; i < R.numRows(); i++ )
	    {
	      if ( std::fabs( R(i,i) ) < 0.001 * std::fabs( R(0,0) ) ) break;
	      rnk += 1;
	    }

	  // Now first we must permute the rows by e
	  IntMatrix p_sub( Teuchos::View, &p[lu_row], num_pts - lu_row, 
			   num_pts - lu_row, 1 );
	  util::permute_matrix_rows( p_sub, evec );
      
	  // And correct by permuting l as well
 
	  RealMatrix l_sub( Teuchos::View, l, num_pts - lu_row, lu_row,lu_row,0);
	  if ( ( l_sub.numRows() > 0 ) && ( l_sub.numCols() > 0 ) )
	    util::permute_matrix_rows( l_sub, evec );

	  // The matrix r gives us inner product information for all rows below 
	  // these in W
	  for ( int i = 0; i < rnk; i++ )
	    {
	      for ( int j = 0; j < num_pts - lu_row; j++ )
		l(lu_row+j,lu_row+i) = R(i,j);
	    }

	  // Now we must find inner products of all the other rows above 
	  // these in W
	  RealMatrix W_sub( Teuchos::View, W, lu_row, W.numCols(), 0, 0 ),
	    Q_sub( Teuchos::View, Q, Q.numRows(), rnk, 0, 0 ),
	    u_sub( Teuchos::View, u, lu_row, rnk, 0, lu_row );

	  if ( ( u_sub.numRows() > 0 ) && ( u_sub.numCols() > 0 ) )
	    u_sub.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, W_sub, 
			    Q_sub, 0.0 );

	  if ( v_index+(current_dim*rnk) > v.length() )
	    v.resize( v.length() + std::max( 1000, current_dim*rnk ) );
	  // The matrix q must be saved in order to characterize basis
	  int cnt = v_index;
	  for ( int j = 0; j < rnk; j++ )
	    {
	      for ( int i = 0; i < Q.numRows(); i++ )
		{
		  v[cnt] = Q(i,j);
		  cnt++;
		}
	    }
	  v_index += ( current_dim * rnk );

	  // Update degree markers, and node and degree count
	  for ( int i = lu_row; i < lu_row+rnk; i++ )
	    k[i] = k_counter;
	  lu_row += rnk;
	  k_counter++;
	}
    }
  // Chop off parts of unnecessarily allocated vector v
  v.resize( v_index );
  
  // Return the indices used by least interpolation. This may be different
  // to the basis_indices specified on entry.
  basis_indices = internal_basis_indices;

  // Make matrix H
  get_least_polynomial_coefficients( v, k, basis_indices, num_vars, num_pts,
				     H );
}


void RegressOrthogPolyApproximation::
get_least_polynomial_coefficients( RealVector &v, IntVector &k,
				   UShort2DArray &basis_indices, int num_vars,
				   int num_pts, RealMatrix &H )
{
  int num_basis_indices = basis_indices.size();
  H.shape( num_pts, num_basis_indices );
  int v_counter = 0, previous_dimension = 0, current_size = 0;
  for ( int i = 0; i < num_pts; i++ )
    {
      if ( ( i == 0 ) || ( k[i] != k[i-1] ) )
	{
	  current_size = 0;
	  for ( int j = 0; j < num_basis_indices; j++ )
	    {
	      if ( l1_norm( basis_indices[j] ) == k[i] )
		current_size++;
	    }
	}
      for ( int j = 0; j < current_size; j++ )
	{
	  H(i,previous_dimension+j) = v[v_counter+j];
	}
      v_counter += current_size;
      if ( ( i+1 < num_pts ) && ( k[i] != k[i+1] ) )
	previous_dimension += current_size;
    }
}


void RegressOrthogPolyApproximation::
update_sparse(Real* dense_coeffs, size_t num_dense_terms)
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  SizetSet& sparse_ind = sparseIndIter->second;

  // just one pass through to define sparse_indices
  sparse_ind.clear();
  update_sparse_indices(dense_coeffs, num_dense_terms, sparse_ind);

  // update exp_coeffs
  update_sparse_coeffs(dense_coeffs, expCoeffsIter->second, sparse_ind);

  // update the sparse Sobol' indices
  update_sparse_sobol(sparse_ind, data_rep->multi_index(),
		      data_rep->sobolIndexMap);
}


void RegressOrthogPolyApproximation::
update_sparse_indices(Real* dense_coeffs, size_t num_dense_terms,
		      SizetSet& sparse_indices)
{
  // always retain leading coefficient (mean)
  if (sparse_indices.empty())
    sparse_indices.insert(0);
  // update sparse_indices based on nonzero coeffs
  for (size_t i=1; i<num_dense_terms; ++i)
    if (std::abs(dense_coeffs[i]) > DBL_EPSILON)
      sparse_indices.insert(i); // discards duplicates (coeffs, coeffgrads)
}


void RegressOrthogPolyApproximation::
update_sparse_coeffs(Real* dense_coeffs, RealVector& exp_coeffs,
		     const SizetSet& sparse_indices)
{
  // build sparse exp_coeffs
  size_t num_exp_terms = sparse_indices.size();
  exp_coeffs.sizeUninitialized(num_exp_terms);
  size_t i; SizetSet::const_iterator cit;
  for (i=0, cit=sparse_indices.begin(); i<num_exp_terms; ++i, ++cit)
    exp_coeffs[i] = dense_coeffs[*cit];
}


void RegressOrthogPolyApproximation::
update_sparse_coeff_grads(Real* dense_coeffs, int row,
			  RealMatrix& exp_coeff_grads,
			  const SizetSet& sparse_indices)
{
  // build sparse exp_coeff_grads
  size_t num_exp_terms = sparse_indices.size();
  if (exp_coeff_grads.numCols() != num_exp_terms)
    exp_coeff_grads.reshape(modSurrData.num_derivative_variables(),
			    num_exp_terms);
  int j; SizetSet::const_iterator cit;
  for (j=0, cit=sparse_indices.begin(); j<num_exp_terms; ++j, ++cit)
    exp_coeff_grads(row, j) = dense_coeffs[*cit];
}


void RegressOrthogPolyApproximation::
update_sparse_sobol(const SizetSet& sparse_indices,
		    const UShort2DArray& shared_multi_index,
		    const BitArrayULongMap& shared_sobol_map)
{
  // define the Sobol' indices based on the sparse multiIndex
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  if ( !data_rep->expConfigOptions.vbdFlag ||
        data_rep->expConfigOptions.vbdOrderLimit == 1 )
    return;

  sparseSobolIndexMap.clear();
  if (sparse_indices.empty()) {
    size_t sobol_len = shared_sobol_map.size();
    if (sobolIndices.length() != sobol_len)
      sobolIndices.sizeUninitialized(sobol_len);
    return;
  }

  // update sparseSobolIndexMap from shared sobolIndexMap.  Since
  // sparseSobolIndexMap is sorted by the sobolIndexMap map-values
  // (indices within a dense sobolIndices), we know that the ordering
  // of the sparse interactions will be consistent, e.g.:
  // sparseSobolIndexMap keys:   0, 2, 4, 5, 9, 11 (from sobolIndexMap)
  // sparseSobolIndexMap values: 0, 1, 2, 3, 4, 5  (new sobolIndices sequence)
  size_t j, num_v = sharedDataRep->numVars;
  StSCIter sit; BAULMCIter cit; BitArray set(num_v);
  for (sit=sparse_indices.begin(); sit!=sparse_indices.end(); ++sit) {
    const UShortArray& sparse_mi = shared_multi_index[*sit];
    for (j=0; j<num_v; ++j)
      if (sparse_mi[j]) set.set(j);   //   activate bit j
      else              set.reset(j); // deactivate bit j
    // define map from shared index to sparse index
    cit = shared_sobol_map.find(set);
    if (cit == shared_sobol_map.end()) {
      if (set.count() <= data_rep->expConfigOptions.vbdOrderLimit) {
	PCerr << "Error: sobolIndexMap lookup failure in RegressOrthogPoly"
	      << "Approximation::update_sparse_sobol() for multi-index\n"
	      << sparse_mi << std::endl;
	abort_handler(-1);
      }
      // else contributions to this set are not tracked
    }
    else // {key,val} = {shared,sparse} index: init sparse to 0 (updated below)
      sparseSobolIndexMap[cit->second] = 0;
  }
  // now that keys are complete, assign new sequence for sparse Sobol indices
  unsigned long sobol_len = 0; ULULMIter mit;
  for (mit=sparseSobolIndexMap.begin(); mit!=sparseSobolIndexMap.end(); ++mit)
    mit->second = sobol_len++;
  // now size sobolIndices
  if (sobolIndices.length() != sobol_len)
    sobolIndices.sizeUninitialized(sobol_len);
}


/** This version clears out all terms not recovered in CS and may
    leave multi-index gaps. */
void RegressOrthogPolyApproximation::
sparse_restriction(UShort2DArray& multi_index, SizetSet& sparse_indices)
{
  if (sparse_indices.empty())
    return; // no further restriction possible

  UShort2DArray dense_mi(multi_index); // copy

  // resize arrays into packed format and copy entries from dense data
  size_t i, num_exp_terms = sparse_indices.size(); SizetSet::const_iterator cit;
  multi_index.resize(num_exp_terms);
  for (i=0, cit=sparse_indices.begin(); i<num_exp_terms; ++i, ++cit)
    multi_index[i] = dense_mi[*cit];
  sparse_indices.clear();
}


/** This version clears out only the dominated multi-index terms not
    recovered in CS to avoid multi-index gaps; it provides a frontier
    with all supporting terms. */
void RegressOrthogPolyApproximation::
frontier_restriction(UShort2DArray& multi_index, SizetSet& sparse_indices)
{
  if (sparse_indices.empty()) // no further restriction possible
    return;

  // first define pareto_mi from the recovered mi terms
  // Note: a complete frontier_mi cannot be defined from recovered sparse terms.
  // Therefore, we employ pareto_mi here for purposes of restriction (+ sparse
  // re-indexing). See advance_multi_index_front() for later use of frontier_mi.
  size_t i, j, num_exp_terms = sparse_indices.size(),
    num_mi = multi_index.size(), new_si;
  UShort2DArray pareto_mi;
  SizetSet::const_iterator si_it;
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  for (i=0, si_it=sparse_indices.begin(); i<num_exp_terms; ++i, ++si_it)
    update_pareto_set(multi_index[*si_it], pareto_mi);
  size_t num_pareto_mi = pareto_mi.size();

  UShort2DArray orig_multi_index(multi_index);  // copy
  SizetSet orig_sparse_indices(sparse_indices); // copy
  multi_index.clear(); sparse_indices.clear();

  // prune out only the multi_index terms that dominate pareto_mi, while
  // updating sparse indices for this pruned multi_index.  We restrict back
  // to the Pareto frontier while ignoring any gaps due to sparse recovery.
  // > this simplifies bookkeeping and forward neighbor searches
  // > JDJ: this is consistent with the observation that the algorithm fails
  //        when there isn't a simple tree structure to the recovered coeffs
  bool i_dominated, j_dominated;
  for (i=0, si_it=orig_sparse_indices.begin(), new_si=0; i<num_mi; ++i) {
    const UShortArray& mi_i = orig_multi_index[i];
    for (j=0; j<num_pareto_mi; ++j) {
      // manage tie: i is "challenger" and is dominated if == frontier member;
      // j is "incumbent" and is not dominated if equal to challenger
      assess_dominance(mi_i, pareto_mi[j], i_dominated, j_dominated);
      if (i_dominated) break;
    }
#ifdef DEBUG
    PCout << "orig mi[" << i << "] dominated = " << i_dominated << std::endl;
#endif // DEBUG
    if (i_dominated) {  // within pareto set
      multi_index.push_back(mi_i);
      if (*si_it == i) { // a recovered term
	sparse_indices.insert(new_si); // reindex within pruned multi_index
	++si_it;         // advance to next original sparse index
      }
      ++new_si;
    }
  }
#ifdef DEBUG
  PCout << "RegressOPA::frontier_restriction(): original multi_index =\n"
	<< orig_multi_index << "original sparse_indices =\n"
	<< orig_sparse_indices << "pareto set =\n" << pareto_mi
	<< "new multi_index =\n" << multi_index << "new sparse_indices =\n"
	<< sparse_indices << std::endl;
#endif // DEBUG
}


void RegressOrthogPolyApproximation::
advance_multi_index(const UShort2DArray& multi_index,
		    UShortArraySetArray& mi_advancements)
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;

  // create a set of advancements from this reference frontier
  size_t i, num_mi,
    num_advance = data_rep->regressConfigOptions.numAdvancements;
  mi_advancements.resize(num_advance);
  add_admissible_forward_neighbors(multi_index, mi_advancements[0]);
  if (num_advance > 1) {
    UShort2DArray combined_mi = multi_index; // copy
    for (i=1; i<num_advance; ++i) {
      append_multi_index(mi_advancements[i-1], combined_mi);
      add_admissible_forward_neighbors(combined_mi, mi_advancements[i]);
    }
  }
}


void RegressOrthogPolyApproximation::
advance_multi_index_front(const UShort2DArray& multi_index,
			  UShortArraySetArray& mi_advancements)
{
  // Compute frontier_mi from the restricted (but complete) multi_index.  This
  // provides a reduced multi_index for purposes of fwd neighbor expansion.
  UShortArraySet mi_frontier;
  update_pareto_frontier(multi_index, mi_frontier);

  // create a set of advancements from this reference frontier
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  size_t i, num_mi,
    num_advance = data_rep->regressConfigOptions.numAdvancements;
  mi_advancements.resize(num_advance);

  add_admissible_forward_neighbors(mi_frontier, mi_advancements[0]);
  for (i=1; i<num_advance; ++i) {
    update_pareto_frontier(mi_advancements[i-1], mi_frontier);
    add_admissible_forward_neighbors(mi_frontier, mi_advancements[i]);
  }
  /*
  for (i=0; i<num_advance; ++i) {
    UShortArraySet& mi_adv_i = mi_advancements[i];
    add_admissible_forward_neighbors(mi_frontier, mi_adv_i);
    if (i+1 < num_advance)
      update_pareto_frontier(mi_adv_i, mi_frontier);
  }
  */
}


void RegressOrthogPolyApproximation::
add_admissible_forward_neighbors(const UShort2DArray& reference_mi,
				 UShortArraySet& fwd_neighbors)
{
  // In current usage, this version does not use a frontier; i.e., reference_mi
  // reflects the full expansion and may contain gaps.  However, both frontier
  // and gappy cases must additionally check for presence of an admissible fwd
  // neighbor within the reference_mi.  Since there is no functional difference,
  // sort ref_mi into a set for find efficiency and employ the frontier version.
  UShortArraySet reference_mi_set(reference_mi.begin(), reference_mi.end());
  add_admissible_forward_neighbors(reference_mi_set, fwd_neighbors);

#ifdef DEBUG
  PCout << "RegressOPA::add_admissible_forward_neighbors(): reference_mi =\n"
	<< reference_mi << "fwd_neighbors =\n" << fwd_neighbors << std::endl;
#endif // DEBUG
}


void RegressOrthogPolyApproximation::
add_admissible_forward_neighbors(const UShortArraySet& reference_mi,
				 UShortArraySet& fwd_neighbors)
{
  // This function is similar to SparseGridDriver::add_active_neighbors()

  UShortArraySet::const_iterator ref_cit, find_cit;
  size_t s, i, j, num_v = sharedDataRep->numVars;
  fwd_neighbors.clear();
  for (ref_cit=reference_mi.begin(); ref_cit!=reference_mi.end(); ++ref_cit) {
    UShortArray neighbor = *ref_cit; // mutable copy
    for (i=0; i<num_v; ++i) {
      // i^{th} candidate for set A (active) computed from forward neighbor:
      // increment by 1 in dimension i
      unsigned short& neighbor_i = neighbor[i];
      ++neighbor_i;
      if (reference_mi.find(neighbor) == reference_mi.end()) {
	// test all backwards neighbors for membership in set O (old)
	bool backward_old = true;
	for (j=0; j<num_v; ++j) {
	  unsigned short& neighbor_j = neighbor[j];
	  if (neighbor_j) { // if 0, then admissible by default
	    --neighbor_j;
	    find_cit = reference_mi.find(neighbor);
	    ++neighbor_j; // restore
	    if (find_cit == reference_mi.end())
	      { backward_old = false; break; }
	  }
	}
	if (backward_old)
	  fwd_neighbors.insert(neighbor);// std::set will discard any duplicates
      }
      --neighbor_i; // restore
    }
  }
#ifdef DEBUG
  PCout << "RegressOPA::add_admissible_forward_neighbors(): reference_mi =\n"
	<< reference_mi << "fwd_neighbors =\n" << fwd_neighbors << std::endl;
#endif // DEBUG
}


void RegressOrthogPolyApproximation::gridSearchFunction( RealMatrix &opts,
						  int M, int N, 
						  int num_function_samples )
{
  // Setup a grid based search
  bool is_under_determined = M < N;

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;

  // Define the 1D grids for under and over-determined LARS, LASSO, OMP, BP and 
  // LS
  RealVectorArray opts1D( 9 );
  opts1D[0].size( 1 ); // solver type
  opts1D[0][0] = CSOpts.solver;
  opts1D[1].size( 1 ); // Solver tolerance. 
  opts1D[1][0] = CSOpts.solverTolerance;
  opts1D[2] = data_rep->regressConfigOptions.noiseTols; // epsilon.
  opts1D[3].size( 1 ); // delta
  opts1D[3] = CSOpts.delta;
  opts1D[4].size( 1 ); // max_number of non_zeros
  opts1D[4] = CSOpts.maxNumIterations; 
  opts1D[5].size( 1 );  // standardizeInputs
  opts1D[5] = false;
  opts1D[6].size( 1 );  // storeHistory
  opts1D[6] = true;  
  opts1D[7].size( 1 );  // Verbosity. Warnings on
  opts1D[7] =  std::max(0, data_rep->expConfigOptions.outputLevel - 1);
  opts1D[8].size( 1 );  // num function samples
  opts1D[8] = num_function_samples;
      
  // Form the multi-dimensional grid
  util::cartesian_product( opts1D, opts, 1 );
};

void RegressOrthogPolyApproximation::
estimate_compressed_sensing_options_via_cross_validation( RealMatrix &vandermonde_matrix, RealMatrix &rhs, std::vector<CompressedSensingOptions> &best_cs_opts, RealVector &best_predictor_indicators, RealMatrixArray &predictor_options_history, RealMatrixArray &predictor_indicators_history, RealMatrixArray &predictor_partition_indicators_history, size_t num_data_pts_fn ){
};


void RegressOrthogPolyApproximation::compute_component_sobol()
{
  if (sparseIndIter == sparseIndices.end() || sparseIndIter->second.empty())
    { OrthogPolyApproximation::compute_component_sobol(); return; }

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
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multi_index();
  RealVector& exp_coeffs = expCoeffsIter->second;
  const SizetSet& sparse_ind = sparseIndIter->second;
  const BitArrayULongMap& index_map = data_rep->sobolIndexMap;
  size_t i, j, num_v = sharedDataRep->numVars; StSCIter cit;
  BitArray set(num_v);
  Real p_var, sum_p_var = 0.;
  for (i=1, cit=++sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit) {
    const UShortArray& mi_i = mi[*cit];
    p_var = exp_coeffs(i) * exp_coeffs(i) * data_rep->norm_squared(mi_i);
    sum_p_var += p_var;

    // determine the bit set corresponding to this expansion term
    for (j=0; j<num_v; ++j)
      if (mi_i[j]) set.set(j);   //   activate bit j
      else         set.reset(j); // deactivate bit j

    // lookup the bit set within sobolIndexMap --> increment the correct
    // Sobol' index with the variance contribution from this expansion term.
    BAULMCIter cit = index_map.find(set);
    if (cit != index_map.end()) { // may not be found if vbdOrderLimit
      // sparseSobolIndexMap definition is bypassed when vbdOrderLimit == 1
      unsigned long sp_index = (data_rep->expConfigOptions.vbdOrderLimit == 1) ?
	cit->second : sparseSobolIndexMap[cit->second];
      sobolIndices[sp_index] += p_var; // divide by sum_p_var below
    }
  }
  if (sum_p_var > SMALL_NUMBER) // don't attribute variance if zero/negligible
    sobolIndices.scale(1./sum_p_var);

#ifdef DEBUG
  PCout << "In RegressOrthogPolyApproximation::compute_component_sobol(), "
	<< "sobolIndices =\n" << sobolIndices;
#endif // DEBUG
}


void RegressOrthogPolyApproximation::compute_total_sobol() 
{
  if (sparseIndIter == sparseIndices.end() || sparseIndIter->second.empty())
    { OrthogPolyApproximation::compute_total_sobol(); return; }

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  size_t j, num_v = sharedDataRep->numVars;
  const UShort2DArray& mi = data_rep->multi_index();
  const SizetSet& sparse_ind = sparseIndIter->second;
  RealVector& exp_coeffs = expCoeffsIter->second;
  totalSobolIndices = 0.;
  if (data_rep->expConfigOptions.vbdOrderLimit) {
    // all component indices may not be available, so compute total indices
    // from scratch by computing and summing variance contributions for each
    // expansion term
    size_t i, num_exp_terms = sparseIndices.size();
    Real p_var, sum_p_var = 0., ratio_i; StSCIter cit;
    for (i=1, cit=++sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit) {
      const UShortArray& mi_i = mi[*cit];
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
    // all component effects have been computed, so add them up:
    // totalSobolIndices parses the bit sets of each of the sobolIndices
    // and adds them to each matching variable bin
    // Note: compact iteration over sparseSobolIndexMap could be done but
    //       requires a value to key mapping to get from uit->first to set.
    for (BAULMCIter cit=index_map.begin(); cit!=index_map.end(); ++cit) {
      ULULMIter uit = sparseSobolIndexMap.find(cit->second);
      if (uit != sparseSobolIndexMap.end()) {
	const BitArray& set = cit->first;
	Real comp_sobol = sobolIndices[uit->second];
	for (j=0; j<num_v; ++j) 
	  if (set[j]) // var j is present in this Sobol' index
	    totalSobolIndices[j] += comp_sobol;
      }
    }
  }

#ifdef DEBUG
  PCout << "In RegressOrthogPolyApproximation::compute_total_sobol(), "
	<< "totalSobolIndices =\n" << totalSobolIndices;
#endif // DEBUG
}


const RealVector& RegressOrthogPolyApproximation::dimension_decay_rates()
{
  if (sparseIndIter == sparseIndices.end() || sparseIndIter->second.empty())
    return OrthogPolyApproximation::dimension_decay_rates();

  size_t i, j, num_exp_terms = sparseIndices.size(),
    num_v = sharedDataRep->numVars;
  if (decayRates.empty())
    decayRates.sizeUninitialized(num_v);

  // define max_orders for each var for sizing LLS matrices/vectors
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multi_index();
  const SizetSet& sparse_ind = sparseIndIter->second;
  RealVector& exp_coeffs = expCoeffsIter->second;
  UShortArray max_orders(num_v, 0); StSCIter cit;
  for (cit=++sparse_ind.begin(); cit!=sparse_ind.end(); ++cit) {
    const UShortArray& mi_i = mi[*cit];
    for (j=0; j<num_v; ++j)
      if (mi_i[j] > max_orders[j])
	max_orders[j] = mi_i[j];
  }

  // size A_vectors and b_vectors and initialize with defaults since
  // sparseIndices could leave gaps
  RealVectorArray A_vectors(num_v), b_vectors(num_v);
  unsigned short order, non_zero, var_index, order_index;
  Real norm;
  for (i=0; i<num_v; ++i) {
    unsigned short max_ord = max_orders[i];
    RealVector& A_i = A_vectors[i]; A_i.sizeUninitialized(max_ord);
    RealVector& b_i = b_vectors[i]; b_i.sizeUninitialized(max_ord);
    BasisPolynomial& poly_i = data_rep->polynomialBasis[i];
    for (j=0; j<max_ord; ++j) {
      order = j + 1;
      A_i[j] = (Real)order;
      // log(norm * 1.e-25) for cut-off coeff value of 1.e-25
      b_i[j] = std::log10(poly_i.norm_squared(order))/2. - 25.; // updated below
    }
  }

  // populate A_vectors and b_vectors
  bool univariate;
  for (cit=++sparse_ind.begin(), i=1; cit!=sparse_ind.end(); ++cit, ++i) {
    const UShortArray& mi_i = mi[*cit];
    univariate = true; non_zero = 0;
    for (j=0; j<num_v; ++j) {
      if (mi_i[j]) {
	++non_zero;
	if (non_zero > 1) { univariate = false; break; }
	else { order = mi_i[j]; var_index = j; order_index = order-1; }
      }
    }
    if (univariate) {
      // find a for y = ax + b with x = term order, y = log(coeff), and
      // b = known intercept for order x = 0
      Real abs_coeff = std::abs(exp_coeffs[i]);
#ifdef DECAY_DEBUG
      PCout << "Univariate contribution: order = " << order << " coeff = "
	    << abs_coeff << " norm = " << std::sqrt(data_rep->
	       polynomialBasis[var_index].norm_squared(order)) << '\n';
#endif // DECAY_DEBUG
      if (abs_coeff > 1.e-25)
	// b = std::log10(abs_coeff * norm), but don't recompute norm
	b_vectors[var_index][order_index] += 25. + std::log10(abs_coeff);
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


RealVector RegressOrthogPolyApproximation::
approximation_coefficients(bool normalized) const
{
  if (sparseIndIter == sparseIndices.end() || sparseIndIter->second.empty())
    return OrthogPolyApproximation::approximation_coefficients(normalized);

  // synchronize the expansion coeffs with the length of the shared multi-index
  // approx_coeffs is an inflation of expansionCoeffs
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multi_index();
  const RealVector& exp_coeffs = expCoeffsIter->second;
  const SizetSet& sparse_ind = sparseIndIter->second;
  RealVector approx_coeffs(mi.size()); // init to 0
  size_t i; StSCIter cit;
  for (i=0, cit=sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit)
    approx_coeffs[*cit] = (normalized) ?
      // basis is divided by norm, so coeff is multiplied by norm
      exp_coeffs[i] * std::sqrt(data_rep->norm_squared(mi[*cit])) :
      exp_coeffs[i];
  return approx_coeffs;
}


void RegressOrthogPolyApproximation::
approximation_coefficients(const RealVector& approx_coeffs, bool normalized)
{
  if (sparseIndIter == sparseIndices.end() || sparseIndIter->second.empty()) {
    OrthogPolyApproximation::
      approximation_coefficients(approx_coeffs,normalized);
    return;
  }

  // won't happen with current import use cases

  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  update_active_iterators(data_rep->activeKey);

  // synchronize the expansion coeffs with the shared multi-index length
  // expansionCoeffs is a deflation of approx_coeffs
  const UShort2DArray& mi = data_rep->multi_index();
  RealVector& exp_coeffs = expCoeffsIter->second;
  const SizetSet& sparse_ind = sparseIndIter->second;
  size_t i, num_sparse = sparseIndices.size(); StSCIter cit;
  if (exp_coeffs.length() != num_sparse)
    exp_coeffs.sizeUninitialized(num_sparse);
  for (i=0, cit=sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit)
    exp_coeffs[i] = (normalized) ?
      approx_coeffs[*cit] / std::sqrt(data_rep->norm_squared(mi[*cit])) :
      approx_coeffs[*cit];

  // allocate arrays in support of external coefficient import (mirrors
  // allocate_arrays() except for redundant size_expansion())
  allocate_total_sobol();
  allocate_component_sobol();
  RealVector& exp_mom = expMomentsIter->second;
  if (exp_mom.length() != 2) exp_mom.sizeUninitialized(2);
}


void RegressOrthogPolyApproximation::
print_coefficients(std::ostream& s, bool normalized)
{
  if (sparseIndIter == sparseIndices.end() || sparseIndIter->second.empty())
    { OrthogPolyApproximation::print_coefficients(s, normalized); return; }

  size_t i, j, num_v = sharedDataRep->numVars;
  StSCIter cit; char tag[10];

  // terms and term identifiers
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multi_index();
  const RealVector& exp_coeffs = expCoeffsIter->second;
  const SizetSet& sparse_ind = sparseIndIter->second;
  for (i=0, cit=sparse_ind.begin(); cit!=sparse_ind.end(); ++i, ++cit) {
    const UShortArray& mi_i = mi[*cit];
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


void RegressOrthogPolyApproximation::
coefficient_labels(std::vector<std::string>& coeff_labels) const
{
  if (sparseIndIter == sparseIndices.end() || sparseIndIter->second.empty())
    { OrthogPolyApproximation::coefficient_labels(coeff_labels); return; }

  size_t i, j, num_v = sharedDataRep->numVars;
  char tag[10];

  const SizetSet& sparse_ind = sparseIndIter->second;
  coeff_labels.reserve(sparse_ind.size());

  // terms and term identifiers
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  const UShort2DArray& mi = data_rep->multi_index();
  for (StSCIter cit=sparse_ind.begin(); cit!=sparse_ind.end(); ++cit) {
    const UShortArray& mi_i = mi[*cit];
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
