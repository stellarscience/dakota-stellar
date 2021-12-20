/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 CombinedSparseGridDriver
//- Description: Implementation code for CombinedSparseGridDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "CombinedSparseGridDriver.hpp"
#include "SharedPolyApproxData.hpp"
#include "sandia_sgmg.hpp"
#include "sandia_sgmga.hpp"
#include "sandia_sgmgg.hpp"
#include "pecos_stat_util.hpp"
#include "pecos_math_util.hpp"

static const char rcsId[]="@(#) $Id: CombinedSparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {

/// initialize static member pointer to active driver instance
CombinedSparseGridDriver* CombinedSparseGridDriver::sgdInstance(NULL);


void CombinedSparseGridDriver::
initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
		const MultivariateDistribution& u_dist,
		const ExpansionConfigOptions& ec_options,
		BasisConfigOptions& bc_options, short growth_rate,
		bool track_colloc, bool track_uniq_prod_wts)
{
  SparseGridDriver::initialize_grid(ssg_level, dim_pref, u_dist, ec_options,
				    bc_options, growth_rate);
  trackCollocDetails     = track_colloc;
  trackUniqueProdWeights = track_uniq_prod_wts;

  // set compute1D{Points,Type1Weights,Type2Weights}
  initialize_rule_pointers();
  // set levelGrowthToOrder
  initialize_growth_pointers();
}


void CombinedSparseGridDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis)
{
  IntegrationDriver::initialize_grid(poly_basis);

  // set compute1D{Points,Type1Weights,Type2Weights}
  initialize_rule_pointers();
  // set levelGrowthToOrder
  initialize_growth_pointers();
}


void CombinedSparseGridDriver::
initialize_grid_parameters(const MultivariateDistribution& mv_dist)
{
  SharedPolyApproxData::
    update_basis_distribution_parameters(mv_dist, polynomialBasis);

  // set a rule-dependent duplicateTol
  initialize_duplicate_tolerance(); // depends on length scale from dist params
}


void CombinedSparseGridDriver::clear_inactive()
{
  SparseGridDriver::clear_inactive();

  std::map<UShortArray, UShort2DArray>::iterator sm_it
    = smolyakMultiIndex.begin();
  std::map<UShortArray, IntArray>::iterator  sc_it = smolyakCoeffs.begin();
  std::map<UShortArray, UShort3DArray>::iterator ck_it = collocKey.begin();
  std::map<UShortArray, Sizet2DArray>::iterator  ci_it = collocIndices.begin();
  std::map<UShortArray, RealVector>::iterator t1_it = type1WeightSets.begin();
  std::map<UShortArray, RealMatrix>::iterator t2_it = type2WeightSets.begin();

  while (sm_it != smolyakMultiIndex.end())
    if (sm_it == smolMIIter) { // preserve active
      ++sm_it; ++sc_it; ++ck_it; ++ci_it;
      if (trackUniqueProdWeights)
	{ ++t1_it; if (computeType2Weights) ++t2_it; }
    }
    else { // clear inactive: postfix increments manage iterator invalidations
      smolyakMultiIndex.erase(sm_it++);   smolyakCoeffs.erase(sc_it++);
      collocKey.erase(ck_it++);           collocIndices.erase(ci_it++);
      if (trackUniqueProdWeights) {
	type1WeightSets.erase(t1_it++);
	if (computeType2Weights)          type2WeightSets.erase(t2_it++);
      }
    }
}


void CombinedSparseGridDriver::initialize_duplicate_tolerance()
{
  bool parameterized_basis = false, numerical_basis = false;
  for (size_t i=0; i<numVars; ++i) {
    short rule = collocRules[i];
    if (rule == GOLUB_WELSCH)
      { numerical_basis = true; break; }
    else if (rule == GEN_GAUSS_LAGUERRE || rule == GAUSS_JACOBI)
      parameterized_basis = true;
  }

  // Allow a looser duplication tolerance for numerically-generated rules than
  // for lookup tables.  sgmg,sgmga use an absolute tolerance which is fine for
  // standardized probability distributions (scaled to O(1)), but exhibits scale
  // dependence for numerically-generated quadrature rules.

  // rules from eigensolves or parameterized solves:
  if (numerical_basis || parameterized_basis) duplicateTol = 1.e-14;
  // rules mostly from lookup tables:
  else duplicateTol = 1.e-15; // ~= 4.5 * DBL_EPSILON

  // For numerically-generated rules, convert duplicateTol to a relative tol by
  // using a length scale.  sandia_rules.cpp uses "(dist <= tol)", so we could
  // modify sandia_rules (in many places) to use "(dist / length_scale <= tol)"
  // or just modify duplicateTol such that new_tol = tol * length_scale.  The
  // length_scale() function returns max(mean,stdev) in u-space (stdev provides
  // a backup for small/zero mean).
  if (numerical_basis) {
    Real length_scale = 0.;
    for (size_t i=0; i<numVars; ++i)
      length_scale += std::pow(polynomialBasis[i].length_scale(), 2);
                    //std::pow(std::max(ranVarMeansU[i], ranVarStdDevsU[i]), 2);
    if (length_scale > DBL_MIN)
      duplicateTol *= std::sqrt(length_scale);
  }
}


void CombinedSparseGridDriver::initialize_rule_pointers()
{
  size_t i, j;
  // compute1DPoints needed for grid_size() and for sgmg/sgmga
  compute1DPoints.resize(numVars);
  for (i=0; i<numVars; ++i)
    compute1DPoints[i] = basis_collocation_points;
  // compute1D{Type1,Type2}Weights only needed for sgmg/sgmga
  if (!refineControl) {
    compute1DType1Weights.resize(numVars);
    for (i=0; i<numVars; i++)
      compute1DType1Weights[i] = basis_type1_collocation_weights;
    /*
    if (computeType2Weights) {
      compute1DType2Weights.resize(numVars);
      for (i=0; i<numVars; ++i) {
	std::vector<CollocFnPtr>& comp_1d_t2_wts_i = compute1DType2Weights[i];
	comp_1d_t2_wts_i.resize(numVars);
	for (j=0; j<numVars; ++j)
	  comp_1d_t2_wts_i[j] = (j==i) ? basis_type2_collocation_weights :
	                                 basis_type1_collocation_weights;
      }
    }
    */
  }
}


void CombinedSparseGridDriver::initialize_growth_pointers()
{
  levelGrowthToOrder.resize(numVars);

  // if INTEGRATION with restricted growth, sync on integrand precision.
  // if INTERPOLATION with restricted growth, sync on number of points
  // (or Lagrange interpolant order = #pts - 1).
  if (driverMode == INTERPOLATION_MODE)
    for (size_t i=0; i<numVars; ++i)
      switch (collocRules[i]) {
      // nested rules with exponential growth:
      case GENZ_KEISTER:
	levelGrowthToOrder[i] = level_to_order_exp_hgk_interp;    break;
      case CLENSHAW_CURTIS: case NEWTON_COTES:
	levelGrowthToOrder[i] = level_to_order_exp_closed_interp; break;
      case GAUSS_PATTERSON: case FEJER2:
	levelGrowthToOrder[i] = level_to_order_exp_open_interp;   break;
      // non-nested or weakly-nested Gauss rules with linear growth
      case GAUSS_HERMITE: case GAUSS_LEGENDRE: // symmetric Gauss linear growth
	levelGrowthToOrder[i] = webbur::level_to_order_linear_wn; break;
      default: // asymmetric Gauss linear growth
	levelGrowthToOrder[i] = webbur::level_to_order_linear_nn; break;
      }
  else // INTEGRATION_MODE or DEFAULT_MODE
    for (size_t i=0; i<numVars; ++i)
      switch (collocRules[i]) {
      // nested rules with exponential growth:
      case GAUSS_PATTERSON:
	levelGrowthToOrder[i] = webbur::level_to_order_exp_gp;    break;
      case GENZ_KEISTER:
	levelGrowthToOrder[i] = webbur::level_to_order_exp_hgk;   break;
      case CLENSHAW_CURTIS: case NEWTON_COTES:
	levelGrowthToOrder[i] = webbur::level_to_order_exp_cc;    break;
      case FEJER2:
	levelGrowthToOrder[i] = webbur::level_to_order_exp_f2;    break;
      // non-nested or weakly-nested Gauss rules with linear growth
      case GAUSS_HERMITE: case GAUSS_LEGENDRE: // symmetric Gauss linear growth
	levelGrowthToOrder[i] = webbur::level_to_order_linear_wn; break;
      default: // asymmetric Gauss linear growth
	levelGrowthToOrder[i] = webbur::level_to_order_linear_nn; break;
      }
}


void CombinedSparseGridDriver::
assign_smolyak_arrays(UShort2DArray& sm_mi, IntArray& sm_coeffs)
{
  // Populate smolyakMultiIndex and smolyakCoeffs.  Identifies
  // use of polynomialBasis[variable][index] based on index 0:num_levels-1.
  // w = q - N = dimension-independent level.  For isotropic,
  //   w + 1 <= |i| <= w + N for i starts at 1 (used for index set defn.)
  //   w - N + 1 <= |j| <= w for j = i - 1 starts at 0 (used for generation)
  // For anisotropic, a weighted linear index set constraint is used.

  size_t i;
  unsigned short ssg_lev = ssgLevIter->second;
  const RealVector& aniso_wts = anisoWtsIter->second;
  if (aniso_wts.empty()) { // initialize multi_index
    UShortArray levels(numVars, ssg_lev);
    SharedPolyApproxData::total_order_multi_index(levels, sm_mi, numVars-1);
    size_t num_terms = sm_mi.size();
    // initialize sm_coeffs
    sm_coeffs.resize(num_terms);
    for (i=0; i<num_terms; i++) {
      int wpNmi = ssg_lev - l1_norm(sm_mi[i]); // w+N-|i| = w-|j|
      sm_coeffs[i] = (int)std::pow(-1., wpNmi)
	* (int)std::floor(BasisPolynomial::n_choose_k(numVars - 1, wpNmi)+.5);
    }
  }
  else { // utilize webbur::sgmga_vcn_{ordered,coef}
    sm_mi.clear();
    sm_coeffs.clear();
    // Utilize webbur::sandia_sgmga_vcn_{ordered,coef} for 0-based index sets
    // (w*alpha_min-|alpha| < |alpha . j| <= w*alpha_min).
    // With scaling alpha_min = 1: w-|alpha| < |alpha . j| <= w.
    // In the isotropic case, reduces to w-N < |j| <= w, which is the same as
    // w-N+1 <= |j| <= w.
    IntArray x(numVars), x_max(numVars); //x_max = ssg_lev;
    UShortArray index_set(numVars);
    Real wt_i, wt_sum = 0., q_max = ssg_lev;
    for (i=0; i<numVars; ++i) {
      wt_i = aniso_wts[i];
      wt_sum += wt_i;
      // minimum nonzero weight is scaled to 1, so just catch special case of 0
      x_max[i] = (wt_i > 1.e-10) ? (int)std::ceil(q_max/wt_i) : 0;
    }
    Real q_min = ssg_lev - wt_sum;
#ifdef DEBUG
    PCout << "q_min = " << q_min << " q_max = " << q_max;
#endif // DEBUG

    bool more = false;
    Real *aniso_wt_vals = aniso_wts.values();
    int  *x0 = &x[0], *xm0 = &x_max[0], coeff;
    webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wt_vals, xm0, x0,
				     q_min, q_max, &more);
    while (more) {
      coeff
	= (int)webbur::sandia_sgmga_vcn_coef(numVars, aniso_wt_vals, x0, q_max);
      if (coeff) {
	sm_coeffs.push_back(coeff);
	for (i=0; i<numVars; ++i)
	  index_set[i] = (unsigned short)x[i];
	sm_mi.push_back(index_set);
      }
      webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wt_vals, xm0, x0, q_min,
				       q_max, &more);
    }
  }

#ifdef DEBUG
  size_t num_terms = sm_coeffs.size();
  PCout << "\nnum Smolyak terms = " << num_terms << '\n';
  for (i=0; i<num_terms; i++)
    PCout << "multi_index[" << i << "]:\n" << sm_mi[i]
	  << "coeffs[" << i << "] = " << sm_coeffs[i] << "\n\n";
#endif // DEBUG
}


void CombinedSparseGridDriver::
update_smolyak_coefficients(size_t start_index, const UShort2DArray& sm_mi,
			    IntArray& sm_coeffs)
{
  size_t num_sets = sm_mi.size();
  if (start_index >= num_sets) return;

  if (sm_coeffs.size() != num_sets)
    sm_coeffs.resize(num_sets);
  size_t j, cntr = 0, len1 = num_sets-1;
  int i, m = numVars, *s1 = new int [numVars*len1], *c1 = new int [len1],
    *s2 = new int [numVars];
  // initialize s1 and c1
  for (i=0; i<start_index; ++i) {
    c1[i] = sm_coeffs[i];
    for (j=0; j<numVars; ++j, ++cntr) // no copy_data() since ushort -> int
      s1[cntr] = sm_mi[i][j]; // sgmgg packs by variable groups
  }
  // for each s2, update sm_coeffs
  for (i=start_index; i<num_sets; ++i) {
    for (j=0; j<numVars; ++j) // no copy_data() since ushort -> int
      s2[j] = sm_mi[i][j];
    webbur::sandia_sgmgg_coef_inc2(m, i, s1, c1, s2, &sm_coeffs[0]);
#ifdef DEBUG
    PCout << "update_smolyak_coefficients(): updated Smolyak coeffs =\n"
	  << sm_coeffs << '\n';
#endif // DEBUG
    if (i<num_sets-1) { // if not last i, update s1 and c1 state for next pass
      for (j=0; j<=i; ++j) // coeffs updated to len i+1; max len = num_sets-1
	c1[j] = sm_coeffs[j];
      for (j=0; j<numVars; ++j, ++cntr)
	s1[cntr] = s2[j]; // max len = (num_sets-1)*numVars
    }
  }
  delete [] s1; delete [] c1; delete [] s2;
}


void CombinedSparseGridDriver::
assign_collocation_key(const UShort2DArray& sm_mi, UShort3DArray& colloc_key)
{
  // define mapping from collocation pts to set of 1d interpolation indices
  size_t i, num_sm_mi = sm_mi.size();
  colloc_key.resize(num_sm_mi);
  UShortArray quad_order(numVars); //, collocation_indices(numVars);
  for (i=0; i<num_sm_mi; ++i) {
    level_to_order(sm_mi[i], quad_order);
    SharedPolyApproxData::
      tensor_product_multi_index(quad_order, colloc_key[i], false);
  }
}


void CombinedSparseGridDriver::
assign_collocation_indices(const UShort3DArray& colloc_key,
			   const IntArray& unique_index_map,
			   Sizet2DArray& colloc_ind, size_t start_index)
{
  // define mapping from collocation pts to set of 1d interpolation indices
  size_t i, j, num_tp_pts, cntr = 0, num_sm_indices = colloc_key.size();
  colloc_ind.resize(num_sm_indices);
  // unique_index_map covers both reference and increment from start_index
  for (i=0; i<start_index; ++i)
    cntr += colloc_key[i].size();
  for (i=start_index; i<num_sm_indices; ++i) {
    num_tp_pts = colloc_key[i].size();
    SizetArray& indices_i = colloc_ind[i];
    indices_i.resize(num_tp_pts);
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      indices_i[j] = unique_index_map[cntr];
#ifdef DEBUG
      PCout << "collocKey[" << i << "][" << j << "]:\n" << colloc_key[i][j]
	    << "collocIndices[" << i << "][" << j << "] = " << indices_i[j]
	    << '\n';
#endif // DEBUG
    }
  }
}


int CombinedSparseGridDriver::grid_size()
{
  int& num_colloc_pts = numPtsIter->second;
  if (num_colloc_pts == 0) { // special value indicated update required
    sgdInstance = this; // sgdInstance required within compute1DPoints below
    unsigned short ssg_lev =   ssgLevIter->second;
    RealVector&  aniso_wts = anisoWtsIter->second;
    num_colloc_pts = (aniso_wts.empty()) ?
      webbur::sgmg_size(numVars, ssg_lev, &compute1DPoints[0], duplicateTol,
	growthRate, &levelGrowthToOrder[0]) :
      webbur::sandia_sgmga_size(numVars, aniso_wts.values(), ssg_lev,
	&compute1DPoints[0], duplicateTol, growthRate, &levelGrowthToOrder[0]);
  }
  return num_colloc_pts;
}


void CombinedSparseGridDriver::
reinterpolated_tensor_grid(const UShortArray& lev_index,
			   const SizetList& reinterp_indices)
{
  std::map<UShortArray, size_t>::iterator map_it = reinterpMap.find(lev_index);
  if (map_it == reinterpMap.end()) {

    // update arrays in place
    activeReinterpIndex = reinterpLevelIndices.size();
    UShortArray usa; UShort2DArray us2a; RealMatrix rm;
    reinterpLevelIndices.push_back(usa);
    reinterpQuadOrders.push_back(usa);
    reinterpVarSets.push_back(rm);
    reinterpCollocKeys.push_back(us2a);
    UShortArray& reinterp_lev_index  = reinterpLevelIndices.back();
    UShortArray& reinterp_quad_order = reinterpQuadOrders.back();
    reinterp_quad_order.resize(numVars); reinterp_lev_index.resize(numVars);

    // adjust level for reinterpolation of covariance (goal = 2x the interpolant
    // order).  For nested rules, this may not double the integrand exactness,
    // but it is unclear how this affects the accuracy of the interpolant.
    unsigned short l, m, target;
    SizetList::const_iterator cit = reinterp_indices.begin();
    for (size_t i=0; i<numVars; ++i) {
      if (cit != reinterp_indices.end() && i == *cit) { // reinterpolated index
	l = lev_index[i]; level_to_order(i, l, m);
	target = 2*m - 1; // target m doubles the interp order (t = 2(m-1)+1)
	while (m < target)
	  { ++l; level_to_order(i, l, m); }
	reinterp_lev_index[i] = l; reinterp_quad_order[i] = m;
	// advance to the next reinterp index
	++cit;
      }
      else { // not a reinterpolated index --> no change from reference
	reinterp_lev_index[i] = lev_index[i];
	level_to_order(i, reinterp_lev_index[i], reinterp_quad_order[i]);
      }
    }

    // compute the reinterpolation grid
    compute_tensor_grid(reinterp_quad_order, reinterp_lev_index,
			reinterp_indices, reinterpVarSets.back(),
			reinterpCollocKeys.back());

    // update reiterpMap bookkeeping
    reinterpMap[lev_index] = activeReinterpIndex;
  }
  else
    activeReinterpIndex = map_it->second;
}


const UShortArray& CombinedSparseGridDriver::maximal_grid()
{
  std::map<UShortArray, RealVector>::const_iterator
    w_cit = type1WeightSets.begin(), max_cit = w_cit;
  size_t num_wts, max_wts = w_cit->second.length(); ++w_cit;
  for (; w_cit!=type1WeightSets.end(); ++w_cit) {
    num_wts = w_cit->second.length();
    if (num_wts > max_wts)
      { max_wts = num_wts; max_cit = w_cit; }
  }
  //maximalKey = max_cit->first;
  //return maximalKey;
  return max_cit->first;
}


void CombinedSparseGridDriver::compute_grid()
{
  assign_smolyak_arrays();

  // For efficiency reasons, incremental sparse grid definition uses
  // different point orderings than sgmg/sgmga.  Therefore, the
  // reference grid computations are kept completely separate.

  // ------------------------------------
  // Compute collocation points
  // ------------------------------------
  grid_size(); // ensure active numCollocPts is up to date
  IntArray& unique_index_map = uniqIndMapIter->second;
  compute_unique_points_weights(ssgLevIter->second, anisoWtsIter->second,
				numPtsIter->second, unique_index_map,
				varSetsIter->second, t1WtIter->second,
				t2WtIter->second);

  if (trackCollocDetails) {
    UShort3DArray& colloc_key = collocKeyIter->second;
    assign_collocation_key(smolMIIter->second, colloc_key); // define collocKey
    assign_collocation_indices(colloc_key, unique_index_map,
			       collocIndIter->second);  // define collocIndices
    assign_1d_collocation_points_weights(); // define 1-D point/weight sets
  }

#ifdef DEBUG
  PCout << "CombinedSparseGridDriver::compute_grid() results:\n"
	<< "unique index mapping:\n" << unique_index_map << "\nvariableSets:\n";
  write_data(PCout, varSetsIter->second, false, true, true);
  if (trackUniqueProdWeights) {
    PCout << "\ntype1WeightSets:\n" << t1WtIter->second;
    if (computeType2Weights) {
      PCout << "\ntype2WeightSets:\n";
      write_data(PCout, t2WtIter->second, false, true, true);
    }
  }
#endif
}


void CombinedSparseGridDriver::combine_grid()
{
  // combinedSmolyakMultiIndex is not a Pareto set and will typically have
  // dominated terms with coeffs < 0.  Special care is therefore required
  // with Pareto logic: we start from a valid reference grid, append all
  // non-dominated multi-index "candidates", update the Smolyak coefficients,
  // and then prune terms for which these coefficients are zero.

  // start from first grid (often maximal)
  std::map<UShortArray, UShort2DArray>::const_iterator sm_cit
    = smolyakMultiIndex.begin();
  combinedSmolyakMultiIndex = sm_cit->second;  ++sm_cit;
  UShort2DArray combined_pareto;
  update_pareto_set(combinedSmolyakMultiIndex, combined_pareto);

  // loop over all other level grids, searching for Pareto-optimal indices
  size_t i, num_mi;  bool batch_update;
  for (; sm_cit != smolyakMultiIndex.end(); ++sm_cit) {
    const UShort2DArray& sm_mi = sm_cit->second;  num_mi = sm_mi.size();

    // One pass with batch Pareto updating: ensure coeff<0 sets get included
    batch_update = false;
    for (i=0; i<num_mi; ++i)
      if (new_dominates_any(sm_mi[i], combined_pareto))
	{ combinedSmolyakMultiIndex.push_back(sm_mi[i]); batch_update = true; }
    if (batch_update) update_pareto_set(sm_mi, combined_pareto);

    /*
    // One pass w/ continual Pareto updating: dangerous -> loss of coeff<0 sets
    for (i=0; i<num_mi; ++i)
      if (update_pareto_set(sm_mi[i], combined_pareto))
	combinedSmolyakMultiIndex.push_back(sm_mi[i]);

    // Two pass: dangerous -> sets needed for coeff<0 could be lost
    UShort2DArray pareto_l;
    update_pareto_set(sm_mi, pareto_l);  num_mi = pareto_l.size();
    for (i=0; i<num_mi; ++i)
      if (update_pareto_set(pareto_l[i], combined_pareto))
	combinedSmolyakMultiIndex.push_back(pareto_l[i]);
    */
  }

  // initialize combinedSmolyakCoeffs from first grid and then update for
  // any Pareto-optimal index set additions
  // Note: sandia_sgmgg_coef_inc2() requires that candidates are advancements
  combinedSmolyakCoeffs = smolyakCoeffs.begin()->second;
  update_smolyak_coefficients(combinedSmolyakCoeffs.size(),
			      combinedSmolyakMultiIndex, combinedSmolyakCoeffs);
  // prune inactive index sets
  prune_inactive(combinedSmolyakMultiIndex, combinedSmolyakCoeffs);

  // recompute combinedCollocKey from scratch
  assign_collocation_key(combinedSmolyakMultiIndex, combinedCollocKey);
  // Define combined points and weights to support expectation() calls
  compute_unique_points_weights(combinedSmolyakMultiIndex,
				combinedSmolyakCoeffs,  combinedCollocKey,
				combinedUniqueIndexMap, false, combinedVarSets,
				combinedT1WeightSets,   combinedT2WeightSets);
  // colloc indices are only regenerated for promotions in combined_to_active()
}


void CombinedSparseGridDriver::combined_to_active(bool clear_combined)
{
  // Replace active arrays with combined arrays

  // Update type2 wts even if inactive, so that 2D array sizes are correct
  // Note: inactive weight sets to be removed by clear_inactive()

  if (clear_combined) {
    std::swap(smolMIIter->second,     combinedSmolyakMultiIndex);
    //std::swap(Iter->second, combinedSmolyakMultiIndexMap); // no corresponding
    std::swap(smolCoeffsIter->second, combinedSmolyakCoeffs);
    std::swap(collocKeyIter->second,  combinedCollocKey);
    std::swap(uniqIndMapIter->second, combinedUniqueIndexMap);
    std::swap(varSetsIter->second,    combinedVarSets);
    std::swap(t1WtIter->second,       combinedT1WeightSets);
    std::swap(t2WtIter->second,       combinedT2WeightSets);

    combinedSmolyakMultiIndex.clear();
    //combinedSmolyakMultiIndexMap.clear();
    combinedSmolyakCoeffs.clear();
    combinedCollocKey.clear();
    combinedUniqueIndexMap.clear();
    combinedVarSets.shapeUninitialized(0,0);
    combinedT1WeightSets.sizeUninitialized(0);
    combinedT2WeightSets.shapeUninitialized(0,0);
  }
  else {
    smolMIIter->second     = combinedSmolyakMultiIndex;
    //Iter->second = combinedSmolyakMultiIndexMap;// no corresponding active
    smolCoeffsIter->second = combinedSmolyakCoeffs;
    collocKeyIter->second  = combinedCollocKey;
    uniqIndMapIter->second = combinedUniqueIndexMap;
    varSetsIter->second    = combinedVarSets;
    t1WtIter->second       = combinedT1WeightSets;
    t2WtIter->second       = combinedT2WeightSets;
  }

  // collocation indices are invalidated by expansion combination since the
  // corresponding combined grids involve overlays of data that no longer
  // reflect individual evaluations (to match restoration of collocIndices,
  // new modSurrData must be defined in {NodalInterp,Orthog}PolyApproximation)
  assign_collocation_indices(collocKeyIter->second, uniqIndMapIter->second,
			     collocIndIter->second);
}


void CombinedSparseGridDriver::
compute_unique_points_weights(unsigned short ssg_lev,
			      const RealVector& aniso_wts, int num_colloc_pts,
			      IntArray& unique_index_map, RealMatrix& var_sets,
			      RealVector& t1_wts, RealMatrix& t2_wts)
{
  // ----------------------------------------------
  // Get collocation points and integration weights
  // ----------------------------------------------
  var_sets.shapeUninitialized(numVars, num_colloc_pts);
  if (trackUniqueProdWeights) {
    t1_wts.sizeUninitialized(num_colloc_pts);
    if (computeType2Weights)
      t2_wts.shapeUninitialized(numVars, num_colloc_pts);
  }
  int* sparse_order = new int [num_colloc_pts*numVars];
  int* sparse_index = new int [num_colloc_pts*numVars];
  sgdInstance = this; // sgdInstance required within compute1D fn pointers
  if (aniso_wts.empty()) { // isotropic sparse grid
    int num_total_pts = webbur::sgmg_size_total(numVars, ssg_lev,
      growthRate, &levelGrowthToOrder[0]);
    unique_index_map.resize(num_total_pts);
    webbur::sgmg_unique_index(numVars, ssg_lev, &compute1DPoints[0],
      duplicateTol, num_colloc_pts, num_total_pts, growthRate,
      &levelGrowthToOrder[0], &unique_index_map[0]);
    webbur::sgmg_index(numVars, ssg_lev, num_colloc_pts, num_total_pts,
      &unique_index_map[0], growthRate, &levelGrowthToOrder[0],
      sparse_order, sparse_index);
    webbur::sgmg_point(numVars, ssg_lev, &compute1DPoints[0], num_colloc_pts,
      sparse_order, sparse_index, growthRate, &levelGrowthToOrder[0],
      var_sets.values());
    if (trackUniqueProdWeights) {
      webbur::sgmg_weight(numVars, ssg_lev, &compute1DType1Weights[0],
	num_colloc_pts, num_total_pts, &unique_index_map[0], growthRate,
	&levelGrowthToOrder[0], t1_wts.values());
      if (computeType2Weights) {
	std::vector<CollocFnPtr> comp_1d_t2_wts = compute1DType1Weights;//copy
	RealVector t2_wt_set(num_colloc_pts);
	for (int i=0; i<numVars; ++i) {
	  comp_1d_t2_wts[i] = basis_type2_collocation_weights;//change ith ptr
	  webbur::sgmg_weight(numVars, ssg_lev, &comp_1d_t2_wts[0],
	    num_colloc_pts, num_total_pts, &unique_index_map[0], growthRate,
	    &levelGrowthToOrder[0], t2_wt_set.values());
	  copy_row(t2_wt_set, t2_wts, i);
	  comp_1d_t2_wts[i] = basis_type1_collocation_weights; // restore ptr
	}
      }
    }
  }
  else {
    int num_total_pts = webbur::sandia_sgmga_size_total(numVars,
      aniso_wts.values(), ssg_lev, growthRate, &levelGrowthToOrder[0]);
    unique_index_map.resize(num_total_pts);
    webbur::sandia_sgmga_unique_index(numVars, aniso_wts.values(), ssg_lev,
      &compute1DPoints[0], duplicateTol, num_colloc_pts, num_total_pts,
      growthRate, &levelGrowthToOrder[0], &unique_index_map[0]);
    webbur::sandia_sgmga_index(numVars, aniso_wts.values(), ssg_lev,
      num_colloc_pts, num_total_pts, &unique_index_map[0], growthRate,
      &levelGrowthToOrder[0], sparse_order, sparse_index);
    webbur::sandia_sgmga_point(numVars, aniso_wts.values(), ssg_lev,
      &compute1DPoints[0], num_colloc_pts, sparse_order, sparse_index,
      growthRate, &levelGrowthToOrder[0], var_sets.values());
    if (trackUniqueProdWeights) {
      webbur::sandia_sgmga_weight(numVars, aniso_wts.values(), ssg_lev,
        &compute1DType1Weights[0], num_colloc_pts, num_total_pts,
	&unique_index_map[0], growthRate, &levelGrowthToOrder[0],
	t1_wts.values());
      if (computeType2Weights) {
	std::vector<CollocFnPtr> comp_1d_t2_wts = compute1DType1Weights;//copy
	RealVector t2_wt_set(num_colloc_pts);
	for (int i=0; i<numVars; ++i) {
	  comp_1d_t2_wts[i] = basis_type2_collocation_weights;//change ith ptr
	  webbur::sandia_sgmga_weight(numVars, aniso_wts.values(), ssg_lev,
	    &comp_1d_t2_wts[0], num_colloc_pts, num_total_pts,
	    &unique_index_map[0], growthRate, &levelGrowthToOrder[0],
	    t2_wt_set.values());
	  copy_row(t2_wt_set, t2_wts, i);
	  comp_1d_t2_wts[i] = basis_type1_collocation_weights; // restore ptr
	}
      }
    }
  }
  delete [] sparse_order;
  delete [] sparse_index;
}


void CombinedSparseGridDriver::
compute_unique_points_weights(const UShort2DArray& sm_mi,
			      const IntArray& sm_coeffs,
			      const UShort3DArray& colloc_key,
			      Sizet2DArray& colloc_ind, int& num_colloc_pts,
			      RealMatrix& a1_pts, RealVector& a1_t1w,
			      RealMatrix& a1_t2w, RealVector& zv,
			      RealVector& r1v, IntArray& sind1, BitArray& isu1,
			      IntArray& uind1, IntArray& uset1, int& num_u1,
			      IntArray& unique_index_map,
			      bool update_1d_pts_wts, RealMatrix& var_sets,
			      RealVector& t1_wts, RealMatrix& t2_wts)
{
  // define a1 pts/wts
  compute_tensor_points_weights(sm_mi, colloc_key, 0, sm_mi.size(),
				update_1d_pts_wts, a1_pts, a1_t1w, a1_t2w);
  // ----
  // INC1
  // ----
  int m = numVars, n1 = a1_pts.numCols(), seed = 1234567;
  zv.sizeUninitialized(m);  r1v.sizeUninitialized(n1);
  sind1.resize(n1);  uind1.resize(n1);
  uset1.resize(n1); // numUnique1 if count_inc1 used
  bool* is_unique1 = new bool[n1];

  webbur::point_radial_tol_unique_index_inc1(m, n1, a1_pts.values(),
    duplicateTol, &seed, zv.values(), r1v.values(), &sind1[0], is_unique1,
    &num_u1, &uset1[0], &uind1[0]);

  copy_data(is_unique1, n1, isu1);
  delete [] is_unique1;

#ifdef DEBUG
  PCout << "Reference unique: numUnique1 = " << num_u1 << "\na1 =\n";
  write_data(PCout, a1_pts, false, true, true);
  PCout << "               r1   indx1 unique1   undx1   xdnu1:\n";
  for (size_t i=0; i<num_u1; ++i)
    PCout << std::setw(17) << r1v[i]   << std::setw(8) << sind1[i]
	  << std::setw(8)  << isu1[i]  << std::setw(8) << uset1[i]
	  << std::setw(8)  << uind1[i] << '\n';
  for (size_t i=num_u1; i<n1; ++i)
    PCout << std::setw(17) << r1v[i]  << std::setw(8)  << sind1[i]
	  << std::setw(8)  << isu1[i] << std::setw(16) << uind1[i] << '\n';
  PCout << std::endl;
#endif // DEBUG

  num_colloc_pts = num_u1;
  assign_unique_indices(isu1, uind1, uset1, unique_index_map);
  assign_collocation_indices(colloc_key, unique_index_map, colloc_ind);
  assign_sparse_points(colloc_ind, 0, isu1, 0, a1_pts, var_sets);
  if (trackUniqueProdWeights)
    assign_sparse_weights(colloc_key, colloc_ind, num_colloc_pts, sm_coeffs,
			  a1_t1w, a1_t2w, t1_wts, t2_wts);
}


void CombinedSparseGridDriver::
compute_tensor_points_weights(const UShort2DArray& sm_mi,
			      const UShort3DArray& colloc_key,
			      size_t start_index, size_t num_indices,
			      bool update_1d_pts_wts, RealMatrix& pts,
			      RealVector& t1_wts, RealMatrix& t2_wts)
{
  // Requirements: updated sm_mi,colloc_key for [start,start+num_indices].
  // 1D Pts/Wts will be updated as indicated by update_1d_pts_wts

  size_t i, j, k, l, cntr, num_tp_pts, num_colloc_pts = 0,
    end = start_index + num_indices;
  // define num_colloc_pts
  for (i=start_index; i<end; ++i)
    num_colloc_pts += colloc_key[i].size();
  // define pts/wts: wts are raw product weights; Smolyak combinatorial
  // coefficient applied in compute_grid()/compute_trial_grid()
  pts.shapeUninitialized(numVars, num_colloc_pts);
  t1_wts.sizeUninitialized(num_colloc_pts);
  if (computeType2Weights)
    t2_wts.shapeUninitialized(numVars, num_colloc_pts);
  for (i=start_index, cntr=0; i<end; ++i) {
    const UShortArray& sm_index = sm_mi[i];
    if (update_1d_pts_wts) { // update collocPts1D, {type1,type2}CollocWts1D
      UShortArray quad_order(numVars);
      level_to_order(sm_index, quad_order);
      update_1d_collocation_points_weights(quad_order, sm_index);
    }
    num_tp_pts = colloc_key[i].size();
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      const UShortArray& key_ij = colloc_key[i][j];
      Real* pt    =    pts[cntr]; // column vector
      Real& t1_wt = t1_wts[cntr]; t1_wt = 1.;
      for (k=0; k<numVars; ++k) {
	pt[k]  =      collocPts1D[sm_index[k]][k][key_ij[k]];
	t1_wt *= type1CollocWts1D[sm_index[k]][k][key_ij[k]];
      }
      if (computeType2Weights) {
	Real* t2_wt = t2_wts[cntr]; // column vector
	for (k=0; k<numVars; ++k) {
	  Real& t2_wt_k = t2_wt[k]; t2_wt_k = 1.;
	  for (l=0; l<numVars; ++l)
	    t2_wt_k *= (l==k) ? type2CollocWts1D[sm_index[l]][l][key_ij[l]] :
	                        type1CollocWts1D[sm_index[l]][l][key_ij[l]];
	}
      }
    }
  }
#ifdef DEBUG
  PCout << "Tensor product weights =\ntype1:\n" << t1_wts;
  if (computeType2Weights)
    { PCout << "type2:\n"; write_data(PCout, t2_wts, false, true, true); }
#endif // DEBUG
}


void CombinedSparseGridDriver::
assign_unique_indices(const BitArray& isu1, const IntArray& xdnu1,
		      const IntArray& undx1, IntArray& unique_index_map)
{
  //const BitArray& isu1  = isUniq1Iter->second;
  //const IntArray& xdnu1 = uniqInd1Iter->second;
  //const IntArray& undx1 = uniqSet1Iter->second;
  //IntArray& unique_index_map = uniqIndMapIter->second;

  size_t i, n1 = xdnu1.size();
  int xdnu_j, a1_index, new_cntr = 0;

  unique_index_map.resize(n1);

  // first pass assigns unique indices
  for (i=0; i<n1; ++i)
    if (isu1[i])
      unique_index_map[i] = new_cntr++;
  // second pass refers back to unique indices and can be a forward reference
  // (dictating two passes)
  for (i=0; i<n1; ++i)
    if (!isu1[i]) {
      // XDNU1[N1] in point_radial_tol_unique_index_inc1() [sandia_rules.cpp]:
      //   the index, in UNDX1, of the tolerably unique point that
      //   "represents" this point.
      xdnu_j = xdnu1[i];
      // UNDX1[UNIQUE_NUM1] in point_radial_tol_unique_index_inc1():
      //   the index, in A1, of the tolerably unique points.
      a1_index = undx1[xdnu_j]; // appears to reproduce i
      unique_index_map[i] = unique_index_map[a1_index];
    }

#ifdef DEBUG
  PCout << "Reference map:\n" << unique_index_map;
#endif // DEBUG
}


void CombinedSparseGridDriver::
assign_sparse_points(const Sizet2DArray& colloc_ind, size_t start_index,
		     const BitArray& raw_is_unique,
		     size_t curr_unique_pts, // 0 (ref) or num_u1 (incr)
		     const RealMatrix& raw_pts, RealMatrix& unique_pts)
{
  // update sizes
  size_t new_unique_pts = raw_is_unique.count();
  unique_pts.reshape(numVars, new_unique_pts + curr_unique_pts);

  size_t i, j, cntr = 0, uniq_index, num_sm_mi = colloc_ind.size(), num_tp_pts;
  for (i=start_index; i<num_sm_mi; ++i) {
    const SizetArray& colloc_ind_i = colloc_ind[i];
    num_tp_pts = colloc_ind_i.size();
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      if (raw_is_unique[cntr]) {
	uniq_index = colloc_ind_i[j];
	copy_data(raw_pts[cntr], numVars, unique_pts[uniq_index]);
      }
    }
  }
}


void CombinedSparseGridDriver::
assign_sparse_weights(const UShort3DArray& colloc_key,
		      const Sizet2DArray& colloc_ind, int num_colloc_pts,
		      const IntArray& sm_coeffs, const RealVector& a1_t1_wts,
		      const RealMatrix& a1_t2_wts, RealVector& unique_t1_wts,
		      RealMatrix& unique_t2_wts)
{
  // update sizes
  unique_t1_wts.size(num_colloc_pts); // init to 0
  if (computeType2Weights)
    unique_t2_wts.shape(numVars, num_colloc_pts); // init to 0

  int uniq_index, delta_coeff, sm_coeff;
  // add contributions for new index sets
  add_sparse_weights(0, colloc_key, colloc_ind, sm_coeffs, a1_t1_wts,
		     a1_t2_wts, unique_t1_wts, unique_t2_wts);

#ifdef DEBUG
  PCout << "reference type1 weight sets:\n" << unique_t1_wts;
  if (computeType2Weights) {
    PCout << "reference type2 weight sets:\n";
    write_data(PCout, unique_t2_wts, false, true, true);
  }
#endif // DEBUG
}


void CombinedSparseGridDriver::
add_sparse_weights(size_t start_index, const UShort3DArray& colloc_key,
		   const Sizet2DArray& colloc_ind, const IntArray& sm_coeffs,
		   const RealVector& raw_t1w, const RealMatrix& raw_t2w,
		   RealVector& unique_t1w, RealMatrix& unique_t2w)
{
  // add contributions for new index sets
  size_t i, j, k, num_sm_mi = colloc_key.size(), uniq_index, num_tp_pts, cntr;
  for (i=start_index, cntr=0; i<num_sm_mi; ++i) {
    int sm_coeff = sm_coeffs[i];
    if (sm_coeff) {
      num_tp_pts = colloc_key[i].size();
      const SizetArray& colloc_ind_i = colloc_ind[i];
      for (j=0; j<num_tp_pts; ++j, ++cntr) {
	uniq_index = colloc_ind_i[j];
	// assign tensor weights to unique weights
	unique_t1w[uniq_index] += sm_coeff * raw_t1w[cntr];
	if (computeType2Weights) {
	  Real*  uniq_t2w_j = unique_t2w[uniq_index];
	  const Real* t2w_j =    raw_t2w[cntr];
	  for (k=0; k<numVars; ++k)
	    uniq_t2w_j[k] += sm_coeff * t2w_j[k];
	}
      }
    }
    else
      cntr += colloc_key[i].size();
  }
}

} // namespace Pecos
