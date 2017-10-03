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
#include "DistributionParams.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: CombinedSparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {

/// initialize static member pointer to active driver instance
CombinedSparseGridDriver* CombinedSparseGridDriver::sgdInstance(NULL);


void CombinedSparseGridDriver::
initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
		const ShortArray& u_types,
		const ExpansionConfigOptions& ec_options,
		BasisConfigOptions& bc_options, short growth_rate,
		bool track_colloc, bool track_uniq_prod_wts)
{
  SparseGridDriver::initialize_grid(ssg_level, dim_pref, u_types, ec_options,
				    bc_options, growth_rate);
  trackCollocDetails     = track_colloc;
  trackUniqueProdWeights = track_uniq_prod_wts;

  // set a rule-dependent duplicateTol
  initialize_duplicate_tolerance();
  // set compute1D{Points,Type1Weights,Type2Weights}
  initialize_rule_pointers();
  // set levelGrowthToOrder
  initialize_growth_pointers();
}


void CombinedSparseGridDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis)
{
  IntegrationDriver::initialize_grid(poly_basis);

  // set a rule-dependent duplicateTol
  initialize_duplicate_tolerance();
  // set compute1D{Points,Type1Weights,Type2Weights}
  initialize_rule_pointers();
  // set levelGrowthToOrder
  initialize_growth_pointers();
}


void CombinedSparseGridDriver::store_grid(size_t index)
{
  size_t stored_len = storedType1WeightSets.size();
  if (index == _NPOS || index == stored_len) { // append
    storedLevMultiIndex.push_back(smolyakMultiIndex);
    storedLevCoeffs.push_back(smolyakCoeffs);
    storedCollocKey.push_back(collocKey);
    storedCollocIndices.push_back(collocIndices);
    storedType1WeightSets.push_back(type1WeightSets);
    storedType2WeightSets.push_back(type2WeightSets);
  }
  else if (index < stored_len) { // replace
    storedLevMultiIndex[index]   = smolyakMultiIndex;
    storedLevCoeffs[index]       = smolyakCoeffs;
    storedCollocKey[index]       = collocKey;
    storedCollocIndices[index]   = collocIndices;
    storedType1WeightSets[index] = type1WeightSets;
    storedType2WeightSets[index] = type2WeightSets;
  }
  else {
    PCerr << "Error: bad index (" << index << ") passed in "
	  << "CombinedSparseGridDriver::store_grid()" << std::endl;
    abort_handler(-1);
  }
}


void CombinedSparseGridDriver::restore_grid(size_t index)
{
  size_t stored_len = storedType1WeightSets.size();
  if (index == _NPOS) {
    smolyakMultiIndex = storedLevMultiIndex.back();
    smolyakCoeffs   = storedLevCoeffs.back();
    collocKey       = storedCollocKey.back();
    collocIndices   = storedCollocIndices.back();
    type1WeightSets = storedType1WeightSets.back();
    type2WeightSets = storedType2WeightSets.back();
  }
  else if (index < stored_len) {
    smolyakMultiIndex = storedLevMultiIndex[index];
    smolyakCoeffs   = storedLevCoeffs[index];
    collocKey       = storedCollocKey[index];
    collocIndices   = storedCollocIndices[index];
    type1WeightSets = storedType1WeightSets[index];
    type2WeightSets = storedType2WeightSets[index];
  }
  else {
    PCerr << "Error: bad index (" << index << ") passed in "
	  << "CombinedSparseGridDriver::restore_grid()" << std::endl;
    abort_handler(-1);
  }
}


void CombinedSparseGridDriver::remove_stored_grid(size_t index)
{
  size_t stored_len = storedType1WeightSets.size();
  if (index == _NPOS || index == stored_len) {
    storedLevMultiIndex.pop_back();
    storedLevCoeffs.pop_back();
    storedCollocKey.pop_back();
    storedCollocIndices.pop_back();
    storedType1WeightSets.pop_back();
    storedType2WeightSets.pop_back();
  }
  else if (index < stored_len) {
    UShort3DArray::iterator u3it = storedLevMultiIndex.begin();
    std::advance(u3it, index); storedLevMultiIndex.erase(u3it);
    Int2DArray::iterator i2it = storedLevCoeffs.begin();
    std::advance(i2it, index); storedLevCoeffs.erase(i2it);
    UShort4DArray::iterator u4it = storedCollocKey.begin();
    std::advance(u4it, index); storedCollocKey.erase(u4it);
    Sizet3DArray::iterator s3it = storedCollocIndices.begin();
    std::advance(s3it, index); storedCollocIndices.erase(s3it);
    RealVectorArray::iterator vit = storedType1WeightSets.begin();
    std::advance(vit, index); storedType1WeightSets.erase(vit);
    RealMatrixArray::iterator mit = storedType2WeightSets.begin();
    std::advance(mit, index); storedType2WeightSets.erase(mit);
  }
}


void CombinedSparseGridDriver::clear_stored()
{
  storedLevMultiIndex.clear(); storedLevCoeffs.clear();
  storedCollocKey.clear();     storedCollocIndices.clear();

  storedType1WeightSets.clear();
  if (computeType2Weights) storedType2WeightSets.clear();
}


size_t CombinedSparseGridDriver::maximal_grid() const
{
  size_t i, num_stored = storedType1WeightSets.size(),
    max_index = _NPOS, max_wts = type1WeightSets.length();
  for (i=0; i<num_stored; ++i)
    if (storedType1WeightSets[i].length() > max_wts)
      { max_index = i; max_wts = storedType1WeightSets[i].length(); }
  return max_index;
}


void CombinedSparseGridDriver::swap_grid(size_t index)
{
  std::swap(storedLevMultiIndex[index], smolyakMultiIndex);
  std::swap(storedLevCoeffs[index],     smolyakCoeffs);
  std::swap(storedCollocKey[index],     collocKey);
  std::swap(storedCollocIndices[index], collocIndices);

  RealVector tmp_vec(type1WeightSets);
  type1WeightSets = storedType1WeightSets[index];
  storedType1WeightSets[index] = tmp_vec;

  if (computeType2Weights) {
    RealMatrix tmp_mat(type2WeightSets);
    type2WeightSets = storedType2WeightSets[index];
    storedType2WeightSets[index] = tmp_mat;
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
  if (refineControl != DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
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
assign_smolyak_arrays(UShort2DArray& multi_index, IntArray& coeffs)
{
  // Populate smolyakMultiIndex and smolyakCoeffs.  Identifies
  // use of polynomialBasis[variable][index] based on index 0:num_levels-1.
  // w = q - N = dimension-independent level.  For isotropic,
  //   w + 1 <= |i| <= w + N for i starts at 1 (used for index set defn.)
  //   w - N + 1 <= |j| <= w for j = i - 1 starts at 0 (used for generation)
  // For anisotropic, a weighted linear index set constraint is used.

  size_t i;
  if (dimIsotropic) { // initialize multi_index
    UShortArray levels(numVars, ssgLevel);
    SharedPolyApproxData::total_order_multi_index(levels, multi_index,
						  numVars-1);
    size_t num_terms = multi_index.size();
    // initialize coeffs
    coeffs.resize(num_terms);
    for (i=0; i<num_terms; i++) {
      int wpNmi = ssgLevel - l1_norm(multi_index[i]); // w+N-|i| = w-|j|
      coeffs[i] = (int)std::pow(-1., wpNmi)
	* (int)std::floor(BasisPolynomial::n_choose_k(numVars - 1, wpNmi)+.5);
    }
  }
  else { // utilize webbur::sgmga_vcn_{ordered,coef}
    multi_index.clear();
    coeffs.clear();
    // Utilize webbur::sandia_sgmga_vcn_{ordered,coef} for 0-based index sets
    // (w*alpha_min-|alpha| < |alpha . j| <= w*alpha_min).
    // With scaling alpha_min = 1: w-|alpha| < |alpha . j| <= w.
    // In the isotropic case, reduces to w-N < |j| <= w, which is the same as
    // w-N+1 <= |j| <= w.
    IntArray x(numVars), x_max(numVars); //x_max = ssgLevel;
    UShortArray index_set(numVars);
    Real wt_sum = 0., q_max = ssgLevel;
    for (i=0; i<numVars; ++i) {
      const Real& wt_i = anisoLevelWts[i];
      wt_sum += wt_i;
      // minimum nonzero weight is scaled to 1, so just catch special case of 0
      x_max[i] = (wt_i > 1.e-10) ? (int)std::ceil(q_max/wt_i) : 0;
    }
    Real q_min = ssgLevel - wt_sum;
#ifdef DEBUG
    PCout << "q_min = " << q_min << " q_max = " << q_max;
#endif // DEBUG

    bool more = false;
    Real *aniso_wts = anisoLevelWts.values();
    int  *x0 = &x[0], *xm0 = &x_max[0], coeff;
    webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wts, xm0, x0,
				     q_min, q_max, &more);
    while (more) {
      coeff = (int)webbur::sandia_sgmga_vcn_coef(numVars, aniso_wts, x0, q_max);
      if (coeff) {
	coeffs.push_back(coeff);
	for (i=0; i<numVars; ++i)
	  index_set[i] = (unsigned short)x[i];
	multi_index.push_back(index_set);
      }
      webbur::sandia_sgmga_vcn_ordered(numVars, aniso_wts, xm0, x0,
				       q_min, q_max, &more);
    }
  }

#ifdef DEBUG
  size_t num_terms = coeffs.size();
  PCout << "\nnum Smolyak terms = " << num_terms << '\n';
  for (i=0; i<num_terms; i++)
    PCout << "multi_index[" << i << "]:\n" << multi_index[i]
	  << "coeffs[" << i << "] = " << coeffs[i] << "\n\n";
#endif // DEBUG
}


void CombinedSparseGridDriver::
update_smolyak_coefficients(size_t start_index,
			    const UShort2DArray& multi_index, IntArray& coeffs)
{
  size_t j, cntr = 0, num_sets = multi_index.size(), len1 = num_sets-1;
  int i, m = numVars;
  if (coeffs.size() != num_sets)
    coeffs.resize(num_sets);
  int *s1 = new int [numVars*len1], *c1 = new int [len1],
      *s2 = new int [numVars];
  // initialize s1 and c1
  for (i=0; i<start_index; ++i) {
    c1[i] = coeffs[i];
    for (j=0; j<numVars; ++j, ++cntr) // no copy_data() since ushort -> int
      s1[cntr] = multi_index[i][j]; // sgmgg packs by variable groups
  }
  // for each s2, update coeffs
  for (i=start_index; i<num_sets; ++i) {
    for (j=0; j<numVars; ++j) // no copy_data() since ushort -> int
      s2[j] = multi_index[i][j];
    webbur::sandia_sgmgg_coef_inc2(m, i, s1, c1, s2, &coeffs[0]);
#ifdef DEBUG
    PCout << "update_smolyak_coefficients(): updated Smolyak coeffs =\n"
	  << coeffs << '\n';
#endif // DEBUG
    if (i<num_sets-1) { // if not last i, update s1 and c1 state for next pass
      for (j=0; j<=i; ++j) // coeffs updated to len i+1; max len = num_sets-1
	c1[j] = coeffs[j];
      for (j=0; j<numVars; ++j, ++cntr)
	s1[cntr] = s2[j]; // max len = (num_sets-1)*numVars
    }
  }
  delete [] s1; delete [] c1; delete [] s2;
}


void CombinedSparseGridDriver::assign_collocation_key()
{
  // define mapping from 1:numCollocPts to set of 1d interpolation indices
  size_t i, num_smolyak_indices = smolyakMultiIndex.size();
  collocKey.resize(num_smolyak_indices);
  UShortArray quad_order(numVars); //, collocation_indices(numVars);
  for (i=0; i<num_smolyak_indices; ++i) {
    level_to_order(smolyakMultiIndex[i], quad_order);
    SharedPolyApproxData::tensor_product_multi_index(quad_order, collocKey[i],
						     false);
  }
}


void CombinedSparseGridDriver::update_collocation_key(size_t start_index)
{
  UShortArray quad_order(numVars);
  size_t i, num_sm_mi = smolyakMultiIndex.size();
  collocKey.resize(num_sm_mi);
  for (i=start_index; i<num_sm_mi; ++i) {
    level_to_order(smolyakMultiIndex[i], quad_order);
    SharedPolyApproxData::tensor_product_multi_index(quad_order, collocKey[i],
						     false);
  }
}


void CombinedSparseGridDriver::assign_collocation_indices()
{
  // define mapping from 1:numCollocPts to set of 1d interpolation indices
  size_t i, j, num_tp_pts, cntr = 0,
    num_smolyak_indices = collocKey.size();
  collocIndices.resize(num_smolyak_indices);
  for (i=0; i<num_smolyak_indices; ++i) {
    num_tp_pts = collocKey[i].size();
    SizetArray& indices_i = collocIndices[i];
    indices_i.resize(num_tp_pts);
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      indices_i[j] = uniqueIndexMapping[cntr];
#ifdef DEBUG
      PCout << "collocKey[" << i << "][" << j << "]:\n" << collocKey[i][j]
	    << "collocIndices[" << i << "][" << j << "] = " << indices_i[j]
	    << '\n';
#endif // DEBUG
    }
  }
}


int CombinedSparseGridDriver::grid_size()
{
  if (updateGridSize) {
    sgdInstance = this; // sgdInstance required within compute1DPoints below
    numCollocPts = (dimIsotropic) ?
      webbur::sgmg_size(numVars, ssgLevel, &compute1DPoints[0], duplicateTol,
	growthRate, &levelGrowthToOrder[0]) :
      webbur::sandia_sgmga_size(numVars, anisoLevelWts.values(), ssgLevel,
	&compute1DPoints[0], duplicateTol, growthRate, &levelGrowthToOrder[0]);
    updateGridSize = false;
  }
  return numCollocPts;
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


void CombinedSparseGridDriver::compute_grid(RealMatrix& var_sets)
{
  assign_smolyak_arrays();

  // For efficiency reasons, incremental sparse grid definition uses
  // different point orderings than sgmg/sgmga.  Therefore, the
  // reference grid computations are kept completely separate.

  if (refineControl == DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
    // compute reference grid only
    assign_collocation_key();               // compute collocKey
    assign_1d_collocation_points_weights(); // define 1-D point/weight sets
    reference_unique(var_sets);             // define reference grid
  }
  else { // compute reference and any refined grids
    // ------------------------------------
    // Compute number of collocation points
    // ------------------------------------
    grid_size(); // ensure numCollocPts is up to date

    // ----------------------------------------------
    // Get collocation points and integration weights
    // ----------------------------------------------
    var_sets.shapeUninitialized(numVars, numCollocPts);
    if (trackUniqueProdWeights) {
      type1WeightSets.sizeUninitialized(numCollocPts);
      if (computeType2Weights)
	type2WeightSets.shapeUninitialized(numVars, numCollocPts);
    }
    int* sparse_order = new int [numCollocPts*numVars];
    int* sparse_index = new int [numCollocPts*numVars];
    sgdInstance = this; // sgdInstance required within compute1D fn pointers
    if (dimIsotropic) {
      int num_total_pts = webbur::sgmg_size_total(numVars, ssgLevel,
	growthRate, &levelGrowthToOrder[0]);
      uniqueIndexMapping.resize(num_total_pts);
      webbur::sgmg_unique_index(numVars, ssgLevel, &compute1DPoints[0],
	duplicateTol, numCollocPts, num_total_pts, growthRate,
	&levelGrowthToOrder[0], &uniqueIndexMapping[0]);
      webbur::sgmg_index(numVars, ssgLevel, numCollocPts, num_total_pts,
	&uniqueIndexMapping[0], growthRate, &levelGrowthToOrder[0],
	sparse_order, sparse_index);
      webbur::sgmg_point(numVars, ssgLevel, &compute1DPoints[0], numCollocPts,
	sparse_order, sparse_index, growthRate, &levelGrowthToOrder[0],
	var_sets.values());
      if (trackUniqueProdWeights) {
	webbur::sgmg_weight(numVars, ssgLevel, &compute1DType1Weights[0],
	  numCollocPts, num_total_pts, &uniqueIndexMapping[0], growthRate,
	  &levelGrowthToOrder[0], type1WeightSets.values());
	if (computeType2Weights) {
	  std::vector<CollocFnPtr> comp_1d_t2_wts = compute1DType1Weights;//copy
	  RealVector t2_wt_set(numCollocPts);
	  for (int i=0; i<numVars; ++i) {
	    comp_1d_t2_wts[i] = basis_type2_collocation_weights;//change ith ptr
	    webbur::sgmg_weight(numVars, ssgLevel, &comp_1d_t2_wts[0],
	      numCollocPts, num_total_pts, &uniqueIndexMapping[0], growthRate,
	      &levelGrowthToOrder[0], t2_wt_set.values());
	    copy_row(t2_wt_set, type2WeightSets, i);
	    comp_1d_t2_wts[i] = basis_type1_collocation_weights; // restore ptr
	  }
	}
      }
    }
    else {
      int num_total_pts = webbur::sandia_sgmga_size_total(numVars,
	anisoLevelWts.values(), ssgLevel, growthRate, &levelGrowthToOrder[0]);
      uniqueIndexMapping.resize(num_total_pts);
      webbur::sandia_sgmga_unique_index(numVars, anisoLevelWts.values(),
	ssgLevel, &compute1DPoints[0], duplicateTol, numCollocPts,num_total_pts,
	growthRate, &levelGrowthToOrder[0], &uniqueIndexMapping[0]);
      webbur::sandia_sgmga_index(numVars, anisoLevelWts.values(), ssgLevel,
        numCollocPts, num_total_pts, &uniqueIndexMapping[0], growthRate,
	&levelGrowthToOrder[0], sparse_order, sparse_index);
      webbur::sandia_sgmga_point(numVars, anisoLevelWts.values(), ssgLevel,
        &compute1DPoints[0], numCollocPts, sparse_order, sparse_index,
	growthRate, &levelGrowthToOrder[0], var_sets.values());
      if (trackUniqueProdWeights) {
	webbur::sandia_sgmga_weight(numVars, anisoLevelWts.values(), ssgLevel,
          &compute1DType1Weights[0], numCollocPts, num_total_pts,
	  &uniqueIndexMapping[0], growthRate, &levelGrowthToOrder[0],
	  type1WeightSets.values());
	if (computeType2Weights) {
	  std::vector<CollocFnPtr> comp_1d_t2_wts = compute1DType1Weights;//copy
	  RealVector t2_wt_set(numCollocPts);
	  for (int i=0; i<numVars; ++i) {
	    comp_1d_t2_wts[i] = basis_type2_collocation_weights;//change ith ptr
	    webbur::sandia_sgmga_weight(numVars, anisoLevelWts.values(),
	      ssgLevel, &comp_1d_t2_wts[0], numCollocPts, num_total_pts,
	      &uniqueIndexMapping[0], growthRate, &levelGrowthToOrder[0],
	      t2_wt_set.values());
	    copy_row(t2_wt_set, type2WeightSets, i);
	    comp_1d_t2_wts[i] = basis_type1_collocation_weights; // restore ptr
	  }
	}
      }
    }
    delete [] sparse_order;
    delete [] sparse_index;

    if (trackCollocDetails) {
      assign_collocation_key();               // compute collocKey
      assign_collocation_indices();           // compute collocIndices
      assign_1d_collocation_points_weights(); // define 1-D point/weight sets
    }
  }

#ifdef DEBUG
  PCout << "CombinedSparseGridDriver::compute_grid() results:\n"
	<< "uniqueIndexMapping:\n" << uniqueIndexMapping << "\nvar_sets:\n";
  write_data(PCout, var_sets, false, true, true);
  if (trackUniqueProdWeights) {
    PCout << "\ntype1WeightSets:\n"; write_data(PCout, type1WeightSets);
    if (computeType2Weights) {
      PCout << "\ntype2WeightSets:\n";
      write_data(PCout, type2WeightSets, false, true, true);
    }
  }
#endif
}


void CombinedSparseGridDriver::compute_trial_grid(RealMatrix& var_sets)
{
  // compute trial variable/weight sets and update collocKey
  UShortArray quad_order(numVars);
  level_to_order(trialSet, quad_order);
  UShort2DArray new_key;
  collocKey.push_back(new_key); // empty array updated in place
  compute_tensor_grid(quad_order, trialSet, a2Points, a2Type1Weights,
		      a2Type2Weights, collocKey.back());

  // track trial sets that have been evaluated (do here since
  // push_trial_set() used for both new trials and restorations)
  computedTrialSets.insert(trialSet);

  // update collocIndices, uniqueIndexMapping, and var_sets,
  // but don't recompute a2 data
  increment_unique(false, true, var_sets);

#ifdef DEBUG
  PCout << "compute_trial_grid() increment:\nunique variable sets:\n"
	<< var_sets;
#endif // DEBUG
}


void CombinedSparseGridDriver::compute_grid_increment(RealMatrix& var_sets)
{
  // TO DO: employ increment_unique for each set in incrementSets

  // build from scratch for now
  compute_grid(var_sets);
}


void CombinedSparseGridDriver::initialize_sets()
{
  // provide a helpful error message in the case where refineControl is not
  // set for generalized adaptation
  // > this traps the error where the reference grid was computed
  //   inconsistently (see logic in compute_grid())
  if (refineControl != DIMENSION_ADAPTIVE_CONTROL_GENERALIZED) {
    PCerr << "Error: CombinedSparseGridDriver::initialize_sets() called for "
	  << "inconsistent refinement control setting (" << refineControl
	  << ")." << std::endl;
    abort_handler(-1);
  }

  // define set O (old) from smolyakMultiIndex and smolyakCoeffs:
  //oldMultiIndex = smolyakMultiIndex;
  oldMultiIndex.clear();
  oldMultiIndex.insert(smolyakMultiIndex.begin(), smolyakMultiIndex.end());
  update_reference();

  // computedTrialSets no longer cleared in finalize_sets(), so do on init
  computedTrialSets.clear();

  // compute initial set A (active) by applying add_active_neighbors()
  // to the frontier of smolyakMultiIndex:
  size_t i, num_old_sets = smolyakCoeffs.size();
  // anisotropic test on coeff==1 is necessary but not sufficient for presence
  // on index set frontier, requiring an additional logic test within
  // add_active_neighbors().  For anisotropic, the weighted norm of the index
  // set may differ from the level --> need to compute Pareto set.
  for (i=0; i<num_old_sets; ++i)
    if ( smolyakCoeffs[i] == 1 && ( !dimIsotropic || // imperfect for aniso
	 ( dimIsotropic && l1_norm(smolyakMultiIndex[i]) == ssgLevel ) ) )
      add_active_neighbors(smolyakMultiIndex[i], dimIsotropic);

#ifdef DEBUG
  PCout << "CombinedSparseGridDriver::initialize_sets():\nold sets:\n"
	<< oldMultiIndex << "active sets:\n" << activeMultiIndex << std::endl;
#endif // DEBUG
}


void CombinedSparseGridDriver::push_trial_set(const UShortArray& set)
{
  trialSet = set;
  size_t last_index = smolyakMultiIndex.size();
  smolyakMultiIndex.push_back(set);

  // update smolyakCoeffs from smolyakMultiIndex
  update_smolyak_coefficients(last_index);

  // collocKey, collocIndices, and uniqueIndexMapping updated within
  // either restore_set() or compute_trial_grid()
}


void CombinedSparseGridDriver::restore_set()
{
  // SparseGridDriver currently retains no memory, so updates are recomputed

  size_t last_index = smolyakMultiIndex.size() - 1;
  // update collocKey
  update_collocation_key(last_index);
  // compute a2; update collocIndices & uniqueIndexMapping; don't update pts/wts
  RealMatrix dummy_set;
  increment_unique(true, false, dummy_set);
}


void CombinedSparseGridDriver::pop_trial_set()
{
  smolyakMultiIndex.pop_back();
  collocKey.pop_back(); collocIndices.pop_back();
  smolyakCoeffs = smolyakCoeffsRef;

  numCollocPts -= numUnique2; // subtract number of trial points
  uniqueIndexMapping.resize(numCollocPts); // prune trial set from end
}


void CombinedSparseGridDriver::
finalize_sets(bool output_sets, bool converged_within_tol)
{
  // For final answer, push all evaluated sets into old and clear active.
  // Multiple trial insertion approach must be compatible with bookkeeping
  // elsewhere (e.g., Dakota::Approximation), i.e., inc2/inc3 set insertions
  // occur one at a time without mixing.

  size_t start_index = smolyakMultiIndex.size();
  // don't insert activeMultiIndex, as this may include sets which have not
  // been evaluated (due to final update_sets() call); use computedTrialSets
  smolyakMultiIndex.insert(smolyakMultiIndex.end(), computedTrialSets.begin(),
			   computedTrialSets.end());
  activeMultiIndex.clear();
  // defer since needed for SharedPolyApproxData::finalization_index()
  //computedTrialSets.clear();

  // update smolyakCoeffs from smolyakMultiIndex
  update_smolyak_coefficients(start_index);
  // update collocKey
  update_collocation_key(start_index);
  // update a2 data, uniqueIndexMapping, collocIndices, numCollocPts
  finalize_unique(start_index);// assure no mixing of discrete a2's
  //merge_unique(); // a1 reference update not needed, no addtnl increments
  //update_reference();

  if (output_sets) {
    size_t i, j, num_sm_mi = smolyakMultiIndex.size();
    if (converged_within_tol) {
      PCout << "Above tolerance index sets:\n";
      size_t last = start_index - 1;
      for (i=0; i<last; ++i)
	print_index_set(PCout, smolyakMultiIndex[i]);
      PCout << "Below tolerance index sets:\n";
      for (i=last; i<num_sm_mi; ++i)
	print_index_set(PCout, smolyakMultiIndex[i]);
    }
    else {
      PCout << "Final index sets:\n";
      for (i=0; i<num_sm_mi; ++i)
	print_index_set(PCout, smolyakMultiIndex[i]);
    }
  }
}


void CombinedSparseGridDriver::reference_unique(RealMatrix& var_sets)
{
  // define a1 pts/wts
  size_t num_sm_mi = smolyakMultiIndex.size();
  compute_tensor_points_weights(0, num_sm_mi, a1Points, a1Type1Weights,
				a1Type2Weights);

  // ----
  // INC1
  // ----
  int m = numVars, n1 = a1Points.numCols(), seed = 1234567;
  zVec.sizeUninitialized(m);  r1Vec.sizeUninitialized(n1);
  sortIndex1.resize(n1);      uniqueIndex1.resize(n1);
  uniqueSet1.resize(n1); // numUnique1 if count_inc1 used
  bool* is_unique1 = new bool[n1];

  webbur::point_radial_tol_unique_index_inc1(m, n1, a1Points.values(),
    duplicateTol, &seed, zVec.values(), r1Vec.values(), &sortIndex1[0],
    is_unique1, &numUnique1, &uniqueSet1[0], &uniqueIndex1[0]);

  copy_data(is_unique1, n1, isUnique1);
  delete [] is_unique1;

#ifdef DEBUG
  PCout << "Reference unique: numUnique1 = " << numUnique1 << "\na1 =\n"
	<< a1Points << "\n               r1   indx1 unique1   undx1   xdnu1:\n";
  for (size_t i=0; i<n1; ++i)
    std::cout << std::setw(17) << r1Vec[i]     << std::setw(8) << sortIndex1[i]
	      << std::setw(8)  << isUnique1[i] << std::setw(8) << uniqueSet1[i]
	      << std::setw(8)  << uniqueIndex1[i] << '\n';
  PCout << std::endl;
#endif // DEBUG

  uniqueIndexMapping = uniqueIndex1; // copy
  assign_tensor_collocation_indices(0, uniqueIndex1);
  numCollocPts = numUnique1;
  update_sparse_points(0, 0, a1Points, isUnique1, uniqueIndex1, var_sets);
  if (trackUniqueProdWeights) {
    type1WeightSets = 0.; if (computeType2Weights) type2WeightSets = 0.;
    update_sparse_weights(0, a1Type1Weights, a1Type2Weights, uniqueIndex1,
			  type1WeightSets, type2WeightSets);
#ifdef DEBUG
    PCout << "\nreference_unique() reference type1WeightSets:\n";
    write_data(PCout, type1WeightSets);
    if (computeType2Weights) {
      PCout << "\nreference_unique() reference type2WeightSets:\n";
      write_data(PCout, type2WeightSets, false, true, true);
    }
#endif // DEBUG
  }
}


void CombinedSparseGridDriver::
increment_unique(bool compute_a2, bool update_sets, RealMatrix& var_sets)
{
  // increment_unique processes the trailing Smolyak index set
  size_t last_index = smolyakMultiIndex.size() - 1;

  // define a1 pts/wts
  if (compute_a2) // else already computed (e.g., within compute_trial_grid())
    compute_tensor_points_weights(last_index, 1, a2Points, a2Type1Weights,
				  a2Type2Weights);

  // ----
  // INC2
  // ----
  int m = numVars, n1 = a1Points.numCols(), n2 = a2Points.numCols();
  r2Vec.sizeUninitialized(n2); sortIndex2.resize(n2);
  uniqueSet2.resize(n2); // numUnique2 if count_inc2 used
  uniqueIndex2.resize(n2);
  bool *is_unique1 = new bool[n1], *is_unique2 = new bool[n2];
  copy_data(isUnique1, is_unique1, n1);

  webbur::point_radial_tol_unique_index_inc2(m, n1, a1Points.values(),
    n2, a2Points.values(), duplicateTol, zVec.values(), r1Vec.values(),
    &sortIndex1[0], is_unique1, numUnique1, &uniqueSet1[0],
    &uniqueIndex1[0], r2Vec.values(), &sortIndex2[0], is_unique2,
    &numUnique2, &uniqueSet2[0], &uniqueIndex2[0]);

  copy_data(is_unique2, n2, isUnique2);
  delete [] is_unique1;
  delete [] is_unique2;

#ifdef DEBUG
  PCout << "Increment unique: numUnique2 = " << numUnique2 << "\na2 =\n"
	<< a2Points << "\n               r2   indx2 unique2   undx2   xdnu2:\n";
  for (size_t i=0; i<n2; ++i)
    std::cout << std::setw(17) << r2Vec[i]     << std::setw(8) << sortIndex2[i]
	      << std::setw(8)  << isUnique2[i] << std::setw(8) << uniqueSet2[i]
	      << std::setw(8)  << uniqueIndex2[i] << '\n';
  PCout << std::endl;
#endif // DEBUG

  uniqueIndexMapping.insert(uniqueIndexMapping.end(), uniqueIndex2.begin(),
			    uniqueIndex2.end());
  assign_tensor_collocation_indices(last_index, uniqueIndex2);
  numCollocPts = numUnique1 + numUnique2;
  if (update_sets) // update unique var_sets
    update_sparse_points(last_index, numUnique1, a2Points, isUnique2,
			 uniqueIndex2, var_sets);
  // update type{1,2}WeightSets
  if (trackUniqueProdWeights) {
    type1WeightSets = type1WeightSetsRef;// to be augmented by last_index data
    if (computeType2Weights)
      type2WeightSets = type2WeightSetsRef; // to be augmented
    update_sparse_weights(last_index, a2Type1Weights, a2Type2Weights,
			  uniqueIndex2, type1WeightSets, type2WeightSets);
#ifdef DEBUG
    PCout << "\nupdated type1 weight sets:\n" << type1WeightSets
	  << "\nupdated type2 weight sets:\n" << type2WeightSets;
#endif // DEBUG
  }
}


void CombinedSparseGridDriver::merge_unique()
{
  int m = numVars, n1 = a1Points.numCols(), n2 = a2Points.numCols(),
    n1n2 = n1+n2, n3, num_unique3;
  RealVector r3_vec(n1n2, false);
  RealMatrix a3_pts(m, n1n2, false);
  IntArray sort_index3(n1n2), unique_set3(n1n2), unique_index3(n1n2);
  bool *is_unique1 = new bool[n1], *is_unique2 = new bool[n2],
       *is_unique3 = new bool[n1n2];
  copy_data(isUnique1, is_unique1, n1);
  copy_data(isUnique2, is_unique2, n2);

  // ----
  // INC3
  // ----
  webbur::point_radial_tol_unique_index_inc3(m, n1, a1Points.values(),
    r1Vec.values(), &sortIndex1[0], is_unique1, numUnique1, &uniqueSet1[0],
    &uniqueIndex1[0], n2, a2Points.values(), r2Vec.values(), &sortIndex2[0],
    is_unique2, numUnique2, &uniqueSet2[0], &uniqueIndex2[0], &n3,
    a3_pts.values(), r3_vec.values(), &sort_index3[0], is_unique3,
    &num_unique3, &unique_set3[0], &unique_index3[0]);

#ifdef DEBUG
  PCout << "Merge unique: num_unique3 = " << num_unique3 << "\na3 =\n" << a3_pts
	<< "\n               r3   indx3 unique3   undx3   xdnu3:\n";
  for (size_t i=0; i<n1n2; ++i)
    std::cout << std::setw(17)    << r3_vec[i] << std::setw(8) << sort_index3[i]
	      << std::setw(8) << is_unique3[i] << std::setw(8) << unique_set3[i]
	      << std::setw(8) << unique_index3[i] << '\n';
  PCout << std::endl;
#endif // DEBUG

  // update reference points/weights (originally defined by _inc1)
  a1Points = a3_pts;
  if (trackUniqueProdWeights) {
    a1Type1Weights.resize(n1n2);
    if (computeType2Weights) a1Type2Weights.reshape(numVars, n1n2);
    size_t i, j;
    for (i=0; i<n2; ++i) {
      a1Type1Weights[n1+i] = a2Type1Weights[i];
      if (computeType2Weights)
	copy_data(a2Type2Weights[i], numVars, a1Type2Weights[n1+i]);
    }
  }
  // update reference indices, counts, radii
  r1Vec        = r3_vec;
  sortIndex1   = sort_index3;
  numUnique1   = num_unique3;
  uniqueSet1   = unique_set3;
  uniqueIndex1 = unique_index3;
  copy_data(is_unique3, n1n2, isUnique1);
  delete [] is_unique1; delete [] is_unique2; delete [] is_unique3;
  // update uniqueIndexMapping, collocIndices, numCollocPts
  uniqueIndexMapping = unique_index3;
  //assign_tensor_collocation_indices(0, unique_index3);
  numCollocPts = num_unique3;
}


void CombinedSparseGridDriver::finalize_unique(size_t start_index)
{
  // This fn supports multiple indices and ensures no order mixing among sets
  // by using inc2/inc3 in careful succession.

  // *** TO DO ***: This doesn't address issue of potential point replication
  // changes between initial trial set status and finalization.  Need an
  // improved mechanism for point restore/finalize in Dakota::Approximation.
  // Could add a virtual fn to interrogate collocation_indices() from the 
  // Approximation level.  Perhaps run some performance tests first to verify
  // that this condition is possible (or does structure of admissible indices
  // prevent replication in trial sets that is not first detected in old sets).

  size_t i, j, num_sm_mi = smolyakMultiIndex.size();
  int m = numVars, n1, n2, n1n2, n3, num_unique3, all_n2 = 0;
  RealVector all_a2t1_wts, r3_vec; RealMatrix a3_pts, all_a2t2_wts;
  IntArray all_unique_index2, sort_index3, unique_set3, unique_index3;
  bool *is_unique1, *is_unique2, *is_unique3;

  for (i=start_index; i<num_sm_mi; ++i) {

    compute_tensor_points_weights(i, 1, a2Points, a2Type1Weights,
				  a2Type2Weights);
    n1 = a1Points.numCols(); n2 = a2Points.numCols();
    all_a2t1_wts.resize(all_n2+n2);
    if (computeType2Weights) all_a2t2_wts.reshape(numVars, all_n2+n2);
    for (j=0; j<n2; ++j) {
      all_a2t1_wts[all_n2+j] = a2Type1Weights[j];
      if (computeType2Weights)
	copy_data(a2Type2Weights[j], numVars, all_a2t2_wts[all_n2+j]);
    }
    all_n2 += n2;

    // INC2
    r2Vec.sizeUninitialized(n2); sortIndex2.resize(n2);
    uniqueSet2.resize(n2);       uniqueIndex2.resize(n2);
    is_unique1 = new bool[n1];   is_unique2 = new bool[n2];
    copy_data(isUnique1, is_unique1, n1);
    webbur::point_radial_tol_unique_index_inc2(m, n1, a1Points.values(),
      n2, a2Points.values(), duplicateTol, zVec.values(), r1Vec.values(),
      &sortIndex1[0], is_unique1, numUnique1, &uniqueSet1[0],
      &uniqueIndex1[0], r2Vec.values(), &sortIndex2[0], is_unique2,
      &numUnique2, &uniqueSet2[0], &uniqueIndex2[0]);
#ifdef DEBUG
    PCout << "Finalize unique: numUnique2 = " << numUnique2 << "\na2 =\n"
	  << a2Points<<"\n               r2   indx2 unique2   undx2   xdnu2:\n";
    for (j=0; j<n2; ++j)
      std::cout << std::setw(17)    << r2Vec[j] << std::setw(8) << sortIndex2[j]
		<< std::setw(8) << isUnique2[j] << std::setw(8) << uniqueSet2[j]
		<< std::setw(8) << uniqueIndex2[j] << '\n';
    PCout << std::endl;
#endif // DEBUG

    all_unique_index2.insert(all_unique_index2.end(), uniqueIndex2.begin(),
			     uniqueIndex2.end());
    numCollocPts += numUnique2;

    if (i < num_sm_mi - 1) {
      // INC3
      n1n2 = n1+n2;                       r3_vec.sizeUninitialized(n1n2);
      a3_pts.shapeUninitialized(m, n1n2); sort_index3.resize(n1n2);
      unique_set3.resize(n1n2);           unique_index3.resize(n1n2);
      is_unique3 = new bool[n1n2];
      webbur::point_radial_tol_unique_index_inc3(m, n1, a1Points.values(),
        r1Vec.values(), &sortIndex1[0], is_unique1, numUnique1, &uniqueSet1[0],
        &uniqueIndex1[0], n2, a2Points.values(), r2Vec.values(), &sortIndex2[0],
        is_unique2, numUnique2, &uniqueSet2[0], &uniqueIndex2[0], &n3,
        a3_pts.values(), r3_vec.values(), &sort_index3[0], is_unique3,
        &num_unique3, &unique_set3[0], &unique_index3[0]);
#ifdef DEBUG
      PCout << "Finalize unique: num_unique3 = " << num_unique3 << "\na3 =\n"
	    << a3_pts<<"\n               r3   indx3 unique3   undx3   xdnu3:\n";
      for (j=0; j<n1n2; ++j)
	std::cout << std::setw(17) << r3_vec[j] << std::setw(8)
		  << sort_index3[j] << std::setw(8) << is_unique3[j]
		  << std::setw(8) << unique_set3[j] << std::setw(8)
		  << unique_index3[j] << '\n';
      PCout << std::endl;
#endif // DEBUG

      // update reference points, indices, counts, radii
      a1Points     = a3_pts;
      r1Vec        = r3_vec;
      sortIndex1   = sort_index3;
      numUnique1   = num_unique3;
      uniqueSet1   = unique_set3;
      uniqueIndex1 = unique_index3;
      copy_data(is_unique3, n1n2, isUnique1);
      delete [] is_unique3;
    }

    delete [] is_unique1; delete [] is_unique2;
  }

  uniqueIndexMapping.insert(uniqueIndexMapping.end(), all_unique_index2.begin(),
			    all_unique_index2.end());
  assign_tensor_collocation_indices(start_index, all_unique_index2);
  if (trackUniqueProdWeights) {
    type1WeightSets = type1WeightSetsRef; // to be augmented
    if (computeType2Weights)
      type2WeightSets = type2WeightSetsRef; // to be augmented
    update_sparse_weights(start_index, all_a2t1_wts, all_a2t2_wts,
			  all_unique_index2, type1WeightSets, type2WeightSets);
#ifdef DEBUG
    PCout << "type1WeightSets =\n"; write_data(PCout, type1WeightSets);
#endif // DEBUG
  }
}


void CombinedSparseGridDriver::
update_sparse_points(size_t start_index, int new_index_offset,
		     const RealMatrix& tensor_pts, const BitArray& is_unique,
		     const IntArray& unique_index, RealMatrix& new_sparse_pts)
{
  size_t i, j, cntr, num_sm_mi = smolyakMultiIndex.size(), num_tp_pts,
    num_pts = is_unique.size(), num_unique_pts = 0;
  for (i=0; i<num_pts; ++i)
    if (is_unique[i])
      ++num_unique_pts;

  // update sizes
  new_sparse_pts.shapeUninitialized(numVars, num_unique_pts);

  int index;
  // add contributions for new index sets
  for (i=start_index, cntr=0; i<num_sm_mi; ++i) {
    num_tp_pts = collocKey[i].size();
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      if (is_unique[cntr]) {
	index = unique_index[cntr] - new_index_offset;
	copy_data(tensor_pts[cntr], numVars, new_sparse_pts[index]);
      }
    }
  }
}


void CombinedSparseGridDriver::
update_sparse_weights(size_t start_index, const RealVector& tensor_t1_wts,
		      const RealMatrix& tensor_t2_wts,
		      const IntArray& unique_index, RealVector& updated_t1_wts,
		      RealMatrix& updated_t2_wts)
{
  size_t i, j, k, cntr, num_sm_mi = smolyakMultiIndex.size(), num_tp_pts;

  // update sizes
  updated_t1_wts.resize(numCollocPts); // new entries initialized to 0
  if (computeType2Weights)
    updated_t2_wts.reshape(numVars, numCollocPts); // new entries init to 0

  int index, delta_coeff, sm_coeff;
  // back out changes in Smolyak coeff for existing index sets
  for (i=0, cntr=0; i<start_index; ++i) {
    delta_coeff = smolyakCoeffs[i] - smolyakCoeffsRef[i];
    if (delta_coeff) {
      num_tp_pts = collocKey[i].size();
      for (j=0; j<num_tp_pts; ++j, ++cntr) {
	index = uniqueIndex1[cntr];
	updated_t1_wts[index] += delta_coeff * a1Type1Weights[cntr];
	if (computeType2Weights) {
	  Real*       up_t2_wts_j = updated_t2_wts[index];
	  const Real* a1_t2_wts_j = a1Type2Weights[cntr];
	  for (k=0; k<numVars; ++k)
	    up_t2_wts_j[k] += delta_coeff * a1_t2_wts_j[k];
	}
      }
    }
    else
      cntr += collocKey[i].size();
  }
  // add contributions for new index sets
  for (i=start_index, cntr=0; i<num_sm_mi; ++i) {
    sm_coeff = smolyakCoeffs[i];
    if (sm_coeff) {
      num_tp_pts = collocKey[i].size();
      for (j=0; j<num_tp_pts; ++j, ++cntr) {
	index = unique_index[cntr];
	updated_t1_wts[index] += sm_coeff * tensor_t1_wts[cntr];
	if (computeType2Weights) {
	  Real*       up_t2_wts_j = updated_t2_wts[index];
	  const Real* te_t2_wts_j = tensor_t2_wts[cntr];
	  for (k=0; k<numVars; ++k)
	    up_t2_wts_j[k] += sm_coeff * te_t2_wts_j[k];
	}
      }
    }
    else
      cntr += collocKey[i].size();
  }
}


void CombinedSparseGridDriver::
compute_tensor_points_weights(size_t start_index, size_t num_indices,
			      RealMatrix& pts, RealVector& t1_wts,
			      RealMatrix& t2_wts)
{
  size_t i, j, k, l, cntr, num_tp_pts, num_colloc_pts = 0,
    end = start_index + num_indices;
  // define num_colloc_pts
  for (i=start_index; i<end; ++i)
    num_colloc_pts += collocKey[i].size();
  // define pts/wts: wts are raw product weights; Smolyak combinatorial
  // coefficient applied in compute_grid()/compute_trial_grid()
  pts.shapeUninitialized(numVars, num_colloc_pts);
  t1_wts.sizeUninitialized(num_colloc_pts);
  if (computeType2Weights)
    t2_wts.shapeUninitialized(numVars, num_colloc_pts);
  for (i=start_index, cntr=0; i<end; ++i) {
    const UShortArray& sm_index = smolyakMultiIndex[i];
    num_tp_pts = collocKey[i].size();
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      const UShortArray& key_ij = collocKey[i][j];
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
    PCout << "Tensor product weights =\ntype1:\n"; write_data(PCout, t1_wts);
    PCout << "type2:\n"; write_data(PCout, t2_wts, false, true, true);
#endif // DEBUG
}


void CombinedSparseGridDriver::
assign_tensor_collocation_indices(size_t start_index, 
				  const IntArray& unique_index)
{
  size_t i, j, cntr, num_tp_pts, num_sm_mi = smolyakMultiIndex.size();
  if (collocIndices.size() < num_sm_mi)
    collocIndices.resize(num_sm_mi);
  for (i=start_index, cntr=0; i<num_sm_mi; ++i) {
    num_tp_pts = collocKey[i].size();
    SizetArray& indices_i = collocIndices[i];
    indices_i.resize(num_tp_pts);
    for (j=0; j<num_tp_pts; ++j, ++cntr)
      indices_i[j] = unique_index[cntr];
  }
}

} // namespace Pecos
