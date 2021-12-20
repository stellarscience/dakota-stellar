/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 TensorProductDriver
//- Description: Implementation code for TensorProductDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "TensorProductDriver.hpp"
#include "PolynomialApproximation.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: TensorProductDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void TensorProductDriver::
initialize_grid(const MultivariateDistribution& u_dist,
		const ExpansionConfigOptions& ec_options,
		const BasisConfigOptions& bc_options)
{
  IntegrationDriver::initialize_grid(u_dist, ec_options, bc_options);
  quadOrder.resize(numVars); levelIndIter->second.resize(numVars);
}


void TensorProductDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis)
{
  IntegrationDriver::initialize_grid(poly_basis);
  quadOrder.resize(numVars); levelIndIter->second.resize(numVars);
}


void TensorProductDriver::precompute_rules()
{
  for (size_t i=0; i<numVars; ++i)
    polynomialBasis[i].precompute_rules(quadOrder[i]);
}


void TensorProductDriver::clear_inactive()
{
  std::map<ActiveKey, UShortArray>::iterator   li_it = levelIndex.begin();
  std::map<ActiveKey, UShort2DArray>::iterator ck_it = collocKey.begin();
  std::map<ActiveKey, RealVector>::iterator    t1_it = type1WeightSets.begin();
  std::map<ActiveKey, RealMatrix>::iterator    t2_it = type2WeightSets.begin();
  while (li_it != levelIndex.end())
    if (li_it == levelIndIter) // preserve active
      { ++li_it, ++ck_it, ++t1_it, ++t2_it; }
    else { // clear inactive: postfix increments manage iterator invalidations
      levelIndex.erase(li_it++);      collocKey.erase(ck_it++);
      type1WeightSets.erase(t1_it++); type2WeightSets.erase(t2_it++);
    }
}


const ActiveKey& TensorProductDriver::maximal_grid()
{
  std::map<ActiveKey, RealVector>::const_iterator
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


void TensorProductDriver::combine_grid()
{
  std::map<ActiveKey, UShortArray>::const_iterator
    li_cit = levelIndex.begin();
  combinedLevelIndex = li_cit->second; ++li_cit;
  for (; li_cit!=levelIndex.end(); ++li_cit) {
    const UShortArray& li = li_cit->second;
    for (size_t v=0; v<numVars; ++v)
      if (li[v] > combinedLevelIndex[v])
	combinedLevelIndex[v] = li[v];
  }

  UShortArray comb_order;
  level_to_order(combinedLevelIndex, comb_order);
  compute_tensor_grid(comb_order, combinedLevelIndex, combinedVarSets,
		      combinedT1WeightSets, combinedT2WeightSets,
		      combinedCollocKey);
}


void TensorProductDriver::combined_to_active(bool clear_combined)
{
  // Replace active arrays with combined arrays

  // Note: inactive weight sets to be removed by clear_inactive()

  if (clear_combined) {
    std::swap(levelIndIter->second,   combinedLevelIndex);
    std::swap(collocKeyIter->second,  combinedCollocKey);
    std::swap(varSetsIter->second,    combinedVarSets);
    std::swap(t1WtIter->second,       combinedT1WeightSets);
    std::swap(t2WtIter->second,       combinedT2WeightSets);

    combinedLevelIndex.clear();
    combinedCollocKey.clear();
    combinedVarSets.shapeUninitialized(0,0);
    combinedT1WeightSets.sizeUninitialized(0);
    combinedT2WeightSets.shapeUninitialized(0,0);
  }
  else {
    levelIndIter->second   = combinedLevelIndex;
    collocKeyIter->second  = combinedCollocKey;
    varSetsIter->second    = combinedVarSets;
    t1WtIter->second       = combinedT1WeightSets;
    t2WtIter->second       = combinedT2WeightSets;
  }

  level_to_order();
}


void TensorProductDriver::
enforce_constraints(const UShortArray& ref_quad_order)
{
  // enforce constraints: ref_quad_order -> quadOrder
  size_t i, len = ref_quad_order.size();
  if (quadOrder.size()            != len)            quadOrder.resize(len);
  if (levelIndIter->second.size() != len) levelIndIter->second.resize(len);
  unsigned short nested_order;
  for (i=0; i<len; ++i) {
    // synchronize on number of points: Lagrange poly order = #pts - 1
    if (driverMode == INTERPOLATION_MODE)
      quadrature_goal_to_nested_quadrature_order(i, ref_quad_order[i],
						 nested_order);
    else // {INTEGRATION,DEFAULT}_MODE: ref_quad_order is non-nested so use
         // non-nested Gauss integrand goal = 2m-1
      integrand_goal_to_nested_quadrature_order(i, 2 * ref_quad_order[i] - 1,
						nested_order);

    // update quadOrder / levelIndex
    if (nested_order == USHRT_MAX) { // required order not available
      PCerr << "Error: order goal could not be attained in TensorProductDriver"
	    << "::enforce_constraints()" << std::endl;
      abort_handler(-1);
    }
    else
      quadrature_order(nested_order, i); // sets quadOrder and levelIndex
  }
}


/** This function selects the smallest nested rule order that meets the
    integrand precision of a corresponding Gauss rule.  It is similar to
    the moderate exponential growth option in sparse grids. */
void TensorProductDriver::
integrand_goal_to_nested_quadrature_order(size_t i,
					  unsigned short integrand_goal,
					  unsigned short& nested_quad_order)
{
  switch (collocRules[i]) {
  case CLENSHAW_CURTIS: case NEWTON_COTES: { // closed rules
    nested_quad_order = 1; // in case while loop falls through
    unsigned short /*level = 0,*/ lev_pow = 1, integrand_actual = 1;
    while (integrand_actual < integrand_goal) {
      lev_pow *= 2; //++level;
      nested_quad_order = lev_pow + 1; //std::pow(2.,level) + 1;
      // integrand exactness for CC; also used for NC
      integrand_actual = (nested_quad_order % 2) ?
	nested_quad_order : nested_quad_order - 1;
    }
    break;
  }
  case FEJER2: { // open rule (integrand exactness for CC also used for F2)
    nested_quad_order = 1; // in case while loop falls through
    unsigned short /*level = 0,*/ lev_pow = 2, integrand_actual = 1;
    while (integrand_actual < integrand_goal) {
      lev_pow *= 2; //++level;
      nested_quad_order = lev_pow - 1; //std::pow(2.,level+1) - 1;
      integrand_actual = (nested_quad_order % 2) ?
	nested_quad_order : nested_quad_order - 1;
    }
    break;
  }
  case GAUSS_PATTERSON: { // open rule
    nested_quad_order = 1; // in case while loop falls through
    unsigned short /*level = 0,*/ lev_pow = 2, integrand_actual = 1,
      previous = nested_quad_order;
    while (integrand_actual < integrand_goal) {
      lev_pow *= 2; //++level;
      // exponential growth
      nested_quad_order = lev_pow - 1; //std::pow(2.,level+1) - 1;
      integrand_actual = 2*nested_quad_order - previous;//2m-1 - constraints + 1
      previous = nested_quad_order;
    }
    break;
  }
  case GENZ_KEISTER: { // open rule with lookup
    unsigned short level = 0, max_level = 5;
    while (level <= max_level && precGenzKeister[level] < integrand_goal)
      ++level;
    nested_quad_order = (level > max_level) ?
      USHRT_MAX : orderGenzKeister[level]; // pass error state up a level
    /*
    nested_quad_order = 1;
    unsigned short integrand_goal = 2*ref_quad_order - 1, level = 0,
      integrand_actual = 1, previous = nested_quad_order, i_rule = GENZ_KEISTER,
      g_rule = FULL_EXPONENTIAL; // map l->o without restriction
    while (integrand_actual < integrand_goal) {
      ++level;
      webbur::level_growth_to_order_new(1, &level, &i_rule, &g_rule,
	                                nested_quad_order);
      integrand_actual = 2*nested_quad_order - previous;//2m-1 - constraints + 1
      previous = nested_quad_order;
    }
    */
    break;
  }
  default: { // Gauss rules
    nested_quad_order = 1;
    unsigned short integrand_actual = 1;
    while (integrand_actual < integrand_goal) {
      // allow even quad order (won't happen for current goal definitions)
      ++nested_quad_order;
      //nested_quad_order += 2; // moderate linear growth: odd for weakly nested
      integrand_actual = 2*nested_quad_order - 1;
    }
    break;
  }
  }
}


/** This function selects the smallest nested rule order that meets the
    quadrature order goal. */
void TensorProductDriver::
quadrature_goal_to_nested_quadrature_order(size_t i, unsigned short quad_goal,
					   unsigned short& nested_quad_order)
{
  switch (collocRules[i]) {
  case CLENSHAW_CURTIS: case NEWTON_COTES: { // closed nested rules
    nested_quad_order = 1; unsigned short /*level = 0,*/ lev_pow = 1;
    while (nested_quad_order < quad_goal) {
      lev_pow *= 2; //++level;
      nested_quad_order = lev_pow + 1; //std::pow(2.,level) + 1;
    }
    break;
  }
  case FEJER2: case GAUSS_PATTERSON:{ // open nested rules
    nested_quad_order = 1; unsigned short /*level = 0,*/ lev_pow = 2;
    while (nested_quad_order < quad_goal) {
      lev_pow *= 2; //++level;
      nested_quad_order = lev_pow - 1; //std::pow(2.,level+1) - 1;
    }
    break;
  }
  case GENZ_KEISTER: { // open nested rule with lookup
    nested_quad_order = 1; unsigned short level = 0, max_level = 5;
    while (level <= max_level && orderGenzKeister[level] < quad_goal)
      ++level;
    nested_quad_order = (level > max_level) ?
      USHRT_MAX : orderGenzKeister[level]; // pass error state up a level
    break;
  }
  default: // open weakly/non-nested Gauss rules
    nested_quad_order = quad_goal; break;
  }
}


void TensorProductDriver::
reinterpolated_tensor_grid(const UShortArray& lev_index,
			   const SizetList& reinterp_indices)
{
  if (lev_index != levelIndIter->second) {
    PCerr << "Error: inconsistent level index in TensorProductDriver::"
	  << "reinterpolated_tensor_grid()." << std::endl;
    abort_handler(-1);
  }

  std::map<UShortArray, size_t>::iterator map_it = reinterpMap.find(lev_index);
  if (map_it == reinterpMap.end()) {

    if (reinterpLevelIndices.empty()) {
      reinterpLevelIndices.resize(1); reinterpQuadOrders.resize(1);
      reinterpVarSets.resize(1);      reinterpCollocKeys.resize(1);
    }
    UShortArray& reinterp_lev_index  = reinterpLevelIndices.back();
    UShortArray& reinterp_quad_order = reinterpQuadOrders.back();
    if (reinterp_lev_index.size() != numVars) {
      reinterp_lev_index.resize(numVars);
      reinterp_quad_order.resize(numVars);
    }

    SizetList::const_iterator cit = reinterp_indices.begin();
    for (size_t i=0; i<numVars; ++i) {
      if (cit != reinterp_indices.end() && i == *cit) { // reinterpolated index
	switch (collocRules[i]) {
	// Note: stick with standard nonlinear progression, even though
	//       CC/NC/F2 could admit other options.
	case CLENSHAW_CURTIS: case NEWTON_COTES: // 1, 3, 5, 9 = 2^i+1
	  reinterp_quad_order[i] = (quadOrder[i] == 1) ? 3 : 2*quadOrder[i] - 1;
	  break;
	case GAUSS_PATTERSON: case FEJER2: // 1, 3, 7, 15 = 2^{i+1}-1
	  reinterp_quad_order[i] = 2*quadOrder[i] + 1;
	  break;
	case GENZ_KEISTER: {
	  size_t index = find_index(orderGenzKeister, quadOrder[i]);
	  if (index < orderGenzKeister.size() - 1) // also manages _NPOS
	    reinterp_quad_order[i] = orderGenzKeister[++index];
	  else {
	    PCerr << "Error: Genz-Keister lookup failure in TensorProductDriver"
		  << "::reinterpolated_tensor_grid()." << std::endl;
	    abort_handler(-1);
	  }
	  break;
	}
	default: // Gauss rules
	  // interp order = m-1; doubled interp order = 2m-2 = (2m-1)-1;
	  // so a reinterp quad order of 2m-1 doubles the interpolant order:
	  reinterp_quad_order[i] = 2*quadOrder[i] - 1;
	}
	reinterp_lev_index[i] = reinterp_quad_order[i] - 1;
	// advance to the next reinterp index
	++cit;
      }
      else { // not a reinterpolated index --> no change from reference
	reinterp_quad_order[i] = quadOrder[i];
	reinterp_lev_index[i]  = lev_index[i];
      }
    }

    // compute the reinterpolation grid
    compute_tensor_grid(reinterp_quad_order, reinterp_lev_index,
			reinterp_indices, reinterpVarSets.back(),
			reinterpCollocKeys.back());

    // update reiterpMap bookkeeping: only 1 index needs to be tracked for TPQ
    reinterpMap.clear();
    reinterpMap[lev_index] = activeReinterpIndex = 0;
  }
  else
    activeReinterpIndex = map_it->second;
}


void TensorProductDriver::compute_grid()
{
#ifdef DEBUG
  // -----------------------------------
  // Output number of collocation points
  // -----------------------------------
  PCout << "Total number of tensor-product quadrature points: "
	<< grid_size() << '\n';
#endif // DEBUG

  // -------------------------------------------------------------------
  // Get collocation points and integration weights and update 1D arrays
  // -------------------------------------------------------------------
  compute_tensor_grid(quadOrder, levelIndIter->second, varSetsIter->second,
		      t1WtIter->second, t2WtIter->second,
		      collocKeyIter->second);
}

} // namespace Pecos
