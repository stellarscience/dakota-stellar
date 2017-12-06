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
initialize_grid(const ShortArray& u_types,
		const ExpansionConfigOptions& ec_options,
		const BasisConfigOptions& bc_options)
{
  IntegrationDriver::initialize_grid(u_types, ec_options, bc_options);
  quadOrder.resize(numVars); levelIndex.resize(numVars);
}


void TensorProductDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis)
{
  IntegrationDriver::initialize_grid(poly_basis);
  quadOrder.resize(numVars); levelIndex.resize(numVars);
}


void TensorProductDriver::precompute_rules()
{
  for (size_t i=0; i<numVars; ++i)
    polynomialBasis[i].precompute_rules(quadOrder[i]);
}


void TensorProductDriver::store_grid(size_t index)
{
  size_t stored_len = storedType1WeightSets.size();
  if (index == _NPOS || index == stored_len) { // append
    storedCollocKey.push_back(collocKey);
    storedLevelIndex.push_back(levelIndex);
    storedType1WeightSets.push_back(type1WeightSets);
    storedType2WeightSets.push_back(type2WeightSets);
  }
  else if (index < stored_len) { // replace
    storedCollocKey[index] = collocKey; storedLevelIndex[index] = levelIndex;
    storedType1WeightSets[index] = type1WeightSets;
    storedType2WeightSets[index] = type2WeightSets;
  }
  else {
    PCerr << "Error: bad index (" << index << ") passed in TensorProductDriver"
	  << "::store_grid()" << std::endl;
    abort_handler(-1);
  }
}


void TensorProductDriver::restore_grid(size_t index)
{
  size_t stored_len = storedType1WeightSets.size();
  if (index == _NPOS) {
    collocKey = storedCollocKey.back(); levelIndex = storedLevelIndex.back();
    type1WeightSets = storedType1WeightSets.back();
    type2WeightSets = storedType2WeightSets.back();
  }
  else if (index < stored_len) {
    collocKey = storedCollocKey[index]; levelIndex = storedLevelIndex[index];
    type1WeightSets = storedType1WeightSets[index];
    type2WeightSets = storedType2WeightSets[index];
  }
  else {
    PCerr << "Error: bad index (" << index << ") passed in TensorProductDriver"
	  << "::restore_grid()" << std::endl;
    abort_handler(-1);
  }
}


void TensorProductDriver::remove_stored_grid(size_t index)
{
  size_t stored_len = storedType1WeightSets.size();
  if (index == _NPOS || index == stored_len) {
    storedCollocKey.pop_back(); storedLevelIndex.pop_back();
    storedType1WeightSets.pop_back();
    storedType2WeightSets.pop_back();
  }
  else if (index < stored_len) {
    UShort3DArray::iterator u3it = storedCollocKey.begin();
    std::advance(u3it, index); storedCollocKey.erase(u3it);
    UShort2DArray::iterator u2it = storedLevelIndex.begin();
    std::advance(u2it, index); storedLevelIndex.erase(u2it);
    RealVectorArray::iterator vit = storedType1WeightSets.begin();
    std::advance(vit, index); storedType1WeightSets.erase(vit);
    RealMatrixArray::iterator mit = storedType2WeightSets.begin();
    std::advance(mit, index); storedType2WeightSets.erase(mit);
  }
}


void TensorProductDriver::clear_stored()
{
  storedLevelIndex.clear();      storedCollocKey.clear();
  storedType1WeightSets.clear(); storedType2WeightSets.clear();
}


size_t TensorProductDriver::maximal_grid() const
{
  size_t i, num_stored = storedType1WeightSets.size(),
    max_index = _NPOS, max_wts = type1WeightSets.length();
  for (i=0; i<num_stored; ++i)
    if (storedType1WeightSets[i].length() > max_wts)
      { max_index = i; max_wts = storedType1WeightSets[i].length(); }
  return max_index;
}


void TensorProductDriver::swap_grid(size_t index)
{
  std::swap(storedCollocKey[index],  collocKey);
  std::swap(storedLevelIndex[index], levelIndex);
  update_quadrature_order_from_level_index();

  RealVector tmp_vec(type1WeightSets);
  type1WeightSets = storedType1WeightSets[index];
  storedType1WeightSets[index] = tmp_vec;

  RealMatrix tmp_mat(type2WeightSets);
  type2WeightSets = storedType2WeightSets[index];
  storedType2WeightSets[index] = tmp_mat;
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
    nested_quad_order = orderGenzKeister[level];
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
    while (level <= max_level && nested_quad_order < quad_goal) {
      ++level;
      nested_quad_order = orderGenzKeister[level];
    }
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
  if (lev_index != levelIndex) {
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
	reinterp_lev_index[i]  = levelIndex[i];
      }
    }

    // compute the reinterpolation grid
    compute_tensor_grid(reinterp_quad_order, reinterp_lev_index,
			reinterp_indices, reinterpVarSets.back(),
			reinterpCollocKeys.back());

    // update reiterpMap bookkeeping: only 1 index needs to be tracked for TPQ
    reinterpMap.clear();
    reinterpMap[levelIndex] = activeReinterpIndex = 0;
  }
  else
    activeReinterpIndex = map_it->second;
}


void TensorProductDriver::compute_grid(RealMatrix& variable_sets)
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
  compute_tensor_grid(quadOrder, levelIndex, variable_sets, type1WeightSets,
		      type2WeightSets, collocKey);
}

} // namespace Pecos
