/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:       SharedProjectOrthogPolyApproxData
//- Description: Implementation code for SharedProjectOrthogPolyApproxData class
//-               
//- Owner:       Mike Eldred

#include "SharedProjectOrthogPolyApproxData.hpp"
#include "TensorProductDriver.hpp"
#include "CombinedSparseGridDriver.hpp"
#include "CubatureDriver.hpp"
#include "pecos_global_defs.hpp"

//#define DEBUG

namespace Pecos {


void SharedProjectOrthogPolyApproxData::allocate_data()
{
  // no combination by default, even if stored{MultiIndex,ExpCoeffs,
  // ExpCoeffGrads} are defined.  Redefined by attribute passed in
  // combine_coefficients(short).
  storedExpCombineType = NO_COMBINE; // reset to initial state (if needed)

  // update_exp_form controls when to update (refinement) and when not to
  // update (subIterator execution) an expansion's multiIndex definition.
  // Simple logic of updating if previous number of points != current number
  // is not robust enough for anisotropic updates --> track using Prev arrays.
  switch (expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: {
    TensorProductDriver* tpq_driver = (TensorProductDriver*)driverRep;
    const UShortArray& quad_order = tpq_driver->quadrature_order();
    bool update_exp_form = (quad_order != quadOrderPrev);
    // *** TO DO: capture updates to parameterized/numerical polynomials?

    if (update_exp_form) {
      UShortArray int_order(numVars);
      quadrature_order_to_integrand_order(driverRep, quad_order, int_order);
      integrand_order_to_expansion_order(int_order, approxOrder);
      tensor_product_multi_index(approxOrder, multiIndex); // include upper bnd

      // precomputation performed by tpqDriver prior to allocate_data()
      //precompute_maximal_rules(approxOrder);

      allocate_component_sobol(multiIndex);
      quadOrderPrev = quad_order;
    }

#ifdef DEBUG
    for (i=0; i<numVars; ++i) {
      OrthogonalPolynomial* poly_rep
	= (OrthogonalPolynomial*)polynomialBasis[i].polynomial_rep();
      for (j=1; j<=quad_order[i]; ++j)
	poly_rep->gauss_check(j);
    }
#endif // DEBUG

    PCout << "Orthogonal polynomial approximation order = { ";
    for (size_t i=0; i<numVars; ++i) PCout << approxOrder[i] << ' ';
    PCout << "} using tensor-product expansion of " << multiIndex.size()
	  << " terms\n";
    break;
  }
  case CUBATURE: {
    CubatureDriver* cub_driver = (CubatureDriver*)driverRep;
    //unsigned short cub_int_order = cub_driver->integrand_order();
    //bool update_exp_form = (cub_int_order != cubIntOrderPrev);

    //if (update_exp_form) {
      UShortArray integrand_order(numVars, cub_driver->integrand_order());
      integrand_order_to_expansion_order(integrand_order, approxOrder);
      total_order_multi_index(approxOrder, multiIndex);

      // See special logic in CubatureDriver::compute_grid() for GOLUB_WELSCH
      //precompute_maximal_rules(approxOrder);

      allocate_component_sobol(multiIndex);
      //cubIntOrderPrev = cub_int_order; // update reference point
    //}

    PCout << "Orthogonal polynomial approximation order = { ";
    for (size_t i=0; i<numVars; ++i)
      PCout << approxOrder[i] << ' ';
    PCout << "} using total-order expansion of " << multiIndex.size()
	  << " terms\n";
    break;
  }
  case COMBINED_SPARSE_GRID: {
    CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
    unsigned short    ssg_level = csg_driver->level();
    const RealVector& aniso_wts = csg_driver->anisotropic_weights();
    bool update_exp_form
      = (ssg_level != ssgLevelPrev || aniso_wts != ssgAnisoWtsPrev ||
	 expConfigOptions.refinementControl ==
	 DIMENSION_ADAPTIVE_CONTROL_GENERALIZED);
    // *** TO DO: capture updates to parameterized/numerical polynomials?

    if (update_exp_form) {
      sparse_grid_multi_index(csg_driver, multiIndex);

      // precomputation performed by ssgDriver prior to allocate_data()
      //precompute_maximal_rules(multiIndex);

      allocate_component_sobol(multiIndex);
      ssgLevelPrev = ssg_level; ssgAnisoWtsPrev = aniso_wts;
    }
    PCout << "Orthogonal polynomial approximation level = " << ssg_level
	  << " using tensor integration and tensor sum expansion of "
	  << multiIndex.size() << " terms\n"; break;
    break;
  }
  default: // SAMPLING
    SharedOrthogPolyApproxData::allocate_data();
    break;
  }
}


void SharedProjectOrthogPolyApproxData::increment_data()
{
  if (expConfigOptions.expCoeffsSolnApproach != COMBINED_SPARSE_GRID) {
    PCerr << "Error: unsupported grid definition in SharedProjectOrthogPoly"
	  << "ApproxData::increment_data()" << std::endl;
    abort_handler(-1);
  }

  // increment tpMultiIndex{,Map,MapRef} arrays, update tpMultiIndex,
  // update multiIndex and append bookkeeping
  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  increment_trial_set(csg_driver, multiIndex);
  // update Sobol' array sizes to pick up new interaction terms
  increment_component_sobol();

  // cleanup
  //if (!reEntrantFlag) {
  //  csg_driver->clear_smolyak_arrays();
  //  csg_driver->clear_collocation_arrays();
  //  tpMultiIndex.clear(); tpMultiIndexMap.clear();
  //}
}


void SharedProjectOrthogPolyApproxData::increment_component_sobol()
{
  if (expConfigOptions.vbdFlag && expConfigOptions.vbdOrderLimit != 1) {
    reset_sobol_index_map_values();
    multi_index_to_sobol_index_map(tpMultiIndex.back());
    assign_sobol_index_map_values();
  }
}


void SharedProjectOrthogPolyApproxData::decrement_data()
{
  if (expConfigOptions.expCoeffsSolnApproach != COMBINED_SPARSE_GRID) {
    PCerr << "Error: unsupported grid definition in SharedProjectOrthogPoly"
	  << "ApproxData::decrement_data()" << std::endl;
    abort_handler(-1);
  }

  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  decrement_trial_set(csg_driver->trial_set(), multiIndex);
}


void SharedProjectOrthogPolyApproxData::pre_push_data()
{
  if (expConfigOptions.expCoeffsSolnApproach != COMBINED_SPARSE_GRID) {
    PCerr << "Error: unsupported grid definition in SharedProjectOrthogPoly"
	  << "ApproxData::pre_push_data()" << std::endl;
    abort_handler(-1);
  }

  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  pre_push_trial_set(csg_driver->trial_set(), multiIndex);
}


void SharedProjectOrthogPolyApproxData::post_push_data()
{
  CombinedSparseGridDriver* csg_driver = (CombinedSparseGridDriver*)driverRep;
  post_push_trial_set(csg_driver->trial_set(), multiIndex);
}


void SharedProjectOrthogPolyApproxData::pre_finalize_data()
{
  if (expConfigOptions.expCoeffsSolnApproach != COMBINED_SPARSE_GRID) {
    PCerr << "Error: unsupported grid definition in SharedProjectOrthogPoly"
	  << "ApproxData::finalize_data()" << std::endl;
    abort_handler(-1);
  }

  // update multiIndex
  std::deque<UShort2DArray>::iterator iit = poppedTPMultiIndex.begin();
  std::deque<SizetArray>::iterator    mit = poppedTPMultiIndexMap.begin();
  std::deque<size_t>::iterator        rit = poppedTPMultiIndexMapRef.begin();
  for (; iit!=poppedTPMultiIndex.end(); ++iit, ++mit, ++rit)
    append_multi_index(*iit, *mit, *rit, multiIndex);
  // move previous expansion data to current expansion
  tpMultiIndex.insert(tpMultiIndex.end(), poppedTPMultiIndex.begin(),
    poppedTPMultiIndex.end());
  tpMultiIndexMap.insert(tpMultiIndexMap.end(), poppedTPMultiIndexMap.begin(),
    poppedTPMultiIndexMap.end());
  tpMultiIndexMapRef.insert(tpMultiIndexMapRef.end(),
    poppedTPMultiIndexMapRef.begin(), poppedTPMultiIndexMapRef.end());
}


void SharedProjectOrthogPolyApproxData::post_finalize_data()
{
  poppedLevMultiIndex.clear();   poppedTPMultiIndex.clear();
  poppedTPMultiIndexMap.clear(); poppedTPMultiIndexMapRef.clear();
}


size_t SharedProjectOrthogPolyApproxData::pre_combine_data(short combine_type)
{
  // based on incoming combine_type, combine the data stored previously
  // by store_coefficients()

  storedExpCombineType = combine_type;

  switch (combine_type) {
  case ADD_COMBINE: {
    // Note: would like to preserve tensor indexing (at least for QUADRATURE
    // case) so that Horner's rule performance opt could be used within
    // tensor_product_value()).  However, a tensor result in the overlay
    // will not occur unless one expansion order dominates the other (partial
    // domination results in sum of tensor expansions as for sparse grids).
    // Therefore, stick with the general-purpose expansion overlay and exclude
    // tensor_product_value() usage for combined coefficient sets.

    // base class version is sufficient; no specialization based on exp form
    return SharedOrthogPolyApproxData::pre_combine_data(combine_type);
    break;
  }
  case MULT_COMBINE: {
    // compute form of product expansion
    switch (expConfigOptions.expCoeffsSolnApproach) {
    case QUADRATURE: { // product of two tensor-product expansions
      size_t max_index = driverRep->maximal_grid();
      if (max_index != _NPOS) swap_data(max_index);

      size_t i, j, num_stored = storedApproxOrder.size();
      for (i=0; i<num_stored; ++i)
	for (j=0; j<numVars; ++j)
	  approxOrder[j] += storedApproxOrder[i][j];

      UShort2DArray multi_index_prod;
      tensor_product_multi_index(approxOrder, combinedMultiIndex);
      allocate_component_sobol(combinedMultiIndex);
      return max_index; break;
    }
    case COMBINED_SPARSE_GRID: { // product of two sums of tensor-product exp.
      size_t max_index = driverRep->maximal_grid();
      if (max_index != _NPOS) swap_data(max_index);

      // filter out dominated Smolyak multi-indices that don't contribute
      // to the definition of the product expansion
      size_t s, i, v, index, num_stored = storedMultiIndex.size();
      UShort2DArray curr_pareto; UShort3DArray stored_pareto(num_stored);
      CombinedSparseGridDriver* csg_driver
	= (CombinedSparseGridDriver*)driverRep;
      update_pareto_set(csg_driver->smolyak_multi_index(), curr_pareto);
      for (s=0; s<num_stored; ++s)
	update_pareto_set(csg_driver->stored_smolyak_multi_index(s),
			  stored_pareto[s]);
      // define a combined multi-index that can enumerate each pareto term;
      // since these are indices and not orders, we exclude the upper bound.
      UShortArray pareto_terms(num_stored+1); UShort2DArray pareto_indices;
      for (s=0; s<num_stored; ++s)
	pareto_terms[s] = stored_pareto[s].size();
      pareto_terms[num_stored] = curr_pareto.size();
      tensor_product_multi_index(pareto_terms, pareto_indices, false);
      // now enumerate all combinations of the non-dominated Smolyak index sets
      size_t num_combinations = pareto_indices.size();
      UShortArray exp_order_prod, stored_order;
      UShort2DArray tp_multi_index_prod;
      for (i=0; i<num_combinations; ++i) {
	UShortArray& pareto_ind_i = pareto_indices[i];
	index = pareto_ind_i[num_stored];
	sparse_grid_level_to_expansion_order(csg_driver, curr_pareto[index],
					     exp_order_prod);
	for (s=0; s<num_stored; ++s) {
	  index = pareto_ind_i[s];
	  sparse_grid_level_to_expansion_order(csg_driver,
					       stored_pareto[s][index],
					       stored_order);
	  for (v=0; v<numVars; ++v)
	    exp_order_prod[v] += stored_order[v];
	}
	// overlay each product expansion from the tensor-product combinations
	tensor_product_multi_index(exp_order_prod, tp_multi_index_prod);
	append_multi_index(tp_multi_index_prod, combinedMultiIndex);
      }	  
      allocate_component_sobol(combinedMultiIndex);
      return max_index; break;
    }
    default:
      // base class version supports product of two total-order expansions
      return SharedOrthogPolyApproxData::pre_combine_data(combine_type);
      break;
    }
    break;
  }
  case ADD_MULT_COMBINE:
    // base class manages this placeholder
    return SharedOrthogPolyApproxData::pre_combine_data(combine_type);
    break;
  }
}


void SharedProjectOrthogPolyApproxData::
sparse_grid_multi_index(CombinedSparseGridDriver* csg_driver,
			UShort2DArray& multi_index)
{
  const UShort2DArray& sm_multi_index = csg_driver->smolyak_multi_index();
  size_t i, num_smolyak_indices = sm_multi_index.size();

  // assemble a complete list of individual polynomial coverage
  // defined from the linear combination of mixed tensor products
  multi_index.clear();
  tpMultiIndex.resize(num_smolyak_indices);
  tpMultiIndexMap.resize(num_smolyak_indices);
  tpMultiIndexMapRef.resize(num_smolyak_indices);
  UShortArray exp_order(numVars);
  for (i=0; i<num_smolyak_indices; ++i) {
    // regenerate i-th exp_order as collocKey[i] cannot be used in general case
    // (i.e., for nested rules GP, CC, F2, or GK).  Rather, collocKey[i] is to
    // be used only as the key to the collocation pts.
    sparse_grid_level_to_expansion_order(csg_driver, sm_multi_index[i],
					 exp_order);
    tensor_product_multi_index(exp_order, tpMultiIndex[i]);
    append_multi_index(tpMultiIndex[i], multi_index, tpMultiIndexMap[i],
		       tpMultiIndexMapRef[i]);
#ifdef DEBUG
    PCout << "level =\n" << sm_multi_index[i] << "expansion_order =\n"
	  << exp_order << "tp_multi_index =\n" << tpMultiIndex[i]
	  << "multi_index =\n" << multi_index << '\n';
#endif // DEBUG
  }

  /*
  case SPARSE_INT_TENSOR_SUM_EXP: case SPARSE_INT_RESTR_TENSOR_SUM_EXP: {
    // assemble a complete list of individual polynomial coverage
    // defined from the linear combination of mixed tensor products
    multi_index.clear();
    UShort2DArray tp_multi_index;
    UShortArray int_order(numVars), exp_order(numVars);
    // Note: restricted rule growth within the sparse grid point set is separate
    // from restricted definition of the expansion terms.  The former makes
    // integrand precision more uniform by delaying exponential sequences, but
    // may still contain some nonuniformity.  The latter may enforce expansion
    // uniformity that would not otherwise be present based on integrand
    // precision alone, in order to reduce the possibility of a response-basis
    // product landing in the concave interior of the integrand resolution.
    short exp_growth = (sparseGridExpansion == SPARSE_INT_RESTR_TENSOR_SUM_EXP)
      ? MODERATE_RESTRICTED_GROWTH : UNRESTRICTED_GROWTH;
    for (i=0; i<num_smolyak_indices; ++i) {
      sparse_grid_level_to_expansion_order(csg_driver, sm_multi_index[i],
                                           exp_order,  exp_growth);
      tensor_product_multi_index(exp_order, tp_multi_index);
      append_multi_index(tp_multi_index, multi_index);
#ifdef DEBUG
      PCout << "level =\n" << sm_multi_index[i] << "integrand order =\n"
	    << int_order << "expansion order =\n" << exp_order << '\n';
	  //<< "tp_multi_index =\n" << tp_multi_index
	  //<< "multi_index =\n" << multi_index << '\n';
#endif // DEBUG
    }
    break;
  }
  case SPARSE_INT_TOTAL_ORD_EXP: {
    // back out approxOrder & use total_order_multi_index()
    UShortArray quad_order(numVars), integrand_order(numVars);
    UShort2DArray pareto(1), total_pareto;
    for (i=0; i<num_smolyak_indices; ++i) {
      csg_driver->level_to_order(sm_multi_index[i], quad_order);
      quadrature_order_to_integrand_order(driverRep, quad_order,
                                          integrand_order);
      // maintain an n-dimensional Pareto front of nondominated multi-indices
      pareto[0] = integrand_order;
      update_pareto_set(pareto, total_pareto);
#ifdef DEBUG
      PCout << "level =\n" << sm_multi_index[i] << "\nquad_order =\n"
	    << quad_order << "\nintegrand_order =\n" << integrand_order << '\n';
#endif // DEBUG
    }
#ifdef DEBUG
    PCout << "total_pareto =\n" << total_pareto << '\n';
#endif // DEBUG

    // first pass: compute max isotropic integrand that fits within Pareto front
    unsigned short order = 0;
    integrand_order.assign(numVars, order);
    bool total_order_dominated = true;
    while (total_order_dominated) {
      // calculate all nondominated polynomials for a total-order expansion
      pareto.clear();
      total_order_multi_index(integrand_order, pareto, 0);
      total_order_dominated = assess_dominance(pareto, total_pareto);
#ifdef DEBUG
      PCout << "integrand_order =\n" << integrand_order << "pareto =\n"
	    << pareto << "total_order_dominated = " << total_order_dominated
	    << '\n';
#endif // DEBUG
      // could increment/decrement by 2's due to expansion_order conversion,
      // but the actual resolvable integrand order is typically odd.
      if (total_order_dominated)
	++order; // advance to next test level
      else
	--order; // exiting loop: rewind to last successful
      integrand_order.assign(numVars, order);
    }
#ifdef DEBUG
    PCout << "Isotropic integrand_order =\n" << integrand_order << '\n';
#endif // DEBUG
    integrand_order_to_expansion_order(integrand_order, approxOrder);
    total_order_multi_index(approxOrder, multi_index);
    break;
  }
  case SPARSE_INT_HEUR_TOTAL_ORD_EXP: // early heuristic
    heuristic_sparse_grid_level_to_expansion_order(csg_driver->level(),
						   approxOrder);
    total_order_multi_index(approxOrder, multi_index);
    break;
  }
  */
}


/* This approach reduces memory requirements but must perform additional
   calculation to regenerate the tp_multi_index instances (previously
   generated in sparse_grid_multi_index()).  Currently, these tp_multi_index
   instances are stored in tpMultiIndex for later use in compute_coefficients().
void SharedProjectOrthogPolyApproxData::
map_tensor_product_multi_index(UShort2DArray& tp_multi_index, size_t tp_index)
{
  const SizetArray& tp_mi_map = tpMultiIndexMap[tp_index];
  size_t i, num_tp_terms = tp_mi_map.size();
  tp_multi_index.resize(num_tp_terms);
  for (i=0; i<num_tp_terms; ++i)
    tp_multi_index[i] = multiIndex[tp_mi_map[i]];
}
*/


Real SharedProjectOrthogPolyApproxData::
tensor_product_value(const RealVector& x, const RealVector& tp_coeffs,
		     const UShortArray& approx_order,
		     const UShort2DArray& tp_mi, RealVector& accumulator)
{
  unsigned short ao_0 = approx_order[0], ao_j, mi_i0, mi_ij;
  size_t i, j, num_tp_coeffs = tp_coeffs.length();
  BasisPolynomial& poly_0 = polynomialBasis[0]; Real x0 = x[0];
  for (i=0; i<num_tp_coeffs; ++i) {
    const UShortArray& tp_mi_i = tp_mi[i]; mi_i0 = tp_mi_i[0];
    if (ao_0)
      accumulator[0] += (mi_i0) ? tp_coeffs[i] * poly_0.type1_value(x0, mi_i0)
	                        : tp_coeffs[i];
    else
      accumulator[0]  = tp_coeffs[i];
    if (mi_i0 == ao_0) {
      // accumulate sums over variables with max key value
      for (j=1; j<numVars; ++j) {
	mi_ij = tp_mi_i[j]; ao_j = approx_order[j];
	if (ao_j)
	  accumulator[j] += (mi_ij) ? accumulator[j-1] *
	    polynomialBasis[j].type1_value(x[j], mi_ij) : accumulator[j-1];
	else
	  accumulator[j]  = accumulator[j-1];
	accumulator[j-1] = 0.;
	if (mi_ij != ao_j)
	  break;
      }
    }
  }
  Real tp_val = accumulator[numVars-1];
  accumulator[numVars-1] = 0.;
  return tp_val;
}


/*
Real SharedProjectOrthogPolyApproxData::
tensor_product_value(const RealVector& x, const RealVector& tp_coeffs,
		     const UShortArray& approx_order,
		     const UShort2DArray& tp_mi, RealVector& accumulator)
{
  //PCout << "test\n";
  unsigned short ao_0 = approx_order[0], ao_j, mi_i0, mi_ij;
  size_t i, j, num_tp_coeffs = tp_coeffs.length();
  BasisPolynomial& poly_0 = polynomialBasis[0]; Real x0 = x[0];
  Teuchos::SerialDenseVector<unsigned short,unsigned short> 
    max_order_1d( numVars );
  std::vector< std::set<unsigned short> > orders_1d( numVars );
  for (i=0; i<num_tp_coeffs; ++i) {
    const UShortArray& tp_mi_i = tp_mi[i];
    for (j=0; j<numVars; ++j) {
      max_order_1d[j] = std::max( max_order_1d[j], tp_mi_i[j] );
      orders_1d[j].insert( tp_mi_i[j] );
    }
  }
  std::vector< RealVector > bases_1d( numVars );
  std::set<unsigned short>::iterator it;
  for (j=0; j<numVars; ++j) {
    bases_1d[j].size( max_order_1d[j] );
    for ( it = orders_1d[j].begin(); it != orders_1d[j].end(); ++it )
      bases_1d[j][*it] = polynomialBasis[j].type1_value( x[j], *it );
  }

  for (i=0; i<num_tp_coeffs; ++i) {
    const UShortArray& tp_mi_i = tp_mi[i]; mi_i0 = tp_mi_i[0];
    if (ao_0)
    //accumulator[0] += (mi_i0) ? tp_coeffs[i] * poly_0.type1_value(x0, mi_i0)
      //: tp_coeffs[i];
      accumulator[0] += (mi_i0) ? tp_coeffs[i] * bases_1d[0][mi_i0] :
	tp_coeffs[i];
    else
      accumulator[0]  = tp_coeffs[i];
    if (mi_i0 == ao_0) {
      // accumulate sums over variables with max key value
      for (j=1; j<numVars; ++j) {
	mi_ij = tp_mi_i[j]; ao_j = approx_order[j];
	if (ao_j)
	  accumulator[j] += (mi_ij) ? accumulator[j-1] *
	    //polynomialBasis[j].type1_value(x[j], mi_ij) : accumulator[j-1];
	    bases_1d[j][mi_ij] : accumulator[j-1];
	else
	  accumulator[j]  = accumulator[j-1];
	accumulator[j-1] = 0.;
	if (mi_ij != ao_j)
	  break;
      }
    }
  }
  Real tp_val = accumulator[numVars-1];
  accumulator[numVars-1] = 0.;
  return tp_val;
}
*/

} // namespace Pecos
