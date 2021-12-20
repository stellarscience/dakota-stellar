/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 IncrementalSparseGridDriver
//- Description: Implementation code for IncrementalSparseGridDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "IncrementalSparseGridDriver.hpp"
#include "SharedPolyApproxData.hpp"
#include "sandia_sgmgg.hpp"
#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: IncrementalSparseGridDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void IncrementalSparseGridDriver::
initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
		const MultivariateDistribution& u_dist,
		const ExpansionConfigOptions& ec_options,
		BasisConfigOptions& bc_options, short growth_rate,
		bool track_uniq_prod_wts)
{
  SparseGridDriver::initialize_grid(ssg_level, dim_pref, u_dist, ec_options,
				    bc_options, growth_rate);
  trackUniqueProdWeights = track_uniq_prod_wts;

  // set compute1D{Points,Type1Weights,Type2Weights}
  //initialize_rule_pointers();
  // set levelGrowthToOrder
  //initialize_growth_pointers();
}


void IncrementalSparseGridDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis)
{
  IntegrationDriver::initialize_grid(poly_basis);

  // set compute1D{Points,Type1Weights,Type2Weights}
  //initialize_rule_pointers();
  // set levelGrowthToOrder
  //initialize_growth_pointers();
}


void IncrementalSparseGridDriver::
update_smolyak_arrays(UShort2DArray& sm_mi, IntArray& sm_coeffs)
{
  // if update indicator is off, Smolyak arrays already updated in grid_size()
  if (numPtsIter->second)
    return;

  // compute new Smolyak multi-index and coefficients, but don't overwrite old
  // (must avoid shifts due to index lower bound)
  // TO DO: remove lower bound in assign_smolyak_arrays()?  --> still requires
  // care since initial grid could be truncated and want to preserve this.
  UShort2DArray new_sm_mi; IntArray new_sm_coeffs;
  assign_smolyak_arrays(new_sm_mi, new_sm_coeffs);


  /*
  // Simple / robust increment but expensive (repeated linear-time look-ups)
  //
  // Old smolyak indices must be preserved but coeffs may be zeroed out.
  // New smolyak indices must be appended, irregardless of level ordering,
  // so that delta_coeff calculations are synched.
  if (new_sm_mi.size() >= sm_mi.size()) {
    size_t i, num_old = sm_mi.size(), num_new = new_sm_mi.size(), old_index;
    sm_coeffs.assign(num_old, 0); // zero out prior to updates from new active
    for (i=0; i<num_new; ++i) {
      UShortArray& new_sm_mi_i = new_sm_mi[i];
      old_index = find_index(sm_mi, new_sm_mi_i);
      if (old_index == _NPOS) { // not found: augment sm_mi and sm_coeffs
	sm_mi.push_back(new_sm_mi_i);
	sm_coeffs.push_back(new_sm_coeffs[i]);
      }
      else // found: retain and update coeff
	sm_coeffs[old_index] = new_sm_coeffs[i];
    }
  }
  else {
    // ...
  }
  */


  // Some assumptions to accelerate increment:


  // 1. See update_smolyak_arrays_aniso() for implementation that accelerates
  //    with assumptions but is still general enough for anisotropic refinement


  // 2. Stronger assumptions for uniform refinement
  // > no interior differences (due to anisotropic refinement) --> only need
  //   to manage differences at front and end
  UShort2DArray::iterator sm_it = sm_mi.begin();  size_t old_index = 0;
  const UShortArray& first_new = new_sm_mi[0];
  // Check for old truncated from beginning of new (preserve in updated)
  while (sm_it != sm_mi.end() && *sm_it != first_new)
    { ++sm_it; ++old_index; }
  // truncate to unmatched leading sets
  sm_mi.resize(old_index);  sm_coeffs.resize(old_index);
  sm_coeffs.assign(old_index, 0); // zero out prior to append of new active
  // augment with active new
  sm_mi.insert(sm_mi.end(), new_sm_mi.begin(), new_sm_mi.end());
  sm_coeffs.insert(sm_coeffs.end(), new_sm_coeffs.begin(), new_sm_coeffs.end());


  /* 3. Another efficient option: unroll and tailor assign_smolyak_arrays():
  unsigned short ssg_lev = ssgLevIter->second;
  UShortArray levels;
  if (isotropic())
    levels.assign(numVars, ssg_lev);
  else
    ;
  SharedPolyApproxData::
    total_order_multi_index(levels, new_sm_mi, ssg_lev-1);//don't know prev lev!
  sm_mi.insert(sm_mi.end(), new_sm_mi.begin(), new_sm_mi.end());
  */
}


void IncrementalSparseGridDriver::
update_smolyak_arrays_aniso(UShort2DArray& sm_mi, IntArray& sm_coeffs)
{
  // if update indicator is off, Smolyak arrays already updated in grid_size()
  if (numPtsIter->second)
    return;

  // compute new Smolyak multi-index and coefficients, but don't overwrite old
  // (must avoid shifts due to index lower bound)
  // TO DO: remove lower bound in assign_smolyak_arrays()?  --> still requires
  // care since initial grid could be truncated and want to preserve this.
  UShort2DArray new_sm_mi; IntArray new_sm_coeffs;
  assign_smolyak_arrays(new_sm_mi, new_sm_coeffs);

  size_t i, num_old = sm_mi.size(), num_new = new_sm_mi.size();
  sm_coeffs.assign(num_old, 0); // zero out old prior to updates from new active
  size_t start = 0, end = num_old;
  for (i=0; i<num_new; ++i) {
    UShortArray& new_sm_mi_i = new_sm_mi[i];
    size_t index = find_index(sm_mi, new_sm_mi_i, start, end);// restrict search
    if (index == _NPOS) {
      sm_mi.push_back(new_sm_mi_i);
      sm_coeffs.push_back(new_sm_coeffs[i]);
    }
    else {
      sm_coeffs[index] = new_sm_coeffs[i];
      // trim contiguous search range when straightforward to do so
      // (this minimal performance optimization is easily defeated)
      if      (index == start) ++start;
      else if (index == end-1) --end;
    }
  }

  /*
  // Assumptions to accelerate increment (relative to simplest option above):
  // > consistent ordering for index sets present in old & new
  // > leading sm_mi index sets that are unmatched at front of new_sm_mi
  //   were truncated from new due to numVars-1 lower bound
  //   --> retain and set coeff to zero.
  // > unmatched internal/trailing indices from new_sm_mi are augmentations
  //   (see bounds enforcement in anisotropic adaptation)
  UShort2DArray::iterator            sm_it =     sm_mi.begin();
  UShort2DArray::const_iterator new_sm_cit = new_sm_mi.begin();
  size_t i, num_old = sm_mi.size(), num_new = new_sm_mi.size(), old_index = 0;
  sm_coeffs.assign(num_old, 0); // zero out old prior to updates from new active
  // Check for old truncated from beginning of new (preserve in updated)
  while (sm_it != sm_mi.end() && *sm_it != *new_sm_cit)
    { ++sm_it; ++old_index; } //sm_coeffs[old_index] = 0;
  // Scan new_sm_mi to update and append sm_{mi,coeffs}
  for (i=0; i<num_new; ++i) {
    UShortArray& new_sm_mi_i = new_sm_mi[i];
    if (old_index < num_old && sm_mi[old_index] == new_sm_mi_i) {
      sm_coeffs[old_index] = new_sm_coeffs[i];
      ++old_index;
    }
    else {
      sm_mi.push_back(new_sm_mi_i);
      sm_coeffs.push_back(new_sm_coeffs[i]);
    }
  }
  // Retired 4/19/21: these simplifying assumptions can be violated, e.g. level
  // 1 isotropic mi {0,0},{1,0},{0,1} w/ coeffs {-1,1,1} incremented to level 2
  // anisotropic mi {0,0},{2,0},{0,1} w/ same coeffs
  */
}


void IncrementalSparseGridDriver::update_collocation_key()
{
  const UShort2DArray& sm_mi = smolMIIter->second;
  UShort3DArray& colloc_key = collocKeyIter->second;
  UShortArray quad_order(numVars);
  size_t i, start_index = colloc_key.size(), num_sm_mi = sm_mi.size();
  colloc_key.resize(num_sm_mi);
  for (i=start_index; i<num_sm_mi; ++i) {
    level_to_order(sm_mi[i], quad_order);
    SharedPolyApproxData::
      tensor_product_multi_index(quad_order, colloc_key[i], false);
  }
}


int IncrementalSparseGridDriver::grid_size()
{
  int& num_colloc_pts = numPtsIter->second;
  if (num_colloc_pts == 0) { // special value indicated update required
    update_smolyak_arrays();
    update_collocation_key();

    RealMatrix a1_pts, a1_t2_wts;  RealVector a1_t1_wts;
    const UShort2DArray& sm_mi = smolMIIter->second;
    compute_tensor_points_weights(sm_mi, collocKeyIter->second, 0, sm_mi.size(),
				  true, a1_pts, a1_t1_wts, a1_t2_wts);

    int m = numVars, n1 = a1_pts.numCols(), seed = 1234567;
    RealVector zv(m, false), r1v(n1, false);  IntArray sind1(n1);
    bool* is_unique1 = new bool[n1];
    // Use _count_inc1 whereas reference_unique uses _index_inc1
    webbur::point_radial_tol_unique_count_inc1(m, n1, a1_pts.values(),
      duplicateTol, &seed, zv.values(), r1v.values(), &sind1[0], is_unique1,
      &num_colloc_pts);
    delete [] is_unique1;
  }
  return num_colloc_pts;
}


void IncrementalSparseGridDriver::compute_grid()
{
  // Note: incremental and combined sparse grid definitions use different point
  // orderings --> IncrementalSparseGridDriver uses a separate implementation
  // to define its reference grid, rather than inheriting from CombinedSGDriver

  update_smolyak_arrays();  // smolyak{MultiIndex,Coeffs}
  update_collocation_key(); // collocKey
  reference_unique();       // compute pts,wts for the reference grid
  update_reference();       // update reference arrays

#ifdef DEBUG
  PCout << "IncrementalSparseGridDriver::compute_grid() results:\n"
	<< "uniqueIndexMapping:\n" << uniqIndMapIter->second << "\nvar_sets:\n";
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


void IncrementalSparseGridDriver::compute_trial_grid(RealMatrix& var_sets)
{
  /* Old approach now embedded within increment_unique():

  // compute trial variable/weight sets and update collocKey
  UShortArray quad_order(numVars);
  const UShortArray& tr_set = trial_set();
  level_to_order(tr_set, quad_order);
  UShort2DArray new_key;
  UShort3DArray& colloc_key = collocKeyIter->second;
  colloc_key.push_back(new_key); // empty array updated in place
  compute_tensor_grid(quad_order, tr_set, a2PIter->second, a2T1WIter->second,
                      a2T2WIter->second, colloc_key.back());
  */

  // trial set already appended, update collocation key
  size_t last_index = collocKeyIter->second.size();
  update_collocation_key(); // needed for compute_tensor_points_weights()
  // compute a2 pts/wts; update collocIndices, uniqueIndexMapping
  increment_unique(last_index);
  // update var_sets with increment (not aggregate) of unique points from a2
  increment_sparse_points(collocIndIter->second, last_index,
			  isUniq2Iter->second, numUniq1Iter->second,
			  a2PIter->second, var_sets);

#ifdef DEBUG
  PCout << "compute_trial_grid():\nunique variable sets:\n" << var_sets;
#endif // DEBUG

  // track trial sets that have been evaluated (do here since
  // increment_smolyak_multi_index() used for both new trials and restorations)
  //computedTrialSets[activeKey].push_back(trial_set());
}


void IncrementalSparseGridDriver::compute_increment(RealMatrix& var_sets)
{
  update_smolyak_arrays();  // update smolyak{MultiIndex,Coeffs}
  update_collocation_key(); // synchronize collocKey
  // update a2 for multiple trial sets
  size_t start_index = smolyakCoeffsRef[activeKey].size();
  increment_unique(start_index);
  // update var_sets with increment (not aggregate) of unique points from a2
  increment_sparse_points(collocIndIter->second, start_index,
			  isUniq2Iter->second, numUniq1Iter->second,
			  a2PIter->second, var_sets);
}


void IncrementalSparseGridDriver::push_increment()
{
  update_smolyak_arrays();  // update smolyak{MultiIndex,Coeffs}
  update_collocation_key(); // synchronize collocKey
  // update a2 for multiple trial sets
  size_t start_index = smolyakCoeffsRef[activeKey].size();
  increment_unique(start_index, false);
}


void IncrementalSparseGridDriver::pop_increment()
{
  IntArray& sm_coeffs_ref = smolyakCoeffsRef[activeKey];
  size_t ref_size = sm_coeffs_ref.size();
  smolMIIter->second.resize(ref_size);
  smolCoeffsIter->second = sm_coeffs_ref;
  collocKeyIter->second.resize(ref_size);
  collocIndIter->second.resize(ref_size);

  numPtsIter->second = numUniq1Iter->second;             // unique ref points
  // pruning of uniqueIndexMapping is not strictly required (it is updated on
  // demand prior to updating collocation indices), but good for completeness
  uniqIndMapIter->second.resize(a1PIter->second.numCols()); // all ref points

  // restore reference points and weights
  varSetsIter->second.reshape(numVars, numPtsIter->second);// = a1PIter->second;
  pop_weights();
}


void IncrementalSparseGridDriver::initialize_sets()
{
  // define set O (old) from smolyakMultiIndex and smolyakCoeffs:
  unsigned short     ssg_lev = ssgLevIter->second;
  const UShort2DArray& sm_mi = smolMIIter->second;
  const IntArray&  sm_coeffs = smolCoeffsIter->second;
  UShortArraySet&     old_mi = oldMultiIndex[activeKey];
  //old_mi = sm_mi;
  old_mi.clear(); old_mi.insert(sm_mi.begin(), sm_mi.end());

  //poppedTrialSets[activeKey].clear();  // cleared in finalize_sets()
  //activeMultiIndex[activeKey].clear(); // cleared in finalize_sets()

  // compute initial set A (active) by applying add_active_neighbors()
  // to the frontier of smolyakMultiIndex:
  size_t i, num_old_sets = sm_coeffs.size();
  // anisotropic test on coeff==1 is necessary but not sufficient for presence
  // on index set frontier, requiring an additional logic test within
  // add_active_neighbors().  For anisotropic, the weighted norm of the index
  // set may differ from the level --> need to compute Pareto set.
  bool dim_iso = isotropic();
  for (i=0; i<num_old_sets; ++i)
    if ( sm_coeffs[i] == 1 && ( !dim_iso || // imperfect for aniso
	 ( dim_iso && l1_norm(sm_mi[i]) == ssg_lev ) ) )
      add_active_neighbors(sm_mi[i], dim_iso);

#ifdef DEBUG
  PCout << "IncrementalSparseGridDriver::initialize_sets():\n  active key:\n"
	<< activeKey << "\n  sm_mi:\n" << sm_mi << "\n  sm_coeffs:\n"
	<< sm_coeffs << "\n  ssg level = " << ssg_lev << "\n  active sets:\n"
	<< activeMultiIndex[activeKey] << std::endl;
#endif // DEBUG
}


void IncrementalSparseGridDriver::
increment_smolyak_multi_index(const UShortArray& set)
{
  UShort2DArray& sm_mi = smolMIIter->second;
  size_t last_index = sm_mi.size();
  sm_mi.push_back(set);

  // update smolyakCoeffs from smolyakMultiIndex
  update_smolyak_coefficients(last_index);

  // collocKey, collocIndices, and uniqueIndexMapping updated within either
  // push_set() or compute_trial_grid()
}


void IncrementalSparseGridDriver::push_set()
{
  // SparseGridDriver currently retains no memory, so updates are recomputed

  // store pushIndex for use by other classes
  UShortArrayDeque& pop_trials = poppedTrialSets[activeKey];
  size_t p_index = find_index(pop_trials, trial_set());
  if (p_index != _NPOS) pop_trials.erase(pop_trials.begin() + p_index);
  pushIndex[activeKey] = p_index;
 
  // synchronize collocKey with smolyakMultiIndex
  update_collocation_key();
  // compute a2; update collocIndices, uniqueIndexMapping
  // no new var_sets and 1D updates have already been performed
  increment_unique(smolMIIter->second.size()-1, false); // no 1D pts/wts update
}


void IncrementalSparseGridDriver::pop_set()
{
  UShort2DArray& sm_mi = smolMIIter->second;
  poppedTrialSets[activeKey].push_back(sm_mi.back());
  pushIndex[activeKey] = _NPOS;

  // restore reference grid state
  sm_mi.pop_back();
  collocKeyIter->second.pop_back();
  collocIndIter->second.pop_back();
  smolCoeffsIter->second = smolyakCoeffsRef[activeKey];

  numPtsIter->second = numUniq1Iter->second;             // unique ref points
  // pruning of uniqueIndexMapping is not strictly required (it is updated on
  // demand prior to updating collocation indices), but good for completeness
  uniqIndMapIter->second.resize(a1PIter->second.numCols()); // all ref points

  // restore reference points and weights
  varSetsIter->second.reshape(numVars, numPtsIter->second);// = a1PIter->second;
  pop_weights();
}


void IncrementalSparseGridDriver::
finalize_sets(bool output_sets, bool converged_within_tol, bool reverted)
{
  // For final answer, push all evaluated sets into old and clear active.
  // Multiple trial insertion approach must be compatible with bookkeeping
  // elsewhere (e.g., Dakota::Approximation), i.e., inc2/inc3 set insertions
  // occur one at a time without mixing.

  UShort2DArray&         sm_mi =         smolMIIter->second;
  UShortArrayDeque& pop_trials = poppedTrialSets[activeKey];
  size_t i, start_index = sm_mi.size();
  // don't insert activeMultiIndex, as this may include sets which have not
  // been evaluated (due to final update_sets() call) -- use poppedTrialSets.
  sm_mi.insert(sm_mi.end(), pop_trials.begin(), pop_trials.end());

  /*
  // finalizeIndex allows external clients to synchronize with
  // poppedTrialSets index ordering
  size_t num_pop_tr = pop_trials.size();
  SizetArray& f_indices = finalizeIndex[activeKey];
  // in latest design, f_indices mapping is trivial
  f_indices.resize(num_pop_tr);
  for (i=0; i<num_pop_tr; ++i)
    f_indices[i] = i;//find_index(pop_trials, pop_mi[i]);
  */
  activeMultiIndex[activeKey].clear();  pop_trials.clear();

  // update smolyakCoeffs from smolyakMultiIndex
  update_smolyak_coefficients(start_index);
  // synchronize collocKey with smolyakMultiIndex
  update_collocation_key();
  // generate final grid, uniqueIndexMapping, collocIndices, numCollocPts
  finalize_unique(start_index);

  if (output_sets) {
    size_t i, j, num_sm_mi = sm_mi.size();
    if (converged_within_tol) {
      // if not reverted, trial set was below tolerance at convergence
      size_t start_below = (reverted) ? start_index : start_index - 1;
      PCout << "Above tolerance index sets:\n";
      for (i=0; i<start_below; ++i)
	print_index_set(PCout, sm_mi[i]);
      PCout << "Below tolerance index sets:\n";
      for (i=start_below; i<num_sm_mi; ++i)
	print_index_set(PCout, sm_mi[i]);
    }
    else {
      PCout << "Final index sets:\n";
      for (i=0; i<num_sm_mi; ++i)
	print_index_set(PCout, sm_mi[i]);
    }
  }
}


void IncrementalSparseGridDriver::
increment_unique_points_weights(size_t start_index, const UShort2DArray& sm_mi,
  const IntArray& sm_coeffs, const IntArray& sm_coeffs_ref,
  const UShort3DArray& colloc_key, Sizet2DArray& colloc_ind,
  int& num_colloc_pts, RealMatrix& a1_pts, RealVector& a1_t1w,
  RealMatrix& a1_t2w,  RealMatrix& a2_pts, RealVector& a2_t1w,
  RealMatrix& a2_t2w, RealVector& zv, RealVector& r1v, RealVector& r2v,
  IntArray& sind1, BitArray& isu1, IntArray& uind1, IntArray& uset1,
  int& num_u1, IntArray& sind2, BitArray& isu2, IntArray& uind2,
  IntArray& uset2, int& num_u2, IntArray& unique_index_map,
  bool update_1d_pts_wts, RealMatrix& pts, RealVector& t1_wts,
  RealMatrix& t2_wts)
{
  size_t i, j, num_sm_mi = sm_mi.size();
  int m = numVars, n1 = a1_pts.numCols(), tp_n2, n2 = 0;
  RealVector tp_t1w; RealMatrix tp_pts, tp_t2w;
  bool *is_unique1, *is_unique2;

  // compute the points/weights for each tensor grid, updating 1D if needed
  for (i=start_index; i<num_sm_mi; ++i) {
    compute_tensor_points_weights(sm_mi, colloc_key, i, 1, update_1d_pts_wts,
				  tp_pts, tp_t1w, tp_t2w);
    tp_n2 = tp_pts.numCols();
    a2_pts.reshape(numVars, n2+tp_n2);  a2_t1w.resize(n2+tp_n2);
    if (computeType2Weights) a2_t2w.reshape(numVars, n2+tp_n2);
    for (j=0; j<tp_n2; ++j) {
      copy_data(tp_pts[j], numVars, a2_pts[n2+j]);
      a2_t1w[n2+j] = tp_t1w[j];
      if (computeType2Weights)
	copy_data(tp_t2w[j], numVars, a2_t2w[n2+j]);
    }
    n2 += tp_n2;
  }

  if (!n2) return;

  // ----
  // INC2
  // ----
  // INC2 detects a2 duplication with respect to a1 reference grid as well as
  // internal duplicates (multiple TP grids in a2 have internal duplication).
  // This allows for a single update spanning multiple Smolyak index sets.
  r2v.sizeUninitialized(n2);  sind2.resize(n2);
  uset2.resize(n2);           uind2.resize(n2);
  is_unique1 = new bool[n1];  copy_data(isu1, is_unique1, n1);
  is_unique2 = new bool[n2];  // bridges inc2 to inc3: isUnique2 not needed
  webbur::point_radial_tol_unique_index_inc2(m, n1, a1_pts.values(), n2,
    a2_pts.values(), duplicateTol, zv.values(), r1v.values(), &sind1[0],
    is_unique1,  num_u1, &uset1[0], &uind1[0],  r2v.values(), &sind2[0],
    is_unique2, &num_u2, &uset2[0], &uind2[0]);
#ifdef DEBUG
  PCout << "Increment unique: numUnique2 = " << num_u2 << "\na2 =\n";
  write_data(PCout, a2_pts, false, true, true);
  PCout << "               r2   indx2 unique2   undx2   xdnu2:\n";
  for (j=0; j<num_u2; ++j)
    PCout << std::setw(17) << r2v[j]        << std::setw(8) << sind2[j]
	  << std::setw(8)  << is_unique2[j] << std::setw(8) << uset2[j]
	  << std::setw(8)  << uind2[j]      << '\n';
  for (j=num_u2; j<n2; ++j)
    PCout << std::setw(17) << r2v[j]       << std::setw(8)  << sind2[j]
	  << std::setw(8) << is_unique2[j] << std::setw(16) << uind2[j] << '\n';
  PCout << std::endl;
#endif // DEBUG
  copy_data(is_unique2, n2, isu2);
  delete [] is_unique1; delete [] is_unique2;

  num_colloc_pts = num_u1 + num_u2;
  update_unique_indices(start_index, num_u1, uind1, uset1, isu2, uind2, uset2,
			unique_index_map);
  assign_collocation_indices(colloc_key, unique_index_map, colloc_ind,
			     start_index);
  assign_sparse_points(colloc_ind, start_index, isu2, num_u1, a2_pts, pts);
  if (trackUniqueProdWeights)
    update_sparse_weights(start_index, colloc_key, colloc_ind, num_colloc_pts,
			  sm_coeffs, sm_coeffs_ref, a1_t1w, a1_t2w, a2_t1w,
			  a2_t2w, t1_wts, t2_wts);
}


void IncrementalSparseGridDriver::
merge_unique_points_weights(const UShort2DArray& sm_mi,
  const IntArray& sm_coeffs, const IntArray& sm_coeffs_ref,
  const UShort3DArray& colloc_key, Sizet2DArray& colloc_ind,
  int& num_colloc_pts, RealMatrix& a1_pts, RealVector& a1_t1w,
  RealMatrix& a1_t2w, RealMatrix& a2_pts, RealVector& a2_t1w,
  RealMatrix& a2_t2w, RealVector& r1v, RealVector& r2v, IntArray& sind1,
  BitArray& isu1, IntArray& uind1, IntArray& uset1, int& num_u1,
  IntArray& sind2, BitArray& isu2, IntArray& uind2, IntArray& uset2,
  int& num_u2, IntArray& unique_index_map, RealMatrix& pts,
  RealVector& t1_wts, RealMatrix& t2_wts)
{
  int i, m = numVars, n1 = a1_pts.numCols(), n2 = a2_pts.numCols(),
    n1n2 = n1+n2, n3;
  RealVector r3v(n1n2, false);
  RealMatrix a3_pts(m, n1n2, false);
  IntArray sind3(n1n2), uset3(n1n2), uind3(n1n2);
  bool *is_unique1 = new bool[n1], *is_unique2 = new bool[n2],
       *is_unique3 = new bool[n1n2];
  copy_data(isu1, is_unique1, n1);
  copy_data(isu2, is_unique2, n2);

  // ----
  // INC3
  // ----
  webbur::point_radial_tol_unique_index_inc3(m, n1, a1_pts.values(),
    r1v.values(), &sind1[0], is_unique1, num_u1, &uset1[0], &uind1[0], n2,
    a2_pts.values(), r2v.values(), &sind2[0], is_unique2, num_u2, &uset2[0],
    &uind2[0], &n3, a3_pts.values(), r3v.values(), &sind3[0], is_unique3,
    &num_colloc_pts, &uset3[0], &uind3[0]);

#ifdef DEBUG
  PCout << "Merge unique: num_unique3 = " << num_colloc_pts << "\na3 =\n";
  write_data(PCout, a3_pts, false, true, true);
  PCout << "               r3   indx3 unique3   undx3   xdnu3:\n";
  for (size_t i=0; i<num_colloc_pts; ++i)
    PCout << std::setw(17) << r3v[i]        << std::setw(8) << sind3[i]
	  << std::setw(8)  << is_unique3[i] << std::setw(8) << uset3[i]
	  << std::setw(8)  << uind3[i] << '\n';
  for (size_t i=num_colloc_pts; i<n1n2; ++i)
    PCout << std::setw(17) << r3v[i]       << std::setw(8)  << sind3[i]
	  << std::setw(8) << is_unique3[i] << std::setw(16) << uind3[i] << '\n';
  PCout << std::endl;
#endif // DEBUG

  // Need to increment again as pop operations need to restore previous state
  // after a non-permanent increment
  /* MSE, 3/19/2019: seems redundant with same steps in increment_unique():
     > all cases of merge_unique follow either increment or push operations
     > only need to redo this is if merge ordering was modified by inc3
       (if needed, check _inc3() ordering by activating DEBUG block above)
  size_t start_index = sm_coeffs_ref.size();
  update_unique_indices(start_index, num_u1, uind1, uset1, isu2, uind2, uset2,
			unique_index_map);
  assign_collocation_indices(colloc_key, unique_index_map, colloc_ind,
			     start_index);
  assign_sparse_points(colloc_ind, start_index, isu2, num_u1, a2_pts, pts);
  if (trackUniqueProdWeights)
    update_sparse_weights(start_index, colloc_key, colloc_ind, num_colloc_pts,
			  sm_coeffs, sm_coeffs_ref, a1_t1w, a1_t2w, a2_t1w,
			  a2_t2w, t1_wts, t2_wts);
  */
  // Promote a3 to a1: update a1 reference points/weights
  //a1_pts = a3_pts; // equivalent, but potentially more copy overhead
  a1_pts.reshape(numVars, n1n2);
  for (i=n1; i<n1n2; ++i)
    copy_data(a3_pts[i], numVars, a1_pts[i]);
  if (trackUniqueProdWeights) {
    a1_t1w.resize(n1n2);
    if (computeType2Weights) a1_t2w.reshape(numVars, n1n2);
    for (i=0; i<n2; ++i) {
      a1_t1w[n1+i] = a2_t1w[i];
      if (computeType2Weights)
	copy_data(a2_t2w[i], numVars, a1_t2w[n1+i]);
    }
  }
  // Promote a3 to a1: update reference indices, counts, radii
  num_u1 = num_colloc_pts;
  r1v = r3v;  sind1 = sind3;  uset1 = uset3;  uind1 = uind3;
  copy_data(is_unique3, n1n2, isu1);
  delete [] is_unique1; delete [] is_unique2; delete [] is_unique3;
}


void IncrementalSparseGridDriver::
update_unique_indices(size_t start_index, int num_uniq1, const IntArray& xdnu1,
		      const IntArray& undx1, const BitArray& isu2,
		      const IntArray& xdnu2, const IntArray& undx2,
		      IntArray& unique_index_map)
{
  int xdnu_j, a1_a2_index, a1_index, new_cntr = num_uniq1;

  size_t i, n1 = xdnu1.size(), n2 = xdnu2.size();
  unique_index_map.resize(n1+n2);

  // first pass assigns unique indices
  for (i=0; i<n2; ++i)
    if (isu2[i])
      unique_index_map[n1+i] = new_cntr++;
  // second pass refers back to unique indices and can be a forward reference
  // (dictating two passes)
  for (i=0; i<n2; ++i)
    if (!isu2[i]) {
      // XDNU2[N2] in point_radial_tol_unique_index_inc2() [sandia_rules.cpp]:
      //   If the value represents an index in UNDX2, this can be inferred by
      //   the fact that its value is >= UNIQUE_NUM1.  To reference UNDX2, the
      //   value should then be decremented by UNIQUE_NUM1.
      xdnu_j = xdnu2[i];
      // UNDX2[UNIQUE_NUM2] in point_radial_tol_unique_index_inc2():
      //   The index in A2 of the tolerably unique points, incremented by N1
      // Note: xdnu --> undx --> recovers original point ordering
      if (xdnu_j >= num_uniq1) {
	a1_a2_index = undx2[xdnu_j - num_uniq1]; // - n1 for a2_index
	unique_index_map[n1+i] = unique_index_map[a1_a2_index];
      }
      else {
	a1_index = undx1[xdnu_j];
	unique_index_map[n1+i] = unique_index_map[a1_index];
      }
    }

#ifdef DEBUG
  PCout << "Incremented map:\n" << unique_index_map;
#endif // DEBUG
}


void IncrementalSparseGridDriver::
increment_sparse_points(const Sizet2DArray& colloc_ind, size_t start_index,
			const BitArray& raw_is_unique,
			size_t colloc_index_offset,// 0 (ref) or num_u1 (incr)
			const RealMatrix& raw_incr_pts,
			RealMatrix&    unique_incr_pts)
{
  // update sizes
  size_t num_unique_pts = raw_is_unique.count();
  unique_incr_pts.shapeUninitialized(numVars, num_unique_pts);

  size_t i, j, cntr = 0, uniq_index, num_sm_mi = colloc_ind.size(), num_tp_pts;
  for (i=start_index; i<num_sm_mi; ++i) {
    const SizetArray& colloc_ind_i = colloc_ind[i];
    num_tp_pts = colloc_ind_i.size();
    for (j=0; j<num_tp_pts; ++j, ++cntr) {
      if (raw_is_unique[cntr]) {
	uniq_index = colloc_ind_i[j] - colloc_index_offset;
	copy_data(raw_incr_pts[cntr], numVars, unique_incr_pts[uniq_index]);
      }
    }
  }
}


void IncrementalSparseGridDriver::
update_sparse_weights(size_t start_index, const UShort3DArray& colloc_key,
		      const Sizet2DArray& colloc_ind, int num_colloc_pts,
		      const IntArray& sm_coeffs, const IntArray& sm_coeffs_ref,
		      const RealVector& a1_t1_wts, const RealMatrix& a1_t2_wts,
		      const RealVector& a2_t1_wts, const RealMatrix& a2_t2_wts,
		      RealVector& unique_t1_wts, RealMatrix& unique_t2_wts)
{
  unique_t1_wts = type1WeightSetsRef[activeKey]; // to be augmented
  if (computeType2Weights)
    unique_t2_wts = type2WeightSetsRef[activeKey]; // to be augmented

  // update sizes
  int delta_coeff, sm_coeff;
  unique_t1_wts.resize(num_colloc_pts); // new entries initialized to 0
  if (computeType2Weights)
    unique_t2_wts.reshape(numVars, num_colloc_pts); // new entries init to 0

  // back out changes in Smolyak coeff for existing index sets
  size_t i, j, k, cntr, num_sm_mi = colloc_key.size(), num_tp_pts, uniq_index;
  for (i=0, cntr=0; i<start_index; ++i) {
    delta_coeff = sm_coeffs[i] - sm_coeffs_ref[i];
    if (delta_coeff) {
      num_tp_pts = colloc_key[i].size();
      const SizetArray& colloc_ind_i = colloc_ind[i];
      for (j=0; j<num_tp_pts; ++j, ++cntr) {
	uniq_index = colloc_ind_i[j];
	// assign tensor weights to unique weights
	unique_t1_wts[uniq_index] += delta_coeff * a1_t1_wts[cntr];
	if (computeType2Weights) {
	  Real*     uniq_t2_wts_j = unique_t2_wts[uniq_index];
	  const Real* a1_t2_wts_j = a1_t2_wts[cntr];
	  for (k=0; k<numVars; ++k)
	    uniq_t2_wts_j[k] += delta_coeff * a1_t2_wts_j[k];
	}
      }
    }
    else
      cntr += colloc_key[i].size();
  }
  // add contributions for new index sets
  add_sparse_weights(start_index, colloc_key, colloc_ind, sm_coeffs,
		     a2_t1_wts, a2_t2_wts, unique_t1_wts, unique_t2_wts);

#ifdef DEBUG
  PCout << "updated type1 weight sets =\n" << unique_t1_wts;
  if (computeType2Weights) {
    PCout << "updated type2 weight sets =\n";
    write_data(PCout, unique_t2_wts, false, true, true);
  }
#endif // DEBUG
}


void IncrementalSparseGridDriver::clear_inactive()
{
  CombinedSparseGridDriver::clear_inactive();

  std::map<ActiveKey, int>::iterator nu1_it = numUnique1.begin();
  std::map<ActiveKey, int>::iterator nu2_it = numUnique2.begin();
  std::map<ActiveKey, RealVector>::iterator z_it = zVec.begin();
  std::map<ActiveKey, RealVector>::iterator r1_it = r1Vec.begin();
  std::map<ActiveKey, RealVector>::iterator r2_it = r2Vec.begin();
  std::map<ActiveKey, RealMatrix>::iterator a1p_it = a1Points.begin();
  std::map<ActiveKey, RealVector>::iterator a11w_it = a1Type1Weights.begin();
  std::map<ActiveKey, RealMatrix>::iterator a12w_it = a1Type2Weights.begin();
  std::map<ActiveKey, RealMatrix>::iterator a2p_it = a2Points.begin();
  std::map<ActiveKey, RealVector>::iterator a21w_it = a2Type1Weights.begin();
  std::map<ActiveKey, RealMatrix>::iterator a22w_it = a2Type2Weights.begin();
  std::map<ActiveKey, IntArray>::iterator si1_it = sortIndex1.begin();
  std::map<ActiveKey, IntArray>::iterator si2_it = sortIndex2.begin();
  std::map<ActiveKey, IntArray>::iterator us1_it = uniqueSet1.begin();
  std::map<ActiveKey, IntArray>::iterator us2_it = uniqueSet2.begin();
  std::map<ActiveKey, IntArray>::iterator ui1_it = uniqueIndex1.begin();
  std::map<ActiveKey, IntArray>::iterator ui2_it = uniqueIndex2.begin();
  std::map<ActiveKey, BitArray>::iterator iu1_it = isUnique1.begin();
  std::map<ActiveKey, BitArray>::iterator iu2_it = isUnique2.begin();
  std::map<ActiveKey, IntArray>::iterator uim_it = uniqueIndexMapping.begin();

  std::map<ActiveKey, IntArray>::iterator scr_it = smolyakCoeffsRef.begin();
  std::map<ActiveKey, RealVector>::iterator t1r_it
    = type1WeightSetsRef.begin();
  std::map<ActiveKey, RealMatrix>::iterator t2r_it
    = type2WeightSetsRef.begin();

  while (a1p_it != a1Points.end())
    if (a1p_it == a1PIter) { // preserve active
      ++nu1_it; ++nu2_it; ++z_it; ++r1_it; ++r2_it; ++a1p_it; ++a11w_it;
      ++a12w_it; ++a2p_it; ++a21w_it; ++a22w_it; ++si1_it; ++si2_it; ++us1_it;
      ++us2_it; ++ui1_it; ++ui2_it; ++iu1_it; ++iu2_it; ++uim_it;
      ++scr_it; //++pmi_it;
      if (trackUniqueProdWeights)
	{ ++t1r_it; if (computeType2Weights) ++t2r_it; }
    }
    else { // clear inactive: postfix increments manage iterator invalidations
      numUnique1.erase(nu1_it++);         numUnique2.erase(nu2_it++);
      zVec.erase(z_it++); r1Vec.erase(r1_it++); r2Vec.erase(r2_it++);
      a1Points.erase(a1p_it++);           a1Type1Weights.erase(a11w_it++);
      a1Type2Weights.erase(a12w_it++);    a2Points.erase(a2p_it++);
      a2Type1Weights.erase(a21w_it++);    a2Type2Weights.erase(a22w_it++);
      sortIndex1.erase(si1_it++);         sortIndex2.erase(si2_it++);
      uniqueSet1.erase(us1_it++);         uniqueSet2.erase(us2_it++);
      uniqueIndex1.erase(ui1_it++);       uniqueIndex2.erase(ui2_it++);
      isUnique1.erase(iu1_it++);          isUnique2.erase(iu2_it++);
      uniqueIndexMapping.erase(uim_it++); smolyakCoeffsRef.erase(scr_it++);
      if (trackUniqueProdWeights) {
	type1WeightSetsRef.erase(t1r_it++);
	if (computeType2Weights) type2WeightSetsRef.erase(t2r_it++);
      }
    }
}

} // namespace Pecos
