/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 HierarchSparseGridDriver
//- Description: Wrapper class for C++ code from packages/quadrature/sparse_grid
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef HIERARCH_SPARSE_GRID_DRIVER_HPP
#define HIERARCH_SPARSE_GRID_DRIVER_HPP

#include "SparseGridDriver.hpp"

namespace Pecos {


/// Derived integration driver class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis, Newton-Cotes, and Gaussian quadrature rules
    within Smolyak sparse grids. */

class HierarchSparseGridDriver: public SparseGridDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  HierarchSparseGridDriver();
  /// constructor
  HierarchSparseGridDriver(unsigned short ssg_level,
			   const RealVector& dim_pref = RealVector(),
			   short growth_rate = MODERATE_RESTRICTED_GROWTH,
			   short refine_control = NO_CONTROL);
  /// destructor
  ~HierarchSparseGridDriver();

  //
  //- Heading: Virtual function redefinitions
  //

  void compute_grid(RealMatrix& var_sets);
  int grid_size();

  void store_grid(size_t index = _NPOS);
  void restore_grid(size_t index = _NPOS);
  void remove_stored_grid(size_t index = _NPOS);
  void clear_stored();

  size_t maximal_grid() const;
  void swap_grid(size_t index);

  void initialize_sets();
  void push_trial_set(const UShortArray& set);
  void restore_set();
  void compute_trial_grid(RealMatrix& var_sets);
  void pop_trial_set();
  //void merge_set();
  void finalize_sets(bool output_sets, bool converged_within_tol);

  const UShortArray& trial_set() const;
  int unique_trial_points() const;

  void compute_grid_increment(RealMatrix& var_sets);

  void print_smolyak_multi_index() const;

  // concatenate type1WeightSets for use in abstract integration functions
  //const RealVector& type1_weight_sets(); // const;
  // concatenate type2WeightSets for use in abstract integration functions
  //const RealMatrix& type2_weight_sets(); // const;

  UShortUShortPair level_to_delta_pair(size_t i, unsigned short lev_i);
  unsigned short level_to_delta_size(size_t i, unsigned short lev_i);
  void level_to_delta_key(size_t i, unsigned short lev_i,
			  UShortArray& delta_key_i);

  /// convert a Smolyak index set into hierarchical quadrature keys
  void levels_to_delta_keys(const UShortArray& levels,
			    UShort2DArray& delta_keys);
  /// convert a Smolyak index set into the sizes of hierarchical
  /// quadrature increments
  void levels_to_delta_sizes(const UShortArray& levels,
			     UShortArray& delta_sizes);

  //
  //- Heading: Member functions
  //

  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
    const ShortArray& u_types, const ExpansionConfigOptions& ec_options,
    BasisConfigOptions& bc_options,
    short growth_rate = MODERATE_RESTRICTED_GROWTH,
    bool track_colloc_indices = true);

  /// return incrementSets
  const UShortArray& increment_sets() const;

  /// return smolyakMultiIndex
  const UShort3DArray& smolyak_multi_index() const;

  /// set trackCollocIndices
  void track_collocation_indices(bool track_colloc_indices);
  /// get trackCollocIndices
  bool track_collocation_indices() const;

  /// return collocKey
  const UShort4DArray& collocation_key() const;
  /// return collocIndices
  const Sizet3DArray& collocation_indices() const;

  /// return storedLevMultiIndex
  const UShort3DArray& stored_smolyak_multi_index(size_t index) const;
  /// return storedCollocKey
  const UShort4DArray& stored_collocation_key(size_t index) const;
  // return storedCollocIndices
  //const Sizet3DArray& stored_collocation_indices(size_t index) const;

  /// discriminate portions of the level-set hierarchy that are
  /// reference sets from those in the current increment
  void partition_keys(UShort2DArray& reference_set_range,
		      UShort2DArray& increment_set_range) const;
  /// discriminate portions of the level-set-point hierarchy that are
  /// in the reference grid from those in the current increment
  void partition_keys(UShort3DArray& reference_pt_range,
		      UShort3DArray& increment_pt_range) const;

  /// return type1WeightSets for use in hierarchical integration functions
  const RealVector2DArray& type1_weight_set_arrays() const;
  /// return type2WeightSets for use in hierarchical integration functions
  const RealMatrix2DArray& type2_weight_set_arrays() const;

private:

  //
  //- Heading: Convenience functions
  //

  void update_smolyak_multi_index(bool clear_sm_mi = false);
  void assign_collocation_key();
  void update_collocation_key();
  void assign_collocation_indices();
  void update_collocation_indices();

  /// kernel routine used for trial set and full sparse grid computations
  void compute_points_weights(RealMatrix& pts, RealVector& t1_wts,
			      RealMatrix& t2_wts, const UShortArray& sm_index,
			      const UShort2DArray& colloc_key);
  /// compute points and weights for a trial set
  void compute_points_weights(RealMatrix& pts, RealVector& t1_wts,
			      RealMatrix& t2_wts);
  /// compute points and weights for all levels of the (initial) sparse grid
  void compute_points_weights(RealMatrix& pts, RealVector2DArray& t1_wts,
			      RealMatrix2DArray& t2_wts);

  //
  //- Heading: Data
  //

  /// flag for use of fully nested 1D rules, allowing formulation using
  /// collocation point increments
  bool nestedGrid;

  /// interpolation depth by index set by numVars array for identifying
  /// the index to use within the polynomialBasis for a particular variable
  /** The index sets correspond to j (0-based) for use as indices, which
      are offset from the i indices (1-based) normally used in the Smolyak
      expressions.  The indices correspond to levels, one within each
      anisotropic tensor-product integration of a Smolyak recursion. */
  UShort3DArray smolyakMultiIndex;

  /// level of trial evaluation set from push_trial_set(); trial set
  /// corresponds to smolyakMultiIndex[trialLevel].back()
  unsigned short trialLevel;
  /// identifies the trailing index set increments within smolyakMultiIndex
  /// due to an isotropic/anistropic grid refinement
  UShortArray incrementSets;

  /// due to the hierarchical structure, collocation indices only need
  /// to be defined in special cases (e.g., generalized sparse grids
  /// for which index sets can appear in different orders).
  bool trackCollocIndices;

  /// levels-by-index sets-by-numDeltaPts-by-numVars array for identifying
  /// the 1-D point indices for sets of tensor-product collocation points
  UShort4DArray collocKey;
  /// levels-by-index sets-by-numTensorProductPts array for linking the
  /// set of tensor products to the unique collocation points evaluated
  Sizet3DArray collocIndices;

  /// stored driver states: copies of smolyakMultiIndex
  UShort4DArray storedLevMultiIndex;
  /// stored driver states: copies of collocKey
  UShort5DArray storedCollocKey;
  // stored driver states: copies of collocIndices
  //Sizet4DArray storedCollocIndices;

  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each point in the sparse grid
  RealVector2DArray type1WeightSets;
  /// the set of type2 weights (for integration of gradient interpolants)
  /// for each derivative component and for each point in the sparse grid
  RealMatrix2DArray type2WeightSets;

  // concatenation of type1WeightSets RealVector2DArray into a RealVector
  //RealVector concatT1WeightSets;
  // concatenation of type2WeightSets RealMatrix2DArray into a RealMatrix
  //RealMatrix concatT2WeightSets;

  /// stored driver state: copy of type1WeightSets
  RealVector3DArray storedType1WeightSets;
  /// stored driver state: copy of type2WeightSets
  RealMatrix3DArray storedType2WeightSets;

  /// type 1 weight sets popped during decrement for later restoration
  /// to type1WeightSets
  std::map<UShortArray, RealVector> poppedT1WtSets;
  /// type 2 weight sets popped during decrement for later restoration
  /// to type2WeightSets
  std::map<UShortArray, RealMatrix> poppedT2WtSets;
};


inline HierarchSparseGridDriver::HierarchSparseGridDriver():
  SparseGridDriver(), nestedGrid(true), trackCollocIndices(true)
{ }


inline HierarchSparseGridDriver::
HierarchSparseGridDriver(unsigned short ssg_level, const RealVector& dim_pref,
			 short growth_rate, short refine_control):
  SparseGridDriver(ssg_level, dim_pref, growth_rate, refine_control),
  nestedGrid(true), trackCollocIndices(true)
{ }


inline HierarchSparseGridDriver::~HierarchSparseGridDriver()
{ }


inline const UShortArray& HierarchSparseGridDriver::trial_set() const
{ return smolyakMultiIndex[trialLevel].back(); }


inline int HierarchSparseGridDriver::unique_trial_points() const
{ return collocKey[trialLevel].back().size(); }


inline const UShortArray& HierarchSparseGridDriver::increment_sets() const
{ return incrementSets; }


inline void HierarchSparseGridDriver::print_smolyak_multi_index() const
{
  size_t i, j, k, num_lev = smolyakMultiIndex.size(), cntr = 1;
  for (i=0; i<num_lev; ++i) {
    const UShort2DArray& sm_mi_i = smolyakMultiIndex[i];
    size_t num_sets = sm_mi_i.size();
    for (j=0; j<num_sets; ++j, ++cntr) {
      PCout << "Smolyak index set " << cntr << ':';
      print_index_set(PCout, sm_mi_i[j]);
    }
  }
}


inline const UShort3DArray& HierarchSparseGridDriver::
smolyak_multi_index() const
{ return smolyakMultiIndex; }


inline void HierarchSparseGridDriver::
track_collocation_indices(bool track_colloc_indices)
{ trackCollocIndices = track_colloc_indices; }


inline bool HierarchSparseGridDriver::track_collocation_indices() const
{ return trackCollocIndices; }


inline const UShort4DArray& HierarchSparseGridDriver::collocation_key() const
{ return collocKey; }


inline const Sizet3DArray& HierarchSparseGridDriver::collocation_indices() const
{ return collocIndices; }


inline const UShort3DArray& HierarchSparseGridDriver::
stored_smolyak_multi_index(size_t index) const
{ return storedLevMultiIndex[index]; }


inline const UShort4DArray& HierarchSparseGridDriver::
stored_collocation_key(size_t index) const
{ return storedCollocKey[index]; }


//inline const Sizet3DArray& HierarchSparseGridDriver::
//stored_collocation_indices(size_t index) const
//{ return storedCollocIndices[index]; }


/*
inline const RealVector& HierarchSparseGridDriver::type1_weight_sets() // const
{
  if (concatT1WeightSets.length() != numCollocPts) {
    concatT1WeightSets.sizeUninitialized(numCollocPts);
    size_t lev, set, pt, cntr = 0, num_levels = type1WeightSets.size(),
      num_sets, num_tp_pts;
    for (lev=0; lev<num_levels; ++lev) {
      num_sets = type1WeightSets[lev].size();
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = type1WeightSets[lev][set].length();
	const SizetArray& colloc_index = collocIndices[lev][set];
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr)
	  concatT1WeightSets[colloc_index[cntr]]
	    = type1WeightSets[lev][set][pt];
      }
    }
  }
  return concatT1WeightSets;
}


inline const RealMatrix& HierarchSparseGridDriver::type2_weight_sets() // const
{
  if (concatT2WeightSets.numCols() != numCollocPts) {
    concatT2WeightSets.shapeUninitialized(numVars, numCollocPts);
    size_t lev, set, pt, v, cntr = 0, num_levels = type2WeightSets.size(),
      num_sets, num_tp_pts;
    for (lev=0; lev<num_levels; ++lev) {
      num_sets = type2WeightSets[lev].size();
      for (set=0; set<num_sets; ++set) {
	num_tp_pts = type2WeightSets[lev][set].numCols();
	const SizetArray& colloc_index = collocIndices[lev][set];
	for (pt=0; pt<num_tp_pts; ++pt, ++cntr) {
	  Real* concat_t2_wts = concatT2WeightSets[colloc_index[cntr]];
	  const Real*  t2_wts = type2WeightSets[lev][set][pt];
	  for (v=0; v<numVars; ++v)
	    concat_t2_wts[v] = t2_wts[v];
	}
      }
    }
  }
  return concatT2WeightSets;
}
*/


inline const RealVector2DArray& HierarchSparseGridDriver::
type1_weight_set_arrays() const
{ return type1WeightSets; }


inline const RealMatrix2DArray& HierarchSparseGridDriver::
type2_weight_set_arrays() const
{ return type2WeightSets; }


inline void HierarchSparseGridDriver::
levels_to_delta_sizes(const UShortArray& levels, UShortArray& delta_sizes)
{
  size_t i, num_lev = levels.size();
  if (delta_sizes.size() != num_lev)
    delta_sizes.resize(num_lev);
  for (i=0; i<num_lev; ++i)
    delta_sizes[i] = level_to_delta_size(i, levels[i]);
}


inline void HierarchSparseGridDriver::
levels_to_delta_keys(const UShortArray& levels, UShort2DArray& delta_keys)
{
  size_t i, num_lev = levels.size();
  if (delta_keys.size() != num_lev)
    delta_keys.resize(num_lev);
  for (i=0; i<num_lev; ++i)
    level_to_delta_key(i, levels[i], delta_keys[i]);
}

} // namespace Pecos

#endif
