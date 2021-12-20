/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 SparseGridDriver
//- Description: Wrapper class for C++ code from packages/quadrature/sparse_grid
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef COMBINED_SPARSE_GRID_DRIVER_HPP
#define COMBINED_SPARSE_GRID_DRIVER_HPP

#include "SparseGridDriver.hpp"

namespace Pecos {

/// pointer to a collocation point or weight evaluation function, matching
/// the GWPointer prototype required by Pecos/packages/VPISparseGrid
typedef void ( *CollocFnPtr ) ( int order, int index, double* data );
/// pointer to a level-growth-to-order mapping function, matching the
/// GWPointer2 prototype required by Pecos/packages/VPISparseGrid
typedef int ( *LevGrwOrdFnPtr ) ( int level, int growth );

class ActiveKey;


/// Derived integration driver class that generates N-dimensional
/// Smolyak sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

/** This class is used by Dakota::NonDSparseGrid, but could also be
    used for general numerical integration of moments.  It employs 1-D
    Clenshaw-Curtis, Newton-Cotes, and Gaussian quadrature rules
    within Smolyak sparse grids. */

class CombinedSparseGridDriver: public SparseGridDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  CombinedSparseGridDriver();
  /// constructor
  CombinedSparseGridDriver(unsigned short ssg_level,
			   const RealVector& dim_pref = RealVector(),
			   short growth_rate = MODERATE_RESTRICTED_GROWTH,
			   short refine_control = NO_CONTROL);
  /// destructor
  ~CombinedSparseGridDriver();

  //
  //- Heading: Virtual function redefinitions
  //

  /// update {smolMI,smolCoeffs,collocKey,collocInd}Iter from activeKey
  void update_active_iterators();

  void compute_grid();
  void compute_grid(RealMatrix& var_sets);
  int grid_size();
  void combine_grid();
  void combined_to_active(bool clear_combined);

  const RealMatrix& variable_sets() const;
  const RealVector& type1_weight_sets() const;
  const RealMatrix& type2_weight_sets() const;
  const RealMatrix& variable_sets(const ActiveKey& key) const;
  const RealVector& type1_weight_sets(const ActiveKey& key) const;
  const RealMatrix& type2_weight_sets(const ActiveKey& key) const;

  const RealMatrix& combined_variable_sets() const;
  const RealVector& combined_type1_weight_sets() const;
  const RealMatrix& combined_type2_weight_sets() const;

  void reinterpolated_tensor_grid(const UShortArray& lev_index,
				  const SizetList& reinterp_indices);

  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis);
  void initialize_grid_parameters(const MultivariateDistribution& mv_dist);

  void clear_inactive();
  void clear_keys();

  const ActiveKey& maximal_grid();

  void print_smolyak_multi_index() const;

  //
  //- Heading: Member functions
  //

  /// initialize all sparse grid settings except for distribution params
  void initialize_grid(unsigned short ssg_level, const RealVector& dim_pref,
    const MultivariateDistribution& u_dist,
    const ExpansionConfigOptions& ec_options, BasisConfigOptions& bc_options,
    short growth_rate = MODERATE_RESTRICTED_GROWTH, bool track_colloc = false,
    bool track_uniq_prod_wts = true);

  /// overloaded form initializes smolyakMultiIndex and smolyakCoeffs
  void assign_smolyak_arrays();
  /// initialize collocKey from smolyakMultiIndex
  void assign_collocation_key();
  /// initialize collocIndices from collocKey and unique_index_map
  void assign_collocation_indices(const UShort3DArray& colloc_key,
				  const IntArray& unique_index_map,
				  Sizet2DArray& colloc_indices,
				  size_t start_index = 0);

  /// set duplicateTol based on the content of collocRules: table lookups will
  /// generally be more precise/repeatable than numerically-generated rules
  void initialize_duplicate_tolerance();

  /// return smolyakMultiIndex[activeKey]
  const UShort2DArray& smolyak_multi_index() const;
  /// return smolyakMultiIndex[key]
  const UShort2DArray& smolyak_multi_index(const ActiveKey& key) const;
  /// return smolyakMultiIndex
  const std::map<ActiveKey, UShort2DArray>& smolyak_multi_index_map() const;

  /// return smolyakCoeffs[activeKey]
  const IntArray& smolyak_coefficients() const;
  /// return smolyakCoeffs[key]
  const IntArray& smolyak_coefficients(const ActiveKey& key) const;

  /// set trackCollocDetails
  void track_collocation_details(bool track_colloc);
  /// get trackCollocDetails
  bool track_collocation_details() const;
  /// set trackUniqueProdWeights
  void track_unique_product_weights(bool track_uniq_prod_wts);
  /// get trackUniqueProdWeights
  bool track_unique_product_weights() const;

  /// return collocKey[activeKey]
  const UShort3DArray& collocation_key() const;
  /// return collocKey[key]
  const UShort3DArray& collocation_key(const ActiveKey& key) const;
  /// return collocIndices[activeKey]
  const Sizet2DArray& collocation_indices() const;
  /// return collocIndices[key]
  const Sizet2DArray& collocation_indices(const ActiveKey& key) const;

  // return duplicateTol
  //Real duplicate_tolerance() const;

  /// return combinedSmolyakMultiIndex
  const UShort2DArray& combined_smolyak_multi_index() const;
  // return combinedSmolyakMultiIndexMap
  //const Sizet2DArray& combined_smolyak_multi_index_map() const;
  /// return combinedCollocKey
  const UShort3DArray& combined_collocation_key() const;

protected:

  //
  //- Heading: Member functions
  //

  /// initialize Smolyak multi-index (index sets defining the set of tensor
  /// products) and Smolyak combinatorial coefficients using an isotropic or
  /// anisotropic index set constraint.  For anisotropic, webbur::sgmga_vcn_*
  /// functions are used to compute index sets satisfying the anisotropic
  /// index set constraint, along with their corresponding coefficients.
  void assign_smolyak_arrays(UShort2DArray& multi_index, IntArray& coeffs);
  /// initialize a collocation key from a smolyak multi-index
  void assign_collocation_key(const UShort2DArray& sm_mi,
			      UShort3DArray& colloc_key);

  /// overloaded form updates smolyakCoeffs from smolyakMultiIndex
  void update_smolyak_coefficients(size_t start_index);
  /// update the coeffs array based on new trailing index sets within
  /// multi_index for incrementally generated generalized sparse grids
  void update_smolyak_coefficients(size_t start_index,
				   const UShort2DArray& sm_mi,
				   IntArray& sm_coeffs);

  /// remove Smolyak multi-indices with zero coefficients
  void prune_inactive(UShort2DArray& sm_mi, IntArray& sm_coeffs);

  /// compute points and type1,2 integration weights for a sparse grid
  /// defined by level and (optionally) anisotropic weights
  void compute_unique_points_weights(unsigned short ssg_lev,
				     const RealVector& aniso_wts,
				     int num_colloc_pts,
				     IntArray& unique_index_map,
				     RealMatrix& var_sets, RealVector& t1_wts,
				     RealMatrix& t2_wts);
  /// compute points and type1,2 integration weights for a sparse grid defined
  /// by an arbitrary multi-index (more general than level + aniso weights)
  void compute_unique_points_weights(const UShort2DArray& sm_mi,
				     const IntArray& sm_coeffs,
				     const UShort3DArray& colloc_key,
				     IntArray& unique_index_map,
				     bool update_1d_pts_wts,
				     RealMatrix& var_sets, RealVector& t1_wts,
				     RealMatrix& t2_wts);
  /// modular helper for public reference_unique(RealMatrix&)
  void compute_unique_points_weights(const UShort2DArray& sm_mi,
    const IntArray& sm_coeffs, const UShort3DArray& colloc_key,
    Sizet2DArray& colloc_ind, int& num_colloc_pts, RealMatrix& a1_pts,
    RealVector& a1_t1w, RealMatrix& a1_t2w, RealVector& zv, RealVector& r1v,
    IntArray& sind1, BitArray& isu1, IntArray& uind1, IntArray& uset1,
    int& num_u1, IntArray& unique_index_map, bool update_1d_pts_wts,
    RealMatrix& var_sets, RealVector& t1_wts, RealMatrix& t2_wts);

  /// aggregate point and weight sets across one or more tensor products
  void compute_tensor_points_weights(const UShort2DArray& sm_mi,
				     const UShort3DArray& colloc_key,
				     size_t start_index, size_t num_indices,
				     bool update_1d_pts_wts, RealMatrix& pts,
				     RealVector& t1_wts, RealMatrix& t2_wts);

  /// define the reference collocation indices
  void assign_unique_indices(const BitArray& isu1, const IntArray& xdnu1,
			     const IntArray& undx1, IntArray& unique_index_map);

  /// convenience function for updating sparse (unique) points from a set of
  /// aggregated (non-unique) tensor points
  void assign_sparse_points(const Sizet2DArray& colloc_ind, size_t start_index,
			    const BitArray& raw_is_unique,
			    size_t curr_unique_pts, const RealMatrix& raw_pts,
			    RealMatrix& unique_pts);

  /// convenience function for assigning sparse weights from a set of
  /// tensor weights
  void assign_sparse_weights(const UShort3DArray& colloc_key,
			     const Sizet2DArray& colloc_ind, int num_colloc_pts,
			     const IntArray& sm_coeffs,
			     const RealVector& a1_t1_wts,
			     const RealMatrix& a1_t2_wts,
			     RealVector& unique_t1_wts,
			     RealMatrix& unique_t2_wts);
  /// convenience function for updating sparse weights from overlaying
  /// a set of tensor weights
  void add_sparse_weights(size_t start_index, const UShort3DArray& colloc_key,
			  const Sizet2DArray& colloc_ind,
			  const IntArray& sm_coeffs, const RealVector& raw_t1w,
			  const RealMatrix& raw_t2w, RealVector& unique_t1w,
			  RealMatrix& unique_t2w);

  //
  //- Heading: Data
  //

  /// numSmolyakIndices-by-numVars array for identifying the index to use
  /// within the polynomialBasis for a particular variable
  /** The index sets correspond to j (0-based) for use as indices, which
      are offset from the i indices (1-based) normally used in the Smolyak
      expressions.  The indices correspond to levels, one for each
      anisotropic tensor-product grid within a Smolyak recursion. */
  std::map<ActiveKey, UShort2DArray> smolyakMultiIndex;
  /// iterator for active entry within smolyakMultiIndex
  std::map<ActiveKey, UShort2DArray>::iterator smolMIIter;

  /// array of Smolyak combinatorial coefficients, one for each tensor
  /// product index set; order is synchronized with smolyakMultiIndex
  std::map<ActiveKey, IntArray> smolyakCoeffs;
  /// iterator for active entry within smolyakCoeffs
  std::map<ActiveKey, IntArray>::iterator smolCoeffsIter;

  /// numSmolyakIndices-by-numTensorProductPts-by-numVars array for identifying
  /// the 1-D point indices for sets of tensor-product collocation points
  std::map<ActiveKey, UShort3DArray> collocKey;
  /// iterator for active entry within collocKey
  std::map<ActiveKey, UShort3DArray>::iterator collocKeyIter;

  /// numSmolyakIndices-by-numTensorProductPts array for linking the set of
  /// tensor products to the unique collocation points evaluated
  std::map<ActiveKey, Sizet2DArray> collocIndices;
  /// iterator for active entry within collocIndices
  std::map<ActiveKey, Sizet2DArray>::iterator collocIndIter;

  /// flag indicating need to track {type1,type2}WeightSets (product weights for
  /// each unique grid point) as opposed to relying on collections of 1D weights
  bool trackUniqueProdWeights;
  /// duplication tolerance used in sgmga routines
  Real duplicateTol;
  // maps indices and bases from sgmga_index() to collocation point index
  //IntArraySizetMap ssgIndexMap;

  /// mapping from points in set of tensor grids to unique sparse grid
  /// collocation points (unrolled array of collocation indices)
  std::map<ActiveKey, IntArray> uniqueIndexMapping;
  /// active entry within uniqueIndexMapping
  std::map<ActiveKey, IntArray>::iterator uniqIndMapIter;

  /// the set of unique collocation points in the sparse grid
  std::map<ActiveKey, RealMatrix> variableSets;
  /// iterator for active entry within variableSets
  std::map<ActiveKey, RealMatrix>::iterator varSetsIter;
  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each unique point in the sparse grid
  std::map<ActiveKey, RealVector> type1WeightSets;
  /// iterator for active entry within type1WeightSets
  std::map<ActiveKey, RealVector>::iterator t1WtIter;
  /// the set of type2 weights (for integration of gradient interpolants) for
  /// each derivative component and for each unique point in the sparse grid
  std::map<ActiveKey, RealMatrix> type2WeightSets;
  /// iterator for active entry within type2WeightSets
  std::map<ActiveKey, RealMatrix>::iterator t2WtIter;

  /// multi-index for maximal grid that is the result of combining a set
  /// of level expansions
  UShort2DArray combinedSmolyakMultiIndex;
  // mapping of terms when aggregating CombinedSparseGridDriver::
  // smolyakMultiIndex into combinedSmolyakMultiIndex in pre_combine_data()
  //Sizet2DArray combinedSmolyakMultiIndexMap;
  /// Smolyak coefficients corresponding to combinedSmolyakMultiIndex
  IntArray combinedSmolyakCoeffs;
  /// collocation key for maximal grid that is the result of combining a
  /// set of level expansions (defined from combinedSmolyakMultiIndex)
  UShort3DArray combinedCollocKey;
  /// mapping from combined sparse grid points to unique collocation points
  IntArray combinedUniqueIndexMap;

  /// variable sets for maximal grid defined by overlaying level grids
  /** Could also be managed within SurrogateData, but would require data
      sharing per PolynomialApproximation instance. */
  RealMatrix combinedVarSets;
  /// combination of CombinedSparseGridDriver::type1WeightSets, consistent
  /// with combination of level expansions
  RealVector combinedT1WeightSets;
  /// combination of CombinedSparseGridDriver::type2WeightSets, consistent
  /// with combination of level expansions
  RealMatrix combinedT2WeightSets;

private:

  //
  //- Heading: Convenience functions
  //

  /// function passed by pointer for computing collocation points for
  /// polynomialBasis[index]
  static void basis_collocation_points(int order, int index, double* data);
  /// function passed by pointer for computing type 1 collocation
  /// weights for polynomialBasis[index]
  static void basis_type1_collocation_weights(int order,int index,double* data);
  /// function passed by pointer for computing type 2 collocation
  /// weights for polynomialBasis[index]
  static void basis_type2_collocation_weights(int order,int index,double* data);

  /// initialize compute1D{Points,Type1Weights,Type2Weights} function pointer
  /// arrays for use within webbur::sgmg() and webbur::sgmga() routines
  void initialize_rule_pointers();
  /// initialize levelGrowthToOrder function pointer arrays for use within
  /// webbur::sgmg() and webbur::sgmga() routines
  void initialize_growth_pointers();

  //
  //- Heading: Data
  //

  /// pointer to instance of this class for use in static member functions
  static CombinedSparseGridDriver* sgdInstance;

  /// flag controls conditional population of collocKey, collocIndices,
  /// collocPts1D and type{1,2}CollocWts1D
  bool trackCollocDetails;

  /// array of pointers to collocation point evaluation functions
  std::vector<CollocFnPtr> compute1DPoints;
  /// array of pointers to type1 collocation weight evaluation functions
  std::vector<CollocFnPtr> compute1DType1Weights;
  // 2D array of pointers to type2 collocation weight evaluation functions
  //std::vector<std::vector<CollocFnPtr> > compute1DType2Weights;

  /// array of pointers to webbur::level_to_growth functions
  std::vector<LevGrwOrdFnPtr> levelGrowthToOrder;

  // store the key identified in the last call to maximal_grid()
  //ActiveKey maximalKey;
};


inline CombinedSparseGridDriver::CombinedSparseGridDriver():
  SparseGridDriver(), smolMIIter(smolyakMultiIndex.end()),
  trackCollocDetails(false), trackUniqueProdWeights(false), duplicateTol(1.e-15)
{
  // careful with invoking virtual fn from ctor
  // recursive lookups are prevented by leading if-checks
  CombinedSparseGridDriver::update_active_iterators();
}


inline CombinedSparseGridDriver::
CombinedSparseGridDriver(unsigned short ssg_level, const RealVector& dim_pref,
			 short growth_rate, short refine_control):
  SparseGridDriver(ssg_level, dim_pref, growth_rate, refine_control),
  smolMIIter(smolyakMultiIndex.end()), trackCollocDetails(false),
  trackUniqueProdWeights(false), duplicateTol(1.e-15)
{
  // careful with invoking virtual fn from ctor
  // recursive lookups are prevented by leading if-checks
  CombinedSparseGridDriver::update_active_iterators();
}


inline CombinedSparseGridDriver::~CombinedSparseGridDriver()
{ }


inline void CombinedSparseGridDriver::update_active_iterators()
{
  // Test for change
  if (smolMIIter != smolyakMultiIndex.end() && smolMIIter->first == activeKey)
    return;

  smolMIIter     = smolyakMultiIndex.find(activeKey);
  smolCoeffsIter = smolyakCoeffs.find(activeKey);
  collocKeyIter  = collocKey.find(activeKey);
  collocIndIter  = collocIndices.find(activeKey);
  uniqIndMapIter = uniqueIndexMapping.find(activeKey);
  varSetsIter    = variableSets.find(activeKey);
  t1WtIter       = type1WeightSets.find(activeKey);
  t2WtIter       = type2WeightSets.find(activeKey);

  /* So long as we only create new keys and avoid modifying existing ones,
     this deep copy is not needed.
  ActiveKey active_copy; // share 1 deep copy of current active key
  if (smolMIIter     == smolyakMultiIndex.end()  ||
      smolCoeffsIter == smolyakCoeffs.end()      ||
      collocKeyIter  == collocKey.end()          ||
      collocIndIter  == collocIndices.end()      ||
      uniqIndMapIter == uniqueIndexMapping.end() ||
      varSetsIter    == variableSets.end()       ||
      t1WtIter       == type1WeightSets.end()    ||
      t2WtIter       == type2WeightSets.end())
    active_copy = activeKey.copy();
  */

  if (smolMIIter == smolyakMultiIndex.end()) {
    std::pair<ActiveKey, UShort2DArray>
      u2a_pair(activeKey/*active_copy*/, UShort2DArray());
    smolMIIter = smolyakMultiIndex.insert(u2a_pair).first;
  }
  if (smolCoeffsIter == smolyakCoeffs.end()) {
    std::pair<ActiveKey, IntArray>
      ia_pair(activeKey/*active_copy*/, IntArray());
    smolCoeffsIter = smolyakCoeffs.insert(ia_pair).first;
  }
  if (collocKeyIter == collocKey.end()) {
    std::pair<ActiveKey, UShort3DArray>
      u3a_pair(activeKey/*active_copy*/, UShort3DArray());
    collocKeyIter = collocKey.insert(u3a_pair).first;
  }
  if (collocIndIter == collocIndices.end()) {
    std::pair<ActiveKey, Sizet2DArray>
      s2a_pair(activeKey/*active_copy*/, Sizet2DArray());
    collocIndIter = collocIndices.insert(s2a_pair).first;
  }
  if (uniqIndMapIter == uniqueIndexMapping.end()) {
    std::pair<ActiveKey, IntArray>
      ua_pair(activeKey/*active_copy*/, IntArray());
    uniqIndMapIter = uniqueIndexMapping.insert(ua_pair).first;
  }
  if (varSetsIter == variableSets.end()) {
    std::pair<ActiveKey, RealMatrix>
      rm_pair(activeKey/*active_copy*/, RealMatrix());
    varSetsIter = variableSets.insert(rm_pair).first;
  }
  if (t1WtIter == type1WeightSets.end()) {
    std::pair<ActiveKey, RealVector>
      rv_pair(activeKey/*active_copy*/, RealVector());
    t1WtIter = type1WeightSets.insert(rv_pair).first;
  }
  if (t2WtIter == type2WeightSets.end()) {
    std::pair<ActiveKey, RealMatrix>
      rm_pair(activeKey/*active_copy*/, RealMatrix());
    t2WtIter = type2WeightSets.insert(rm_pair).first;
  }

  SparseGridDriver::update_active_iterators();
}


inline void CombinedSparseGridDriver::clear_keys()
{
  SparseGridDriver::clear_keys();

  smolyakMultiIndex.clear();   smolMIIter = smolyakMultiIndex.end();
  smolyakCoeffs.clear();       smolCoeffsIter = smolyakCoeffs.end();

  collocKey.clear();           collocKeyIter = collocKey.end();
  collocIndices.clear();       collocIndIter = collocIndices.end();
  uniqueIndexMapping.clear();  uniqIndMapIter = uniqueIndexMapping.end();

  variableSets.clear();        varSetsIter = variableSets.end();
  type1WeightSets.clear();     t1WtIter = type1WeightSets.end();
  type2WeightSets.clear();     t2WtIter = type2WeightSets.end();
}


inline const UShort2DArray& CombinedSparseGridDriver::
smolyak_multi_index() const
{ return smolMIIter->second; }


inline const UShort2DArray& CombinedSparseGridDriver::
smolyak_multi_index(const ActiveKey& key) const
{
  std::map<ActiveKey, UShort2DArray>::const_iterator cit
    = smolyakMultiIndex.find(key);
  if (cit == smolyakMultiIndex.end()) {
    PCerr << "Error: key not found in CombinedSparseGridDriver::"
	  << "smolyak_multi_index()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const std::map<ActiveKey, UShort2DArray>& CombinedSparseGridDriver::
smolyak_multi_index_map() const
{ return smolyakMultiIndex; }


inline const IntArray& CombinedSparseGridDriver::smolyak_coefficients() const
{ return smolCoeffsIter->second; }


inline const IntArray& CombinedSparseGridDriver::
smolyak_coefficients(const ActiveKey& key) const
{
  std::map<ActiveKey, IntArray>::const_iterator cit = smolyakCoeffs.find(key);
  if (cit == smolyakCoeffs.end()) {
    PCerr << "Error: key not found in CombinedSparseGridDriver::"
	  << "smolyak_coefficients()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline void CombinedSparseGridDriver::
track_collocation_details(bool track_colloc)
{ trackCollocDetails = track_colloc; }


inline bool CombinedSparseGridDriver::track_collocation_details() const
{ return trackCollocDetails; }


inline void CombinedSparseGridDriver::
track_unique_product_weights(bool track_uniq_prod_wts)
{ trackUniqueProdWeights = track_uniq_prod_wts; }


inline bool CombinedSparseGridDriver::track_unique_product_weights() const
{ return trackUniqueProdWeights; }


inline const UShort3DArray& CombinedSparseGridDriver::collocation_key() const
{ return collocKeyIter->second; }


inline const UShort3DArray& CombinedSparseGridDriver::
collocation_key(const ActiveKey& key) const
{
  std::map<ActiveKey, UShort3DArray>::const_iterator cit
    = collocKey.find(key);
  if (cit == collocKey.end()) {
    PCerr << "Error: key not found in CombinedSparseGridDriver::"
	  << "collocation_key()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const Sizet2DArray& CombinedSparseGridDriver::collocation_indices() const
{ return collocIndIter->second; }


inline const Sizet2DArray& CombinedSparseGridDriver::
collocation_indices(const ActiveKey& key) const
{
  std::map<ActiveKey, Sizet2DArray>::const_iterator cit
    = collocIndices.find(key);
  if (cit == collocIndices.end()) {
    PCerr << "Error: key not found in CombinedSparseGridDriver::"
	  << "collocation_indices()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


//inline Real CombinedSparseGridDriver::duplicate_tolerance() const
//{ return duplicateTol; }


inline void CombinedSparseGridDriver::print_smolyak_multi_index() const
{
  const UShort2DArray& sm_mi = smolMIIter->second;
  const IntArray&  sm_coeffs = smolCoeffsIter->second;
  size_t i, sm_mi_len = sm_mi.size(), cntr = 0;
  for (i=0; i<sm_mi_len; ++i)
    if (sm_coeffs[i]) {
      PCout << "Smolyak index set " << ++cntr << " (coeff = "
	    << sm_coeffs[i] << "):";
      print_index_set(PCout, sm_mi[i]);
    }
}


inline void CombinedSparseGridDriver::assign_smolyak_arrays()
{ assign_smolyak_arrays(smolMIIter->second, smolCoeffsIter->second); }


inline void CombinedSparseGridDriver::
update_smolyak_coefficients(size_t start_index)
{
  update_smolyak_coefficients(start_index, smolMIIter->second,
			      smolCoeffsIter->second);
}


inline void CombinedSparseGridDriver::
prune_inactive(UShort2DArray& sm_mi, IntArray& sm_coeffs)
{
  size_t i, num_coeffs = sm_coeffs.size(), cntr;
  for (i=0, cntr=0; i<num_coeffs; ++i)
    if (sm_coeffs[i])
      ++cntr;
  if (cntr == num_coeffs) return; // all coeffs already non-zero

  IntArray nonzero_sm_coeffs(cntr);  UShort2DArray nonzero_sm_mi(cntr);
  for (i=0, cntr=0; i<num_coeffs; ++i)
    if (sm_coeffs[i]) {
      nonzero_sm_mi[cntr]     = sm_mi[i];
      nonzero_sm_coeffs[cntr] = sm_coeffs[i];
      ++cntr;
    }
  std::swap(sm_mi, nonzero_sm_mi);  std::swap(sm_coeffs, nonzero_sm_coeffs);
}


inline void CombinedSparseGridDriver::assign_collocation_key()
{ assign_collocation_key(smolMIIter->second, collocKeyIter->second); }


inline void CombinedSparseGridDriver::
compute_unique_points_weights(const UShort2DArray& sm_mi,
			      const IntArray& sm_coeffs,
			      const UShort3DArray& colloc_key,
			      IntArray& unique_index_map,
			      bool update_1d_pts_wts, RealMatrix& var_sets,
			      RealVector& t1_wts, RealMatrix& t2_wts)
{
  RealMatrix a1_pts, a1_t2w;  RealVector a1_t1w, zv, r1v;
  Sizet2DArray colloc_ind;    int num_colloc_pts, num_u1;
  BitArray isu1;              IntArray sind1, uind1, uset1;
  compute_unique_points_weights(sm_mi, sm_coeffs, colloc_key, colloc_ind,
				num_colloc_pts, a1_pts, a1_t1w, a1_t2w, zv,
				r1v, sind1, isu1, uind1, uset1, num_u1,
				unique_index_map, update_1d_pts_wts,
				var_sets, t1_wts, t2_wts);
}


inline void CombinedSparseGridDriver::compute_grid(RealMatrix& var_sets)
{
  compute_grid();
  var_sets = varSetsIter->second; // copy
}


inline const RealMatrix& CombinedSparseGridDriver::variable_sets() const
{ return varSetsIter->second; }


inline const RealMatrix& CombinedSparseGridDriver::
variable_sets(const ActiveKey& key) const
{
  std::map<ActiveKey, RealMatrix>::const_iterator cit
    = variableSets.find(key);
  if (cit == variableSets.end()) {
    PCerr << "Error: key not found in CombinedSparseGridDriver::"
	  << "variable_sets()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const RealVector& CombinedSparseGridDriver::type1_weight_sets() const
{ return t1WtIter->second; }


inline const RealVector& CombinedSparseGridDriver::
type1_weight_sets(const ActiveKey& key) const
{
  std::map<ActiveKey, RealVector>::const_iterator cit
    = type1WeightSets.find(key);
  if (cit == type1WeightSets.end()) {
    PCerr << "Error: key not found in CombinedSparseGridDriver::"
	  << "type1_weight_sets()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const RealMatrix& CombinedSparseGridDriver::type2_weight_sets() const
{ return t2WtIter->second; }


inline const RealMatrix& CombinedSparseGridDriver::
type2_weight_sets(const ActiveKey& key) const
{
  std::map<ActiveKey, RealMatrix>::const_iterator cit
    = type2WeightSets.find(key);
  if (cit == type2WeightSets.end()) {
    PCerr << "Error: key not found in CombinedSparseGridDriver::"
	  << "type2_weight_sets()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const UShort2DArray& CombinedSparseGridDriver::
combined_smolyak_multi_index() const
{ return combinedSmolyakMultiIndex; }


//inline const Sizet2DArray& CombinedSparseGridDriver::
//combined_smolyak_multi_index_map() const
//{ return combinedSmolyakMultiIndexMap; }


inline const UShort3DArray& CombinedSparseGridDriver::
combined_collocation_key() const
{ return combinedCollocKey; }


inline const RealMatrix& CombinedSparseGridDriver::
combined_variable_sets() const
{ return combinedVarSets; }


inline const RealVector& CombinedSparseGridDriver::
combined_type1_weight_sets() const
{ return combinedT1WeightSets; }//return type1_weight_sets(maximalKey);


inline const RealMatrix& CombinedSparseGridDriver::
combined_type2_weight_sets() const
{ return combinedT2WeightSets; }//return type2_weight_sets(maximalKey);


inline void CombinedSparseGridDriver::
basis_collocation_points(int order, int index, double* data)
{
  const RealArray& colloc_pts
    = sgdInstance->polynomialBasis[index].collocation_points(order);
  std::copy(colloc_pts.begin(), colloc_pts.begin()+order, data);
}


inline void CombinedSparseGridDriver::
basis_type1_collocation_weights(int order, int index, double* data)
{
  const RealArray& colloc_wts
    = sgdInstance->polynomialBasis[index].type1_collocation_weights(order);
  std::copy(colloc_wts.begin(), colloc_wts.begin()+order, data);
}


inline void CombinedSparseGridDriver::
basis_type2_collocation_weights(int order, int index, double* data)
{
  const RealArray& colloc_wts
    = sgdInstance->polynomialBasis[index].type2_collocation_weights(order);
  std::copy(colloc_wts.begin(), colloc_wts.begin()+order, data);
}

} // namespace Pecos

#endif
