/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 TensorProductDriver
//- Description: 
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef TENSOR_PRODUCT_DRIVER_HPP
#define TENSOR_PRODUCT_DRIVER_HPP

#include "IntegrationDriver.hpp"
#include "ActiveKey.hpp"

namespace Pecos {


/// generates N-dimensional tensor-product quadrature grids for
/// numerical evaluation of expectation integrals over independent
/// standard random variables.

/** This class is used by Dakota::NonDQuadrature, but could also be
    used for general numerical integration of moments. */

class TensorProductDriver: public IntegrationDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  TensorProductDriver();
  /// constructor
  TensorProductDriver(const UShortArray& quad_order);
  /// constructor
  TensorProductDriver(const UShortArray& quad_order, const ActiveKey& key);
  /// destructor
  ~TensorProductDriver();

  //
  //- Heading: Virtual function redefinitions
  //

  void active_key(const ActiveKey& key);
  void clear_keys();
  void clear_inactive();

  void compute_grid();
  void compute_grid(RealMatrix& var_sets);
  int  grid_size();
  void reinterpolated_tensor_grid(const UShortArray& lev_index,
				  const SizetList& reinterp_indices);
  const ActiveKey& maximal_grid();
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

  //
  //- Heading: Member functions
  //

  /// set quadOrder
  void quadrature_order(const UShortArray& quad_order);
  /// set ith entry in quadOrder
  void quadrature_order(unsigned short order, size_t i);
  /// return quadOrder
  const UShortArray& quadrature_order();
  /// return ith entry in quadOrder
  unsigned short quadrature_order(size_t i);

  /// update quadOrder and levelIndex from ref_quad_order while
  /// satisfying nested rule constraints
  void reference_quadrature_order(const UShortArray& ref_quad_order,
				  bool constraints);
  /// return active entry in refQuadOrder
  const UShortArray& reference_quadrature_order();

  /// determine the lowest quadrature order that provides integrand
  /// exactness at least as great as the specified goal, while
  /// satisfying any nestedness constraints
  void integrand_goal_to_nested_quadrature_order(size_t i,
    unsigned short integrand_goal, unsigned short& nested_quad_order);
  /// determine the lowest quadrature order that provides at least as many
  /// points as the specified goal, while satisfying any nestedness constraints
  void quadrature_goal_to_nested_quadrature_order(size_t i,
    unsigned short quad_goal, unsigned short& nested_quad_order);

  /// return levelIndex[activeKey]
  const UShortArray& level_index() const;
  /// return levelIndex[key]
  const UShortArray& level_index(const ActiveKey& key) const;
  /// return collocKey[activeKey]
  const UShort2DArray& collocation_key() const;
  /// return collocKey[key]
  const UShort2DArray& collocation_key(const ActiveKey& key) const;

  /// return combinedLevelIndex
  const UShortArray& combined_level_index() const;
  /// return combinedCollocKey
  const UShort2DArray& combined_collocation_key() const;

  /// stand-alone initializer of tensor grid settings (except for
  /// distribution params)
  void initialize_grid(const MultivariateDistribution& u_dist,
		       const ExpansionConfigOptions& ec_options,
		       const BasisConfigOptions& bc_options);
  /// helper initializer of tensor grid settings (except distribution params)
  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis);

  /// precompute quadrature rules to the maximum current order for each basis
  /// polynomial (efficiency optimization when rules are expensive to compute)
  void precompute_rules();

private:

  //
  //- Heading: Convenience functions
  //

  /// update {levelInd,collocKey}Iter based on activeKey
  void update_active_iterators();

  /// map from ref_quad_order to quadOrder and levelIndex, enforcing any
  /// (nested) constraints
  void enforce_constraints(const UShortArray& ref_quad_order);

  /// update lev_index from q_ord
  void order_to_level(const UShortArray& q_ord,	UShortArray& lev_index);
  /// update levelIndex from quadOrder
  void order_to_level();
  /// update levelIndex[i] from quadOrder[i]
  void order_to_level(size_t i);
  /// update level from order
  void order_to_level(size_t i, unsigned short order, unsigned short& level);

  /// update q_ord from lev_index
  void level_to_order(const UShortArray& lev_index, UShortArray& q_ord);
  /// update quadOrder from levelIndex
  void level_to_order();
  /// update quadOrder[i] from levelIndex[i]
  void level_to_order(size_t i);
  /// update order from level
  void level_to_order(size_t i, unsigned short level, unsigned short& order);

  //
  //- Heading: Data
  //

  /// the isotropic/anisotropic quadrature order
  UShortArray quadOrder;

  /// reference quadrature order: provides a reference for generating a
  /// quadOrder that satisfies applicable constraints (e.g., nestedness, etc.)
  std::map<ActiveKey, UShortArray> refQuadOrder;

  /// quadrature order offset by one for use as 0-based indices
  std::map<ActiveKey, UShortArray> levelIndex;
  /// iterator to active entry within levelIndex
  std::map<ActiveKey, UShortArray>::iterator levelIndIter;

  /// num points-by-numVars array for identifying the 1-D point
  /// indices for sets of tensor-product collocation points
  std::map<ActiveKey, UShort2DArray> collocKey;
  /// iterator to active entry within levelIndex
  std::map<ActiveKey, UShort2DArray>::iterator collocKeyIter;

  /// the set of unique collocation points in the sparse grid
  std::map<ActiveKey, RealMatrix> variableSets;
  /// iterator for active entry within variableSets
  std::map<ActiveKey, RealMatrix>::iterator varSetsIter;
  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each point in the tensor grid
  std::map<ActiveKey, RealVector> type1WeightSets;
  /// iterator for active entry within type1WeightSets
  std::map<ActiveKey, RealVector>::iterator t1WtIter;
  /// the set of type2 weights (for integration of gradient interpolants)
  /// for each derivative component and for each point in the tensor grid
  std::map<ActiveKey, RealMatrix> type2WeightSets;
  /// iterator for active entry within type2WeightSets
  std::map<ActiveKey, RealMatrix>::iterator t2WtIter;

  /// multi-index for maximal grid that is the result of combining a set
  /// of level expansions
  UShortArray combinedLevelIndex;
  /// collocation key for maximal grid that is the result of combining a
  /// set of level expansions (defined from combinedSmolyakMultiIndex)
  UShort2DArray combinedCollocKey;

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

  /// database key indicating the currently active integration configuration.
  /// the key is a multi-index managing multiple modeling dimensions such as
  /// model form, discretization level, etc.
  ActiveKey activeKey;
  // store the key identified in the last call to maximal_grid()
  //ActiveKey maximalKey;
};


inline TensorProductDriver::TensorProductDriver():
  IntegrationDriver(BaseConstructor()), levelIndIter(levelIndex.end())
{ TensorProductDriver::update_active_iterators(); }// default activeKey is empty


inline TensorProductDriver::TensorProductDriver(const UShortArray& quad_order):
  IntegrationDriver(BaseConstructor()), levelIndIter(levelIndex.end())
{
  TensorProductDriver::update_active_iterators(); // default activeKey is empty
  quadrature_order(quad_order);
}


inline TensorProductDriver::
TensorProductDriver(const UShortArray& quad_order, const ActiveKey& key):
  IntegrationDriver(BaseConstructor()), levelIndIter(levelIndex.end()),
  activeKey(key) // shallow copy
{
  TensorProductDriver::update_active_iterators(); // activeKey set init above
  quadrature_order(quad_order);
}


inline TensorProductDriver::~TensorProductDriver()
{ }


inline void TensorProductDriver::active_key(const ActiveKey& key)
{
  if (activeKey != key) {
    activeKey = key; // shallow copy
    update_active_iterators();
  }
}


inline void TensorProductDriver::update_active_iterators()
{
  // Test for change
  if (levelIndIter != levelIndex.end() && levelIndIter->first == activeKey)
    return;

  levelIndIter  = levelIndex.find(activeKey);
  collocKeyIter = collocKey.find(activeKey);
  varSetsIter   = variableSets.find(activeKey);
  t1WtIter      = type1WeightSets.find(activeKey);
  t2WtIter      = type2WeightSets.find(activeKey);

  /* So long as we only create new keys and avoid modifying existing ones,
     this deep copy is not needed.
  ActiveKey active_copy; // share 1 deep copy of current active key
  if (levelIndIter == levelIndex.end()  || collocKeyIter == collocKey.end() ||
      varSetsIter == variableSets.end() || t1WtIter == type1WeightSets.end() ||
      t2WtIter == type2WeightSets.end())
    active_copy = activeKey.copy();
  */

  if (levelIndIter == levelIndex.end()) {
    std::pair<ActiveKey, UShortArray>
      ua_pair(activeKey/*active_copy*/, UShortArray());
    levelIndIter = levelIndex.insert(ua_pair).first;
  }
  level_to_order(); // empty for new levelIndex

  if (collocKeyIter == collocKey.end()) {
    std::pair<ActiveKey, UShort2DArray>
      u2a_pair(activeKey/*active_copy*/, UShort2DArray());
    collocKeyIter = collocKey.insert(u2a_pair).first;
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
}


inline void TensorProductDriver::clear_keys()
{
  activeKey.clear();
  levelIndex.clear();       levelIndIter =      levelIndex.end();
  collocKey.clear();       collocKeyIter =       collocKey.end();
  variableSets.clear();      varSetsIter =    variableSets.end();
  type1WeightSets.clear();      t1WtIter = type1WeightSets.end();
  type2WeightSets.clear();      t2WtIter = type2WeightSets.end();

  // this database is shared among all keys (not currently used by TPQ)
  //clear_1d_collocation_points_weights();
}


inline void TensorProductDriver::
order_to_level(const UShortArray& q_ord, UShortArray& lev_index)
{
  size_t i, len = q_ord.size();
  if (lev_index.size() != len) lev_index.resize(len);
  for (i=0; i<len; ++i)
    lev_index[i] = q_ord[i] - 1;
}


inline void TensorProductDriver::order_to_level()
{ order_to_level(quadOrder, levelIndIter->second); }


inline void TensorProductDriver::order_to_level(size_t i)
{ levelIndIter->second[i] = quadOrder[i] - 1; }


inline void TensorProductDriver::
order_to_level(size_t i, unsigned short order, unsigned short& level)
{ level = order - 1; }


inline void TensorProductDriver::
level_to_order(const UShortArray& lev_index, UShortArray& q_ord)
{
  size_t i, len = lev_index.size();
  if (q_ord.size() != len) q_ord.resize(len);
  for (i=0; i<len; ++i)
    q_ord[i] = lev_index[i] + 1;
}


inline void TensorProductDriver::level_to_order()
{ level_to_order(levelIndIter->second, quadOrder); }


inline void TensorProductDriver::level_to_order(size_t i)
{ quadOrder[i] = levelIndIter->second[i] + 1; }


inline void TensorProductDriver::
level_to_order(size_t i, unsigned short level, unsigned short& order)
{ order = level + 1; }


inline void TensorProductDriver::quadrature_order(const UShortArray& quad_order)
{ quadOrder = quad_order;  order_to_level(); }


inline void TensorProductDriver::
quadrature_order(unsigned short order, size_t i)
{ quadOrder[i] = order;  order_to_level(i); }


inline const UShortArray& TensorProductDriver::quadrature_order()
{ level_to_order();  return quadOrder; }


inline unsigned short TensorProductDriver::quadrature_order(size_t i)
{ level_to_order(i); return quadOrder[i]; }


inline const UShortArray& TensorProductDriver::reference_quadrature_order()
{
  std::map<ActiveKey, UShortArray>::const_iterator cit
    = refQuadOrder.find(activeKey);
  return (cit == refQuadOrder.end()) ? quadrature_order() : cit->second;
}


inline void TensorProductDriver::
reference_quadrature_order(const UShortArray& ref_quad_order, bool constraints)
{
  if (constraints) {
    refQuadOrder[activeKey] = ref_quad_order; // cache reference
    enforce_constraints(ref_quad_order); // map ref to quadOrder
  }
  else // no constraints to enforce -> no need for refQuadOrder
    quadrature_order(ref_quad_order);
}


inline void TensorProductDriver::compute_grid(RealMatrix& var_sets)
{
  compute_grid();
  var_sets = varSetsIter->second; // copy
}


inline const RealMatrix& TensorProductDriver::variable_sets() const
{ return varSetsIter->second; }


inline const RealMatrix& TensorProductDriver::
variable_sets(const ActiveKey& key) const
{
  std::map<ActiveKey, RealMatrix>::const_iterator cit
    = variableSets.find(key);
  if (cit == variableSets.end()) {
    PCerr << "Error: key not found in TensorProductDriver::variable_sets()."
	  << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const RealVector& TensorProductDriver::type1_weight_sets() const
{ return t1WtIter->second; }


inline const RealVector& TensorProductDriver::
type1_weight_sets(const ActiveKey& key) const
{
  std::map<ActiveKey, RealVector>::const_iterator cit
    = type1WeightSets.find(key);
  if (cit == type1WeightSets.end()) {
    PCerr << "Error: key not found in TensorProductDriver::type1_weight_sets()."
	  << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const RealMatrix& TensorProductDriver::type2_weight_sets() const
{ return t2WtIter->second; }


inline const RealMatrix& TensorProductDriver::
type2_weight_sets(const ActiveKey& key) const
{
  std::map<ActiveKey, RealMatrix>::const_iterator cit
    = type2WeightSets.find(key);
  if (cit == type2WeightSets.end()) {
    PCerr << "Error: key not found in TensorProductDriver::type2_weight_sets()."
	  << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const UShortArray& TensorProductDriver::level_index() const
{ return levelIndIter->second; }


inline const UShortArray& TensorProductDriver::
level_index(const ActiveKey& key) const
{
  std::map<ActiveKey, UShortArray>::const_iterator cit = levelIndex.find(key);
  if (cit == levelIndex.end()) {
    PCerr << "Error: key not found in TensorProductDriver::level_index()."
	  << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const UShort2DArray& TensorProductDriver::collocation_key() const
{ return collocKeyIter->second; }


inline const UShort2DArray& TensorProductDriver::
collocation_key(const ActiveKey& key) const
{
  std::map<ActiveKey, UShort2DArray>::const_iterator cit
    = collocKey.find(key);
  if (cit == collocKey.end()) {
    PCerr << "Error: key not found in TensorProductDriver::"
	  << "collocation_key()." << std::endl;
    abort_handler(-1);
  }
  return cit->second;
}


inline const RealMatrix& TensorProductDriver::combined_variable_sets() const
{ return combinedVarSets; }      //variable_sets(maximalKey);


inline const RealVector& TensorProductDriver::combined_type1_weight_sets() const
{ return combinedT1WeightSets; } //type1_weight_sets(maximalKey);


inline const RealMatrix& TensorProductDriver::combined_type2_weight_sets() const
{ return combinedT2WeightSets; } //type2_weight_sets(maximalKey);


inline const UShortArray& TensorProductDriver::combined_level_index() const
{ return combinedLevelIndex; }


inline const UShort2DArray& TensorProductDriver::
combined_collocation_key() const
{ return combinedCollocKey; }


inline int TensorProductDriver::grid_size()
{
  int size = 1;
  for (size_t i=0; i<numVars; ++i)
    size *= quadOrder[i];
  return size;
}

} // namespace Pecos

#endif
