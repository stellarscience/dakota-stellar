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

  TensorProductDriver();                              ///< default constructor
  TensorProductDriver(const UShortArray& quad_order); ///< constructor
  ~TensorProductDriver();                             ///< destructor

  //
  //- Heading: Virtual function redefinitions
  //

  void compute_grid(RealMatrix& variable_sets);
  int  grid_size();
  void reinterpolated_tensor_grid(const UShortArray& lev_index,
				  const SizetList& reinterp_indices);
  void store_grid(size_t index = _NPOS);
  void restore_grid(size_t index = _NPOS);
  void remove_stored_grid(size_t index = _NPOS);
  void clear_stored();

  size_t maximal_grid() const;
  void swap_grid(size_t index);

  //
  //- Heading: Member functions
  //

  /// set quadOrder
  void quadrature_order(const UShortArray& quad_order);
  /// set ith entry in quadOrder
  void quadrature_order(unsigned short order, size_t i);
  /// return quadOrder
  const UShortArray& quadrature_order() const;
  /// return ith entry in quadOrder
  unsigned short quadrature_order(size_t i) const;

  /// determine the lowest quadrature order that provides integrand
  /// exactness at least as great as the specified goal, while
  /// satisfying any nestedness constraints
  void integrand_goal_to_nested_quadrature_order(size_t i,
    unsigned short integrand_goal, unsigned short& nested_quad_order);
  /// determine the lowest quadrature order that provides at least as many
  /// points as the specified goal, while satisfying any nestedness constraints
  void quadrature_goal_to_nested_quadrature_order(size_t i,
    unsigned short quad_goal, unsigned short& nested_quad_order);
  /// update quadOrder and levelIndex from ref_quad_order while
  /// satisfying nested rule constraints
  void nested_quadrature_order(const UShortArray& ref_quad_order);

  /// return type1WeightSets
  const RealVector& type1_weight_sets() const;
  /// return type2WeightSets
  const RealMatrix& type2_weight_sets() const;

  /// return levelIndex
  const UShortArray& level_index() const;
  /// return collocKey
  const UShort2DArray& collocation_key() const;

  /// return storedLevelIndex[index]
  const UShortArray& stored_level_index(size_t index) const;
  /// return storedCollocKey[index]
  const UShort2DArray& stored_collocation_key(size_t index) const;

  /// stand-alone initializer of tensor grid settings (except for
  /// distribution params)
  void initialize_grid(const ShortArray& u_types,
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

  /// update levelIndex from quadOrder
  void update_level_index_from_quadrature_order();
  /// update levelIndex[i] from quadOrder[i]
  void update_level_index_from_quadrature_order(size_t i);

  /// update quadOrder from levelIndex
  void update_quadrature_order_from_level_index();

  //
  //- Heading: Data
  //

  /// the isotropic/anisotropic quadrature order
  UShortArray quadOrder;

  /// quadrature order offset by one for use as 0-based indices
  UShortArray levelIndex;
  /// num points-by-numVars array for identifying the 1-D point
  /// indices for sets of tensor-product collocation points
  UShort2DArray collocKey;

  /// stored driver states: copies of levelIndex
  UShort2DArray storedLevelIndex;
  /// stored driver states: copies of collocKey
  UShort3DArray storedCollocKey;

  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each point in the tensor grid
  RealVector type1WeightSets;
  /// the set of type2 weights (for integration of gradient interpolants)
  /// for each derivative component and for each point in the tensor grid
  RealMatrix type2WeightSets;

  /// stored driver states: copies of type1WeightSets
  RealVectorArray storedType1WeightSets;
  /// stored driver states: copies of type2WeightSets
  RealMatrixArray storedType2WeightSets;
};


inline void TensorProductDriver::update_level_index_from_quadrature_order()
{
  size_t i, len = quadOrder.size();
  if (levelIndex.size() != len) levelIndex.resize(len);
  for (i=0; i<len; ++i)
    levelIndex[i] = quadOrder[i] - 1;
}


inline void TensorProductDriver::
update_level_index_from_quadrature_order(size_t i)
{ levelIndex[i] = quadOrder[i] - 1; }


inline void TensorProductDriver::update_quadrature_order_from_level_index()
{
  size_t i, len = levelIndex.size();
  if (quadOrder.size() != len) quadOrder.resize(len);
  for (i=0; i<len; ++i)
    quadOrder[i] = levelIndex[i] + 1;
}


inline void TensorProductDriver::quadrature_order(const UShortArray& quad_order)
{ quadOrder = quad_order; update_level_index_from_quadrature_order(); }


inline void TensorProductDriver::
quadrature_order(unsigned short order, size_t i)
{ quadOrder[i] = order; update_level_index_from_quadrature_order(i); }


inline const UShortArray& TensorProductDriver::quadrature_order() const
{ return quadOrder; }


inline unsigned short TensorProductDriver::quadrature_order(size_t i) const
{ return quadOrder[i]; }


inline void TensorProductDriver::
nested_quadrature_order(const UShortArray& ref_quad_order)
{
  unsigned short nested_order;
  switch (driverMode) {
  case INTERPOLATION_MODE: // synchronize on number of points
                           // (Lagrange interpolant order = #pts - 1)
    for (size_t i=0; i<numVars; ++i) {
      quadrature_goal_to_nested_quadrature_order(i, ref_quad_order[i],
						 nested_order);
      quadrature_order(nested_order, i); // sets quadOrder and levelIndex
    }
    break;
  default: // {INTEGRATION,DEFAULT}_MODE: synchronize on integrand prec 2m-1
    for (size_t i=0; i<numVars; ++i) {
      integrand_goal_to_nested_quadrature_order(i, 2 * ref_quad_order[i] - 1,
						nested_order);
      quadrature_order(nested_order, i); // sets quadOrder and levelIndex
    }
    break;
  }
}


inline const RealVector& TensorProductDriver::type1_weight_sets() const
{ return type1WeightSets; }


inline const RealMatrix& TensorProductDriver::type2_weight_sets() const
{ return type2WeightSets; }


inline const UShortArray& TensorProductDriver::level_index() const
{ return levelIndex; }


inline const UShort2DArray& TensorProductDriver::collocation_key() const
{ return collocKey; }


inline const UShortArray& TensorProductDriver::
stored_level_index(size_t index) const
{ return storedLevelIndex[index]; }


inline const UShort2DArray& TensorProductDriver::
stored_collocation_key(size_t index) const
{ return storedCollocKey[index]; }


inline int TensorProductDriver::grid_size()
{
  int size = 1;
  for (size_t i=0; i<numVars; ++i)
    size *= quadOrder[i];
  return size;
}


inline TensorProductDriver::TensorProductDriver():
  IntegrationDriver(BaseConstructor())
{ }


inline TensorProductDriver::TensorProductDriver(const UShortArray& quad_order):
  IntegrationDriver(BaseConstructor())
{ quadrature_order(quad_order); }


inline TensorProductDriver::~TensorProductDriver()
{ }

} // namespace Pecos

#endif
