/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 IntegrationDriver
//- Description: 
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef INTEGRATION_DRIVER_HPP
#define INTEGRATION_DRIVER_HPP

#include "pecos_data_types.hpp"
#include "BasisPolynomial.hpp"

namespace Pecos {

class MultivariateDistribution;
class ExpansionConfigOptions;
class BasisConfigOptions;
class ActiveKey;


/// base class for generating N-dimensional grids for numerical evaluation
/// of expectation integrals over independent standard random variables.

/** This class enables Dakota::NonD{Quadrature,Cubature,SparseGrid}. */

class IntegrationDriver
{
public:

  //
  //- Heading: Constructors, destructor, and operator=
  //

  /// default constructor
  IntegrationDriver();
  /// standard constructor for envelope
  IntegrationDriver(short driver_type);
  /// copy constructor
  IntegrationDriver(const IntegrationDriver& driver);

  /// destructor
  virtual ~IntegrationDriver();

  /// assignment operator
  IntegrationDriver operator=(const IntegrationDriver& driver);

  //
  //- Heading: Virtual functions
  //

  /// initialize all grid settings (distribution params already
  /// set within poly_basis)
  virtual void initialize_grid(const std::vector<BasisPolynomial>& poly_basis);
  /// set int_rules and growth_rules from u_types and mode booleans
  virtual void initialize_grid(const MultivariateDistribution& u_dist,
			       const ExpansionConfigOptions& ec_options,
			       const BasisConfigOptions& bc_options);
  /// update polynomialBasis with data from MultivariateDistribution
  virtual void initialize_grid_parameters(
			       const MultivariateDistribution& mv_dist);

  /// compute variable and weight sets for the grid
  virtual void compute_grid();
  /// compute variable and weight sets for the grid
  virtual void compute_grid(RealMatrix& var_sets);
  /// compute number of collocation points
  virtual int grid_size();

  /// computes and stores data for reinterpolation of covariance on a
  /// higher-order tensor grid
  virtual void reinterpolated_tensor_grid(const UShortArray& lev_index,
					  const SizetList& reinterp_indices);

  /// set key identifying active data set
  virtual void active_key(const ActiveKey& key);
  /// remove all keyed data sets
  virtual void clear_keys();
  /// clear inactive grid settings following their usage/combination
  virtual void clear_inactive();

  /// reset the grid state (undo refinements)
  virtual void reset();

  /// return the key of the maximal stored grid state
  virtual const ActiveKey& maximal_grid();

  /// combine grid data and points/weights
  virtual void combine_grid();
  /// promote combined grid data and points/weights to active data
  virtual void combined_to_active(bool clear_combined);

  /// return active variableSets for Cubature/TensorProduct/CombinedSparseGrid
  /// (HierarchSparseGridDriver::variableSets is 2DArray)
  virtual const RealMatrix& variable_sets() const;
  /// return active type1WeightSets from Cubature/TensorProduct/
  /// CombinedSparseGrid (HierarchSparseGridDriver::type1WeightSets is 2DArray)
  virtual const RealVector& type1_weight_sets() const;
  /// return active type2WeightSets from Cubature/TensorProduct/
  /// CombinedSparseGrid (HierarchSparseGridDriver::type2WeightSets is 2DArray)
  virtual const RealMatrix& type2_weight_sets() const;
  /// return variableSets[key] for TensorProduct/CombinedSparseGrid
  /// (HierarchSparseGridDriver::variableSets is 2DArray)
  virtual const RealMatrix& variable_sets(const ActiveKey& key) const;
  /// return type1WeightSets[key] from TensorProduct/CombinedSparseGrid
  /// (HierarchSparseGridDriver::type1WeightSets is 2DArray)
  virtual const RealVector& type1_weight_sets(const ActiveKey& key) const;
  /// return type2WeightSets[key] from TensorProduct/CombinedSparseGrid
  /// (HierarchSparseGridDriver::type2WeightSets is 2DArray)
  virtual const RealMatrix& type2_weight_sets(const ActiveKey& key) const;

  /// return combinedVarSets for TensorProduct/CombinedSparseGrid
  virtual const RealMatrix& combined_variable_sets() const;
  /// return combinedT1WeightSets for TensorProduct/CombinedSparseGrid
  virtual const RealVector& combined_type1_weight_sets() const;
  /// return combinedT2WeightSets for TensorProduct/CombinedSparseGrid
  virtual const RealMatrix& combined_type2_weight_sets() const;

  //
  //- Heading: Member functions
  //

  /// assign letter or replace existing letter with a new one
  void assign_rep(std::shared_ptr<IntegrationDriver> driver_rep);

  /// compute variable sets for a tensor-product grid
  void compute_tensor_grid(const UShortArray& quad_order,
			   const UShortArray& lev_index,
			   const SizetList& subset_indices, 
			   RealMatrix& variable_sets,
			   UShort2DArray& colloc_key);

  /// assign collocPts1D and type{1,2}CollocWts1D for level/order
  void assign_1d_collocation_points_weights(const UShortArray& quad_order,
					    const UShortArray& lev_index);
  /// assign collocPts1D and type{1,2}CollocWts1D for level/order
  /// for subset variables
  void assign_1d_collocation_points_weights(const UShortArray& quad_order,
					    const UShortArray& lev_index,
					    const SizetList& subset_indices);

  /// return polynomialBasis
  const std::vector<BasisPolynomial>& polynomial_basis() const;
  /// return polynomialBasis
  std::vector<BasisPolynomial>& polynomial_basis();
  /// set polynomialBasis
  void polynomial_basis(const std::vector<BasisPolynomial>& poly_basis);

  /// return basisParamUpdates
  const BitArray& polynomial_basis_parameter_updates() const;

  /// set driverMode
  void mode(short driver_mode);
  /// get driverMode
  short mode() const;

  /// return collocPts1D
  const Real3DArray& collocation_points_1d()  const;
  /// return type1CollocWts1D
  const Real3DArray& type1_collocation_weights_1d() const;
  /// return type2CollocWts1D
  const Real3DArray& type2_collocation_weights_1d() const;

  /// return collocRules
  const ShortArray& collocation_rules() const;

  /// return reinterpLevelIndices[activeReinterpIndex]
  const UShortArray& reinterpolated_level_index() const;
  /// return reinterpQuadOrders[activeReinterpIndex]
  const UShortArray& reinterpolated_quadrature_order() const;
  /// return reinterpVarSets[activeReinterpIndex]
  const RealMatrix& reinterpolated_variable_sets() const;
  /// return reinterpCollocKeys[activeReinterpIndex]
  const UShort2DArray& reinterpolated_collocation_key() const;

  /// return orderGenzKeister
  const UShortArray& genz_keister_order()     const;
  /// return precGenzKeister
  const UShortArray& genz_keister_precision() const;

  /// returns driverRep for access to derived class member functions
  /// that are not mapped to the base level
  std::shared_ptr<IntegrationDriver> driver_rep() const;

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  IntegrationDriver(BaseConstructor);

  //
  //- Heading: Member functions
  //

  /// compute variable and weight sets for a tensor-product grid
  void compute_tensor_grid(const UShortArray& quad_order,
			   const UShortArray& lev_index,
			   RealMatrix& variable_sets,
			   RealVector& t1_weight_sets,
			   RealMatrix& t2_weight_sets,
			   UShort2DArray& colloc_key);

  /// resize collocPts1D and type{1,2}CollocWts1D
  void resize_1d_collocation_points_weights(const UShortArray& lev_index);
  /// update collocPts1D[lev_index][i] and type{1,2}CollocWts1D[lev_index][i]
  /// using points/weights of order quad_order
  void assign_1d_collocation_points_weights(size_t i, unsigned short quad_order,
					    unsigned short lev_index);
  /// clear collocPts1D and type{1,2}CollocWts1D
  void clear_1d_collocation_points_weights();

  //
  //- Heading: Data
  //

  /// number of variables in the tensor-product grid
  size_t numVars;

  /// enumeration value indicating INTEGRATION_MODE or INTERPOLATION_MODE
  short driverMode;

  /// enumeration codes for integration rule options.  Manages internal
  /// mode switches for 1D polynomial types: e.g., GAUSS_LEGENDRE or
  /// GAUSS_PATTERSON for Legendre, CLENSHAW_CURTIS or FEJER2 for
  /// Chebyshev, GAUSS_HERMITE or GENZ_KEISTER for Hermite.
  ShortArray collocRules;

  /// array of one-dimensional orthogonal polynomials used in
  /// computing Gaussian quadrature points and weights
  std::vector<BasisPolynomial> polynomialBasis;
  /// set of flags indicating parameter updates to polynomialBasis
  BitArray basisParamUpdates;

  /// num_levels_per_var x numVars sets of 1D collocation points
  Real3DArray collocPts1D;
  /// num_levels_per_var x numVars sets of 1D type1 collocation weights
  Real3DArray type1CollocWts1D;
  /// num_levels_per_var x numVars sets of 1D type2 collocation weights
  Real3DArray type2CollocWts1D;

  /// flag indicating usage of compute1DType2Weights to define type2WeightSets
  bool computeType2Weights;

  /// bookkeeping for reinterpolation of covariance: stored level indices
  UShort2DArray reinterpLevelIndices;
  /// bookkeeping for reinterpolation of covariance: stored quadrature orders
  UShort2DArray reinterpQuadOrders;
  /// bookkeeping for reinterpolation of covariance: stored variable sets
  RealMatrixArray reinterpVarSets;
  /// bookkeeping for reinterpolation of covariance: stored collocation keys
  UShort3DArray reinterpCollocKeys;

  /// tracks existing reinterpolation grids to avoid unnecessary recomputation
  std::map<UShortArray, size_t> reinterpMap;
  /// bookkeeping for reinterpolation of covariance: active index into arrays
  size_t activeReinterpIndex;

  /// lookup for set of 1-D Genz-Keister quadrature orders
  static UShortArray orderGenzKeister;
  /// lookup for set of 1-D Genz-Keister integrand precisions
  static UShortArray precGenzKeister;

private:

  //
  //- Heading: Member functions
  //

  /// Used only by the standard envelope constructor to initialize
  /// driverRep to the appropriate derived type.
  std::shared_ptr<IntegrationDriver> get_driver(short driver_type);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  std::shared_ptr<IntegrationDriver> driverRep;
};


inline const std::vector<BasisPolynomial>& 
IntegrationDriver::polynomial_basis() const
{ return (driverRep) ? driverRep->polynomialBasis : polynomialBasis; }


inline std::vector<BasisPolynomial>& IntegrationDriver::polynomial_basis()
{ return (driverRep) ? driverRep->polynomialBasis : polynomialBasis; }


inline void IntegrationDriver::
polynomial_basis(const std::vector<BasisPolynomial>& poly_basis)
{
  if (driverRep) driverRep->polynomialBasis = poly_basis;
  else           polynomialBasis = poly_basis;
}


inline const BitArray& IntegrationDriver::
polynomial_basis_parameter_updates() const
{ return (driverRep) ? driverRep->basisParamUpdates : basisParamUpdates; }


inline void IntegrationDriver::mode(short driver_mode)
{
  if (driverRep) driverRep->driverMode = driver_mode;
  else           driverMode = driver_mode;
}


inline short IntegrationDriver::mode() const
{ return (driverRep) ? driverRep->driverMode : driverMode; }


inline const Real3DArray& IntegrationDriver::collocation_points_1d() const
{ return collocPts1D; }


inline const Real3DArray& IntegrationDriver::
type1_collocation_weights_1d() const
{ return type1CollocWts1D; }


inline const Real3DArray& IntegrationDriver::
type2_collocation_weights_1d() const
{ return type2CollocWts1D; }


/** Don't worry about preserving layout as {update,assign}_1d can manage. */
inline void IntegrationDriver::clear_1d_collocation_points_weights()
{ collocPts1D.clear(); type1CollocWts1D.clear(); type2CollocWts1D.clear(); }


inline void IntegrationDriver::
assign_1d_collocation_points_weights(const UShortArray& quad_order,
				     const UShortArray& lev_index)
{
  // resize arrays
  resize_1d_collocation_points_weights(lev_index);
  // assign values
  for (size_t i=0; i<numVars; ++i)
    assign_1d_collocation_points_weights(i, quad_order[i], lev_index[i]);
}


inline void IntegrationDriver::
assign_1d_collocation_points_weights(const UShortArray& quad_order,
				     const UShortArray& lev_index,
				     const SizetList& subset_indices)
{
  // resize arrays (all variables for simplicity)
  resize_1d_collocation_points_weights(lev_index);
  // assign values for subset variables (for memory efficiency)
  SizetList::const_iterator cit;  size_t i;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    i = *cit;
    assign_1d_collocation_points_weights(i, quad_order[i], lev_index[i]);
  }
}


inline const ShortArray& IntegrationDriver::collocation_rules() const
{ return collocRules; }


inline const UShortArray& IntegrationDriver::reinterpolated_level_index() const
{ return reinterpLevelIndices[activeReinterpIndex]; }


inline const UShortArray& IntegrationDriver::
reinterpolated_quadrature_order() const
{ return reinterpQuadOrders[activeReinterpIndex]; }


inline const RealMatrix& IntegrationDriver::reinterpolated_variable_sets() const
{ return reinterpVarSets[activeReinterpIndex]; }


inline const UShort2DArray& IntegrationDriver::
reinterpolated_collocation_key() const
{ return reinterpCollocKeys[activeReinterpIndex]; }


inline const UShortArray& IntegrationDriver::genz_keister_order() const
{ return (driverRep) ? driverRep->orderGenzKeister : orderGenzKeister; }


inline const UShortArray& IntegrationDriver::genz_keister_precision() const
{ return (driverRep) ? driverRep->precGenzKeister : precGenzKeister; }


inline std::shared_ptr<IntegrationDriver> IntegrationDriver::driver_rep() const
{ return driverRep; }

} // namespace Pecos

#endif
