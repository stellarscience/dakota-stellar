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

class AleatoryDistParams;
class ExpansionConfigOptions;
class BasisConfigOptions;


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
  virtual void initialize_grid(const ShortArray& u_types,
			       const ExpansionConfigOptions& ec_options,
			       const BasisConfigOptions& bc_options);
  /// update polynomialBasis with data from AleatoryDistParams
  virtual void initialize_grid_parameters(const ShortArray& u_types,
					  const AleatoryDistParams& adp);

  /// compute scaled variable and weight sets for the TPQ grid
  virtual void compute_grid(RealMatrix& var_sets);
  /// compute number of collocation points
  virtual int grid_size();

  /// computes and stores data for reinterpolation of covariance on a
  /// higher-order tensor grid
  virtual void reinterpolated_tensor_grid(const UShortArray& lev_index,
					  const SizetList& reinterp_indices);

  /// store configuration settings for the current grid before advancing to the
  /// next settings within a prescribed grid sequence (default is push_back)
  virtual void store_grid(size_t index = _NPOS);
  /// restore configuration settings from a previously stored grid
  virtual void restore_grid(size_t index = _NPOS);
  /// remove configuration settings for a stored grid (default is pop_back)
  virtual void remove_stored_grid(size_t index = _NPOS);
  /// clear stored grid settings following their usage/combination
  virtual void clear_stored();

  /// return the index of the maximal stored grid state (_NPOS if the
  /// current unstored grid state)
  virtual size_t maximal_grid() const;
  /// swap settings between the current grid and the stored grid
  /// identified by index
  virtual void swap_grid(size_t index);

  /// return type1WeightSets from Cubature/TensorProduct/CombinedSparseGrid
  /// or concatenate type1WeightSets in HierarchSparseGrid
  virtual const RealVector& type1_weight_sets() const;
  /// return type2WeightSets from Cubature/TensorProduct/CombinedSparseGrid
  /// or concatenate type2WeightSets in HierarchSparseGrid
  virtual const RealMatrix& type2_weight_sets() const;

  //
  //- Heading: Member functions
  //

  /// assign letter or replace existing letter with a new one
  void assign_rep(IntegrationDriver* driver_rep, bool ref_count_incr);

  /// compute variable sets for a tensor-product grid
  void compute_tensor_grid(const UShortArray& quad_order,
			   const UShortArray& lev_index,
			   const SizetList& subset_indices, 
			   RealMatrix& variable_sets,
			   UShort2DArray& colloc_key);

  /// update collocPts1D and type{1,2}CollocWts1D
  void update_1d_collocation_points_weights(const UShortArray& quad_order,
					    const UShortArray& lev_index);

  /// return polynomialBasis
  const std::vector<BasisPolynomial>& polynomial_basis() const;

  /// set driverMode
  void mode(short driver_mode);
  /// get driverMode
  short mode() const;

  // append to end of type1WeightSets
  //void append_type1_weight_sets(const RealVector& t1_wts);
  // append to end of type2WeightSets
  //void append_type2_weight_sets(const RealMatrix& t2_wts);

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
  IntegrationDriver* driver_rep() const;

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

  /// update collocPts1D and type{1,2}CollocWts1D for subset variables
  void update_1d_collocation_points_weights(const UShortArray& quad_order,
					    const UShortArray& lev_index,
					    const SizetList& subset_indices);
  /// update collocPts1D[lev_index][i] and type{1,2}CollocWts1D[lev_index][i]
  void assign_1d_collocation_points_weights(size_t i, unsigned short quad_order,
					    unsigned short lev_index);

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

  // the set of type1 weights (for integration of value interpolants)
  // associated with each point in the {TPQ,SSG,Cub} grid
  //RealVector type1WeightSets;
  // the set of type2 weights (for integration of gradient interpolants)
  // for each derivative component and for each point in the {TPQ,SSG} grid
  //RealMatrix type2WeightSets;

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
  IntegrationDriver* get_driver(short driver_type);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  IntegrationDriver* driverRep;
  /// number of objects sharing driverRep
  int referenceCount;
};


inline const std::vector<BasisPolynomial>& 
IntegrationDriver::polynomial_basis() const
{ return (driverRep) ? driverRep->polynomialBasis : polynomialBasis; }


inline void IntegrationDriver::mode(short driver_mode)
{
  if (driverRep) driverRep->driverMode = driver_mode;
  else           driverMode = driver_mode;
}


inline short IntegrationDriver::mode() const
{ return (driverRep) ? driverRep->driverMode : driverMode; }


/*
inline void IntegrationDriver::
append_type1_weight_sets(const RealVector& t1_wts)
{
  if (driverRep)
    driverRep->append_type1_weight_sets(t1_wts);
  else {
    size_t i, num_curr_t1_wts = type1WeightSets.length(),
      num_new_t1_wts = t1_wts.length(),
      num_total_t1_wts = num_curr_t1_wts + num_new_t1_wts;
    type1WeightSets.resize(num_total_t1_wts);
    for (i=0; i<num_new_t1_wts; ++i)
      type1WeightSets[num_curr_t1_wts+i] = t1_wts[i];
  }
}


inline void IntegrationDriver::
append_type2_weight_sets(const RealMatrix& t2_wts)
{
  if (driverRep)
    driverRep->append_type2_weight_sets(t2_wts);
  else {
    size_t i, j, num_curr_t2_wts = type2WeightSets.numCols(),
      num_new_t2_wts = t2_wts.numCols(),
      num_total_t2_wts = num_curr_t2_wts + num_new_t2_wts;
    type2WeightSets.reshape(numVars, num_total_t2_wts);
    for (i=0; i<num_new_t2_wts; ++i) {
      Real*      curr_t2_i = type2WeightSets[num_curr_t2_wts+i];
      const Real* new_t2_i = t2_wts[i];
      for (j=0; j<numVars; ++j)
	curr_t2_i[j] = new_t2_i[j];
    }
  }
}
*/


inline const Real3DArray& IntegrationDriver::collocation_points_1d() const
{ return collocPts1D; }


inline const Real3DArray& IntegrationDriver::
type1_collocation_weights_1d() const
{ return type1CollocWts1D; }


inline const Real3DArray& IntegrationDriver::
type2_collocation_weights_1d() const
{ return type2CollocWts1D; }


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


inline IntegrationDriver* IntegrationDriver::driver_rep() const
{ return driverRep; }

} // namespace Pecos

#endif
