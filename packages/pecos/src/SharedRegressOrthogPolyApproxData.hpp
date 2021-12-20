/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        SharedRegressPolyApproxData
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        John Jakeman

#ifndef SHARED_REGRESS_ORTHOG_POLY_APPROX_DATA_HPP
#define SHARED_REGRESS_ORTHOG_POLY_APPROX_DATA_HPP

#include "SharedOrthogPolyApproxData.hpp"
#include "LightweightSparseGridDriver.hpp"
#include "LinearSolverPecosSrc.hpp"

namespace Pecos {


/// Container class for various regression configuration options

/** The RegressionConfigOptions class provides a simple container class
    for regression configuration options related to ... */

class RegressionConfigOptions
{
  //
  //- Heading: Friends
  //

  //friend class PolynomialApproximation;

public:

  /// default constructor
  RegressionConfigOptions();
  /// constructor
  RegressionConfigOptions(bool cv, bool cv_noise_only,
			  unsigned short max_cv_order, bool scaling, int seed,
			  const RealVector& noise_tols, Real l2_penalty,
			  bool normalize_cv, unsigned short init_lev,
			  unsigned short growth_fact,
			  unsigned short num_advance);
  /// copy constructor
  RegressionConfigOptions(const RegressionConfigOptions& rc_options);
  /// destructor
  ~RegressionConfigOptions();

//private:

  /// flag for use of automatic cross-validation for parameter
  /// selection in regression approaches
  bool crossValidation;
  /// flag to restrict cross-validation to only estimate the noise
  /// tolerance in order to manage computational cost
  bool crossValidNoiseOnly;
  /// maximum number of order candidates for cross-validation
  unsigned short maxCVOrderCandidates;
  /// scale range of response data to [0,1] to align with regression tols
  bool respScaling;

  /// random seed for data fold definition within cross validation
  int randomSeed;

  /// noise tolerance(s) for compressed sensing algorithms; vector form
  /// used in cross-validation
  RealVector noiseTols;
  /// the L2 penalty parameter for LASSO (elastic net variant)
  Real l2Penalty;

  /// flag indicating whether to scale the CV error estimates by the size
  /// of the candidate basis expansion within adapted basis selection
  bool normalizeCV;

  /// initial sparse grid level that provides the starting point for basis
  /// adaptation in ADAPTED_BASIS_GENERALIZED mode
  unsigned short initSGLevel;
  /// a scalar growth factor for defining the tpMultiIndex contribution
  /// from a particular trial index set for basis adaptation in
  /// ADAPTED_BASIS_GENERALIZED mode
  unsigned short multiIndexGrowthFactor;

  /// number of front expansions per iteration for basis adaptation in
  /// ADAPTED_BASIS_EXPANDING_FRONT mode
  unsigned short numAdvancements;

  /// flag indicating restriction and front expansion using a multi-index
  /// frontier
  /** This option reduces memory requirements somewhat, but allowing
      multi-index gaps in the general case results in reduced mutual
      coherence and better numerical performance. */
  bool advanceByFrontier;
};


inline RegressionConfigOptions::RegressionConfigOptions():
  crossValidation(false), crossValidNoiseOnly(false),
  maxCVOrderCandidates(USHRT_MAX), respScaling(false), randomSeed(0),
  l2Penalty(0.), normalizeCV(false), initSGLevel(0), multiIndexGrowthFactor(2),
  numAdvancements(3), advanceByFrontier(false)
{ }


inline RegressionConfigOptions::
RegressionConfigOptions(bool cv, bool cv_noise_only,
			unsigned short max_cv_order, bool scaling, int seed,
			const RealVector& noise_tols, Real l2_penalty,
			bool normalize_cv, unsigned short init_lev,
			unsigned short growth_fact, unsigned short num_advance):
  crossValidation(cv), crossValidNoiseOnly(cv_noise_only),
  maxCVOrderCandidates(max_cv_order), respScaling(scaling), randomSeed(seed),
  noiseTols(noise_tols), l2Penalty(l2_penalty), normalizeCV(normalize_cv),
  initSGLevel(init_lev), multiIndexGrowthFactor(growth_fact),
  numAdvancements(num_advance), advanceByFrontier(false)
{ }


inline RegressionConfigOptions::
RegressionConfigOptions(const RegressionConfigOptions& rc_options):
  crossValidation(rc_options.crossValidation),
  crossValidNoiseOnly(rc_options.crossValidNoiseOnly),
  maxCVOrderCandidates(rc_options.maxCVOrderCandidates),
  respScaling(rc_options.respScaling),
  randomSeed(rc_options.randomSeed), noiseTols(rc_options.noiseTols),
  l2Penalty(rc_options.l2Penalty), normalizeCV(rc_options.normalizeCV),
  initSGLevel(rc_options.initSGLevel),
  multiIndexGrowthFactor(rc_options.multiIndexGrowthFactor),
  numAdvancements(rc_options.numAdvancements),
  advanceByFrontier(rc_options.advanceByFrontier)
{ }


inline RegressionConfigOptions::~RegressionConfigOptions()
{ }


/// Derived approximation class for multivariate orthogonal polynomial
/// approximation with coefficient estimation via regression.

/** The SharedRegressOrthogPolyApproxData class provides a global
    approximation based on multivariate orthogonal polynomials, where
    the coefficients are computed using regression approaches such as
    least squares (L2) or compressed sensing (L1).  It is used
    primarily for polynomial chaos expansion aproaches to UQ. */

class SharedRegressOrthogPolyApproxData: public SharedOrthogPolyApproxData
{
  //
  //- Heading: Friends
  //

  friend class RegressOrthogPolyApproximation;

public:

  //
  //- Heading: Constructor and destructor
  //

  /// lightweight constructor
  SharedRegressOrthogPolyApproxData(short basis_type,
				    const UShortArray& approx_order,
				    size_t num_vars);
  /// full constructor
  SharedRegressOrthogPolyApproxData(short basis_type,
				    const UShortArray& approx_order,
				    size_t num_vars,
				    const ExpansionConfigOptions& ec_options,
				    const BasisConfigOptions& bc_options,
				    const RegressionConfigOptions& rc_options);
  /// destructor
  ~SharedRegressOrthogPolyApproxData();

  //
  //- Heading: Member functions
  //

  /// helper function for incrementing that is modular on trial set and
  /// multi-index
  void increment_trial_set(const UShortArray& trial_set,
			   UShort2DArray& aggregated_mi);

  /// define the TP index sets that can be removed from a generalized
  /// sparse grid based on the recovered sparse_indices; used for the
  /// restriction operation within basis adaptation
  bool set_restriction(UShort2DArray& aggregated_mi, SizetSet& sparse_indices,
		       SizetSet& save_tp);

  /// clear adapted basis bookkeeping
  void clear_adapted();

  /// set regressConfigOptions
  void configuration_options(const RegressionConfigOptions& rc_options);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void allocate_data();
  void increment_data();
  void decrement_data();

  void pre_push_data();
  //void post_push_data();

  //void pre_finalize_data();
  //void post_finalize_data();

  bool advancement_available();

  //
  //- Heading: Member functions
  //

  /// append multi-indices from append_mi that do not already appear in
  /// combined_mi (consistent ordering assumed); define append_mi_map
  /// and append_mi_map_ref
  void append_leading_multi_index(const UShort2DArray& append_mi,
				  UShort2DArray& combined_mi,
				  SizetSet& append_mi_map,
				  size_t& append_mi_map_ref);
  /// append multi-indices from append_mi that do not already appear in
  /// combined_mi, updating sparse_indices, exp_coeffs, and exp_coeff_grads
  void append_sparse_multi_index(SizetSet& sparse_indices,
				 const UShort2DArray& append_mi,
				 UShort2DArray& combined_mi,
				 RealVector& exp_coeffs,
				 RealMatrix& exp_coeff_grads);

private:

  //
  //- Heading: Member functions
  //

  /// update shared approxOrder based on settings computed for one of the QoI
  void update_approx_order(unsigned short new_order);
  /// convert active approxOrder to active multiIndex according to basis type
  void approx_order_to_multi_index();

  /// pack polynomial contributions to Psi matrix for regression
  void pack_polynomial_data(const RealVector& c_vars, const UShortArray& mi,
			    bool add_val,  double* pack_val,  size_t& pv_cntr,
			    bool add_grad, double* pack_grad, size_t& pg_cntr);
  /// pack response contributions to RHS for regression
  void pack_response_data(const SurrogateDataResp& sdr,
			  bool add_val,  double* pack_val,  size_t& pv_cntr,
			  bool add_grad, double* pack_grad, size_t& pg_cntr);
  /// pack scaled response contributions to RHS for regression
  void pack_response_data(const SurrogateDataResp& sdr,
			  const RealRealPair& factors,
			  bool add_val,  double* pack_val,  size_t& pv_cntr,
			  bool add_grad, double* pack_grad, size_t& pg_cntr);

  //
  //- Heading: Data
  //

  /// container for regression configuration options
  RegressionConfigOptions regressConfigOptions;

  /// Wrapper class that is used to solve regression problems
  CompressedSensingTool CSTool;

  /// lower matrix factor in factorization
  RealMatrix lowerFactor;
  /// upper matrix factor in factorization
  RealMatrix upperFactor;
  /// pivoting history of block-LU factorization
  RealMatrix pivotHistory;
  /// pivoting vector
  IntVector pivotVect;

  /// sparse grid driver for adapting a CS candidate basis using greedy
  /// adaptation via the generalized sparse grid algorithm; it's state
  /// is reset for each response QoI
  LightweightSparseGridDriver lsgDriver;
};


inline SharedRegressOrthogPolyApproxData::
SharedRegressOrthogPolyApproxData(short basis_type,
				  const UShortArray& approx_order,
				  size_t num_vars):
  SharedOrthogPolyApproxData(basis_type, approx_order, num_vars)
{ }


inline SharedRegressOrthogPolyApproxData::
SharedRegressOrthogPolyApproxData(short basis_type,
				  const UShortArray& approx_order,
				  size_t num_vars,
				  const ExpansionConfigOptions& ec_options,
				  const BasisConfigOptions& bc_options,
				  const RegressionConfigOptions& rc_options):
  SharedOrthogPolyApproxData(basis_type, approx_order, num_vars,
			     ec_options, bc_options),
  regressConfigOptions(rc_options)
{ }


inline SharedRegressOrthogPolyApproxData::~SharedRegressOrthogPolyApproxData()
{ }


inline void SharedRegressOrthogPolyApproxData::
configuration_options(const RegressionConfigOptions& rc_options)
{ regressConfigOptions = rc_options; }


inline bool SharedRegressOrthogPolyApproxData::advancement_available()
{
  // if CV saturation is possible, delegate test down to individual QoI in
  // RegressOrthogPolyApproximation::advancement_available()
  return ( regressConfigOptions.crossValidation &&
	  !regressConfigOptions.crossValidNoiseOnly ) ? false : true;
}


inline void SharedRegressOrthogPolyApproxData::clear_adapted()
{
  switch (expConfigOptions.expBasisType) {
  case ADAPTED_BASIS_GENERALIZED:
    // reset tensor-product bookkeeping and save restorable data
    //poppedLevMultiIndex[activeKey].clear();
    poppedMultiIndex[activeKey].clear();
    poppedMultiIndexMap[activeKey].clear();
    poppedMultiIndexMapRef[activeKey].clear();

    tpMultiIndex[activeKey].clear();       // will be rebuilt in allocate_data()
    tpMultiIndexMap[activeKey].clear();    // will be rebuilt in allocate_data()
    tpMultiIndexMapRef[activeKey].clear(); // will be rebuilt in allocate_data()
    break;
  }
}

} // namespace Pecos

#endif
