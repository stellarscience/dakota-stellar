/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchInterpPolyApproximation
//- Description:  Class for polynomial approximation by hierarchical
//-               interpolation
//-               
//- Owner:        Mike Eldred

#ifndef HIERARCH_INTERP_POLY_APPROXIMATION_HPP
#define HIERARCH_INTERP_POLY_APPROXIMATION_HPP

#include "InterpPolyApproximation.hpp"
#include "SharedHierarchInterpPolyApproxData.hpp"
#include "HierarchSparseGridDriver.hpp"

namespace Pecos {


/// Derived approximation class for hierarchical interpolation polynomials
/// (interpolating values and potentially gradients at collocation points).

/** The HierarchInterpPolyApproximation class provides a polynomial
    approximation based on hierarchical interpolation.  Both local and
    global hierarchical basis functions are available.  It is used
    primarily for stochastic collocation approaches to uncertainty
    quantification. */

class HierarchInterpPolyApproximation: public InterpPolyApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// Default constructor
  HierarchInterpPolyApproximation(const SharedBasisApproxData& shared_data);
  /// destructor
  ~HierarchInterpPolyApproximation();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// size expansionType{1,2}Coeffs and expansionType1CoeffGrads
  void allocate_arrays();

  /// initialize productType{1,2}Coeffs using pointers to other QoI
  void initialize_covariance(PolynomialApproximation* poly_approx_2);
  /// clear covariancePointers
  void clear_covariance_pointers();
  /// initialize product interpolant accumulators (prodType{1,2}Coeffs)
  /// from covariancePointers
  void initialize_products();
  /// check if prodType{1,2}Coeffs are defined (non-empty)
  bool product_interpolants();

  void compute_coefficients();

  /// update the coefficients for the expansion of interpolation polynomials:
  /// increment expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void increment_coefficients();
  /// restore the coefficients to their previous state prior to last increment:
  /// decrement expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void pop_coefficients(bool save_data);
  /// restore the coefficients to a previously incremented state as
  /// identified by the current increment to the Smolyak multi index:
  /// push expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void push_coefficients();
  /// finalize the coefficients by applying all previously evaluated increments:
  /// finalize expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void finalize_coefficients();

  /// update combinedExpT{1Coeffs,2Coeffs,1CoeffGrads}
  void combine_coefficients();

  bool update_active_iterators(const UShortArray& key);
  void combined_to_active(bool clear_combined = true);
  void clear_inactive();

  void integrate_response_moments(size_t num_moments, bool combined_stats);
  void integrate_expansion_moments(size_t num_moments, bool combined_stats);

  Real value(const RealVector& x);
  const RealVector& gradient_basis_variables(const RealVector& x);
  const RealVector& gradient_basis_variables(const RealVector& x,
					     const SizetArray& dvv);
  const RealVector& gradient_nonbasis_variables(const RealVector& x);
  const RealSymMatrix& hessian_basis_variables(const RealVector& x);

  Real stored_value(const RealVector& x, const UShortArray& key);
  const RealVector& stored_gradient_basis_variables(const RealVector& x,
						    const UShortArray& key);
  const RealVector& stored_gradient_basis_variables(const RealVector& x,
						    const SizetArray& dvv,
						    const UShortArray& key);
  const RealVector& stored_gradient_nonbasis_variables(const RealVector& x,
						       const UShortArray& key);
  const RealSymMatrix& stored_hessian_basis_variables(const RealVector& x,
						      const UShortArray& key);

  Real mean();
  Real mean(const RealVector& x);
  Real combined_mean();
  Real combined_mean(const RealVector& x);
  
  const RealVector& mean_gradient();
  const RealVector& mean_gradient(const RealVector& x,
				  const SizetArray& dvv);
  const RealVector& variance_gradient();
  const RealVector& variance_gradient(const RealVector& x,
				      const SizetArray& dvv);

  Real covariance(PolynomialApproximation* poly_approx_2);
  Real covariance(const RealVector& x,
		  PolynomialApproximation* poly_approx_2);
  Real combined_covariance(PolynomialApproximation* poly_approx_2);
  Real combined_covariance(const RealVector& x,
			   PolynomialApproximation* poly_approx_2);

  Real beta(bool cdf_flag, Real z_bar);
  Real beta(const RealVector& x, bool cdf_flag, Real z_bar);
  Real combined_beta(bool cdf_flag, Real z_bar);
  Real combined_beta(const RealVector& x, bool cdf_flag, Real z_bar);

  Real delta_mean();
  Real delta_mean(const RealVector& x);
  Real delta_combined_mean();
  Real delta_combined_mean(const RealVector& x);

  Real delta_std_deviation();
  Real delta_std_deviation(const RealVector& x);
  Real delta_combined_std_deviation();
  Real delta_combined_std_deviation(const RealVector& x);

  Real delta_covariance(PolynomialApproximation* poly_approx_2);
  Real delta_covariance(const RealVector& x,
			PolynomialApproximation* poly_approx_2);
  Real delta_combined_covariance(PolynomialApproximation* poly_approx_2);
  Real delta_combined_covariance(const RealVector& x,
				 PolynomialApproximation* poly_approx_2);

  Real delta_beta(bool cdf_flag, Real z_bar);
  Real delta_beta(const RealVector& x, bool cdf_flag, Real z_bar);
  Real delta_combined_beta(bool cdf_flag, Real z_bar);
  Real delta_combined_beta(const RealVector& x, bool cdf_flag, Real z_bar);

  Real delta_z(bool cdf_flag, Real beta_bar);
  Real delta_z(const RealVector& x, bool cdf_flag, Real beta_bar);
  Real delta_combined_z(bool cdf_flag, Real beta_bar);
  Real delta_combined_z(const RealVector& x, bool cdf_flag, Real beta_bar);

  void compute_total_sobol_indices();
  void compute_partial_variance(const BitArray& set_value);

  void clear_computed_bits();

private:

  //
  //- Heading: Convenience functions
  //

  /// emulate new surrogate data following promotion of combined expansion
  /// to active (simplifies stats computation for FINAL_RESULTS)
  void synthetic_surrogate_data(SurrogateData& surr_data);

  /// reset computedRef* to zero
  void clear_reference_computed_bits();
  /// reset computedDelta* to zero
  void clear_delta_computed_bits();
  /// reset all computed bit states to zero
  void clear_current_computed_bits();

  /// compute the value at a point for a particular interpolation level
  Real value(const RealVector& x, const UShort3DArray& sm_mi,
	     const UShort4DArray& key, const RealVector2DArray& t1_coeffs,
	     const RealMatrix2DArray& t2_coeffs, unsigned short level,
	     const UShort2DArray& set_partition = UShort2DArray());
  /// compute the value at a point for a particular interpolation
  /// level and for a specified subset of the variables
  Real value(const RealVector& x, const UShort3DArray& sm_mi,
	     const UShort4DArray& key, const RealVector2DArray& t1_coeffs,
	     const RealMatrix2DArray& t2_coeffs, unsigned short level,
	     const SizetList& subset_indices,
	     const UShort2DArray& set_partition = UShort2DArray());

  /// compute the approximate gradient with respect to the basis variables
  /// at a particular point for a particular interpolation level
  const RealVector& gradient_basis_variables(const RealVector& x,
    const UShort3DArray& sm_mi, const UShort4DArray& key,
    const RealVector2DArray& t1_coeffs, const RealMatrix2DArray& t2_coeffs,
    unsigned short level, const UShort2DArray& set_partition = UShort2DArray());
  /// compute the approximate gradient with respect to the basis variables
  /// at a particular point for a particular interpolation level
  const RealVector& gradient_basis_variables(const RealVector& x,
    const UShort3DArray& sm_mi, const UShort4DArray& key,
    const RealVector2DArray& t1_coeffs, const RealMatrix2DArray& t2_coeffs,
    unsigned short level, const SizetList& subset_indices,
    const UShort2DArray& set_partition = UShort2DArray());
  /// compute the approximate gradient with respect to the basis variables
  /// for a particular point, interpolation level, and DVV
  const RealVector& gradient_basis_variables(const RealVector& x,
    const UShort3DArray& sm_mi, const UShort4DArray& key,
    const RealVector2DArray& t1_coeffs, const RealMatrix2DArray& t2_coeffs,
    const SizetArray& dvv, unsigned short level,
    const UShort2DArray& set_partition = UShort2DArray());

  /// compute the approximate gradient with respect to the nonbasis
  /// variables at a particular point for a particular interpolation level
  const RealVector& gradient_nonbasis_variables(const RealVector& x,
    const UShort3DArray& sm_mi, const UShort4DArray& key,
    const RealMatrix2DArray& t1_coeff_grads, unsigned short level,
    const UShort2DArray& set_partition = UShort2DArray());

  /// compute the approximate Hessian with respect to the basis variables
  /// at a particular point for a particular interpolation level
  const RealSymMatrix& hessian_basis_variables(const RealVector& x,
    const UShort3DArray& sm_mi,	const UShort4DArray& colloc_key,
    const RealVector2DArray& t1_coeffs, unsigned short level,
    const UShort2DArray& set_partition = UShort2DArray());

  /// update bookkeeping when adding a grid increment relative to the
  /// grid reference
  void increment_reference_to_current();
  /// update bookkeeping when removing a grid increment and returning
  /// to the grid reference
  void decrement_current_to_reference();

  /// compute the reference mean, excluding the current grid
  /// increment, using ref_key
  Real reference_mean(const UShort2DArray& ref_key);
  /// compute the reference mean, excluding the current grid
  /// increment, using ref_key
  Real reference_mean(const RealVector& x, const UShort2DArray& ref_key);
  /// compute the reference mean, excluding the current grid
  /// increment, using ref_key
  Real reference_combined_mean(
    const std::map<UShortArray, UShort2DArray>& ref_key_map);
  /// compute the reference mean, excluding the current grid
  /// increment, using ref_key
  Real reference_combined_mean(const RealVector& x,
    const std::map<UShortArray, UShort2DArray>& ref_key_map);
  /// compute the reference variance, excluding the current grid
  /// increment, using ref_key
  Real reference_variance(const UShort2DArray& ref_key);
  /// compute the reference variance, excluding the current grid
  /// increment, using ref_key
  Real reference_variance(const RealVector& x, const UShort2DArray& ref_key);
  /// compute the reference variance, excluding the current grid
  /// increment, using ref_key
  Real reference_combined_variance(
    const std::map<UShortArray, UShort2DArray>& ref_key_map);
  /// compute the reference variance, excluding the current grid
  /// increment, using ref_key
  Real reference_combined_variance(const RealVector& x,
    const std::map<UShortArray, UShort2DArray>& ref_key_map);

  /// compute the covariance increment due to the current grid increment
  /// using the active coefficients and weights
  Real delta_covariance(const RealVector2DArray& r1_t1_coeffs,
    const RealMatrix2DArray& r1_t2_coeffs,const RealVector2DArray& r2_t1_coeffs,
    const RealMatrix2DArray& r2_t2_coeffs, bool same,
    const RealVector2DArray& r1r2_t1_coeffs,
    const RealMatrix2DArray& r1r2_t2_coeffs, const RealVector2DArray& t1_wts,
    const RealMatrix2DArray& t2_wts, const UShort2DArray& ref_key,
    const UShort2DArray& incr_key);
  /// compute the covariance increment due to the current grid increment
  /// using all of the coefficients and weights
  Real delta_covariance(
    const std::map<UShortArray, RealVector2DArray>& r1_t1c_map,
    const std::map<UShortArray, RealMatrix2DArray>& r1_t2c_map,
    const std::map<UShortArray, RealVector2DArray>& r2_t1c_map,
    const std::map<UShortArray, RealMatrix2DArray>& r2_t2c_map, bool same,
    const RealVector2DArray& r1r2_t1c, const RealMatrix2DArray& r1r2_t2c,
    const std::map<UShortArray, RealVector2DArray>& t1_wts_map,
    const std::map<UShortArray, RealMatrix2DArray>& t2_wts_map,
    const UShortArray& active_key,
    const std::map<UShortArray, UShort2DArray>& ref_key_map,
    const std::map<UShortArray, UShort2DArray>& incr_key_map);
  /// compute the covariance increment at x due to the current grid increment
  Real delta_covariance(const RealVector& x,
    const RealVector2DArray& r1_t1_coeffs,const RealMatrix2DArray& r1_t2_coeffs,
    const RealVector2DArray& r2_t1_coeffs,const RealMatrix2DArray& r2_t2_coeffs,
    bool same, const RealVector2DArray& r1r2_t1_coeffs,
    const RealMatrix2DArray& r1r2_t2_coeffs, const UShort3DArray& sm_mi,
    const UShort4DArray& colloc_key, const UShort2DArray& ref_key,
    const UShort2DArray& incr_key);
  /// compute the covariance increment due to the current grid increment
  /// using all of the coefficients and weights
  Real delta_covariance(const RealVector& x,
    const std::map<UShortArray, RealVector2DArray>& r1_t1c_map,
    const std::map<UShortArray, RealMatrix2DArray>& r1_t2c_map,
    const std::map<UShortArray, RealVector2DArray>& r2_t1c_map,
    const std::map<UShortArray, RealMatrix2DArray>& r2_t2c_map, bool same,
    const RealVector2DArray& r1r2_t1c, const RealMatrix2DArray& r1r2_t2c,
    const std::map<UShortArray, UShort3DArray>& sm_mi_map,
    const std::map<UShortArray, UShort4DArray>& colloc_key_map,
    const UShortArray& active_key,
    const std::map<UShortArray, UShort2DArray>& ref_key_map,
    const std::map<UShortArray, UShort2DArray>& incr_key_map);

  /// compute the mean increment due to the current grid increment
  Real delta_mean(const UShort2DArray& incr_key);
  /// compute the mean increment due to the current grid increment
  Real delta_mean(const RealVector& x, const UShort2DArray& incr_key);
  /// compute the mean increment due to the current grid increment
  Real delta_combined_mean(
    const std::map<UShortArray, UShort2DArray>& incr_key_map);
  /// compute the mean increment due to the current grid increment
  Real delta_combined_mean(const RealVector& x,
    const std::map<UShortArray, UShort2DArray>& incr_key_map);

  /// compute the variance increment due to the current grid increment
  Real delta_variance(const UShort2DArray& ref_key,
		      const UShort2DArray& incr_key);
  /// compute the variance increment due to the current grid increment
  Real delta_variance(const RealVector& x, const UShort2DArray& ref_key,
		      const UShort2DArray& incr_key);
  /// compute the variance increment due to the current grid increment
  Real delta_combined_variance(
    const std::map<UShortArray, UShort2DArray>& ref_key_map,
    const std::map<UShortArray, UShort2DArray>& incr_key_map);
  /// compute the variance increment due to the current grid increment
  Real delta_combined_variance(const RealVector& x,
    const std::map<UShortArray, UShort2DArray>& ref_key_map,
    const std::map<UShortArray, UShort2DArray>& incr_key_map);

  /// compute the standard deviation increment due to the current grid increment
  Real delta_std_deviation(const UShort2DArray& ref_key,
			   const UShort2DArray& incr_key);
  /// compute the standard deviation increment due to the current grid increment
  Real delta_std_deviation(const RealVector& x, const UShort2DArray& ref_key,
			   const UShort2DArray& incr_key);
  /// compute the standard deviation increment due to the current grid increment
  Real delta_combined_std_deviation(
    const std::map<UShortArray, UShort2DArray>& ref_key_map,
    const std::map<UShortArray, UShort2DArray>& incr_key_map);
  /// compute the standard deviation increment due to the current grid increment
  Real delta_combined_std_deviation(const RealVector& x,
    const std::map<UShortArray, UShort2DArray>& ref_key_map,
    const std::map<UShortArray, UShort2DArray>& incr_key_map);

  /// compute the reliability index increment due to the current grid increment
  Real delta_beta(bool cdf_flag, Real z_bar, const UShort2DArray& ref_key,
		  const UShort2DArray& incr_key);
  /// compute the reliability index increment due to the current grid increment
  Real delta_beta(const RealVector& x, bool cdf_flag, Real z_bar,
		  const UShort2DArray& ref_key, const UShort2DArray& incr_key);
  /// compute the reliability index increment due to the current grid increment
  Real delta_combined_beta(bool cdf_flag, Real z_bar,
    const std::map<UShortArray, UShort2DArray>& ref_key_map,
    const std::map<UShortArray, UShort2DArray>& incr_key_map);
  /// compute the reliability index increment due to the current grid increment
  Real delta_combined_beta(const RealVector& x, bool cdf_flag, Real z_bar,
    const std::map<UShortArray, UShort2DArray>& ref_key_map,
    const std::map<UShortArray, UShort2DArray>& incr_key_map);

  /// compute the response level increment due to the current grid increment
  Real delta_z(bool cdf_flag, Real beta_bar, const UShort2DArray& ref_key,
	       const UShort2DArray& incr_key);
  /// compute the response level increment due to the current grid increment
  Real delta_z(const RealVector& x, bool cdf_flag, Real beta_bar,
	       const UShort2DArray& ref_key, const UShort2DArray& incr_key);
  /// compute the response level increment due to the current grid increment
  Real delta_combined_z(bool cdf_flag, Real beta_bar,
    const std::map<UShortArray, UShort2DArray>& ref_key_map,
    const std::map<UShortArray, UShort2DArray>& incr_key_map);
  /// compute the response level increment due to the current grid increment
  Real delta_combined_z(const RealVector& x, bool cdf_flag, Real beta_bar,
    const std::map<UShortArray, UShort2DArray>& ref_key_map,
    const std::map<UShortArray, UShort2DArray>& incr_key_map);

  /// shared logic for handling exceptional cases
  Real beta_map(Real mu, Real var, bool cdf_flag, Real z_bar);
  /// shared logic for handling exceptional cases
  Real delta_beta_map(Real mu0, Real delta_mu, Real var0, Real delta_sigma,
		      bool cdf_flag, Real z_bar);

  /// build the active product interpolant, with or without
  /// corresponding SurrogateData
  void product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
    RealVector2DArray& prod_t1c, RealMatrix2DArray& prod_t2c,
    const UShort2DArray& set_partition = UShort2DArray());
  /// form type 1/2 coefficients for interpolation of R_1 R_2
  void product_interpolant(const SDVArray& sdv_array,
    const SDRArray& sdr_array_1, const SDRArray& sdr_array_2,
    const UShort3DArray& sm_mi, const UShort4DArray& colloc_key,
    const Sizet3DArray& colloc_index, RealVector2DArray& r1r2_t1_coeffs,
    RealMatrix2DArray& r1r2_t2_coeffs,
    const UShort2DArray& set_partition = UShort2DArray());
  /// form type 1/2 coefficients for interpolation of R_1 R_2 when
  /// corresponding SurrogateData is not available
  void product_interpolant(const RealMatrix2DArray& var_sets,
    const UShort3DArray& sm_mi, const UShort4DArray& colloc_key,
    const RealVector2DArray& r1_t1_coeffs,
    const RealMatrix2DArray& r1_t2_coeffs,
    const RealVector2DArray& r2_t1_coeffs,
    const RealMatrix2DArray& r2_t2_coeffs, bool same,
    RealVector2DArray& r1r2_t1_coeffs, RealMatrix2DArray& r1r2_t2_coeffs,
    const UShort2DArray& set_partition = UShort2DArray());

  /// form type 1/2 coefficients for interpolation of delta [ R_1 R_2 ] across
  /// model levels/fidelities when corresponding SurrogateData is available
  void product_difference_interpolant(
    HierarchInterpPolyApproximation* hip_approx_2, RealVector2DArray& prod_t1c,
    RealMatrix2DArray& prod_t2c, const UShortArray& lf_key,
    const UShort2DArray& set_partition = UShort2DArray());
  /// form type 1/2 coefficients for interpolation of delta [ R_1 R_2 ] across
  /// model levels/fidelities when corresponding SurrogateData is available
  void product_difference_interpolant(const SurrogateData& surr_data_1,
    const SurrogateData& surr_data_2, const UShort3DArray& sm_mi,
    const UShort4DArray& colloc_key, const Sizet3DArray& colloc_index,
    RealVector2DArray& prod_t1c, RealMatrix2DArray& prod_t2c,
    const UShortArray& lf_key,
    const UShort2DArray& set_partition = UShort2DArray());

  /// build the active central product interpolant, with or without
  /// corresponding SurrogateData
  void central_product_interpolant(
    HierarchInterpPolyApproximation* hip_approx_2, Real mean_1, Real mean_2,
    RealVector2DArray& cov_t1_coeffs, RealMatrix2DArray& cov_t2_coeffs,
    const UShort2DArray& set_partition = UShort2DArray());
  /// build the central product interpolant map, with or without
  /// corresponding SurrogateData
  void central_product_interpolant(
    HierarchInterpPolyApproximation* hip_approx_2, Real mean_1, Real mean_2,
    std::map<UShortArray, RealVector2DArray>& cov_t1c_map,
    std::map<UShortArray, RealMatrix2DArray>& cov_t2c_map,
    const std::map<UShortArray, UShort2DArray>& set_partition_map);
  /// form type 1/2 coefficients for interpolation of (R_1 - mu_1)(R_2 - mu_2)
  /// using the SurrogateData and collocation indices provided
  void central_product_interpolant(const SDVArray& sdv_array,
    const SDRArray& sdr_array_1, const SDRArray& sdr_array_2, Real mean_1,
    Real mean_2, const UShort3DArray& sm_mi, const UShort4DArray& colloc_key,
    const Sizet3DArray&  colloc_index, RealVector2DArray& cov_t1c,
    RealMatrix2DArray& cov_t2c,
    const UShort2DArray& set_partition = UShort2DArray());
  /// form type 1/2 coefficients for interpolation of (R_1 - mu_1)(R_2 - mu_2)
  /// when corresponding SurrogateData is not available
  void central_product_interpolant(const RealMatrix2DArray& var_sets,
    const UShort3DArray& sm_mi, const UShort4DArray& colloc_key,
    const RealVector2DArray& r1_t1_coeffs,
    const RealMatrix2DArray& r1_t2_coeffs,
    const RealVector2DArray& r2_t1_coeffs,
    const RealMatrix2DArray& r2_t2_coeffs, bool same, Real mean_1,Real mean_2,
    RealVector2DArray& cov_t1_coeffs, RealMatrix2DArray& cov_t2_coeffs,
    const UShort2DArray& set_partition = UShort2DArray());

  /// build the active central product gradient interpolant, with or without
  /// corresponding SurrogateData
  void central_product_gradient_interpolant(
    HierarchInterpPolyApproximation* hip_approx_2, Real mean_1, Real mean_2,
    const RealVector& mean1_grad, const RealVector& mean2_grad,
    RealMatrix2DArray& cov_t1c_grads,
    const UShort2DArray& set_partition = UShort2DArray());
  /// form type1 coefficient gradients for interpolation of 
  /// d/ds [(R_1 - mu_1)(R_2 - mu_2)] when corresponding SurrogateData
  /// is available
  void central_product_gradient_interpolant(const SDVArray& sdv_array,
    const SDRArray& sdr_array_1, const SDRArray& sdr_array_2, Real mean_r1,
    Real mean_r2, const RealVector& mean1_grad, const RealVector& mean2_grad,
    const UShort3DArray& sm_mi, const UShort4DArray& colloc_key,
    const Sizet3DArray& colloc_index, RealMatrix2DArray& cov_t1c_grads,
    const UShort2DArray& set_partition = UShort2DArray());
  /// form type1 coefficient gradients for interpolation of 
  /// d/ds [(R_1 - mu_1)(R_2 - mu_2)] when corresponding surrData
  /// is not available
  void central_product_gradient_interpolant(const RealMatrix2DArray& var_sets,
    const UShort3DArray& sm_mi, const UShort4DArray& colloc_key,
    const RealVector2DArray& r1_t1_coeffs,
    const RealMatrix2DArray& r1_t2_coeffs,
    const RealMatrix2DArray& r1_t1_coeff_grads,
    const RealVector2DArray& r2_t1_coeffs,
    const RealMatrix2DArray& r2_t2_coeffs,
    const RealMatrix2DArray& r2_t1_coeff_grads, bool same, Real mean_1,
    Real mean_2, const RealVector& mean1_grad, const RealVector& mean2_grad, 
    RealMatrix2DArray& cov_t1c_grads,
    const UShort2DArray& set_partition = UShort2DArray());

  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using active weights from the HierarchSparseGridDriver
  Real expectation(const RealVector2DArray& t1_coeffs,
		   const RealMatrix2DArray& t2_coeffs,
		   const UShort2DArray& set_partition = UShort2DArray());
  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using t{1,2}_wts
  Real expectation(const RealVector2DArray& t1_coeffs,
		   const RealMatrix2DArray& t2_coeffs,
		   const RealVector2DArray& t1_wts,
		   const RealMatrix2DArray& t2_wts,
		   const UShort2DArray& set_partition = UShort2DArray());
  // compute the expected value of the interpolant given by t{1,2}_coeffs
  //Real expectation(const RealVector2DArray& t1_coeffs,
  //		   const RealMatrix2DArray& t2_coeffs,
  //		   const UShort3DArray& pt_partition);
  /// compute the expected value of the interpolant given by maps of t{1,2}
  /// coefficients and weights
  Real expectation(const std::map<UShortArray, RealVector2DArray>& t1c_map,
		   const std::map<UShortArray, RealMatrix2DArray>& t2c_map,
		   const std::map<UShortArray, RealVector2DArray>& t1_wts_map,
		   const std::map<UShortArray, RealMatrix2DArray>& t2_wts_map,
		   const std::map<UShortArray, UShort2DArray>&
		     set_partition_map);
  /// compute the expected value of the interpolant given by maps of t{1,2}
  /// coefficients and weights
  Real expectation(
    const std::map<UShortArray, std::map<PolynomialApproximation*,
      RealVector2DArray> >& prod_t1c_map,
    const std::map<UShortArray, std::map<PolynomialApproximation*,
      RealMatrix2DArray> >& prod_t2c_map,
    PolynomialApproximation* poly_approx_2,
    const std::map<UShortArray, RealVector2DArray>& t1_wts_map,
    const std::map<UShortArray, RealMatrix2DArray>& t2_wts_map,
    const std::map<UShortArray, UShort2DArray>& set_partition_map);

  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using active weights from the HierarchSparseGridDriver
  Real expectation(const RealVector& x, const RealVector2DArray& t1_coeffs,
		   const RealMatrix2DArray& t2_coeffs,
		   const UShort2DArray& set_partition = UShort2DArray());
  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using partial weights determined from sm_mi and colloc_key
  Real expectation(const RealVector& x, const RealVector2DArray& t1_coeffs,
		   const RealMatrix2DArray& t2_coeffs,
		   const UShort3DArray& sm_mi, const UShort4DArray& colloc_key,
		   const UShort2DArray& set_partition = UShort2DArray());
  /// compute the expected value of the interpolant given by maps of t{1,2}
  /// coefficients using partial weights determined from Smolyak multi-index
  /// and collocation key maps
  Real expectation(const RealVector& x,
		   const std::map<UShortArray, RealVector2DArray>& t1c_map,
		   const std::map<UShortArray, RealMatrix2DArray>& t2c_map,
		   const std::map<UShortArray, UShort3DArray>& sm_mi_map,
		   const std::map<UShortArray, UShort4DArray>& colloc_key_map,
		   const std::map<UShortArray, UShort2DArray>&
		     set_partition_map);
  /// compute the expected value of the interpolant given by maps of t{1,2}
  /// coefficients using partial weights determined from Smolyak multi-index
  /// and collocation key maps
  Real expectation(const RealVector& x,
    const std::map<UShortArray, std::map<PolynomialApproximation*,
      RealVector2DArray> >& prod_t1c_map,
    const std::map<UShortArray, std::map<PolynomialApproximation*,
      RealMatrix2DArray> >& prod_t2c_map,
    PolynomialApproximation* poly_approx_2,
    const std::map<UShortArray, UShort3DArray>& sm_mi_map,
    const std::map<UShortArray, UShort4DArray>& colloc_key_map,
    const std::map<UShortArray, UShort2DArray>& set_partition_map);

  /// compute the expected value of the gradient interpolant given by
  /// t1_coeff_grads using weights from the HierarchSparseGridDriver
  const RealVector& expectation_gradient(
    const RealMatrix2DArray& t1_coeff_grads);
  /// compute the expected value of the gradient interpolant given by
  /// t1_coeff_grads using t1_wts
  const RealVector& expectation_gradient(
    const RealMatrix2DArray& t1_coeff_grads, const RealVector2DArray& t1_wts);

  /// compute the expectation of t1_coeff_grads for index t1cg_index
  Real expectation_gradient(const RealVector& x,
			    const RealMatrix2DArray& t1_coeff_grads,
			    size_t t1cg_index);
  /// compute the expectation of t1_coeff_grads for index t1cg_index
  Real expectation_gradient(const RealVector& x,
			    const RealMatrix2DArray& t1_coeff_grads,
			    const UShort3DArray& sm_mi,
			    const UShort4DArray& colloc_key, size_t t1cg_index);

  /// compute the expectation of the gradient of {t1,t2}_coeffs for
  /// index deriv_index
  Real expectation_gradient(const RealVector& x,
			    const RealVector2DArray& t1_coeffs,
			    const RealMatrix2DArray& t2_coeffs,
			    size_t deriv_index);
  /// compute the expectation of the gradient of {t1,t2}_coeffs for
  /// index deriv_index
  Real expectation_gradient(const RealVector& x,
			    const RealVector2DArray& t1_coeffs,
			    const RealMatrix2DArray& t2_coeffs,
			    const UShort3DArray& sm_mi,
			    const UShort4DArray& colloc_key,
			    size_t deriv_index);

  /// increment expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  /// for a single index_set
  void increment_coefficients(const UShortArray& index_set);

  /// increment coefficients of product interpolants
  void increment_products(const UShort2DArray& set_partition = UShort2DArray());

  /// move all coefficient sets from poppedExp{T1Coeffs,T2Coeffs,T1CoeffGrads}
  /// to expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void promote_all_popped_coefficients();

  /// helper function for common case where coefficients and modSurrData
  /// are synchronized
  void integrate_response_moments(size_t num_moments,
				  const UShort3DArray& sm_mi,
				  const UShort4DArray& colloc_key,
				  const Sizet3DArray&  colloc_index,
				  const SDVArray& sdv_array,
				  const SDRArray& sdr_array);
  /// helper function for expansion combination case where modSurrData
  /// does not span the aggregate set of combined terms
  void integrate_response_moments(size_t num_moments,
				  const RealMatrix2DArray& var_sets,
				  const UShort3DArray&     sm_mi,
				  const UShort4DArray&     colloc_key,
				  const RealVector2DArray& t1_coeffs,
				  const RealMatrix2DArray& t2_coeffs,
				  const RealVector2DArray& t1_wts,
				  const RealMatrix2DArray& t2_wts);

  /// compute member expansion for Sobol' index integration
  void member_coefficients_weights(const BitArray&    member_bits,
    RealVector2DArray& member_t1_coeffs, RealVector2DArray& member_t1_wts,
    RealMatrix2DArray& member_t2_coeffs, RealMatrix2DArray& member_t2_wts,
    UShort4DArray& member_colloc_key,    Sizet3DArray& member_colloc_index);
  /// form hierarchical interpolant of (h-mean)^2 from member-variable
  /// expansion of h
  void central_product_member_coefficients(const BitArray& m_bits,
    const RealVector2DArray& m_t1_coeffs, const RealMatrix2DArray& m_t2_coeffs,
    const UShort4DArray&    m_colloc_key, const Sizet3DArray&   m_colloc_index,
    Real mean, RealVector2DArray& cprod_m_t1_coeffs,
    RealMatrix2DArray& cprod_m_t2_coeffs);

  //
  //- Heading: Data
  //

  /// performance setting: reuse accumulated productType{1,2}Coeffs for fast
  /// computation of covariance metrics, despite loss of precision from
  /// subtractive cancellation (E[R_i R_j] - Mean_i Mean_j)
  bool speedOverPrecision;

  /// bookkeeping to track computation of reference mean to avoid
  /// unnecessary recomputation
  std::map<UShortArray, short> computedRefMean;
  /// iterator to active entry in computedRefMean
  std::map<UShortArray, short>::iterator compRefMeanIter;
  /// bookkeeping to track computation of mean increment to avoid
  /// unnecessary recomputation
  std::map<UShortArray, short> computedDeltaMean;
  /// iterator to active entry in computedDeltaMean
  std::map<UShortArray, short>::iterator compDeltaMeanIter;
  /// bookkeeping to track computation of reference variance to avoid
  /// unnecessary recomputation
  std::map<UShortArray, short> computedRefVariance;
  /// iterator to active entry in computedRefVariance
  std::map<UShortArray, short>::iterator compRefVarIter;
  /// bookkeeping to track computation of variance increment to avoid
  /// unnecessary recomputation
  std::map<UShortArray, short> computedDeltaVariance;
  /// iterator to active entry in computedDeltaVariance
  std::map<UShortArray, short>::iterator compDeltaVarIter;

  /// track previous evaluation point for all_variables reference mean
  /// to avoid unnecessary recomputation
  std::map<UShortArray, RealVector> xPrevRefMean;
  /// track previous evaluation point for all_variables mean increment
  /// to avoid unnecessary recomputation
  std::map<UShortArray, RealVector> xPrevDeltaMean;
  /// track previous evaluation point for all_variables reference
  /// variance to avoid unnecessary recomputation
  std::map<UShortArray, RealVector> xPrevRefVar;
  /// track previous evaluation point for all_variables variance
  /// increment to avoid unnecessary recomputation
  std::map<UShortArray, RealVector> xPrevDeltaVar;

  /// storage for reference mean and variance
  std::map<UShortArray, RealVector> referenceMoments;
  /// iterator to active entry in referenceMoments
  std::map<UShortArray, RealVector>::iterator refMomentsIter;
  /// storage for mean and variance increments
  std::map<UShortArray, RealVector> deltaMoments;
  /// iterator to active entry in deltaMoments
  std::map<UShortArray, RealVector>::iterator deltaMomentsIter;
  /// storage for reference moment gradients (mean, variance)
  std::map<UShortArray, RealVectorArray> momentRefGradients;

  /// the type1 coefficients of the expansion for interpolating values
  std::map<UShortArray, RealVector2DArray> expansionType1Coeffs;
  /// iterator pointing to active node in expansionType1Coeffs
  std::map<UShortArray, RealVector2DArray>::iterator expT1CoeffsIter;

  /// the type2 coefficients of the expansion for interpolating gradients
  std::map<UShortArray, RealMatrix2DArray> expansionType2Coeffs;
  /// iterator pointing to active node in expansionType2Coeffs
  std::map<UShortArray, RealMatrix2DArray>::iterator expT2CoeffsIter;

  /// the gradients of the type1 expansion coefficients
  /** may be interpreted as either the gradients of the expansion coefficients
      or the coefficients of expansions for the response gradients.  This
      array is used when sensitivities of moments are needed with respect to
      variables that do not appear in the expansion (e.g., with respect to
      design variables for an expansion only over the random variables). */
  std::map<UShortArray, RealMatrix2DArray> expansionType1CoeffGrads;
  /// iterator pointing to active node in expansionType1CoeffGrads
  std::map<UShortArray, RealMatrix2DArray>::iterator expT1CoeffGradsIter;

  /// type 1 expansion coefficients popped during decrement for later
  /// restoration to expansionType1Coeffs
  std::map<UShortArray, RealVectorDequeArray> poppedExpT1Coeffs;
  /// type 2 expansion coefficients popped during decrement for later
  /// restoration to expansionType2Coeffs
  std::map<UShortArray, RealMatrixDequeArray> poppedExpT2Coeffs;
  /// type 1 expansion coefficient gradients popped during decrement
  /// for later restoration to expansionType1CoeffGrads
  std::map<UShortArray, RealMatrixDequeArray> poppedExpT1CoeffGrads;

  /// roll up of expansion type 1 coefficients across all keys
  RealVector2DArray combinedExpT1Coeffs;
  /// roll up of expansion type 2 coefficient gradients across all keys
  RealMatrix2DArray combinedExpT2Coeffs;
  /// roll up of expansion type 1 coefficient gradients across all keys
  RealMatrix2DArray combinedExpT1CoeffGrads;

  /// the type1 coefficients of the expansion for interpolating values
  std::map<UShortArray, std::map<PolynomialApproximation*, RealVector2DArray> >
    productType1Coeffs;
  /// iterator pointing to active node in productType1Coeffs
  std::map<UShortArray, std::map<PolynomialApproximation*,
    RealVector2DArray> >::iterator prodT1CoeffsIter;
  /// the type2 coefficients of the expansion for interpolating values
  std::map<UShortArray, std::map<PolynomialApproximation*, RealMatrix2DArray> >
    productType2Coeffs;
  /// iterator pointing to active node in productType2Coeffs
  std::map<UShortArray, std::map<PolynomialApproximation*,
    RealMatrix2DArray> >::iterator prodT2CoeffsIter;

  /// the type1 coefficients of the expansion for interpolating values
  std::map<UShortArray, std::map<PolynomialApproximation*,
    RealVectorDequeArray> > poppedProdType1Coeffs;
  /// the type2 coefficients of the expansion for interpolating values
  std::map<UShortArray, std::map<PolynomialApproximation*,
    RealMatrixDequeArray> > poppedProdType2Coeffs;

  /// array of pointers to the set of QoI used in covariance calculations
  std::set<PolynomialApproximation*> covariancePointers;
};


inline HierarchInterpPolyApproximation::
HierarchInterpPolyApproximation(const SharedBasisApproxData& shared_data):
  InterpPolyApproximation(shared_data),
  expT1CoeffsIter(expansionType1Coeffs.end()),
  prodT1CoeffsIter(productType1Coeffs.end()),
  prodT2CoeffsIter(productType2Coeffs.end()), speedOverPrecision(false)
{ }


inline HierarchInterpPolyApproximation::~HierarchInterpPolyApproximation()
{ }


inline bool HierarchInterpPolyApproximation::
update_active_iterators(const UShortArray& key)
{
  // Test for change
  if (expT1CoeffsIter != expansionType1Coeffs.end() &&
      expT1CoeffsIter->first == key)
    return false;
  
  expT1CoeffsIter = expansionType1Coeffs.find(key);
  if (expT1CoeffsIter == expansionType1Coeffs.end()) {
    std::pair<UShortArray, RealVector2DArray> rv_pair(key, RealVector2DArray());
    expT1CoeffsIter = expansionType1Coeffs.insert(rv_pair).first;
  }
  expT2CoeffsIter = expansionType2Coeffs.find(key);
  if (expT2CoeffsIter == expansionType2Coeffs.end()) {
    std::pair<UShortArray, RealMatrix2DArray> rm_pair(key, RealMatrix2DArray());
    expT2CoeffsIter = expansionType2Coeffs.insert(rm_pair).first;
  }
  expT1CoeffGradsIter = expansionType1CoeffGrads.find(key);
  if (expT1CoeffGradsIter == expansionType1CoeffGrads.end()) {
    std::pair<UShortArray, RealMatrix2DArray> rm_pair(key, RealMatrix2DArray());
    expT1CoeffGradsIter = expansionType1CoeffGrads.insert(rm_pair).first;
  }

  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  if (data_rep->expConfigOptions.refineControl) {
    prodT1CoeffsIter = productType1Coeffs.find(key);
    if (prodT1CoeffsIter == productType1Coeffs.end()) {
      std::map<PolynomialApproximation*, RealVector2DArray> empty_rvm;
      std::pair<UShortArray, std::map<PolynomialApproximation*,
	RealVector2DArray> > prv_pair(key, empty_rvm);
      prodT1CoeffsIter = productType1Coeffs.insert(prv_pair).first;
    }
    prodT2CoeffsIter = productType2Coeffs.find(key);
    if (prodT2CoeffsIter == productType2Coeffs.end()) {
      std::map<PolynomialApproximation*, RealMatrix2DArray> empty_rmm;
      std::pair<UShortArray, std::map<PolynomialApproximation*,
        RealMatrix2DArray> > prm_pair(key, empty_rmm);
      prodT2CoeffsIter = productType2Coeffs.insert(prm_pair).first;
    }
  }

  compRefMeanIter = computedRefMean.find(key);
  if (compRefMeanIter == computedRefMean.end()) {
    std::pair<UShortArray, short> us_pair(key, 0);
    compRefMeanIter = computedRefMean.insert(us_pair).first;
  }
  compDeltaMeanIter = computedDeltaMean.find(key);
  if (compDeltaMeanIter == computedDeltaMean.end()) {
    std::pair<UShortArray, short> us_pair(key, 0);
    compDeltaMeanIter = computedDeltaMean.insert(us_pair).first;
  }
  compRefVarIter = computedRefVariance.find(key);
  if (compRefVarIter == computedRefVariance.end()) {
    std::pair<UShortArray, short> us_pair(key, 0);
    compRefVarIter = computedRefVariance.insert(us_pair).first;
  }
  compDeltaVarIter = computedDeltaVariance.find(key);
  if (compDeltaVarIter == computedDeltaVariance.end()) {
    std::pair<UShortArray, short> us_pair(key, 0);
    compDeltaVarIter = computedDeltaVariance.insert(us_pair).first;
  }

  refMomentsIter = referenceMoments.find(key);
  if (refMomentsIter == referenceMoments.end()) {
    std::pair<UShortArray, RealVector> rv_pair(key, RealVector());
    refMomentsIter = referenceMoments.insert(rv_pair).first;
  }
  deltaMomentsIter = deltaMoments.find(key);
  if (deltaMomentsIter == deltaMoments.end()) {
    std::pair<UShortArray, RealVector> rv_pair(key, RealVector());
    deltaMomentsIter = deltaMoments.insert(rv_pair).first;
  }

  InterpPolyApproximation::update_active_iterators(key);
  return true;
}


inline bool HierarchInterpPolyApproximation::product_interpolants()
{
  return ( ( prodT1CoeffsIter != productType1Coeffs.end() &&
	    !prodT1CoeffsIter->second.empty() ) ||
	   ( prodT2CoeffsIter != productType2Coeffs.end() &&
	    !prodT2CoeffsIter->second.empty() ) );
}


inline void HierarchInterpPolyApproximation::
initialize_covariance(PolynomialApproximation* poly_approx_2)
{
  // don't know all active (multilevel) keys at initialization time and don't
  // know all pointers within a single approximation's compute_coefficients()
  // --> cache pointers here for use in initializing productType{1,2}Coeffs
  // within compute_coefficients() (which calls initialize_products())
  covariancePointers.insert(poly_approx_2);
}


inline void HierarchInterpPolyApproximation::clear_covariance_pointers()
{ covariancePointers.clear(); }


inline void HierarchInterpPolyApproximation::clear_reference_computed_bits()
{ compRefMeanIter->second = compRefVarIter->second = 0; }// clear reference bits


inline void HierarchInterpPolyApproximation::clear_delta_computed_bits()
{ compDeltaMeanIter->second = compDeltaVarIter->second = 0; }// clear delta bits


inline void HierarchInterpPolyApproximation::clear_current_computed_bits()
{ compMeanIter->second = compVarIter->second = 0; }


inline void HierarchInterpPolyApproximation::clear_computed_bits()
{
  clear_reference_computed_bits();
  clear_delta_computed_bits();
  clear_current_computed_bits();
}


inline void HierarchInterpPolyApproximation::increment_reference_to_current()
{
  // update reference bits
  short computed_mean = compMeanIter->second,
        computed_var  =  compVarIter->second;
  compRefMeanIter->second = computed_mean;
  compRefVarIter->second  = computed_var;

  // update reference data
  if ( (computed_mean & 1) || (computed_var & 1) )
    refMomentsIter->second = numMomentsIter->second;
  if ( (computed_mean & 2) || (computed_var & 2) ) {
    SharedHierarchInterpPolyApproxData* data_rep
      = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
    momentRefGradients[data_rep->activeKey] = momentGradsIter->second;
  }

  clear_current_computed_bits(); clear_delta_computed_bits();
}


inline void HierarchInterpPolyApproximation::decrement_current_to_reference()
{
  // update current bits
  short comp_ref_mean = compRefMeanIter->second,
        comp_ref_var  = compRefVarIter->second;
  compMeanIter->second = comp_ref_mean;
  compVarIter->second  = comp_ref_var;

  // update current data
  if ( (comp_ref_mean & 1) || (comp_ref_var & 1) )
    numMomentsIter->second = refMomentsIter->second;
  if ( (comp_ref_mean & 2) || (comp_ref_var & 2) ) {
    SharedHierarchInterpPolyApproxData* data_rep
      = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
    momentGradsIter->second = momentRefGradients[data_rep->activeKey];
  }

  clear_delta_computed_bits(); // clear delta bits, but retain reference
}


inline void HierarchInterpPolyApproximation::
integrate_response_moments(size_t num_moments, bool combined_stats)
{
  // standard variables mode supports four moments using the collocation rules
  // as integration rules
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();

  // Support combined_stats for completeness
  // > use of combined_to_active() prior to full_stats computation makes
  //   this moot / unused for HIPA
  if (combined_stats)
    integrate_response_moments(num_moments,
      hsg_driver->combined_hierarchical_variable_sets(),
      hsg_driver->combined_smolyak_multi_index(),
      hsg_driver->combined_collocation_key(),
      combinedExpT1Coeffs, combinedExpT2Coeffs,
      hsg_driver->combined_type1_hierarchical_weight_sets(),
      hsg_driver->combined_type2_hierarchical_weight_sets());
  else { // compute response moments for active expansion
    const UShort3DArray&      sm_mi = hsg_driver->smolyak_multi_index();
    const UShort4DArray& colloc_key = hsg_driver->collocation_key();
    const Sizet3DArray&  colloc_ind = hsg_driver->collocation_indices();
    // check for colloc indices that were invalidated by expansion combination
    if (hsg_driver->track_collocation_indices() && colloc_ind.empty())
      integrate_response_moments(num_moments,
	hsg_driver->hierarchical_variable_sets(), sm_mi, colloc_key,
	expT1CoeffsIter->second, expT2CoeffsIter->second,
        hsg_driver->type1_hierarchical_weight_sets(),
        hsg_driver->type2_hierarchical_weight_sets());
    else // colloc_index is valid -> can pull from modSurrData vars/responses
      integrate_response_moments(num_moments, sm_mi, colloc_key, colloc_ind,
        modSurrData.variables_data(), modSurrData.response_data());
  }
}


inline void HierarchInterpPolyApproximation::
integrate_expansion_moments(size_t num_moments, bool combined_stats)
{
  // for now: nested interpolation is exact
  expMomentsIter->second = numMomentsIter->second;

  // a couple different ways to go with this in the future:
  // (1) evaluate hierarchical value(lev) - value(lev-1) with HSGDriver wts
  // (2) evaluate value() with CSGDriver wts
  //  > promote Nodal implementation of this function to base class
  //  > redefine HierarchSparseGridDriver::type1_weight_sets() to generate
  //    from 1D weights array in CSG-style approach (not simple concatenation)
}


inline void HierarchInterpPolyApproximation::
product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
		    RealVector2DArray& prod_t1c, RealMatrix2DArray& prod_t2c,
		    const UShort2DArray& set_partition)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  if (hsg_driver->track_collocation_indices() &&
      hsg_driver->collocation_indices().empty()) { // invalidated by combination
    bool same = (this == hip_approx_2);    
    product_interpolant(hsg_driver->hierarchical_variable_sets(),
      hsg_driver->smolyak_multi_index(), hsg_driver->collocation_key(),
      expT1CoeffsIter->second, expT2CoeffsIter->second,
      hip_approx_2->expT1CoeffsIter->second,
      hip_approx_2->expT2CoeffsIter->second, same,
      prod_t1c, prod_t2c, set_partition);
  }
  // use SurrogateData instance + colloc_indices for forming product interp
  else if (true) // current uses are modSurrData; possible future conditional
    product_interpolant(modSurrData.variables_data(),
      modSurrData.response_data(), hip_approx_2->modSurrData.response_data(),
      hsg_driver->smolyak_multi_index(), hsg_driver->collocation_key(),
      hsg_driver->collocation_indices(), prod_t1c, prod_t2c, set_partition);
  else
    product_interpolant(surrData.variables_data(),
      surrData.response_data(), hip_approx_2->surrData.response_data(),
      hsg_driver->smolyak_multi_index(), hsg_driver->collocation_key(),
      hsg_driver->collocation_indices(), prod_t1c, prod_t2c, set_partition);
}


inline void HierarchInterpPolyApproximation::
product_difference_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
  RealVector2DArray& prod_t1c, RealMatrix2DArray& prod_t2c,
  const UShortArray& lf_key, const UShort2DArray& set_partition)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  product_difference_interpolant(surrData, hip_approx_2->surrData,
    hsg_driver->smolyak_multi_index(), hsg_driver->collocation_key(),
    hsg_driver->collocation_indices(), prod_t1c, prod_t2c, lf_key,
    set_partition);
}


inline void HierarchInterpPolyApproximation::
central_product_interpolant(HierarchInterpPolyApproximation* hip_approx_2,
			    Real mean_1, Real mean_2,
			    RealVector2DArray& cov_t1_coeffs,
			    RealMatrix2DArray& cov_t2_coeffs,
			    const UShort2DArray& set_partition)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  if (hsg_driver->track_collocation_indices() &&
      hsg_driver->collocation_indices().empty()) {// invalidated by combination
    bool same = (this == hip_approx_2);
    central_product_interpolant(hsg_driver->hierarchical_variable_sets(),
      hsg_driver->smolyak_multi_index(), hsg_driver->collocation_key(),
      expT1CoeffsIter->second, expT2CoeffsIter->second,
      hip_approx_2->expT1CoeffsIter->second,
      hip_approx_2->expT2CoeffsIter->second, same, mean_1, mean_2,
      cov_t1_coeffs, cov_t2_coeffs, set_partition);
  }
  // use SurrogateData instance + colloc_indices for forming product interp
  else if (true) // current uses are modSurrData; possible future conditional
    central_product_interpolant(modSurrData.variables_data(),
      modSurrData.response_data(), hip_approx_2->modSurrData.response_data(),
      mean_1, mean_2, hsg_driver->smolyak_multi_index(),
      hsg_driver->collocation_key(), hsg_driver->collocation_indices(),
      cov_t1_coeffs, cov_t2_coeffs, set_partition);
  else
    central_product_interpolant(surrData.variables_data(),
      surrData.response_data(), hip_approx_2->surrData.response_data(),
      mean_1, mean_2, hsg_driver->smolyak_multi_index(),
      hsg_driver->collocation_key(), hsg_driver->collocation_indices(),
      cov_t1_coeffs, cov_t2_coeffs, set_partition);
}


inline void HierarchInterpPolyApproximation::
central_product_gradient_interpolant(
  HierarchInterpPolyApproximation* hip_approx_2, Real mean_1, Real mean_2,
  const RealVector& mean1_grad, const RealVector& mean2_grad,
  RealMatrix2DArray& cov_t1c_grads, const UShort2DArray& set_partition)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  if (hsg_driver->track_collocation_indices() &&
      hsg_driver->collocation_indices().empty()) {// invalidated by combination
    bool same = (this == hip_approx_2);
    central_product_gradient_interpolant(
      hsg_driver->hierarchical_variable_sets(),
      hsg_driver->smolyak_multi_index(), hsg_driver->collocation_key(),
      expT1CoeffsIter->second, expT2CoeffsIter->second,
      expT1CoeffGradsIter->second, hip_approx_2->expT1CoeffsIter->second,
      hip_approx_2->expT2CoeffsIter->second,
      hip_approx_2->expT1CoeffGradsIter->second, same, mean_1, mean_2,
      mean1_grad, mean2_grad, cov_t1c_grads, set_partition);
  }
  // use SurrogateData instance + colloc_indices for forming product interp
  else if (true) // current uses are modSurrData; possible future conditional
    central_product_gradient_interpolant(modSurrData.variables_data(),
      modSurrData.response_data(), hip_approx_2->modSurrData.response_data(),
      mean_1, mean_2, mean1_grad, mean2_grad, hsg_driver->smolyak_multi_index(),
      hsg_driver->collocation_key(), hsg_driver->collocation_indices(),
      cov_t1c_grads, set_partition);
  else
    central_product_gradient_interpolant(surrData.variables_data(),
      surrData.response_data(), hip_approx_2->surrData.response_data(),
      mean_1, mean_2, mean1_grad, mean2_grad, hsg_driver->smolyak_multi_index(),
      hsg_driver->collocation_key(), hsg_driver->collocation_indices(),
      cov_t1c_grads, set_partition);
}


inline Real HierarchInterpPolyApproximation::
beta_map(Real mu, Real var, bool cdf_flag, Real z_bar)
{
  if (var > 0.) {
    Real stdev = std::sqrt(var);
    return (cdf_flag) ? (mu - z_bar)/stdev : (z_bar - mu)/stdev;
  }
  else
    return ( (cdf_flag && mu <= z_bar) || (!cdf_flag && mu > z_bar) ) ?
      Pecos::LARGE_NUMBER : -Pecos::LARGE_NUMBER;
}


inline Real HierarchInterpPolyApproximation::beta(bool cdf_flag, Real z_bar)
{ return beta_map(mean(), variance(), cdf_flag, z_bar); }


inline Real HierarchInterpPolyApproximation::
beta(const RealVector& x, bool cdf_flag, Real z_bar)
{ return beta_map(mean(x), variance(x), cdf_flag, z_bar); }


inline Real HierarchInterpPolyApproximation::
combined_beta(bool cdf_flag, Real z_bar)
{ return beta_map(combined_mean(), combined_variance(), cdf_flag, z_bar); }


inline Real HierarchInterpPolyApproximation::
combined_beta(const RealVector& x, bool cdf_flag, Real z_bar)
{ return beta_map(combined_mean(x), combined_variance(x), cdf_flag, z_bar); }


inline Real HierarchInterpPolyApproximation::delta_std_deviation()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key);

  return delta_std_deviation(ref_key, incr_key);
}


inline Real HierarchInterpPolyApproximation::
delta_std_deviation(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key);

  return delta_std_deviation(x, ref_key, incr_key);
}


inline Real HierarchInterpPolyApproximation::
delta_combined_std_deviation()
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  std::map<UShortArray, UShort2DArray> ref_key_map, incr_key_map;
  data_rep->hsg_driver()->partition_keys(ref_key_map, incr_key_map);
  // Note: these keys are for partition of active expansion and should not be
  //       used for partitioning a combined expansion (see {ref,incr}_key_map)

  return delta_combined_std_deviation(ref_key_map, incr_key_map);
}


inline Real HierarchInterpPolyApproximation::
delta_combined_std_deviation(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  std::map<UShortArray, UShort2DArray> ref_key_map, incr_key_map;
  data_rep->hsg_driver()->partition_keys(ref_key_map, incr_key_map);
  // Note: these keys are for partition of active expansion and should not be
  //       used for partitioning a combined expansion (see {ref,incr}_key_map)

  return delta_combined_std_deviation(x, ref_key_map, incr_key_map);
}


inline Real HierarchInterpPolyApproximation::
delta_beta(bool cdf_flag, Real z_bar)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key);

  return delta_beta(cdf_flag, z_bar, ref_key, incr_key);
}


inline Real HierarchInterpPolyApproximation::
delta_beta(const RealVector& x, bool cdf_flag, Real z_bar)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key);

  return delta_beta(x, cdf_flag, z_bar, ref_key, incr_key);
}


inline Real HierarchInterpPolyApproximation::
delta_combined_beta(bool cdf_flag, Real z_bar)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  std::map<UShortArray, UShort2DArray> ref_key_map, incr_key_map;
  data_rep->hsg_driver()->partition_keys(ref_key_map, incr_key_map);
  // Note: these keys are for partition of active expansion and should not be
  //       used for partitioning a combined expansion (see {ref,incr}_key_map)

  return delta_combined_beta(cdf_flag, z_bar, ref_key_map, incr_key_map);
}


inline Real HierarchInterpPolyApproximation::
delta_combined_beta(const RealVector& x, bool cdf_flag, Real z_bar)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  std::map<UShortArray, UShort2DArray> ref_key_map, incr_key_map;
  data_rep->hsg_driver()->partition_keys(ref_key_map, incr_key_map);
  // Note: these keys are for partition of active expansion and should not be
  //       used for partitioning a combined expansion (see {ref,incr}_key_map)

  return delta_combined_beta(x, cdf_flag, z_bar, ref_key_map, incr_key_map);
}


inline Real HierarchInterpPolyApproximation::
delta_z(bool cdf_flag, Real beta_bar)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key);

  return delta_z(cdf_flag, beta_bar, ref_key, incr_key);
}


inline Real HierarchInterpPolyApproximation::
delta_z(const RealVector& x, bool cdf_flag, Real beta_bar)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  UShort2DArray ref_key, incr_key;
  data_rep->hsg_driver()->partition_keys(ref_key, incr_key);

  return delta_z(x, cdf_flag, beta_bar, ref_key, incr_key);
}


inline Real HierarchInterpPolyApproximation::
delta_combined_z(bool cdf_flag, Real beta_bar)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  std::map<UShortArray, UShort2DArray> ref_key_map, incr_key_map;
  data_rep->hsg_driver()->partition_keys(ref_key_map, incr_key_map);
  // Note: these keys are for partition of active expansion and should not be
  //       used for partitioning a combined expansion (see {ref,incr}_key_map)

  return delta_combined_z(cdf_flag, beta_bar, ref_key_map, incr_key_map);
}


inline Real HierarchInterpPolyApproximation::
delta_combined_z(const RealVector& x, bool cdf_flag, Real beta_bar)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  std::map<UShortArray, UShort2DArray> ref_key_map, incr_key_map;
  data_rep->hsg_driver()->partition_keys(ref_key_map, incr_key_map);
  // Note: these keys are for partition of active expansion and should not be
  //       used for partitioning a combined expansion (see {ref,incr}_key_map)

  return delta_combined_z(x, cdf_flag, beta_bar, ref_key_map, incr_key_map);
}


inline Real HierarchInterpPolyApproximation::
expectation(const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs,
	    const UShort2DArray& set_partition)
{
  // This version defaults to active type1/2 wts from HierarchSparseGridDriver
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  return expectation(t1_coeffs, t2_coeffs,
		     hsg_driver->type1_hierarchical_weight_sets(),
		     hsg_driver->type2_hierarchical_weight_sets(),
		     set_partition);
}


inline Real HierarchInterpPolyApproximation::
expectation(const RealVector& x, const RealVector2DArray& t1_coeffs,
	    const RealMatrix2DArray& t2_coeffs,
	    const UShort2DArray& set_partition)
{
  // This version defaults to active sm_mi/key from HierarchSparseGridDriver
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  return expectation(x, t1_coeffs, t2_coeffs, hsg_driver->smolyak_multi_index(),
		     hsg_driver->collocation_key(), set_partition);
}


inline const RealVector& HierarchInterpPolyApproximation::
expectation_gradient(const RealMatrix2DArray& t1_coeff_grads)
{
  // This version defaults to active type1 wts from HierarchSparseGridDriver
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  return expectation_gradient(t1_coeff_grads,
    data_rep->hsg_driver()->type1_hierarchical_weight_sets());
}


inline Real HierarchInterpPolyApproximation::
expectation_gradient(const RealVector& x,
		     const RealMatrix2DArray& t1_coeff_grads, size_t t1cg_index)
{
  // This version defaults to active sm_mi/key from HierarchSparseGridDriver
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  return expectation_gradient(x, t1_coeff_grads,
			      hsg_driver->smolyak_multi_index(),
			      hsg_driver->collocation_key(), t1cg_index);
}


inline Real HierarchInterpPolyApproximation::
expectation_gradient(const RealVector& x, const RealVector2DArray& t1_coeffs,
		     const RealMatrix2DArray& t2_coeffs, size_t deriv_index)
{
  // This version defaults to active sm_mi/key from HierarchSparseGridDriver
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  return expectation_gradient(x, t1_coeffs, t2_coeffs,
			      hsg_driver->smolyak_multi_index(),
			      hsg_driver->collocation_key(), deriv_index);
}


inline Real HierarchInterpPolyApproximation::value(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return value(x, sm_mi, hsg_driver->collocation_key(), expT1CoeffsIter->second,
	       expT2CoeffsIter->second, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_basis_variables(x, sm_mi, hsg_driver->collocation_key(),
				  expT1CoeffsIter->second,
				  expT2CoeffsIter->second, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_basis_variables(x, sm_mi, hsg_driver->collocation_key(),
				  expT1CoeffsIter->second,
				  expT2CoeffsIter->second, dvv, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_nonbasis_variables(x, sm_mi, hsg_driver->collocation_key(),
				     expT1CoeffGradsIter->second, max_level);
}


inline const RealSymMatrix& HierarchInterpPolyApproximation::
hessian_basis_variables(const RealVector& x)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index();
  unsigned short   max_level = sm_mi.size() - 1;
  return hessian_basis_variables(x, sm_mi, hsg_driver->collocation_key(),
				 expT1CoeffsIter->second, max_level);
}


inline Real HierarchInterpPolyApproximation::
stored_value(const RealVector& x, const UShortArray& key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index(key);
  unsigned short   max_level = sm_mi.size() - 1;
  return value(x, sm_mi, hsg_driver->collocation_key(key),
	       expansionType1Coeffs[key], expansionType2Coeffs[key], max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
stored_gradient_basis_variables(const RealVector& x, const UShortArray& key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index(key);
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_basis_variables(x, sm_mi, hsg_driver->collocation_key(key),
				  expansionType1Coeffs[key],
				  expansionType2Coeffs[key], max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
stored_gradient_basis_variables(const RealVector& x, const SizetArray& dvv,
				const UShortArray& key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index(key);
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_basis_variables(x, sm_mi, hsg_driver->collocation_key(key),
				  expansionType1Coeffs[key],
				  expansionType2Coeffs[key], dvv, max_level);
}


inline const RealVector& HierarchInterpPolyApproximation::
stored_gradient_nonbasis_variables(const RealVector& x, const UShortArray& key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index(key);
  unsigned short   max_level = sm_mi.size() - 1;
  return gradient_nonbasis_variables(x, sm_mi, hsg_driver->collocation_key(key),
				     expansionType1CoeffGrads[key], max_level);
}


inline const RealSymMatrix& HierarchInterpPolyApproximation::
stored_hessian_basis_variables(const RealVector& x, const UShortArray& key)
{
  SharedHierarchInterpPolyApproxData* data_rep
    = (SharedHierarchInterpPolyApproxData*)sharedDataRep;
  HierarchSparseGridDriver* hsg_driver = data_rep->hsg_driver();
  const UShort3DArray& sm_mi = hsg_driver->smolyak_multi_index(key);
  unsigned short   max_level = sm_mi.size() - 1;
  return hessian_basis_variables(x, sm_mi, hsg_driver->collocation_key(key),
				 expansionType1Coeffs[key], max_level);
}

} // namespace Pecos

#endif
