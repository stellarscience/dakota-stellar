/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        NodalInterpPolyApproximation
//- Description:  Class for Nodal Interpolation Polynomial Approximation
//-               
//- Owner:        Mike Eldred

#ifndef NODAL_INTERP_POLY_APPROXIMATION_HPP
#define NODAL_INTERP_POLY_APPROXIMATION_HPP

#include "InterpPolyApproximation.hpp"

namespace Pecos {


/// Derived approximation class for nodal interpolation polynomials
/// (global approximation interpolating function values and
/// potentially gradients at collocation points).

/** The NodalInterpPolyApproximation class provides a global polynomial
    approximation based on either Lagrange or Hermite interpolation
    polynomials using a nodal basis approach.  It is used primarily
    for stochastic collocation approaches to uncertainty quantification. */

class NodalInterpPolyApproximation: public InterpPolyApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  NodalInterpPolyApproximation(const SharedBasisApproxData& shared_data);
  /// destructor
  ~NodalInterpPolyApproximation();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// size expansionType{1,2}Coeffs and expansionType1CoeffGrads
  void allocate_arrays();

  void compute_expansion_coefficients();

  /// update the coefficients for the expansion of interpolation polynomials:
  /// increment expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void increment_coefficients();
  /// restore the coefficients to their previous state prior to last increment:
  /// decrement expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void decrement_coefficients();
  /// restore the coefficients to a previously incremented state as
  /// identified by the current increment to the Smolyak multi index:
  /// push expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void push_coefficients();
  /// finalize the coefficients by applying all previously evaluated increments:
  /// finalize expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  void finalize_coefficients();

  /// store current state within storedExpType{1Coeffs,2Coeffs,1CoeffGrads}
  void store_coefficients(size_t index = _NPOS);
  /// restore previous state from storedExpType{1Coeffs,2Coeffs,1CoeffGrads}
  void restore_coefficients(size_t index = _NPOS);
  /// swap storedExpType{1Coeffs,2Coeffs,1CoeffGrads}[index] with active
  /// current data
  void swap_coefficients(size_t index);
  /// remove a redundant entry from storedExpType{1Coeffs,2Coeffs,1CoeffGrads}
  /// prior to combine_coefficients (default is pop_back)
  void remove_stored_coefficients(size_t index = _NPOS);

  /// augment current interpolant using
  /// storedExpType{1Coeffs,2Coeffs,1CoeffGrads}[index]
  void combine_coefficients(short combine_type, size_t swap_index);

  void integrate_response_moments(size_t num_moments);
  void integrate_expansion_moments(size_t num_moments);

  Real value(const RealVector& x);
  const RealVector& gradient_basis_variables(const RealVector& x);
  const RealVector& gradient_basis_variables(const RealVector& x,
					     const SizetArray& dvv);
  const RealVector& gradient_nonbasis_variables(const RealVector& x);
  const RealSymMatrix& hessian_basis_variables(const RealVector& x);

  Real stored_value(const RealVector& x, size_t index);
  const RealVector& stored_gradient_basis_variables(const RealVector& x,
						    size_t index);
  const RealVector& stored_gradient_nonbasis_variables(const RealVector& x,
						       size_t index);

  Real mean();
  Real mean(const RealVector& x);
  const RealVector& mean_gradient();
  const RealVector& mean_gradient(const RealVector& x, const SizetArray& dvv);

  Real variance();
  Real variance(const RealVector& x);
  const RealVector& variance_gradient();
  const RealVector& variance_gradient(const RealVector& x,
				      const SizetArray& dvv);

  Real covariance(PolynomialApproximation* poly_approx_2);
  Real covariance(const RealVector& x, PolynomialApproximation* poly_approx_2);

  void compute_total_sobol_indices();
  void compute_partial_variance(const BitArray& set_value);

  RealVector approximation_coefficients(bool normalized) const;
  void approximation_coefficients(const RealVector& approx_coeffs,
				  bool normalized);

private:

  //
  //- Heading: Convenience functions
  //

  /// compute the mean of a tensor interpolant on a tensor grid;
  /// contributes to mean(x)
  Real tensor_product_mean(const RealVector& x, const UShortArray& lev_index,
    const UShort2DArray& key, const SizetArray& colloc_index);

  /// compute the gradient of the mean of a tensor interpolant on a
  /// tensor grid; contributes to mean_gradient(x)
  const RealVector& tensor_product_mean_gradient(const RealVector& x,
    const UShortArray& lev_index,   const UShort2DArray& key,
    const SizetArray& colloc_index, const SizetArray& dvv);

  /// compute the covariance of two tensor interpolants on the same
  /// tensor grid using an interpolation of products or product of
  /// interpolants approach; contributes to covariance(x, poly_approx_2)
  Real tensor_product_covariance(const RealVector& x,
    const UShortArray& lev_index,   const UShort2DArray& key,
    const SizetArray& colloc_index, NodalInterpPolyApproximation* nip_approx_2);
  /// compute the covariance of two tensor interpolants on different
  /// tensor grids using a product of interpolants approach;
  /// contributes to covariance(x, poly_approx_2)
  Real tensor_product_covariance(const RealVector& x,
    const UShortArray& lev_index_1, const UShort2DArray& key_1,
    const SizetArray& colloc_index_1, const UShortArray& lev_index_2,
    const UShort2DArray& key_2, const SizetArray& colloc_index_2,
    NodalInterpPolyApproximation* nip_approx_2);

  /// compute the gradient of the variance of a tensor interpolant on
  /// a tensor grid using an interpolation of products or product of
  /// interpolants approach; contributes to variance_gradient(x)
  const RealVector& tensor_product_variance_gradient(const RealVector& x,
    const UShortArray& lev_index, const UShort2DArray& key,
    const SizetArray& colloc_index, const SizetArray& dvv);

  /// compute value of reduced-dimension interpolant
  Real value(const RealVector& x, const RealVectorArray& t1_coeffs,
	     const RealMatrixArray& t2_coeffs, const UShort3DArray& colloc_key,
	     const SizetList& subset_indices);
  /// compute gradient of reduced-dimension interpolant with respect
  /// to basis variables
  const RealVector& gradient_basis_variables(const RealVector& x,
					     const RealVectorArray& t1_coeffs,
					     const RealMatrixArray& t2_coeffs,
					     const UShort3DArray& colloc_key,
					     const SizetList& subset_indices);

  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using weights from the CombinedSparseGridDriver
  Real expectation(const RealVector& t1_coeffs, const RealMatrix& t2_coeffs);
  /// compute the expected value of the interpolant given by t{1,2}_coeffs
  /// using t{1,2}_wts
  Real expectation(const RealVector& t1_coeffs, const RealVector& t1_wts,
		   const RealMatrix& t2_coeffs, const RealMatrix& t2_wts);

  /// computes higher-order grid for tensor reinterpolation of the
  /// covariance fn for non-integrated dimensions in all_variables mode
  void reinterpolated_level(const UShortArray& lev_index);

  /// compute integral for total Sobol' index for variables in a set
  Real member_integral(const BitArray& member_bits, Real mean);
  /// defines member_coeffs and member_wts for a particular membership set
  void member_coefficients_weights(const BitArray& member_bits,
    const UShortArray& quad_order,    const UShortArray& lev_index,
    const UShort2DArray& colloc_key,  const SizetArray& colloc_index,
    RealVector& member_t1_coeffs,     RealVector& member_t1_wts,
    RealMatrix& member_t2_coeffs,     RealMatrix& member_t2_wts,
    UShort2DArray& member_colloc_key, SizetArray& member_colloc_index);
  /// create a unique map key for value() and gradient() calculation reuse
  void update_member_key(const UShortArray& data,
			 const SizetList& member_indices,
			 UShortArray& member_map_key, size_t cntr);

  //
  //- Heading: Data
  //

  /// the type1 coefficients of the expansion for interpolating values
  RealVector expansionType1Coeffs;
  /// the type2 coefficients of the expansion for interpolating gradients
  RealMatrix expansionType2Coeffs;
  /// the gradients of the type1 expansion coefficients
  /** may be interpreted as either the gradients of the expansion
      coefficients or the coefficients of expansions for the response
      gradients.  This array is used when sensitivities of moments are
      needed with respect to variables that do not appear in the
      expansion (e.g., with respect to design variables for an
      expansion only over the random variables). */
  RealMatrix expansionType1CoeffGrads;

  /// copies of expansionType1Coeffs state for subsequent restoration
  RealVectorArray storedExpType1Coeffs;
  /// copies of expansionType2Coeffs state for subsequent restoration
  RealMatrixArray storedExpType2Coeffs;
  /// copies of expansionType1CoeffGrads state for subsequent restoration
  RealMatrixArray storedExpType1CoeffGrads;

  /// the gradient of the mean of a tensor-product interpolant; a
  /// contributor to meanGradient
  RealVector tpMeanGrad;     // TO DO: move to shared data? (2nd pass tuning)
  /// the gradient of the variance of a tensor-product interpolant; a
  /// contributor to varianceGradient
  RealVector tpVarianceGrad; // TO DO: move to shared data? (2nd pass tuning)
};


inline NodalInterpPolyApproximation::
NodalInterpPolyApproximation(const SharedBasisApproxData& shared_data):
  InterpPolyApproximation(shared_data)
{ }


inline NodalInterpPolyApproximation::~NodalInterpPolyApproximation()
{ }


inline void NodalInterpPolyApproximation::increment_coefficients()
{ push_coefficients(); allocate_component_sobol(); }


inline void NodalInterpPolyApproximation::decrement_coefficients()
{
  size_t new_colloc_pts = surrData.points();
  if (expansionCoeffFlag) {
   expansionType1Coeffs.resize(new_colloc_pts);
   SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
   if (data_rep->basisConfigOptions.useDerivs) {
     size_t num_deriv_vars = expansionType2Coeffs.numRows();
     expansionType2Coeffs.reshape(num_deriv_vars, new_colloc_pts);
   }
  }
  if (expansionCoeffGradFlag) {
   size_t num_deriv_vars = expansionType1CoeffGrads.numRows();
   expansionType1CoeffGrads.reshape(num_deriv_vars, new_colloc_pts);
  }
}


inline void NodalInterpPolyApproximation::finalize_coefficients()
{ push_coefficients(); }


inline RealVector NodalInterpPolyApproximation::
approximation_coefficients(bool normalized) const
{
  if (normalized)
    PCerr << "Warning: normalized coefficients not supported in "
	  << "InterpPolyApproximation export." << std::endl;
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  if (data_rep->basisConfigOptions.useDerivs) {
    PCerr << "Error: approximation_coefficients() not supported in "
	  << "InterpPolyApproximation for type2 coefficients." << std::endl;
    return abort_handler_t<const RealVector&>(-1);
  }
  else
    return RealVector(Teuchos::View, expansionType1Coeffs.values(),
		      expansionType1Coeffs.length());
}


inline void NodalInterpPolyApproximation::
approximation_coefficients(const RealVector& approx_coeffs, bool normalized)
{
  if (normalized)
    PCerr << "Warning: normalized coefficients not supported in "
	  << "InterpPolyApproximation import." << std::endl;
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  if (data_rep->basisConfigOptions.useDerivs) {
    PCerr << "Error: approximation_coefficients() not supported in "
	  << "InterpPolyApproximation for type2 coefficients." << std::endl;
    abort_handler(-1);
  }
  else {
    expansionType1Coeffs = approx_coeffs;

    allocate_total_sobol();
    allocate_component_sobol();

    if (numericalMoments.empty()) {
      SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
      size_t num_moments = (data_rep->nonRandomIndices.empty()) ? 4 : 2;
      numericalMoments.sizeUninitialized(num_moments);
    }
  }
}


/** In this case, all expansion variables are random variables and the
    variance of the expansion uses an interpolation of response products. */
inline Real NodalInterpPolyApproximation::variance()
{ return covariance(this); }


/** In this case, a subset of the expansion variables are random
    variables and the variance of the expansion involves integration
    over this subset and evaluation over the subset's complement. */
inline Real NodalInterpPolyApproximation::variance(const RealVector& x)
{ return covariance(x, this); }


inline Real NodalInterpPolyApproximation::
expectation(const RealVector& t1_coeffs, const RealMatrix& t2_coeffs)
{
  // This version defaults to type1/2 weights from CombinedSparseGridDriver
  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  IntegrationDriver* int_driver = data_rep->driverRep;
  return expectation(t1_coeffs, int_driver->type1_weight_sets(),
		     t2_coeffs, int_driver->type2_weight_sets());
}


inline void NodalInterpPolyApproximation::
update_member_key(const UShortArray& data,
		  const SizetList&   member_indices,
		  UShortArray& member_map_key, size_t cntr)
{
  for (SizetList::const_iterator cit=member_indices.begin();
       cit!=member_indices.end(); ++cit, ++cntr)
    member_map_key[cntr] = data[*cit];
}

} // namespace Pecos

#endif
