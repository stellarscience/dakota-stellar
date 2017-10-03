/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        ProjectOrthogPolyApproximation
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        Mike Eldred

#ifndef PROJECT_ORTHOG_POLY_APPROXIMATION_HPP
#define PROJECT_ORTHOG_POLY_APPROXIMATION_HPP

#include "OrthogPolyApproximation.hpp"

namespace Pecos {


/// Derived approximation class for multivariate orthogonal polynomial
/// approximation with coefficient estimation via numerical integration.

/** The ProjectOrthogPolyApproximation class provides a global approximation
    based on multivariate orthogonal polynomials, where the coefficients are
    computed using numerical integration approaches such as quadrature,
    cubature, sparse grids, and random sampling.  It is used primarily for
    polynomial chaos expansion aproaches to UQ. */

class ProjectOrthogPolyApproximation: public OrthogPolyApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  ProjectOrthogPolyApproximation(const SharedBasisApproxData& shared_data);
  /// destructor
  ~ProjectOrthogPolyApproximation();

  //
  //- Heading: Member functions
  //

  //
  //- Cosmin Heading: Virtual function redefinitions and member functions
  //
  void compute_coefficients();
  void push_coefficients();
  void increment_coefficients();
  void decrement_coefficients();
  void finalize_coefficients();

protected:

  //
  //- Heading: Virtual function redefinitions and member functions
  //

  int  min_coefficients() const;

  /// initialize polynomialBasis, multiIndex, et al.
  void allocate_arrays();

  Real value(const RealVector& x);
  Real stored_value(const RealVector& x, size_t index);

  /// compute numerical moments to order 4 and expansion moments to order 2
  void compute_moments();
  /// compute expansion moments in all-variables mode to order 2
  void compute_moments(const RealVector& x);

private:

  //
  //- Heading: Member functions
  //

  /// initialize multi_index using a sparse grid expansion
  void sparse_grid_multi_index(UShort2DArray& multi_index);
  // initialize tp_multi_index from tpMultiIndexMap
  //void map_tensor_product_multi_index(UShort2DArray& tp_multi_index,
  //				        size_t tp_index);

  /// extract tp_data_points from surrData and tp_weights from
  /// driverRep->type1CollocWts1D
  void integration_data(size_t tp_index, SDVArray& tp_data_vars,
			SDRArray& tp_data_resp, RealVector& tp_weights);
  /// computes the chaosCoeffs via numerical integration (expCoeffsSolnApproach
  /// can be QUADRATURE, CUBATURE, or COMBINED_SPARSE_GRID)
  void integrate_expansion(const UShort2DArray& multi_index,
			   const SDVArray& data_vars, const SDRArray& data_resp,
			   const RealVector& wt_sets, RealVector& exp_coeffs,
			   RealMatrix& exp_coeff_grads);

  /// computes the chaosCoeffs via averaging of samples
  /// (expCoeffsSolnApproach is SAMPLING)
  void expectation();

  /// update expansion{Coeffs,CoeffGrads} by adding one or more tensor-product
  /// expansions and updating all Smolyak coefficients
  void append_tensor_expansions(size_t start_index);

  /// update numericalMoments using numerical integration applied
  /// directly to surrData
  void integrate_response_moments(size_t num_moments);

  //
  //- Heading: Data
  //

  /// previous expansionCoeffs (aggregated total, not tensor-product
  /// contributions) prior to append_tensor_expansions()
  RealVector prevExpCoeffs;
  /// previous expansionCoeffGrads (aggregated total, not tensor-product
  /// contributions) prior to append_tensor_expansions()
  RealMatrix prevExpCoeffGrads;

  /// the set of tensor-product contributions to expansionCoeffs
  RealVectorArray tpExpansionCoeffs;
  /// the set of tensor-product contributions to expansionCoeffGrads
  RealMatrixArray tpExpansionCoeffGrads;

  /// popped instances of tpExpansionCoeffs that were computed but not selected
  std::deque<RealVector> poppedTPExpCoeffs;
  /// popped tpExpansionCoeffGrads instances that were computed but not selected
  std::deque<RealMatrix> poppedTPExpCoeffGrads;
};


inline ProjectOrthogPolyApproximation::
ProjectOrthogPolyApproximation(const SharedBasisApproxData& shared_data):
  OrthogPolyApproximation(shared_data)
{ }


inline ProjectOrthogPolyApproximation::~ProjectOrthogPolyApproximation()
{ }


inline void ProjectOrthogPolyApproximation::compute_moments()
{
  // standard variables mode supports two expansion and four numerical moments
  mean(); variance(); // updates expansionMoments[0] and [1]
  //standardize_moments(expansionMoments);

  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  if (data_rep->expConfigOptions.expCoeffsSolnApproach != SAMPLING)
    integrate_response_moments(4);
}


inline void ProjectOrthogPolyApproximation::compute_moments(const RealVector& x)
{
  // all variables mode only supports first two moments
  mean(x); variance(x); // updates expansionMoments[0] and [1]
  //standardize_moments(expansionMoments);

  //SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  //if (data_rep->expConfigOptions.expCoeffsSolnApproach != SAMPLING)
  //  integrate_response_moments(2, x); // TO DO
}

} // namespace Pecos

#endif
