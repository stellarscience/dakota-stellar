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
#include "SharedProjectOrthogPolyApproxData.hpp"

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
  //- Heading: Virtual function redefinitions
  //

  void compute_coefficients();
  void increment_coefficients();
  void pop_coefficients(bool save_data);
  void push_coefficients();
  void finalize_coefficients();

  //
  //- Heading: Member functions
  //

protected:

  //
  //- Heading: Virtual function redefinitions and member functions
  //

  int min_coefficients() const;

  void allocate_arrays();

  Real value(const RealVector& x);
  Real stored_value(const RealVector& x, const UShortArray& key);

  /// compute numerical moments to order 4 and expansion moments to order 2
  void compute_moments(bool full_stats = true, bool combined_stats = false);
  // compute expansion moments in all-variables mode to order 2
  //void compute_moments(const RealVector& x, bool full_stats = true,
  //		       bool combined_stats = false);

  void combined_to_active(bool clear_combined = true);

private:

  //
  //- Heading: Member functions
  //

  /// perform sanity checks prior to numerical integration
  void integration_checks();

  /// initialize multi_index using a sparse grid expansion
  void sparse_grid_multi_index(UShort2DArray& multi_index);
  // initialize tp_multi_index from tpMultiIndexMap
  //void map_tensor_product_multi_index(UShort2DArray& tp_multi_index,
  //				        size_t tp_index);

  /// extract tp_data_points from modSurrData and tp_weights from
  /// driverRep->type1CollocWts1D
  void integration_data(size_t tp_index, SDVArray& tp_data_vars,
			SDRArray& tp_data_resp, RealVector& tp_weights);
  /// computes the chaosCoeffs via numerical integration (expCoeffsSolnApproach
  /// can be QUADRATURE, CUBATURE, or COMBINED_SPARSE_GRID)
  void integrate_expansion(const UShort2DArray& multi_index,
			   const SDVArray& data_vars, const SDRArray& data_resp,
			   const RealVector& wt_sets, RealVector& exp_coeffs,
			   RealMatrix& exp_coeff_grads);

  /// update surr_data with synthetic data (from evaluating the expansion)
  /// following promotion of combined expansion to active (simplifies stats
  /// computation for FINAL_RESULTS)
  void synthetic_surrogate_data(SurrogateData& surr_data);

  /// computes the chaosCoeffs via averaging of samples
  /// (expCoeffsSolnApproach is SAMPLING)
  void expectation();

  /// update expansion{Coeffs,CoeffGrads} by adding one or more tensor-product
  /// expansions and updating all Smolyak coefficients
  void append_tensor_expansions(size_t start_index);

  /// update numericalMoments using numerical integration applied
  /// directly to modSurrData
  void integrate_response_moments(size_t num_moments);//, bool combined_stats);

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
  std::map<UShortArray, RealVectorArray> tpExpansionCoeffs;
  /// the set of tensor-product contributions to expansionCoeffGrads
  std::map<UShortArray, RealMatrixArray> tpExpansionCoeffGrads;

  /// popped instances of either expansionCoeffs or tpExpansionCoeffs,
  /// depending on exp soln approach, that were computed but not selected
  std::map<UShortArray, RealVectorDeque> poppedExpCoeffs;
  /// popped instances of either expansionCoeffGrads or tpExpansionCoeffGrads,
  /// depending on exp soln approach, that were computed but not selected
  std::map<UShortArray, RealMatrixDeque> poppedExpCoeffGrads;
};


inline ProjectOrthogPolyApproximation::
ProjectOrthogPolyApproximation(const SharedBasisApproxData& shared_data):
  OrthogPolyApproximation(shared_data)
{ }


inline ProjectOrthogPolyApproximation::~ProjectOrthogPolyApproximation()
{ }


inline void ProjectOrthogPolyApproximation::
combined_to_active(bool clear_combined)
{
  OrthogPolyApproximation::combined_to_active(clear_combined);

  // Create synthetic modSurrData for the combined-now-active coeffs, for
  // supporting FINAL_RESULTS processing (numerical moments on combined grid)
  // Note: exclude CUBATURE and SAMPLING, which lack combined grids.
  SharedProjectOrthogPolyApproxData* data_rep
    = (SharedProjectOrthogPolyApproxData*)sharedDataRep;
  switch (data_rep->expConfigOptions.expCoeffsSolnApproach) {
  case QUADRATURE: case COMBINED_SPARSE_GRID: case INCREMENTAL_SPARSE_GRID:
    synthetic_surrogate_data(modSurrData); // overwrite data for activeKey
    break;
  }
}


inline void ProjectOrthogPolyApproximation::
compute_moments(bool full_stats, bool combined_stats)
{
  // standard variables mode supports two expansion and four numerical moments

  RealVector& exp_mom = expMomentsIter->second;
  if (exp_mom.length() != 2) exp_mom.resize(2);
  if (combined_stats)
    { combined_mean(); combined_variance(); } // for combinedExpCoeffs
  else {
    mean(); variance(); // updates first two expansionMoments
    //standardize_moments(exp_mom);
  }

  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  // if full stats, augment analytic expansion moments with numerical moments
  // (from quadrature applied to the SurrogateData)
  if (full_stats &&
      // > currently supported by TPQ, SSG, Cubature (Sampling also possible)
      data_rep->expConfigOptions.expCoeffsSolnApproach != SAMPLING) { //&&
      // > combined expansions do not admit a unified set of collocation data
      //   and backfilling direct response data with interpolated surrogate
      //   values violates some of the intent (Other considerations: adds post-
      //   processing expense, but adds support for higher order moments)
    //!data_rep->expConfigOptions.combineType)

    if (combined_stats) {
      // current uses follow combined_to_active(), so don't need this for now
      PCerr << "Error: combined mode unavailable for final stats.  Project"
	    << "OrthogPolyApproximation::compute_moments()\n       currently "
	    << "requires promotion of combined to active." << std::endl;
      abort_handler(-1);
    }
    integrate_response_moments(4);//, combined_stats); // TO DO
  }
  else
    numMomentsIter->second.resize(0);
}


/* OrthogPolyApproximation::compute_moments(const RealVector&) is used for now

inline void ProjectOrthogPolyApproximation::
compute_moments(const RealVector& x, bool full_stats, bool combined_stats)
{
  // all variables mode currently supports two expansion moments

  RealVector& exp_mom = expMomentsIter->second;
  if (exp_mom.length() != 2) exp_mom.resize(2);
  if (combined_stats)
    { combined_mean(x); combined_variance(x); } // for combinedExpCoeffs
  else {
    mean(x); variance(x); // updates first two expansionMoments
    //standardize_moments(exp_mom);
  }

  SharedPolyApproxData* data_rep = (SharedPolyApproxData*)sharedDataRep;
  if (full_stats &&
      data_rep->expConfigOptions.expCoeffsSolnApproach != SAMPLING) {

    // current uses follow combined_to_active(), so don't need this for now
    if (combined_stats) {
      PCerr << "Error: combined_stats unavailable.  ProjectOrthogPoly"
            << "Approximation::compute_moments()\n       currently "
            << "requires promotion of combined to active." << std::endl;
      abort_handler(-1);
    }
    integrate_response_moments(2, x);//, combined_stats); // TO DO
  }
  else
    numMomentsIter->second.resize(0);
}
*/

} // namespace Pecos

#endif
