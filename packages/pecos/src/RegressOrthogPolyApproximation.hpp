/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        RegressOrthogPolyApproximation
//- Description:  Class for Multivariate Orthogonal Polynomial Approximations
//-               
//- Owner:        John Jakeman

#ifndef REGRESS_ORTHOG_POLY_APPROXIMATION_HPP
#define REGRESS_ORTHOG_POLY_APPROXIMATION_HPP

#include "OrthogPolyApproximation.hpp"
#include "LinearSolver.hpp"
#include "FaultTolerance.hpp"
#include "SharedRegressOrthogPolyApproxData.hpp"

namespace Pecos {


/// Derived approximation class for multivariate orthogonal polynomial
/// approximation with coefficient estimation via regression.

/** The RegressOrthogPolyApproximation class provides a global
    approximation based on multivariate orthogonal polynomials, where
    the coefficients are computed using regression approaches such as
    least squares (L2) or compressed sensing (L1).  It is used
    primarily for polynomial chaos expansion aproaches to UQ. */

class RegressOrthogPolyApproximation: public OrthogPolyApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  RegressOrthogPolyApproximation(const SharedBasisApproxData& shared_data);
  /// destructor
  ~RegressOrthogPolyApproximation();

  //
  //- Heading: Member functions
  //

  /// Build the pce vandermonde matrix A and extract the function (and gradient)
  /// data b so that we can solve (possible approximately) Ax=b; also populate
  /// the matrix of points corresponding to the random variable sample set
  void build_linear_system( RealMatrix &A, RealMatrix &B, RealMatrix &points );
  /// Build the pce vandermonde matrix A and extract the function (and gradient)
  /// data b so that we can solve (possible approximately) Ax=b using the
  /// provided multi-index; also populate the matrix of points corresponding
  /// to the random variable sample set
  void build_linear_system( RealMatrix &A, RealMatrix &B, RealMatrix &points,
			    const UShort2DArray& multi_index );

  /// Build the pce vandermonde matrix A and extract the function (and
  /// gradient) data b so that we can solve (possible approximately) Ax=b
  void build_linear_system( RealMatrix &A, RealMatrix &B );
  /// Build the pce vandermonde matrix A and extract the function (and
  /// gradient) data b so that we can solve (possible approximately)
  /// Ax=b using the provided multi-index
  void build_linear_system( RealMatrix &A, RealMatrix &B,
			    const UShort2DArray& multi_index );

  /// Build the pce vandermonde matrix A used in solving the system Ax=b
  void build_linear_system( RealMatrix &A );
  /// Build the pce vandermonde matrix A used in solving the system Ax=b
  /// using the provided multi-index
  void build_linear_system( RealMatrix &A, const UShort2DArray& multi_index );

  /// augment a Vandermonde matrix with additional sample points for a fixed
  /// multi-index
  void augment_linear_system(const RealVectorArray& samples, RealMatrix& A,
			     const UShort2DArray& multi_index);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  int min_coefficients() const;
  void compute_coefficients();

  void increment_coefficients();
  void store_coefficients(size_t index = _NPOS);
  void restore_coefficients(size_t index = _NPOS);
  void swap_coefficients(size_t index);
  void remove_stored_coefficients(size_t index = _NPOS);

  void combine_coefficients(short combine_type, size_t swap_index);

  void allocate_arrays();

  size_t expansion_terms() const;
  const RealVector& dimension_decay_rates();

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

  Real mean(const RealVector& x);
  const RealVector& mean_gradient(const RealVector& x, const SizetArray& dvv);

  Real covariance(PolynomialApproximation* poly_approx_2);
  Real covariance(const RealVector& x, PolynomialApproximation* poly_approx_2);
  const RealVector& variance_gradient();
  const RealVector& variance_gradient(const RealVector& x,
				      const SizetArray& dvv);

  void compute_component_sobol();
  void compute_total_sobol();
  ULongULongMap sparse_sobol_index_map() const;

  RealVector approximation_coefficients(bool normalized) const;
  void approximation_coefficients(const RealVector& approx_coeffs,
				  bool normalized);

  void print_coefficients(std::ostream& s, bool normalized);
  void coefficient_labels(std::vector<std::string>& all_coeff_tags) const;

private:

  //
  //- Heading: Member functions
  //

  /// set the information needed to ensure fault tolerance
  void set_fault_info();

  /// selects the solver for L1 or L2 minimization based on user input
  void select_solver(bool cv_active);

  /// Run the regression method set in select_solver() to compute the
  /// expansion coefficients using L1 or L2 minimization
  void run_regression();

  /// perform an adaptive selection of candidate basis for fixed data,
  /// employing a generalized sparse grid to define the candidate basis
  /// index sets
  void adapt_regression();
  /// from among the active index sets, select the candidate refinement that
  /// provides the greatest reduction in cross-validation error
  Real select_best_active_multi_index();
  /// from among the candidate basis expansions, select the option that
  /// provides the greatest reduction in cross-validation error
  Real select_best_basis_expansion();

  /// Use cross validation to choose solver hyper-parameters when
  /// solving the linear system Ax=b. e.g. if the linear solver has an
  /// epsilon tolerance internally select the best epsilon and return
  /// the corresponding solution
  Real run_cross_validation_solver(const UShort2DArray& multi_index,
				   RealVector& exp_coeffs,
				   SizetSet& sparse_indices);

  /// Use cross validation to find the hyper-parameters of the polynomial
  /// chaos expansion. e.g. find the 'best' total degree basis
  Real run_cross_validation_expansion();

  /// encapsulate usage of CSTool.solve() and bookkeeping of its sparse solution
  void compressed_sensing(RealMatrix &A, RealMatrix &B);

  /// For a specific vandermonde matrix find the compressed sennsing 
  // options that produce the best PCE
  void estimate_compressed_sensing_options_via_cross_validation(
    RealMatrix &vandermonde_matrix, RealMatrix &rhs,
    std::vector<CompressedSensingOptions> &best_cs_opts,
    RealVector &best_predictor_indicators,
    RealMatrixArray &predictor_options_history,
    RealMatrixArray &predictor_indicators_history,
    RealMatrixArray &predictor_partition_indicators_history,
    size_t num_data_pts_fn );

  /// define multiIndex and expansionCoeffs from nonzero dense_coeffs
  void update_sparse(Real* dense_coeffs, size_t num_dense_terms);
  /// augment sparse_indices based on nonzero dense_coeffs
  void update_sparse_indices(Real* dense_coeffs, size_t num_dense_terms,
			     SizetSet& sparse_indices);
  /// define sparse expansionCoeffs from dense_coeffs and sparse_indices
  void update_sparse_coeffs(Real* dense_coeffs, RealVector& exp_coeffs,
			    const SizetSet& sparse_indices);
  /// define a row of sparse expansionCoeffGrads from dense_coeffs and
  /// sparse_indices
  void update_sparse_coeff_grads(Real* dense_coeffs, int row,
				 RealMatrix& exp_coeff_grads,
				 const SizetSet& sparse_indices);
  /// define sparseSobolIndexMap from sparseIndices, shared multi_index,
  /// and shared sobolIndexMap
  void update_sparse_sobol(const SizetSet& sparse_indices,
			   const UShort2DArray& shared_multi_index,
			   const BitArrayULongMap& shared_sobol_map);

  /// Perform restriction from dense arrays + sparse_indices key into
  /// packed arrays without key
  void sparse_restriction(UShort2DArray& multi_index, SizetSet& sparse_indices);
  /// Perform restriction from original multi_index by defining a Pareto
  /// frontier of recovered terms and then pruning multi_index back to a
  /// complete set (no gaps) within this Pareto frontier
  void frontier_restriction(UShort2DArray& multi_index,
			    SizetSet& sparse_indices);

  /// perform SharedOrthogPolyApproxData::numAdvancements expansions of 
  /// multi_index to create the candidates array
  void advance_multi_index(const UShort2DArray& multi_index,
			   UShortArraySetArray& mi_advancements);
  /// perform SharedOrthogPolyApproxData::numAdvancements expansions of 
  /// multi_index to create the candidates array
  void advance_multi_index_front(const UShort2DArray& multi_index,
				 UShortArraySetArray& mi_advancements);
  /// generate a set of admissible forward neighbors from a reference
  /// multi-index (non-frontier version)
  void add_admissible_forward_neighbors(const UShort2DArray& reference_mi,
					UShortArraySet& fwd_neighbors);
  /// generate a set of admissible forward neighbors from a reference
  /// multi-index (frontier version)
  void add_admissible_forward_neighbors(const UShortArraySet& reference_mi,
					UShortArraySet& fwd_neighbors);

  /// define a multi-index frontier from the incoming multi_index.  This differs
  /// from a Pareto frontier in that the definition of dominated is relaxed
  /// (must not be < another term in all dimensions, but can be equal).
  void define_frontier(const UShort2DArray& multi_index,
		       UShortArraySet& combined_pareto);
  /// update a multi-index frontier from an incoming multi_index term.  This
  /// differs from a Pareto frontier in that the definition of dominated is
  /// relaxed (must not be < another term in all dimensions, but can be equal).
  void define_frontier(const UShortArray& mi_i,
		       UShortArraySet& combined_pareto);

  /// define a default definition for sparse_ind: 0 to num_terms-1
  void inflate(SizetSet& sparse_ind, size_t num_terms);

  /// overlay the passed expansion with the aggregate
  /// expansion{Coeffs,CoeffGrads} as managed by the multi_index_map
  void overlay_expansion(const SizetSet& sparse_ind_2,
			 const SizetArray& append_mi_map,
			 const RealVector& exp_coeffs_2,
			 const RealMatrix& exp_grads_2, int coeff_2);
  /// multiply current expansion ("a") with incoming expansion ("b")
  /// and store in product expansion ("c")
  void multiply_expansion(const SizetSet& sparse_ind_b,
			  const UShort2DArray& multi_index_b,
			  const RealVector& exp_coeffs_b,
			  const RealMatrix& exp_grads_b,
			  const UShort2DArray& multi_index_c);

  /**
   * \brief Define the set of options used in the cross validation grid search
   * 
   * \param opts (output) the options to be used in the grid search
   * \param M The number of rows of the vandermonde matrix
   * \param N The number of columns of the vandermonde matrix
   */
  void gridSearchFunction( RealMatrix &opts, int M, int N, 
			   int num_function_samples );

  void least_interpolation( RealMatrix &pts, 
			    RealMatrix &vals );

  void transform_least_interpolant( RealMatrix &L,
				    RealMatrix &U,
				    RealMatrix &H,
				    IntVector &p,
				    RealMatrix &vals );

  void least_factorization( RealMatrix &x,
			    UShort2DArray &basis_indices,
			    RealMatrix &l,
			    RealMatrix &u, 
			    RealMatrix &H, 
			    IntVector &p,
			    IntVector &k );

  void get_least_polynomial_coefficients( RealVector &v, IntVector &k,
					  UShort2DArray &basis_indices,
					  int num_dims, int num_pts,
					  RealMatrix &H );

  //
  //- Heading: Data
  //

  /// order of orthogonal best polynomial expansion found using cross validation
  IntVector bestApproxOrder;

  /// Stuct use to define the options of a compressed sensing solve
  CompressedSensingOptions CSOpts;

  /// store the fault info about the response data
  FaultInfo faultInfo;
  /// tracks use of sparse solvers, indicated the need to employ
  /// sparseIndices and sparseSobolIndexMap
  bool sparseSoln;

  /// tracks sparse terms within multiIndex and expansion{Coeffs,CoeffGrads}
  /// that are retained from an original candidate set
  /** a set is used to manage unique indices among expansionCoeff{s,Grads}.
      Sorting also simplifies covariance calculations, but care must be 
      exercised to retain synchronization with expansionCoeff{s,Grads}
      ordering when merging sparse multi-indices. */
  SizetSet sparseIndices;
  /// copy of sparseIndices stored in store_coefficients() for use in
  /// combine_coefficients()
  SizetSetArray storedSparseIndices;

  /// maps shared index from sobolIndexMap values to sparse index into
  /// sparse sobolIndices
  ULongULongMap sparseSobolIndexMap;

  /// PCE multi-index during the basis adaptation process.  Once complete,
  /// the shared multiIndex and sparseIndices are updated.
  UShort2DArray adaptedMultiIndex;
  /// sparse indices identifying recovered expansion coefficients within
  /// adaptedMultiIndex during the basis adaptation process.  Once complete,
  /// the shared multiIndex and sparseIndices are updated.
  SizetSet adaptedSparseIndices;
  /// the adapted multi-index that corresponds to the best solution identified
  /// Due to frontier/sparse restriction operations, bestAdaptedMultiIndex
  /// cannot be assumed to be a subset of adaptedMultiIndex.
  UShort2DArray bestAdaptedMultiIndex;
  /// the cross validation error reference point for adapting a CS
  /// candidate basis; it's state is reset for each response QoI
  Real cvErrorRef;
};


inline RegressOrthogPolyApproximation::
RegressOrthogPolyApproximation(const SharedBasisApproxData& shared_data):
  OrthogPolyApproximation(shared_data), sparseSoln(false)
{ }


inline RegressOrthogPolyApproximation::~RegressOrthogPolyApproximation()
{ }


inline void RegressOrthogPolyApproximation::build_linear_system( RealMatrix &A )
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  build_linear_system( A, data_rep->multiIndex );
}


inline void RegressOrthogPolyApproximation::
build_linear_system( RealMatrix &A, RealMatrix &B )
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  build_linear_system( A, B, data_rep->multiIndex );
}


inline void RegressOrthogPolyApproximation::
build_linear_system( RealMatrix &A, RealMatrix &B, RealMatrix &points )
{
  SharedRegressOrthogPolyApproxData* data_rep
    = (SharedRegressOrthogPolyApproxData*)sharedDataRep;
  build_linear_system( A, B, points, data_rep->multiIndex );
}


inline void RegressOrthogPolyApproximation::
inflate(SizetSet& sparse_ind, size_t num_terms)
{
  sparse_ind.clear();
  for (size_t i=0; i<num_terms; ++i)
    sparse_ind.insert(i);
}


inline ULongULongMap RegressOrthogPolyApproximation::
sparse_sobol_index_map() const
{ return sparseSobolIndexMap; }


inline size_t RegressOrthogPolyApproximation::expansion_terms() const
{
  return (sparseIndices.empty()) ? OrthogPolyApproximation::expansion_terms()
                                 : sparseIndices.size();
}

} // namespace Pecos

#endif
