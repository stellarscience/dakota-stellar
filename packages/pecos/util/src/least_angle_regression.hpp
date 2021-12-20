/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_LEAST_ANGLE_REGRESSION_HPP
#define PECOS_UTIL_LEAST_ANGLE_REGRESSION_HPP

#include "teuchos_data_types.hpp"
#include <set>

namespace Pecos {
namespace util {

Real find_max_correlation( const RealMatrix &correlation,
			   const std::set<int> &inactive_indices,
			   const RealVector &column_norms,
			   bool normalize, bool non_negative);

void find_indices_to_add( const RealMatrix &correlation,
			  const std::set<int> &inactive_indices,
			  Real max_abs_correlation,
			  std::vector<int> &new_indices,
			  const RealVector &column_norms,
			  bool normalize );

void resize_memory( RealMatrix &chol_factor, int num_covariates, 
		    RealMatrix &solutions, RealMatrix &metrics, 
		    int homotopy_iter, int memory_chunk_size,
		    bool store_history );

int update_cholesky_factor( const RealMatrix &Amatrix, 
			    RealMatrix &A_sparse, 
			    RealMatrix &chol_factor,
			    const std::vector<int> &new_indices,
			    int verbosity, Real delta );

int update_active_index_set( std::vector<int> &active_indices,
			     std::set<int> &inactive_indices,
			     const std::vector<int> &new_indices,
			     int homotopy_iter,
			     int verbosity );

void downdate_cholesky_factor( RealMatrix &chol_factor, 
			       const std::vector<int> &active_indices,
			       int sparse_index_to_drop,
			       RealMatrix &A_sparse);

void downdate_active_index_set( std::vector<int> &active_indices,
				std::set<int> &inactive_indices,
				int sparse_index_to_drop,
				int homotopy_iter,
				int num_covariates, int verbosity );

void compute_equidistant_vector( const RealMatrix &chol_factor,
				 const RealMatrix &correlation,
				 const std::vector<int> &active_indices,
				 const RealMatrix &Amatrix,
				 const RealMatrix &A_sparse,
				 RealMatrix &equiangular_vec,
				 RealMatrix &angles,
				 RealMatrix &w_sparse,
				 Real &normalisation_factor,
				 bool non_negative);

void find_indices_to_drop( const RealVector &solution, 
			   const std::vector<int> &active_indices,
			   const RealMatrix &w_sparse,
			   Real &gamma_tilde,
			   int &sparse_index_to_drop );


Real compute_step_size( Real max_abs_correlation,
			const std::set<int> &inactive_indices,
			const RealMatrix &correlation,
			const RealMatrix &angles,
			int num_covariates, int N,
			Real normalisation_factor,
			bool non_negative);

bool check_exit_conditions( int homotopy_iter,
			    Real residual_norm, int num_covariates,
			    Real residual_tol, int max_num_covariates, 
			    int max_num_iter, int verbosity, 
			    bool remove_index, Real prev_residual_norm );

/**
 * \brief Compute the least angle regression ( and the lasso modification ) 
 * solution to \f[ \arg \! \min \|Ax-b\|^2_2 \quad\text{subject to}\quad
 \|x\|_1\le \tau \f]
 *
 * \param A ( M x N ) matrix of the linear system Ax=b
 *
 * \param b ( M x 1 ) vector of the linear system Ax=b
 *
 * \param epsilon controls the exit condition. If set to zero
 * the method will return the set of solutions for all \f$ \tau \f$
 * with sparsity ranging from 1 to M. If epsilon is non zero
 * then the method will terminate when the residual 
 * \f$ \|Ax-b\|_2<\varepsilon \f$
 *
 * \param results_0 (output) On exit result_0 contains the solution
 * at each iteration of the homotopy algorithm.
 *
 * \param result_1 (output) Contains metrics about the solutions
 * contained in solutions. Specifically for each solution the residual
 * and number of non-zero-terms is stored.
 *
 * \param solver specify whether to compute the least angle regression (LARS)
 * or the lasso solution
 *
 * \param max_num_iterations specify the maximum number of iterations.
 * For LARS this is also the maximum number of non-zeros. For lasso
 * the numbre of non-zeros will likely be less than the number of iterations
 *
 * \param verbosity turn print statements on and off. 
 * 0: off, 1: warnings on,  2: all print statements on.
 *
 * \param non_negative specify whether to enforce non-negative solutions
 *
 * \param store_history specify whether to store the solution at every step
 * of the algorithm or just final step.

 * \param memory_chunk_size LAR internally allocates memory for solutions in 
 * chunks to avoid allocating unnecessary memory for sparsity << A.numCols()
 */
void least_angle_regression( const RealMatrix &A, 
			     const RealVector &b, 
			     RealMatrix &result_0,
			     RealMatrix &result_1,
			     Real epsilon, 
			     int solver,
			     Real delta,
			     int max_num_iterations,
			     int max_num_covariates,
			     int verbosity,
			     bool normalise_choice,
			     bool non_negative,
			     bool store_history,
			     int memory_chunk_size);

}  // namespace util
}  // namespace Pecos

#endif  // include guard
