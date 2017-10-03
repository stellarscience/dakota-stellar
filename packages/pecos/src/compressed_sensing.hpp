#ifndef COMPRESSED_SENSING_HPP
#define COMPRESSED_SENSING_HPP

#include "LinearAlgebra.hpp"
#include "pecos_global_defs.hpp"

namespace Pecos {

//! @name Interior-point optimization methods auxilary functions.
//@{ 

/**
 * \brief Compute the duality gap of the primal-dual interior point 
 * optimization algorithm. This is called by the function 
 * BasisPursuitPrimalDual().
 */
Real BP_surrogate_duality_gap( RealVector &primal_residual,
			       RealVector &fu1, RealVector &fu2, 
			       RealVector &lamu1, RealVector &lamu2, 
			       RealVector &AtV, Real mu, Real pdtol,
			       Real &tau, Real &slackness_norm );

/**
 * \brief Compute the central point on the central path computed by
 * BPDN_log_barrier_interior_point_method using newton's method.
 */
int BPDN_compute_central_point( RealMatrix &A, RealVector &b, RealVector &x, 
				RealVector &u, RealMatrix &AtA, 
				Real epsilon, Real &tau, 
				Real newton_tol, int newton_maxiter,
				Real conjugate_gradients_tol = 1e-8,
				int verbosity = 0 );

//@}

//! @name Interior-point l1 minimization methods.
//@{ 

/**
 * \brief Compute the Basis Pursuit ( BP ) sparse solution to 
 \f[ \arg \! \min \|x\|_1\quad\text{subject to}\quad Ax=b\f] by 
 * solving the associated linear program using the primal-dual interior-point
 * method.
 *
 * \param A ( M x N ) matrix of the linear system AX=B
 *
 * \param B ( M x num_rhs ) matrix of the linear system AX=B
 *
 * \param X ( N x num_rhs ) solutions to AX=B
 * 
 * \param primal_dual_tol controls the accuracy of the internal optimzation 
 * alogrithm
 *
 * \param conjugate_gradients_tol Specfies the tolerance of the 
 * the conjugate gradients method used to solve the newton step. If 
 * conjugate_gradients_tol < 0 cholesky factorization will be used.
 * 
 * \param verbosity turn print statements on and off. 
 */
void BP_primal_dual_interior_point_method( RealMatrix &A, RealVector &b, 
					   RealMatrix &result, 
					   Real primal_dual_tol,
					   Real conjugate_gradients_tol,
					   int verbosity );

/**
 * \brief Compute the Basis Pursuit Denoising sparse solution to 
 \f[ \arg \! \min \|x\|_1\quad\text{subject to}\quad\|Ax-b\|_2\le
 \varepsilon\f] using 
 * the log barrier method to solve the associated quadratic cone problem.
 *
 * \param A ( M x N ) matrix of the linear system AX=B
 *
 * \param B ( M x num_rhs ) matrix of the linear system AX=B
 *
 * \param X ( N x num_rhs ) solutions to AX=B
 *
 * \param epsilon defines the inequality constraint.
 *
 * \param log_barrier_tol controls the accuracy of the internal optimzation 
 * alogrithm   
 *
 * \param conjugate_gradients_tol Specfies the tolerance of the 
 * the conjugate gradients method used to solve the newton step. If 
 * conjugate_gradients_tol < 0 cholesky factorization will be used.
 *
 * \param verbosity turn print statements on and off. 
 * 0: off, 1: warnings on,  2: all print statements on.
 */
void BPDN_log_barrier_interior_point_method( RealMatrix &A, RealVector &b, 
					     RealMatrix &result, 
					     Real epsilon, 
					     Real log_barrier_tol,
					     Real conjugate_gradients_tol,
					     int verbosity );
//@}

//! @name Greedy methods.
//@{

/**
 * \brief Compute a greedy approximation to 
 \f[ \arg \! \min \|x\|_1\quad\text{subject to}\quad\|Ax-b\|_2\le
 \varepsilon\f]
 *
 * If  max_num_iterations > min(M,N) then if the algorithm completes
 * the last solution will be the least squares solution. However sometimes
 * due to numerical roundoff error the algorithm may exit early and 
 * the least squares solution will not be computed.
 *
 * \param A ( M x N ) matrix of the linear system Ax=b
 *
 * \param b ( M x 1 ) vector of the linear system Ax=b
 *
 * \param result_0 (output) On exit result_0 contains the solution
 * at each iteration of the algorithm.
 *
 * \param result_1 (output) Contains metrics about the solutions
 * contained in solutions. Specifically for each solution the residual
 * and number of non-zero-terms is stored.
 *
 * \param epsilon defines the inequality constraint.  If set to zero
 * the method will return all possible solutions with sparsity
 * increasing in unit increments from 1 to M
 *
 * \param max_num_iterations specify the maximum number of non-zeros
 * 
 * \param verbosity turn print statements on and off. 
 * 0: off, 1: warnings on,  2: all print statements on.
 *
 * \param ordering enforce a set of columns to be chosen first
 */
void orthogonal_matching_pursuit( RealMatrix &A, RealVector &b, 
				  RealMatrix &result_0,
				  RealMatrix &result_1,
				  Real epsilon ,
				  int max_num_iterations,
				  int verbosity,
				  IntVector &ordering );

/**
 * Orthogonal matching pursuit which uses a cholesky updating algorithm.
 * This allows one to quickly compute the inverse of A'A and thus
 * compute fast cross validation leave p out errors.
 */
void orthogonal_matching_pursuit_cholesky( RealMatrix &A, RealVector &b, 
					   RealMatrix &result_0,
					   RealMatrix &result_1,
		   std::vector< IntVector > &training_indices,
		   std::vector< IntVector > &validation_indices,
		   std::vector< std::vector < RealVector > > &cv_residuals,
					   RealVector &cv_scores,
					   Real epsilon ,
					   int max_num_iterations,
					   int verbosity );

//@}

//! @name Homotopy methods.
//@{

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
 */
void least_angle_regression( RealMatrix &A, 
			     RealVector &b, 
			     RealMatrix &result_0,
			     RealMatrix &result_1,
			     Real epsilon, 
			     int solver,
			     Real delta,
			     int max_num_iterations,
			     int verbosity );
//@}

int loo_step_lsq_cross_validation( RealMatrix &A, RealVector &b, 
				   IntVector &ordering, RealMatrix &result_0,
				   RealVector &result_1, int verbosity,
				   bool use_tpn );

void cosamp( RealMatrix &A, 
	     RealVector &b, 
	     RealMatrix &result_0,
	     RealMatrix &result_1,
	     int sparsity,
	     int max_iter,
	     int verbosity );

} // namespace Pecos

#endif //COMPRESSED_SENSING_HPP
