/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "least_angle_regression.hpp"
#include "math_tools.hpp"
#include "LinearSystemSolver.hpp" // include regressiontype enums

namespace Pecos {
namespace util {

  Real find_max_correlation( const RealMatrix &correlation,
                             const std::set<int> &inactive_indices,
                             const RealVector &column_norms,
                             bool normalize, bool non_negative )
  {
    Real max_abs_correlation( -1.0 );
    std::set<int>::iterator inactive_index_iter;
    Real abs_correlation = 0.;
    for ( inactive_index_iter = inactive_indices.begin();
          inactive_index_iter != inactive_indices.end();
          inactive_index_iter++ )
      {
        int index = *inactive_index_iter;
        abs_correlation = correlation(index,0);
        if (!non_negative) abs_correlation = std::abs( abs_correlation );
        if ( normalize ) abs_correlation /= column_norms[index];
        max_abs_correlation=std::max( abs_correlation, max_abs_correlation );
      }
    return max_abs_correlation;
  }

  void find_indices_to_add( const RealMatrix &correlation,
                            const std::set<int> &inactive_indices,
                            Real max_abs_correlation,
                            std::vector<int> &new_indices,
                            const RealVector &column_norms,
                            bool normalize )
  {
    Real eps = 3e-16;
    std::set<int>::iterator inactive_index_iter;
    for ( inactive_index_iter = inactive_indices.begin();
          inactive_index_iter != inactive_indices.end();
          inactive_index_iter++ ){
      int index = *inactive_index_iter;
      Real abs_correlation = std::abs( correlation(index,0) );
      if ( normalize ) abs_correlation /= column_norms[index];
      if ( std::abs( abs_correlation - max_abs_correlation ) < eps ){
        new_indices.push_back( index );
      }
    }
  }

  void resize_memory( RealMatrix &chol_factor, int num_covariates,
                      RealMatrix &solutions, RealMatrix &metrics,
                      int homotopy_iter, int memory_chunk_size,
                      bool store_history )
  {
    if ( chol_factor.numRows() <= num_covariates )
      chol_factor.reshape( chol_factor.numRows() + memory_chunk_size,
                           chol_factor.numCols() + memory_chunk_size );

    if ( (store_history) && (solutions.numCols() <= homotopy_iter) ) {
      solutions.reshape( solutions.numRows(),
                         solutions.numCols() + memory_chunk_size );
      metrics.reshape( metrics.numRows(), solutions.numCols() + memory_chunk_size);
    }
  }

  int update_cholesky_factor( const RealMatrix &Amatrix,
                              RealMatrix &A_sparse,
                              RealMatrix &chol_factor,
                              const std::vector<int> &new_indices,
                              int verbosity,
                              Real delta ){

    int num_rows = Amatrix.numRows();
    int prev_num_covariates = A_sparse.numCols();
    int num_new_covariates = 0;
    std::vector<int> new_independent_indices( new_indices.size() );
    // Update the Cholesky factorization of the grammian of the A sparse matrix
    //for ( int i=0; i<(int)new_indices.size(); i++ ){
    // only add 1st variable
    for ( int i=0; i<1; i++ ){
      int index_to_add = new_indices[i];
      RealMatrix A_col( Teuchos::View, Amatrix, num_rows, 1, 0,
                        index_to_add );

      int colinear =
        cholesky_factorization_update_insert_column( A_sparse, chol_factor,
                                                     A_col, prev_num_covariates+
                                                     num_new_covariates,
                                                     delta );
      // Add new covariate to the sparse A matrix
      A_sparse.reshape( num_rows, A_sparse.numCols() + 1);
      for ( int i=0; i<num_rows; i++ )
        A_sparse(i,prev_num_covariates) = Amatrix(i,index_to_add);

      if ( colinear ){
        if ( verbosity > 0  ){
          // Usually occurs when A is a vandermonde matrix and the inputs
          // have been standardized. The first column will be colinear
          // This condition should only occur when
          // num_covariates = max_num_covariates - 1
          std::stringstream msg;
          msg << "Exiting: attempted to add colinear vector\n";
          std::cout << msg.str();
        }
        return 1;
      }
      else{
        new_independent_indices[num_new_covariates] = index_to_add;
        num_new_covariates+=1;
      }
    }
    return 0;
  }


  int update_active_index_set( std::vector<int> &active_indices,
                               std::set<int> &inactive_indices,
                               const std::vector<int> &new_indices,
                               int homotopy_iter,
                               int verbosity ){
    // TODO: this will break down if num new indices > 1. Throw error for now
    // but remove when fixed
    //if ( (int)new_indices.size() > 1 )
    //  throw( std::runtime_error("tried to add more than one covariate.") );

    std::set<int>::iterator inactive_index_iter;
    int index_to_add = new_indices[0];
    active_indices.push_back( index_to_add );
    inactive_index_iter = inactive_indices.find( index_to_add );
    inactive_indices.erase( inactive_index_iter );
    return index_to_add;
  }

  void downdate_cholesky_factor( RealMatrix &chol_factor,
                                 const std::vector<int> &active_indices,
                                 int sparse_index_to_drop,
                                 RealMatrix &A_sparse){
    int num_covariates = A_sparse.numCols();
    cholesky_factorization_update_delete_column( chol_factor,
                                                 sparse_index_to_drop,
                                                 num_covariates );
    delete_column( sparse_index_to_drop, A_sparse );
  }

  void downdate_active_index_set( std::vector<int> &active_indices,
                                  std::set<int> &inactive_indices,
                                  int sparse_index_to_drop,
                                  int homotopy_iter,
                                  int num_covariates, int verbosity ){
    int index_to_drop = active_indices[sparse_index_to_drop];
    inactive_indices.insert( index_to_drop );
    std::vector<int>::iterator it;
    it =  active_indices.begin() + sparse_index_to_drop;
    active_indices.erase( active_indices.begin() + sparse_index_to_drop );
  }

  void compute_equidistant_vector( const RealMatrix &chol_factor,
                                   const RealMatrix &correlation,
                                   const std::vector<int> &active_indices,
                                   const RealMatrix &Amatrix,
                                   const RealMatrix &A_sparse,
                                   RealMatrix &equiangular_vec,
                                   RealMatrix &angles,
                                   RealMatrix &w_sparse,
                                   Real &normalisation_factor,
                                   bool non_negative ){
    Teuchos::BLAS<int, Real> blas;

    // Get the signs of the correlations
    int num_covariates = active_indices.size();
    RealMatrix signs_sparse( num_covariates, 1, false );
    for ( int sparse_index = 0; sparse_index < num_covariates; sparse_index++ ){
      if (non_negative)
        signs_sparse(sparse_index,0) = 1;
      else
        signs_sparse(sparse_index,0) =
          sgn( correlation(active_indices[sparse_index],0) );
    }

    // normalisation_factor = 1 / sqrt( s'*inv(A'*A)*s )
    // inv(A'A)*s  = > solve A'*A*z \ s => chol_factor'chol_factor \ s
    // so solve two problems w = chol_factor' \ s then z = chol_factor \ w
    RealMatrix w, z;
    RealMatrix prev_chol_factor( Teuchos::View, chol_factor, num_covariates,
                                 num_covariates, 0, 0 );
    // warning changed substitution solve to not use defaults,
    // teuchos enums passed here may be incorrect. I have not tested
    substitution_solve( prev_chol_factor, signs_sparse, w, Teuchos::TRANS,
                        Teuchos::UPPER_TRI, Teuchos::NON_UNIT_DIAG );
    substitution_solve( prev_chol_factor, w, z, Teuchos::NO_TRANS,
                        Teuchos::UPPER_TRI, Teuchos::NON_UNIT_DIAG );
    normalisation_factor = 1.0 / std::sqrt( blas.DOT( num_covariates,
                                                      signs_sparse[0],
                                                      1, z[0], 1 ) );

    // Compute unit vector making equal angles, less than 90 degrees ,
    // with the columns of A_sparse
    w_sparse.shapeUninitialized( z.numRows(), z.numCols() );
    w_sparse.assign( z );
    w_sparse *= normalisation_factor;

    // Compute the equiangular vector
    // Assumes that equiangular_vec has the right size
    equiangular_vec.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS,
                              1.0, A_sparse, w_sparse, 0.0 );

    // Compute angles betwen A_j and equiangular_vec
    // Assume that angles has the right size
    angles.multiply( Teuchos::TRANS, Teuchos::NO_TRANS,
                     1.0, Amatrix, equiangular_vec, 0.0 );

  }

  void find_indices_to_drop( const RealVector &solution,
                             const std::vector<int> &active_indices,
                             const RealMatrix &w_sparse,
                             Real &gamma_tilde,
                             int &sparse_index_to_drop ){
    // Find the first occurance of a solution element changing sign
    // (3.5) (Efron 2004)
    int num_covariates = active_indices.size();
    gamma_tilde = std::numeric_limits<Real>::max();
    sparse_index_to_drop = -1;
    for ( int n = 0; n < num_covariates; n++ ){
      Real gamma = -solution[active_indices[n]] / (w_sparse(n,0));
      if ( ( gamma > 0 ) && ( gamma < gamma_tilde ) ){
        sparse_index_to_drop = n;
        gamma_tilde = gamma;
      }
    }
  }


  Real compute_step_size( Real max_abs_correlation,
                          const std::set<int> &inactive_indices,
                          const RealMatrix &correlation,
                          const RealMatrix &angles,
                          int num_covariates, int N,
                          Real normalisation_factor,
                          bool non_negative){
    // Compute the lars step size. (2.13) (Efron 2004)
    // When the active_indices contains all covariates, gamma_hat
    // is not defined. By convention in this situation the algorithm takes
    // gamma_hat = corelationMax, making the solution x equal to the
    // ordinary least squares estimate for the full set of N covariates.
    Real gamma_hat = max_abs_correlation / normalisation_factor;
    Real eps = 2*std::numeric_limits<double>::epsilon();
    if ( num_covariates < N ){
      std::set<int>::iterator inactive_index_iter;
      for ( inactive_index_iter = inactive_indices.begin();
            inactive_index_iter != inactive_indices.end();
            inactive_index_iter++ ){
        int n = *inactive_index_iter;
        Real gamma_hat1 = ( max_abs_correlation - correlation(n,0) ) /
          ( normalisation_factor - angles(n,0) + eps);

        gamma_hat = ( ( gamma_hat1 < gamma_hat ) && ( gamma_hat1 > 0 ) )
          ? gamma_hat1 : gamma_hat;

        if (!non_negative){
          Real gamma_hat2 =  ( max_abs_correlation + correlation(n,0) ) /
            ( normalisation_factor + angles(n,0) + eps);
          gamma_hat = ( ( gamma_hat2 < gamma_hat ) && ( gamma_hat2 > 0 ) )
            ? gamma_hat2 : gamma_hat;
        }
      }
    }
    return gamma_hat;
  }

  bool check_exit_conditions( int homotopy_iter,
                              Real residual_norm, int num_covariates,
                              Real residual_tol, int max_num_covariates,
                              int max_num_iters, int verbosity,
                              bool remove_index, Real prev_residual_norm ){
    bool done = false;

    if ( residual_norm <= residual_tol ){
      if ( verbosity > 1 )
        std::cout << "\nExiting: residual norm lower than tolerance\n";
      done = true;
    }

    if ( homotopy_iter == max_num_iters ){
      if ( verbosity > 1 )
        std::cout << "\nExiting: maximum number of iterations reached\n";
      done = true;
    }

    if ( (!remove_index) && ( num_covariates >= max_num_covariates ) ){
      if ( verbosity > 1 )
        std::cout << "\nExiting: maximum number of covariates reached\n";
      done = true;
    }

    if (prev_residual_norm < residual_norm ){
      if ( verbosity > 1 )
        std::cout << "\nExiting: residual started increasing\n";
      done = true;
    }

    if ( done && num_covariates==0)
      throw(std::runtime_error("Residual tolerance was too large and no basis functions were chosen"));

    return done;

  }

  /*void add_first_covariate( const RealMatrix &Amatrix, RealMatrix &A_sparse,
    const RealMatrix &correlation,
    std::set<int> &inactive_indices,
    std::vector<int> &active_indices, RealMatrix &metrics,
    RealMatrix &chol_factor, RealVector &column_norms,
    bool normalize_choice, Real delta,
    int homotopy_iter, int verbosity ){
    Real max_abs_correlation = find_max_correlation( correlation,
    inactive_indices,
    column_norms,
    normalize_choice );

    std::vector<int> new_indices;
    find_indices_to_add( correlation, inactive_indices, max_abs_correlation,
    new_indices, column_norms, normalize_choice );

    update_cholesky_factor( Amatrix, A_sparse, chol_factor,
    new_indices, verbosity,
    delta );

    update_active_index_set( active_indices, inactive_indices, new_indices,
    metrics, homotopy_iter, verbosity );
    }*/

  void least_angle_regression( const RealMatrix &Amatrix,
                               const RealVector &b,
                               RealMatrix &result_0,
                               RealMatrix &result_1,
                               Real residual_tol,
                               int solver,
                               Real delta,
                               int max_num_iters,
                               int max_num_covariates,
                               int verbosity,
                               bool normalize_choice,
                               bool non_negative, bool store_history,
                               int memory_chunk_size)
  {
    if (delta>0)
      throw(std::runtime_error("delta > 0 not currently supported"));
  
    RealVector column_norms;
    get_column_norms( Amatrix, column_norms );

    Teuchos::BLAS<int, Real> blas;

    int M( Amatrix.numRows() ), N( Amatrix.numCols() );

    if ( M != b.length() )
      throw( std::runtime_error("least_angle_regresssion: Matrix and RHS are inconistent") );

    // Apply defaults
    if ( max_num_covariates == 0 ){
      max_num_covariates = std::min( M, N );
      //if ( delta > std::numeric_limits<Real>::epsilon() )
      //  max_num_covariates = std::min( N, max_num_iters );
    }else{
      max_num_covariates = std::min( max_num_covariates, std::min( M, N ) );
    }

    // Lasso will usually require more iterations than max covariates
    // as variables are added and deleted during the algorithm. However
    // to allow cross validation I restrict the number of iterations
    // max_num_iters to max_num_covariates. If I do not do this then
    // LASSO can use different iteration numbers for different training
    // data but the same cross validation model options

    // I do not think the above comment is valid any more, now that I am
    // doing cross validation on residual norm. Need to double check this
    // though

    memory_chunk_size = std::min( memory_chunk_size, std::min(M, 500));

    // Memory for active indices
    std::vector<int>::iterator vec_iter;
    std::vector<int> active_indices;

    // Store inactive indices
    std::set<int> inactive_indices;
    for ( int n = 0;  n < N; n++) inactive_indices.insert( n );

    // Matrix to store the full rank matrix associated with x_sparse
    RealMatrix A_sparse;

    // Matrix to store cholesky factorization of A'A and initialize to zero
    int size = std::min( N, memory_chunk_size );
    RealMatrix chol_factor( size, size );
    //RealMatrix chol_factor( N, N );

    // Matrix to store the equidistant anglse
    RealMatrix angles( N, 1, false );

    // Matrix to store the equiangular vector
    RealMatrix equiangular_vec( M, 1, false );

    int initial_num_stored_coeff = ((store_history) ? memory_chunk_size : 1);
    // Initialise all entries of x to zero
    result_0.shape( N, initial_num_stored_coeff );
    // Allocate memory to store solution metrics
    result_1.shapeUninitialized( 2, initial_num_stored_coeff );


    // Compute residual
    RealVector residual( Teuchos::Copy, b.values(), b.length() );
    Real residual_norm = residual.normFrobenius();

    // Compute correlation of columns with residual
    RealMatrix Atb( N, 1, false );
    Atb.multiply( Teuchos::TRANS, Teuchos::NO_TRANS,
                  1.0, Amatrix, residual, 0.0 );
    RealMatrix correlation( Teuchos::Copy, Atb, N, 1 );

    // The index that violates the lasso sign condition
    int index_to_drop = -1;

    if ( verbosity > 1 )
      {
        if ( solver == LASSO_REGRESSION )
          std::cout << "LASSO ( delta = " << delta << " )\n";
        else
          std::cout << "LARS ( delta = " << delta << " )\n";
        std::printf("%-4s %-5s %-7s %-8s %-21s %-21s %-21s\n",
                    "Iter","Added","Dropped","Sparsity","Max Correlation",
                    "Residual Norm","l1 Norm x");
      }

    bool done = false;
    bool drop_covariate = false;

    // The lasso approximation of b. Initialize to zero
    RealMatrix b_hat( M, 1 );
    int homotopy_iter( 0 );
    Real prev_residual_norm = std::numeric_limits<double>::max();
    RealVector prev_solution(N,false); prev_solution=0.0;
    while ( !done )
      {
        Real max_abs_correlation = find_max_correlation( correlation,
                                                         inactive_indices,
                                                         column_norms,
                                                         normalize_choice,
                                                         non_negative );
        //
        if (max_abs_correlation <=std::numeric_limits<double>::epsilon()){
          // if non-negative is true then we stop when we the maximum correlation
          // is negative. We also stop if correlation drops below machine precision
          if ( verbosity > 1 )
            std::cout << "max_correlation became to small" << std::endl;
          break;
        }

        std::vector<int> new_indices;
        find_indices_to_add( correlation, inactive_indices, max_abs_correlation,
                             new_indices, column_norms, normalize_choice );

        done = check_exit_conditions( homotopy_iter,
                                      residual_norm, active_indices.size(),
                                      residual_tol, max_num_covariates,
                                      max_num_iters, verbosity,
                                      drop_covariate, prev_residual_norm );
        prev_residual_norm = residual_norm;
        if (done){
          break;
        }

        resize_memory( chol_factor, active_indices.size(), result_0, result_1,
                       homotopy_iter, memory_chunk_size, store_history );

        // where to store solution and metrics in result_0 and result_1
        // respectively. Must be done after resize memory or prev_solution
        // will point to wrong location in memory
        int storage_index = (store_history) ? homotopy_iter : 0;

        if ( !drop_covariate ){
          int colinear = update_cholesky_factor( Amatrix, A_sparse, chol_factor,
                                                 new_indices, verbosity, delta );

          if (colinear) break;

          int index_to_add =
            update_active_index_set( active_indices, inactive_indices, new_indices,
                                     homotopy_iter, verbosity  );
          // store which variable was added to the active index set
          result_1(1,storage_index) = index_to_add;
        }

        RealMatrix w_sparse;
        Real normalisation_factor = -1;
        compute_equidistant_vector( chol_factor, correlation, active_indices,
                                    Amatrix, A_sparse, equiangular_vec,
                                    angles, w_sparse,
                                    normalisation_factor,
                                    non_negative);


        if ((store_history) && (homotopy_iter>0)){
          RealVector prev_solution_view(Teuchos::View, result_0[homotopy_iter-1], N);
          prev_solution=prev_solution_view;
        }else{
          // When a covariate is dropped result will be close to zero
          // set it to zero here. This minimizes numerical errors.
          // Also if I do not do this results will differ
          // depending on value of store_history.
          if (drop_covariate)
            result_0(index_to_drop,0)=0.;

          for (int i=0; i<N; ++i)
            prev_solution[i]=result_0(i,0);
        }


        Real gamma_tilde = std::numeric_limits<Real>::max();
        int violating_sparse_index = -1;
        if ( ( ( solver == LASSO_REGRESSION ) || (non_negative) ) && ( homotopy_iter > 0 ) ){
          find_indices_to_drop( prev_solution, active_indices, w_sparse,
                                gamma_tilde, violating_sparse_index );
        }

        Real gamma_hat = compute_step_size( max_abs_correlation, inactive_indices,
                                            correlation, angles,
                                            active_indices.size(), N,
                                            normalisation_factor,
                                            non_negative);

        Real gamma_min = std::min( gamma_hat, gamma_tilde );

        // Update the solution.
        for ( int n = 0; n < (int)active_indices.size(); n++ ){
          result_0(active_indices[n],storage_index) =
            prev_solution[active_indices[n]]+gamma_min*w_sparse(n,0);
        }

        // Update the rhs and residual.
        // b_new = b_old + gamma_min * equiangular_vec (2.12) (Efron 2004)
        // => r_new = b - b_new = b - ( b_old + gamma_min * equiangular_vec )
        //          = r_old - gamma_min * equiangular_vec
        for ( int m = 0; m < M; m++ ){
          b_hat(m,0) += gamma_min * equiangular_vec(m,0);
          residual[m] -= gamma_min * equiangular_vec(m,0);
        }

        // Update the correlation (2.15) (Efron 2004)
        for ( int n = 0; n < N; n++ )
          correlation(n,0) -= ( gamma_min * angles(n,0) );

        residual_norm = residual.normFrobenius();
        //result_1(0,homotopy_iter) = residual_norm;
        result_1(0,storage_index) = residual_norm;

        if ( verbosity > 1 ){
          RealVector x( Teuchos::View, result_0[storage_index], N );
          if ( !drop_covariate ){
            std::printf("%-4d %-13d %-8d ", homotopy_iter,
                        new_indices[0], (int)active_indices.size());
            std::printf("%-21.15e %-21.15e %-21.15e\n", max_abs_correlation,
                   residual_norm, x.normOne());
          }else{
            std::printf( "%-10d %-7d %-8d ", homotopy_iter,
                         index_to_drop, (int)active_indices.size());
            std::printf("%-21.15e %-21.15e %-21.15e\n",
                        std::abs(correlation(index_to_drop,0)),
                        residual_norm, x.normOne());
          }
        }

        drop_covariate = ( gamma_hat >= gamma_tilde );

        if ( ( drop_covariate ) && ( solver == LASSO_REGRESSION || (non_negative))){
          downdate_cholesky_factor( chol_factor, active_indices,
                                    violating_sparse_index,
                                    A_sparse );

          index_to_drop = active_indices[violating_sparse_index];
          downdate_active_index_set( active_indices, inactive_indices,
                                     violating_sparse_index,
                                     homotopy_iter, active_indices.size(),
                                     verbosity );
          // Store which variable was removed from active index set
          // minus sign indicates variable was removed.
          result_1(1,storage_index) = -index_to_drop;


          // update residual and correlation
          //RealVector x( Teuchos::View, result_0[homotopy_iter], N );
          RealVector x( Teuchos::View, result_0[storage_index], N );
          correlation(index_to_drop,0) = 0.;
          for ( int m = 0; m < M; m++ ){
            residual[m] -= Amatrix(m,index_to_drop) *
              x[index_to_drop];
            correlation(index_to_drop,0) +=
              Amatrix(m,index_to_drop) * residual[m];
          }
        }//else{}
        homotopy_iter++;
      }

    // remove unused memory
    if (store_history){
      result_0.reshape( N, homotopy_iter );
      result_1.reshape( result_1.numRows(), homotopy_iter );
    }
  };

}  // namespace util
}  // namespace Pecos
