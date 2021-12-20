/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_ORTHOGONAL_MATCHING_PURSUIT_HPP
#define PECOS_UTIL_ORTHOGONAL_MATCHING_PURSUIT_HPP

#include "linear_algebra.hpp"
#include "math_tools.hpp"

// BMA TODO: Move this implementation to a .cpp file

namespace Pecos {
namespace util {

  template<typename ScalarType>
  class MatrixVectorMultiplicationOperator{
  private:
    typedef typename Teuchos::SerialDenseMatrix<int,ScalarType> MatrixType;
    typedef typename Teuchos::SerialDenseVector<int,ScalarType> VectorType;

    MatrixType matrix_;

  public:
    MatrixVectorMultiplicationOperator(){};
    ~MatrixVectorMultiplicationOperator(){};

    void set_matrix( const MatrixType &matrix ){
      matrix_ = matrix;
    };

    /** Multiply \c A and \c vector and add them to \e this; 
        \e result = \c beta * \e result + \c alpha*matrix*vector
      
        \param mode - Use the transpose of \c A if mode = 1, 
        else don't use the transpose if mode=0. If A is complex apply 
        the conjugate transpose.
        \param alpha - The scaling factor for \c matrix * \c vector.
        \param vector - SerialDenseMatrix must have numCols()==1
        \param beta - The scaling factor for \e result.
    **/
    virtual void apply( const VectorType &vector, int mode,
                        ScalarType alpha, ScalarType beta,
                        VectorType& result ){

      Teuchos::ETransp trans = Teuchos::NO_TRANS;
      if ( Teuchos::ScalarTraits<ScalarType>::isComplex && mode)
        // if matrix is complex apply conjugate transpose
        trans = Teuchos::CONJ_TRANS;
      else if ( mode )
        //mode==1 => apply transpose
        // if matrix is real apply transpose
        trans = Teuchos::TRANS;

      GEMV<int,ScalarType>(trans, true, alpha, matrix_, vector, beta, result);
    };
  };

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
   *
   * \param store_history specify whether to store the solution at every step
   * of the algorithm or just final step.
   *
   * \param memory_chunk_size OMP internally allocates memory for solutions in 
   * chunks to avoid allocating unnecessary memory for sparsity << A.numCols()
   */

  template<typename ScalarType>
  void orthogonal_matching_pursuit(const Teuchos::SerialDenseMatrix<int,ScalarType> &A, 
                                   const Teuchos::SerialDenseVector<int,ScalarType> &b, 
                                   Teuchos::SerialDenseMatrix<int,ScalarType> &result_0,
                                   RealMatrix &result_1,
                                   Real epsilon, 
                                   int max_nnz,
                                   int verbosity,
                                   IntVector &ordering,
                                   bool normalise_choice,
                                   bool store_history,
                                   int memory_chunk_size) {

    typedef typename Teuchos::SerialDenseVector<int,ScalarType> VectorType;
    typedef typename Teuchos::SerialDenseMatrix<int,ScalarType> MatrixType;

    MatrixVectorMultiplicationOperator<ScalarType> mat_vec_op;
    mat_vec_op.set_matrix(A);

    if ( max_nnz < 1 )
      throw( std::runtime_error("OMP() Ensure nnz>0") );
  
    RealVector column_norms;
    get_column_norms( A, column_norms );
  
    Teuchos::BLAS<int, ScalarType> blas;

    int M( A.numRows() ), N( A.numCols() );

    if ( M != b.length() )
      throw( std::runtime_error("OMP() A and rhs are inconsistent") );

    // Apply default number of iterations
    if ( max_nnz == 0 )
      max_nnz = std::min( M, N );

    // Determine the maximum number of iterations
    int max_num_indices( std::min( M, max_nnz ) );
    max_num_indices = std::min( N, max_num_indices );

    // VectorType to store non-zero solution entries
    VectorType x_sparse;	

    memory_chunk_size = std::min(memory_chunk_size,std::min(M,N));
    int initial_N = (store_history) ? memory_chunk_size : 1;

    // Initialise entries of all solutions to zero
    result_0.shape( N, initial_N );

    // Allocate memory to store solution metrics
    result_1.shapeUninitialized( 2, initial_N );

    // MatrixType to store Q and R factorization
    MatrixType Q(M, memory_chunk_size), R(memory_chunk_size, memory_chunk_size);

    // Compute residual
    VectorType residual( Teuchos::Copy, b.values(), b.length() );

    // Compute correlation of columns with residual
    VectorType Atb( N, false );
    mat_vec_op.apply(residual, 1, 1.0, 0.0, Atb );
    VectorType correlation( Atb );  
  
    // Compute norm of residual
    Real residual_norm(b.normFrobenius());

    // MatrixType to store the full rank matrix associated with x_sparse
    MatrixType A_sparse_memory( M, memory_chunk_size, false );
    VectorType Atb_sparse_memory( N, false );

    if ( verbosity > 1 ){
      std::cout << "Orthogonal Matching Pursuit\n";
      std::cout << "Store history: " << store_history << "\n";
      std::cout << "Residual tolerance: " << epsilon << "\n";
      std::cout << "A: (" << A.numRows() << "," << A.numCols() << ")\n";
      std::printf( "Iter\tAdded\tResidual\tCorrelation\tl1 norm of x\n" );
    }

    int num_active_indices( 0 );
    IntVector active_index_set( max_num_indices );
    bool done = false;
    while ( !done ){
      int active_index = -1;
      Real max_abs_correlation;
      if ( num_active_indices >= ordering.length() ){
        // Find the column that has the largest inner product with
        // the residuals
        max_abs_correlation = -std::numeric_limits<Real>::max();
        for ( int n = 0; n < N; n++ ){
          Real abs_correlation_n = std::abs(correlation[n]);
          if ( normalise_choice )
            abs_correlation_n /= column_norms[n];
          if ( abs_correlation_n > max_abs_correlation ){
            max_abs_correlation = abs_correlation_n;
            active_index = n;
          }
        }
        // For complex variables IMAX computes abs(a+bi) as |a|+|b|
        // but need sqrt(a^2+b^2)
        // Warning IAMX returns the index of the element with the 
        // largest magnitude but IAMAX assumes indexing 
        // 1,..,N not 0,...,N-1
        //active_index = blas.IAMAX( correlation.length, 
        //			   correlation.values(), 1 ) - 1;
        //max_abs_correlation = std::abs( correlation[active_index] );
      }else{
        active_index = ordering[num_active_indices];
        //VectorType active_col( Teuchos::View, A[active_index], A.numRows() );
        VectorType active_col =
          Teuchos::getCol(Teuchos::View, const_cast<RealMatrix &>(A), active_index);

        max_abs_correlation = std::abs( active_col.dot( residual ) );
        /// Warning not sure if dot will do correct thing for complex variables. I think it does
      }
    
      // todo define active_index_set as std::set and use find function
      for ( int i = 0; i < num_active_indices; i++ ){
        if ( active_index_set[i] == active_index ){
          if ( verbosity > 1 ){
            std::cout << "Exiting: New active index " <<  active_index
                      << " has already been added. "
                      << "This has likely occured because all correlations "
                      << "are roughly the same size. "
                      << "This means problem has been solved to roughly "
                      << "machine precision. Check this before continuing. "
                      << "Correlation of active_index is: "
                      << max_abs_correlation << std::endl;
          }
          done = true;
        }
      }
      if ( done ) break;

      if ( Q.numCols() <= num_active_indices ){
        Q.reshape( Q.numRows(), Q.numCols() + memory_chunk_size );
        R.reshape( R.numRows() + memory_chunk_size, 
                   R.numCols() + memory_chunk_size );
        A_sparse_memory.reshape( M, A_sparse_memory.numCols() + 
                                 memory_chunk_size);
        if (store_history){
          result_0.reshape(N, result_0.numCols() + memory_chunk_size);
          result_1.reshape(2, result_1.numCols() + memory_chunk_size);
        }
      }
    
      // Update the QR factorisation.	 
      //VectorType A_col( Teuchos::View, A[active_index], M );
      VectorType A_col =
        Teuchos::getCol(Teuchos::View, const_cast<RealMatrix &>(A), active_index);
      int colinear = qr_factorization_update_insert_column( Q, R, A_col, 
                                                            num_active_indices );
       
      if ( !colinear ){
        active_index_set[num_active_indices] = active_index;
        VectorType Atb_sparse( Teuchos::View, Atb_sparse_memory.values(),
                               num_active_indices+1);
      
        int index( active_index_set[num_active_indices] );
        Atb_sparse[num_active_indices] = Atb[index];

        //Solve R'z = A'b via back substitution
        MatrixType z;
        MatrixType R_new( Teuchos::View, R, num_active_indices+1, 
                          num_active_indices+1, 0, 0 );
      
        substitution_solve<ScalarType>( R_new, Atb_sparse, z, Teuchos::CONJ_TRANS , Teuchos::UPPER_TRI, Teuchos::NON_UNIT_DIAG);
      
        //Solve Rx = z via back substitution to obtain signal
        substitution_solve( R_new, z, x_sparse, Teuchos::NO_TRANS, Teuchos::UPPER_TRI, Teuchos::NON_UNIT_DIAG );

        MatrixType A_sparse( Teuchos::View, A_sparse_memory, 
                             M, num_active_indices+1, 0, 0 );
        for ( int m = 0; m < M; m++ )
          A_sparse(m,num_active_indices) = A(m,active_index);
        VectorType residual( Teuchos::Copy, b.values(), b.length() );
        // residual = b - A_sparse * x_sparse: O(Mk) k < M
        GEMV( Teuchos::NO_TRANS, true,
              -Teuchos::ScalarTraits<ScalarType>::one(), A_sparse, x_sparse,
              Teuchos::ScalarTraits<ScalarType>::one(), residual );
        // Only compute correlation if needed. 
        // Note: Correlation is effectively the Monte Carlo quadrature 
        // pseudo spectral estimate of the coefficient size when 
        // approximating the residual 
        if ( num_active_indices >= ordering.length()-1 ){
          // -1 because ordering has specified its last enforced index.
          // correlation = A' * residual: O(MN) N >= M
          mat_vec_op.apply(residual,1,1.,0.,correlation);
        }

        residual_norm = residual.normFrobenius();
      
        num_active_indices++;   
      }else{
        //New column was co linear so ignore
        if ( verbosity > 0 ){
          std::stringstream msg;
          msg << "No variable added. Column " << active_index;
          msg << " was colinear " << std::endl;
          std::cout << msg.str();
        }
      }
      int storage_index = (store_history) ? num_active_indices-1 : 0;
      for ( int n = 0; n < num_active_indices; n++ ){
        result_0(active_index_set[n],storage_index) = x_sparse[n];
      }
      result_1(0,storage_index) = residual_norm;
      result_1(1,storage_index) = active_index;     
    
      if ( verbosity > 1 )
        std::printf( "%d\t%d\t%1.5e\t%1.5e\t%1.5e\n", num_active_indices, 
                     active_index, residual_norm, max_abs_correlation,
                     x_sparse.normOne() );
    
      if ( residual_norm <= epsilon ){
        if ( verbosity > 1 )
          std::cout << "Exiting: residual norm lower than tolerance\n";
        done = true;
      }
      
      if ( num_active_indices >= max_num_indices ){
        if ( verbosity > 1 )
          std::cout <<  "Exiting: maximum number of covariates reached\n";
        done = true;
      }

      if ( colinear ){
        // Usually occurs when A is a vandermonde matrix and the inputs 
        // have been standardized. The first column will be colinear
        // This condition should only occur when 
        // num_covariates = max_num_covariates - 1
        if ( verbosity > 1 )
          std::cout << "Exiting: attempted to add colinear vector\n";
        done = true;
      }
    }

    // remove unused memory
    if (store_history){
      result_0.reshape( N, num_active_indices );
      result_1.reshape( 2, num_active_indices );
    }
  }


}  // namespace util
}  // namespace Pecos

#endif  // include guard
