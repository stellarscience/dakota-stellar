/**
 * \file LinearAlgebra.hpp
 * \author John D. Jakeman
 * \date 2 January 2012
 * \brief Functions used to solve systems of linear equations.
 */

#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_LAPACK_wrappers.hpp"

#include "pecos_data_types.hpp"

namespace Pecos{

// typedef Teuchos::SerialDenseMatrix<int,int> IntMatrix;
// typedef Teuchos::SerialDenseMatrix<int,Real> RealMatrix;
// typedef Teuchos::SerialDenseVector<int,int> IntVector;
// typedef Teuchos::SerialDenseVector<int,Real> RealVector;
// typedef std::vector<RealMatrix> RealMatrixArray;

/// Reshape a matrix only if the sizes of the matrices differ
template < typename O, typename T >
void reshape( Teuchos::SerialDenseMatrix<O,T> &matrix, O M, O N )
{
  if ( ( matrix.numRows() != M )  || ( matrix.numCols() != N ) )
    matrix.reshape( M, N );
};

/// Reshape a vector only if the sizes of the vectors differ
template < typename O, typename T >
void resize( Teuchos::SerialDenseVector<O,T> &v, O n )
{
  if ( v.length() != n ) v.resize( n );
};

/// Fill a column of a matrix with the contents of a flattened matrix.
/// The matrix is flattened so that one column proceeds the next.
template < typename O, typename T >
void fill_column( int col_number, const Teuchos::SerialDenseMatrix<O,T> &source,
		  Teuchos::SerialDenseMatrix<O,T> &target )
{
  if ( source.numRows() * source.numCols() != target.numRows() )
    {
      std::stringstream msg;
      msg << "fill_column( matrix, matrix ) Matrix shapes are inconsistent.";
      msg << "\nsource is " << source.numRows() << " x " << source.numCols();
      msg << " and target is " << target.numRows() << " x " << target.numCols();
      msg << "\n";
      throw( std::runtime_error( msg.str() ) );
    }

  int iter( 0 );
  for ( int j = 0; j < source.numCols(); j++ )
    {
      for ( int i = 0; i < source.numRows(); i++ )
	{
	  target(iter,col_number) = source(i,j);
	  iter++;
	}
    }
};

/// Append a column vector to a matrix.
template < typename O, typename T >
void append_column( const Teuchos::SerialDenseVector<O,T> &vector, 
		    Teuchos::SerialDenseMatrix<O,T> &matrix )
{
  int N( matrix.numCols() );
  matrix.reshape( vector.length(), N + 1 );
  fill_column( N, vector, matrix );
};

/// Delete a column from a matrix
template < typename O, typename T >
void delete_column( O col_num, 
		    Teuchos::SerialDenseMatrix<O,T> &matrix,
		    bool resize = true )
{
  int M( matrix.numRows() ), N( matrix.numCols() );
  for ( O n = col_num+1; n < N; n++ )
    {
      for ( O m = 0; m < M; m++ )
	matrix(m,n-1) = matrix(m,n);
    }
  if ( resize )
    matrix.reshape( M, N-1 );
};

/// Append a row vector to a matrix.
template < typename O, typename T >
void append_row( const Teuchos::SerialDenseVector<O,T> &vector, 
		 Teuchos::SerialDenseMatrix<O,T> &matrix )
{
  int M( matrix.numRows() );
  matrix.reshape( M + 1, vector.length() );
  for ( int j = 0; j < vector.length(); j++ )
    matrix(M,j) = vector[j];
};

// Delete a row from a matrix
template < typename O, typename T >
void deleteRow( O row_num, 
		Teuchos::SerialDenseMatrix<O,T> &matrix,
		bool resize = true )
{
  int M( matrix.numRows() ), N( matrix.numCols() );
  for ( O n = 0; n < N; n++ )
    {
      for ( O m = row_num+1; m < M; m++ )
	matrix(m-1,n) = matrix(m,n);
    }
  if ( resize )
    matrix.reshape( M-1, N );
};

/// Append a matrix to the bottom of another matrix
template < typename O, typename T >
void row_append( const Teuchos::SerialDenseMatrix<O,T> &source, 
		 Teuchos::SerialDenseMatrix<O,T> &target )
{

  int Ms( source.numRows() ), Ns( source.numCols() ),
    Mt( target.numRows() ), Nt( target.numCols() );

  if ( ( Nt != Ns ) && ( Mt > 0 ) )
    {
      std::stringstream msg;
      msg << "row_append() Matrix shapes are inconsistent.";
      msg << "\nsource is " << Ms << " x " << Ns << " and target is ";
      msg << Mt << " x " << Nt << "\n";
      throw( std::runtime_error( msg.str() ) );
    }

  target.reshape( Mt + Ms, Ns );
  for ( int j = 0; j < Ns; j++ )
    {
      for ( int i = 0; i < Ms; i++ )
	{
	  target(Mt+i,j) = source(i,j);
	}
    }
};

/// Append a matrix to the right hand side of another matrix
template < typename O, typename T >
void column_append( const Teuchos::SerialDenseMatrix<O,T> &source, 
		    Teuchos::SerialDenseMatrix<O,T> &target )
{

  int Ms( source.numRows() ), Ns( source.numCols() ),
    Mt( target.numRows() ), Nt( target.numCols() );

  if ( ( Mt != Ms ) && ( Nt > 0 ) )
    {
      std::stringstream msg;
      msg << "column_append() Matrix shapes are inconsistent.";
      msg << "\nsource is " << Ms << " x " << Ns << " and target is ";
      msg << Mt << " x " << Nt << "\n";
      throw( std::runtime_error( msg.str() ) );
    }
  target.reshape( Ms, Nt + Ns );
  for ( int j = 0; j < Ns; j++ )
    {
      for ( int i = 0; i < Ms; i++ )
	{
	  target(i,Nt+j) = source(i,j);
	}
    }
};

/// Fill a row of a matrix with the contents of a flattened matrix.
/// The matrix is flattened so that one column proceeds the next.
template < typename O, typename T >
void fill_row( int row_number, const Teuchos::SerialDenseMatrix<O,T> &source,
	       Teuchos::SerialDenseMatrix<O,T> &target )
{
  if ( source.numRows() * source.numCols() != target.numCols() )
    {
      std::stringstream msg;
      msg << "fill_row( matrix, matrix ) Matrix shapes are inconsistent.";
      msg << "\nsource is " << source.numRows() << " x " << source.numCols();
      msg << " and target is " << target.numRows() << " x " << target.numCols();
      msg << "\n";
      throw( std::runtime_error( msg.str() ) );
    }

  int iter( 0 );
  for ( int j = 0; j < source.numCols(); j++ )
    {
      for ( int i = 0; i < source.numRows(); i++ )
	{
	  target(row_number,iter) = source(i,j);
	  iter++;
	}
    }
};

/** \brief Fill a column of a matrix with the contents of an array.
 *
 * The user must ensure that the length of source is consistent 
 * with the target. This cannot be used to store a flattened matrix
 * using source.values() if source is a sub view.
 */
template < typename O, typename T >
void fill_column( int col_number, const T* source,
		  Teuchos::SerialDenseMatrix<O,T> &target )
{
  for ( int i = 0; i < target.numRows(); i++ )
    {
      target(i,col_number) = source[i];
    }
};

/** \brief Return the transpose of a matrix.
 *
 * \param matrix (input) the original matrix 
 *
 * \param matrix_transpose (output) the transpose of the original matrix
 */
template < typename O, typename T >
void transpose( Teuchos::SerialDenseMatrix<O,T> &matrix, 
		Teuchos::SerialDenseMatrix<O,T> &matrix_transpose )
{
  int num_rows = matrix.numRows(), num_cols = matrix.numCols();
  matrix_transpose.shapeUninitialized( num_cols, num_rows );
  for ( O j = 0; j < num_rows; j++ )
    {
      for ( O i = 0; i < num_cols; i++ )
	{
	  matrix_transpose(i,j) = matrix(j,i);
	}
    }
};

/*
 * Define a weak ordering of vectors
 * @param v1 a "ScalarVector.hpp"
 * @param v2 a "ScalarVector.hpp"
 * @return true if any entries of v1 < v2
 */
template <typename O, typename T>
bool weakScalarVectorLessThan( const Teuchos::SerialDenseVector<O,T>& v1, 
			       const Teuchos::SerialDenseVector<O,T>& v2 )
{
  int n1 = v1.length(), n2 = v2.length();
  if ( n1 < n2 ) return true;
  if ( n1 > n2 ) return false;
  for (int i = 0; i < n1; i++)
    {
      if ( v1[i] < v2[i] ) return true;
      if ( v1[i] > v2[i] ) return false; 
    }
  
  return false;
};

template <typename O, typename T>
void permute_matrix_rows( Teuchos::SerialDenseMatrix<O,T> &A, IntVector &P )
{
  #ifdef DEBUG
  if ( A.numRows() != P.length() )
    throw( std::runtime_error("permute rows: A and P are inconsistent.") );
  #endif
  Teuchos::SerialDenseMatrix<O,T> tmp(Teuchos::Copy, A, A.numRows(), A.numCols() );
  for ( int j = 0; j < A.numCols(); j++ )
    {
      for ( int i = 0; i < P.length(); i++ )
	A(i,j) = tmp(P[i],j);
    }
};

template <typename O, typename T>
void permute_matrix_columns( Teuchos::SerialDenseMatrix<O,T> &A, IntVector &P )
{
  #ifdef DEBUG
  if ( A.numCols() != P.length() )
    throw( std::runtime_error("permute columns: A and P are inconsistent.") );
  #endif
  Teuchos::SerialDenseMatrix<O,T> tmp(Teuchos::Copy, A, A.numRows(), A.numCols() );
  for ( int j = 0; j < P.length(); j++ )
    {
      for ( int i = 0; i < A.numRows(); i++ )
	A(i,j) = tmp(i,P[j]);
    }
};

template <typename O, typename T>
T trace( Teuchos::SerialDenseMatrix<O,T> &A )
{
  if ( A.numRows() != A.numCols() )
    {
      std::string msg = "trace() A must be square";
      throw( std::runtime_error( msg ) );
    }
  T result = 0;
  for ( int i = 0; i < A.numRows(); i++ )
    result +=  A(i,i);
  return result;
}

/** \brief computes the Cholesky factorization of a real symmetric
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U**T * U,  if UPLO = 'U', or
*     A = L  * L**T,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  \param A (input) DOUBLE PRECISION positive-definite matrix, dimension (N,N)
*
*  \param for_lapack flag specifying whether the cholesky factor is to be used 
*  for a lapack function. If true then if A = L * L**T the upper part is not 
*  referenced. If false the upper part is filled in with zeros. 
*  Similarly for A = U**T * U
*/
int cholesky( RealMatrix &A, RealMatrix &result, Teuchos::EUplo uplo,
	      bool for_lapack );

/**
 * \brief Solves a system of linear equations A*X = B with a symmetric
 *  positive definite matrix A using a precomputed Cholesky factorization.
 * 
 *  Several right hand side vectors b and solution vectors x can be
 *  handled in a single call; they are stored as the columns of the
 *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
 *  matrix X.
 */
int solve_using_cholesky_factor( RealMatrix &L, RealMatrix& B, 
				 RealMatrix& result, Teuchos::EUplo uplo );

/**
 * \brief Solves a system of linear equations A*X = B with a symmetric
 *  positive definite matrix A using the Cholesky factorization. The cholesky 
 * factorization is computed internally.
 * 
 *  Several right hand side vectors b and solution vectors x can be
 *  handled in a single call; they are stored as the columns of the
 *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
 *  matrix X.
 *
 * \param rcond (INPUT/OUTPUT) if rcond < 0 compute the reciprocal of the 
 *                       condition number of A and store the result in rcond on 
 *			 output. if rcond > 0 rcond is not computed. 
 */
int cholesky_solve( RealMatrix& A, RealMatrix& B, RealMatrix& result,  
		    Real &rcond );

/** \brief Solves overdetermined or underdetermined real linear systems
 *  involving an M-by-N matrix A, or its transpose, using a QR or LQ
 *  factorization of A.  It is assumed that A has full rank.
 *
 *  Several right hand side vectors b and solution vectors x can be
 *  handled in a single call; they are stored as the columns of the
 *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
 *  matrix X.
 *
 * \todo at the moment A cannot be rank-defficient. implement using dgelsy.
 */
void qr_solve( RealMatrix &A, RealMatrix &B, RealMatrix &result, 
	       Teuchos::ETransp trans = Teuchos::NO_TRANS );

/** \brief Compute the minimum-norm solution to a real linear least
 * squares problem:
 *     minimize 2-norm(| b - A*x |)
 * using the singular value decomposition (SVD) of A. A is an M-by-N
 * matrix which may be rank-deficient.
 *
 * Several right hand side vectors b and solution vectors x can be
 * handled in a single call; they are stored as the columns of the
 * M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
 * X.
 *
 * The effective rank of A is determined by treating as zero those
 * singular values which are less than RCOND times the largest singular
 * value.
 *
 *  \param A (input) DOUBLE PRECISION matrix, dimension (M,N)
 *          On entry, the M-by-N matrix A.
 *          
 *  \param B (input) DOUBLE PRECISION matrix, dimension (M,NRHS)
 *          On entry, the M-by-NRHS right hand side matrix B.
 *          
 *  \param result_1 (output) DOUBLE PRECISION matrix, dimension (min(M,N))
 *          The singular values of A in decreasing order.
 *          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
 *
 *  \param RCOND (input) DOUBLE PRECISION
 *          RCOND is used to determine the effective rank of A.
 *          Singular values S(i) <= RCOND*S(1) are treated as zero.
 *          If RCOND < 0, machine precision is used instead.
 *
 *  \param RANK (output) INTEGER
 *          The effective rank of A, i.e., the number of singular values
 *          which are greater than RCOND*S(1).
 *
 *  \param result_0 (output) DOUBLE PRECISION matrix, dimension (M,NRHS)
 *          On exit, X is the N-by-NRHS solution
 *          matrix.  If m >= n and RANK = n, the residual
 *          sum-of-squares for the solution in the i-th column is given
 *          by the sum of squares of elements n+1:m in that column.
 *
 *  \todo implement dgelsd which is faster than dgelss
 */
void svd_solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0,
		RealVector &result_1, int &rank, Real rcond = -1 );

/** \brief Solves a triangular system
 *
 *  Solves a triangular system of the form
 *     A * X = B  or  A**T * X = B,
 *
 *  where A is a triangular matrix of order N, and B is an N-by-NRHS
 *  matrix.  A check is made to verify that A is nonsingular.
 *
 *  Arguments
 *  =========
 *
 *  \param A (input) DOUBLE PRECISION array, dimension (LDA,N)
 *          The triangular matrix A.  If uplo = Teuchos::UPPER_TRI, the leading N-by-N
 *          upper triangular part of the array A contains the upper
 *          triangular matrix, and the strictly lower triangular part of
 *          A is not referenced.  If UPLO = LOWER_TRU, the leading N-by-N lower
 *          triangular part of the array A contains the lower triangular
 *          matrix, and the strictly upper triangular part of A is not
 *          referenced.  If DIAG = UNIT_DIAG, the diagonal elements of A are
 *          also not referenced and are assumed to be 1.
 *
 *  \param B (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
 *          On entry, the right hand side matrix B.
 *          On exit, if INFO = 0, the solution matrix X.
 *
 *  \param result (output) DOUBLE PRECISION matrix, dimension (N,NRHS)
 *          On exit, X is the N-by-NRHS solution
 *          matrix.
 *
 *  \param uplo Specifies whether A is upper or lower traingular
 *
 *  \param trans Specifies the form of the system of equations:
 *          = Teuchos::NO_TRANS:  A * X = B  (No transpose)
 *          = TRANS :  A**T * X = B  (Transpose)
 *
 *  \param diag Specifies if the matrix is unit diagonal
 *
 */
void substitution_solve( RealMatrix &A, 
			 RealMatrix &B, 
			 RealMatrix &result, 
			 Teuchos::ETransp trans = Teuchos::NO_TRANS,
			 Teuchos::EUplo uplo = Teuchos::UPPER_TRI,
			 Teuchos::EDiag diag = Teuchos::NON_UNIT_DIAG );

/**
 * \brief Returns the QR factorization of [A a] given a qr factorization of A
 *
 * QR facorization is \f$O(N^3)\f$ but this update is only \f$O(N^2)\f$
 *
 * \param Q (input/output) The ( M x iter ) orthogonal matrix. 
 * On entry contains the original
 * ( M x tier-1 ) orthogonal matrix Q with an additional ( M x 1 ) 
 * possibly unitialized column allocated to allow the storage of the 
 * larger updated matrix Q. On exit Q is the updated orthogonal matrix factor
 * appended to it . 
 * Must have all excess memory initialized to zero.
 *
 * \param R (input/output) The ( iter x iter ) upper triangular matrix. 
 * On entry contains the ( iter-1 x iter-1 ) original matrix R 
 * with memory allocated to allow the storage of the larger updated matrix R
 * Must have all excess memory initialized to zero.
 *
 * \param col (intput) the new column being added to A.
 * 
 * \param iter specifies the size of the original Q and R matricies.
 *
 * \return info = 0 update sucessful. If info = 1, the new column was colinear
 * with the active set.
 */
int qr_factorization_update_insert_column( RealMatrix &Q, RealMatrix &R, 
					   RealMatrix &col, int iter );

/**
 * \brief Return the cholesky factorization of a positive definite grammian 
 matrix [A a]'*[A a] given a cholesky factorization of A'*A.
 *
 * Cholesky facorization is \f$O(N^3)\f$ but this update is only \f$O(N^2)\f$
 *
 * \param A (output) The ( M x iter ) Matrix used to compute the gramian A'A
 *
 * \param U (input/output) The ( iter x iter ) upper triangular matrix. 
 * On entry contains the ( iter-1 x iter-1 ) original matrix U 
 * with memory allocated to allow the storage of the larger updated matrix U.
 * Must have all excess memory initialized to zero.
 *
 * \param col (intput) the new column being added to A.
 * 
 * \param iter specifies the size of the original U matrix.

 * \param delta Regularization parameter. If delta > 0 the function will return 
 * chol( [A a]'*[A A] + delta*I )
 *
 * \return info = 0 update sucessful. If info = 1, the new column was colinear
 * with the active set.
 */
int cholesky_factorization_update_insert_column( RealMatrix &A, RealMatrix &U, 
						 RealMatrix &col, int iter,
						 Real delta = 0 );

/**
 * \brief Compute the givens rotation of a (2x1) vector x and also returned
 * the rotated vector.
 *
 * \param x (intput) (2x1) vector to be rotated
 *
 * \param x_rot (ouput) (2x1) rotated vector
 *
 * \parm givensMatrix (ouput) (2x2) Givens rotation matrix
 */
void givens_rotation( RealVector &x, RealVector &x_rot, 
		      RealMatrix &givensMatrix );

/**
 * \brief Update the cholesky factorization of a positive definite grammian 
 * matrix A'A when a column is deleted from A
 *
 * \param U (input/output) The ( N x N ) upper triangular matrix. 
 * On entry contains the ( N x N ) original matrix U 
 * On exit containt the ( N-1 x N-1 ) new matrix U with all entries in the 
 * Nth row and column set to zero.
 * 
 * \param col_index The index of the column to be deleted from A
 *
 * \param N the number of rows and columns of U
 */
void cholesky_factorization_update_delete_column( RealMatrix &U, 
						  int col_index,
						  int N);

// For qr updating for deleting and including a row go to
// http://www.maths.manchester.ac.uk/~clucas/updating/
// This code is in fortran and must be compiled and wrapped correctly

/**
 * Calculate an approximate solution to Ax=b using the conjugate gradient method
 *
 * \return info = 0 update sucessful. If info = 1, the matrix is not positive
 * definite. info = 2. the matrix contains inf or nan.
 */
int conjugate_gradients_solve( RealMatrix &A, RealVector &b, RealVector &x, 
			       Real &r_norm, 
			       Real cg_tol = 1e-8, int max_iter = 50,
			       int verbosity = 1 );

/** \brief Solve the linear equality-constrained least squares (LSE)
*  problem:
*
*          minimize || b - A*x ||_2   subject to   C*x = d
*
*  where A is an M-by-N matrix, C is a P-by-N matrix, b is a given
*  M-vector, and d is a given P-vector. It is assumed that
*  P <= N <= M+P, and
*
*           rank(C) = P and  rank( (A) ) = N.
*                                ( (C) )
*
*  These conditions ensure that the LSE problem has a unique solution,
*  which is obtained using a generalized RQ factorization of the
*  matrices (C, A) given by
*
*     C = (0 R)*Q,   A = Z*T*Q.
*
* \param A (input) DOUBLE PRECISION M-by-N matrix A.
*
* \param C (input) DOUBLE PRECISION P-by-N matrix B.
*
* \param b (input) DOUBLE PRECISION M vector
*          On entry, b contains the right hand side vector for the
*          least squares part of the LSE problem.
*
* \param d (input) DOUBLE PRECISION P vector
*          On entry, d contains the right hand side vector for the
*          constrained equation.
*
* \param x (output) DOUBLE PRECISION N -vector
*          On exit, X is the solution of the LSE problem.
*/
void equality_constrained_least_squares_solve( RealMatrix &A, 
					       RealVector &b,
					       RealMatrix &C, 
					       RealVector &d,
					       RealMatrix &x, 
					       int verbosity = 0 );
/**
 * \brief Computes the inverse of a real symmetric positive definite
 *  matrix A using the Cholesky factorization A = L*L**T or A = U**T*U
 *
 * \param U (input) DOUBLE PRECISION NxN matrix U. U is the  
 * lower or upper traingular factor L or U from the Cholesky factorization 
 * A = L*L**T or A = U**T*U
 *
 * \param result (output)  DOUBLE PRECISION NxN matrix. On exit result is the
 * inverse of A
 */
void cholesky_inverse(  RealMatrix &L, RealMatrix &result,
			Teuchos::EUplo uplo );

void pivoted_qr_factorization( RealMatrix &A, RealMatrix &result_0, 
			       RealMatrix &result_1, IntVector &result );

/**
 * \brief Compute an LU factorization with parital row pivoting of the
 * m-by-n matrix A. The factorization has the form A = P * L * U
 * where P is the permutation matrix, L is lower triangular with
 * unit diagonal elements and U is upper triangular.

 * /param max_iters the maximum number of pivots to perform. This is 
                    distinguishing feature from pivoted_lu_factorization()
 * /param num_initial_rows The first n rows of A which must be included in the
 * factorization before any other row can be included
 */
void truncated_pivoted_lu_factorization( RealMatrix &A,
					 RealMatrix &result_0, 
					 RealMatrix &result_1,
					 IntVector &values_0,
					 int max_iters,
					 int num_initial_rows );

template<typename O, typename T>
void eye( int N, Teuchos::SerialDenseMatrix<O,T> &result )
{
  result.shape( N, N );
  for ( O i = 0; i < N; i++ )
    result(i,i) = 1.;
};

template<typename O, typename T>
void unit_vector( int n, int k, Teuchos::SerialDenseVector<O,T> &result )
{
#ifdef DEBUG
  if ( k >=n ) 
    {
      std::string msg = "unit_vector() ensure k < n";
      throw( std::runtime_error( msg ) );
    }
#endif
  result.size( n );
  result[k] = 1.;
};

void lu_inverse( RealMatrix &L, RealMatrix &U, IntVector &p, 
		 RealMatrix &result );

void lu_solve( RealMatrix &A, 
	       RealMatrix &B, 
	       RealMatrix &result,
	       bool copy,
	       Teuchos::ETransp trans );


/*extern "C"
{
  void dgeqp3_( const int *M, const int *N, double *A, 
		const int *LDA, int *JPVT, double *TAU, 
		double *WORK, const int *LWORK, int *INFO );
		}*/

}  // namespace Pecos

#endif
