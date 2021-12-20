/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "linear_algebra.hpp"
#include "Teuchos_BLAS_wrappers.hpp"
#include "Teuchos_LAPACK_wrappers.hpp"
#include "Teuchos_SerialSpdDenseSolver.hpp"
#include "teuchos_data_types.hpp"
#define disp( a ) std::cout << #a << ": " << ( a ) << std::endl

// DLASWP and DGEQP3 are in Teuchos_LAPACK_wrappers.hpp
#define DPSTRF_F77 F77_BLAS_MANGLE(dpstrf,DPSTRF)
// unused:
// #define DGETC2_F77 F77_BLAS_MANGLE(dgetc2,DGETC2)
// #define DORGRQ_F77 F77_BLAS_MANGLE(dorgrq,DORGRQ)
extern "C"
{
  void DPSTRF_F77( const char* UPLO, int *N, double *A,
		   int *LDA, int *PIV, int *RANK, double *TOL, double *WORK,
		   int *INFO );

  // void DGETC2_F77( int *N, double *A, int *LDA, int *IPIV, int *JPIV, int *INFO );

  // void DORGRQ_F77( int *M, int *N, int *K, double *A, int* LDA, double* TAU, 
  // 		   double *WORK, const int *LWORK, int *INFO );
}

namespace Pecos {
namespace util {

template<> void GEMV<int,double>(
      Teuchos::ETransp trans, bool conjugate_trans, double alpha,
      const Teuchos::SerialDenseMatrix<int,double> &matrix,
      const Teuchos::SerialDenseVector<int,double> &vector,
      double beta,
      Teuchos::SerialDenseVector<int,double> &result){

  int result_len = matrix.numRows();
  if ((trans==Teuchos::TRANS) || (trans==Teuchos::CONJ_TRANS))
    result_len = matrix.numCols();
  
  if (result.length() != result_len){
    if (beta!=0.)
      throw(std::runtime_error("result inconsistent with matrix but beta!=0"));
    result.sizeUninitialized(result_len);
  }
  
  int inc=1;
  char T=Teuchos::ETranspChar[trans];
  int M=matrix.numRows(), N=matrix.numCols(), stride=matrix.stride();
  DGEMV_F77( &T, &M, &N, &alpha, 
	     matrix.values(), &stride, 
	     vector.values(), &inc, &beta, result.values(), &inc );
}

// BMA: Disabled the complex variant since no current use cases
// When reactivated will need to edit Dakota and Pecos CMakeLists.txt
//#ifdef HAVE_COMPLEX_BLAS
//template<> void GEMV<int,std::complex<double> >(
//      Teuchos::ETransp trans, bool conjugate_trans, std::complex<double> alpha,
//      const Teuchos::SerialDenseMatrix<int,std::complex<double> >&matrix,
//      const Teuchos::SerialDenseVector<int,std::complex<double> >&vector,
//      std::complex<double> beta,
//      Teuchos::SerialDenseVector<int,std::complex<double> > &result){
//
//  int result_len = matrix.numRows();
//  if ((trans==Teuchos::TRANS) || (trans==Teuchos::CONJ_TRANS))
//    result_len = matrix.numCols();
//  
//  if (result.length() != result_len){
//    if (beta!=0.)
//      throw(std::runtime_error("result inconsistent with matrix but beta!=0"));
//    result.sizeUninitialized(result_len);
//  }
//
//  if ((trans==Teuchos::TRANS) && conjugate_trans)
//    trans = Teuchos::CONJ_TRANS;
//  
//  int inc=1;
//  char T=Teuchos::ETranspChar[trans];
//  int M=matrix.numRows(), N=matrix.numCols(), stride=matrix.stride();
//  ZGEMV_F77( &T, &M, &N, &alpha, 
//	     matrix.values(), &stride, 
//	     vector.values(), &inc, &beta, result.values(), &inc );
//}
//#endif

int cholesky( const RealMatrix &A, RealMatrix &result, Teuchos::EUplo uplo, 
              bool for_lapack )
{
  Teuchos::LAPACK<int, Real> la;
  int M = A.numRows();
  result.reshape( M, M );
  result.assign( A );
  
  int info;
  la.POTRF( Teuchos::EUploChar[uplo], M, result.values(), result.stride(), 
	    &info );

  if ( info > 0 ) 
    {
      //std::cout << "cholesky() The matrix A is not positive definite\n";
      return info;
    }
  if ( info < 0 ) 
    {
      std::stringstream msg;
      msg << "cholesky() POTRF failed\n";
      msg << "The " << std::abs( info ) << "-th argument had an ";
      msg << "illegal value";
      throw( std::runtime_error( msg.str() ) );
      throw( std::runtime_error( msg.str() ) );
    }
  
  // lapack returns the lower/upper triangle of the (symmetric) inverse of A
  // so we must use symmetry to fill in the entries above/below the diagonal
  if ( !for_lapack )
    {
      if ( uplo == Teuchos::LOWER_TRI )
	{
	  for ( int j = 1; j < M; j++ )
	    for ( int i = 0; i < j; i++ )
	      result(i,j) = 0.;
	}
      else
	{
	  for ( int j = 1; j < M; j++ )
	    for ( int i = 0; i < j; i++ )
	      result(j,i) = 0.;
	}
    }

  return info;
};

int solve_using_cholesky_factor( const RealMatrix &L, const RealMatrix& B,
                                 RealMatrix& result, Teuchos::EUplo uplo )
{
  Teuchos::LAPACK<int, Real> la;
  int m( L.numRows() ), num_rhs( B.numCols() );
  result.reshape( B.numRows(), num_rhs );
  result.assign( B );

  // Solves the system of linear equations A*X = B with a symmetric
  // positive definite matrix A=LL' using the Cholesky factorization
  int info;
  la.POTRS( Teuchos::EUploChar[uplo], m, num_rhs,
	    L.values(), L.stride(),
	    result.values(), result.stride(), &info );
  result.reshape( L.numRows(), num_rhs );

  return info;
};

int teuchos_cholesky_solve( RealMatrix &A, RealMatrix &B, 
                            RealMatrix& result, Real &rcond ){
  // Copy A so that it is not affected outside the scope of this function
  RealSymMatrix chol_factor( A.numRows() );
  for (int j=0; j<A.numRows(); j++)
    {
      chol_factor(j,j) = A(j,j);
      for (int i=j+1; i<A.numRows();i++)
	{
	  chol_factor(i,j) = A(i,j);
	  chol_factor(j,i) = A(i,j);
	}
    }
  // Copy B so that it is not affected outside the scope of this function
  RealMatrix rhs( Teuchos::Copy, B, B.numRows(), B.numCols() );

  result.shapeUninitialized( rhs.numRows(), rhs.numCols() );
  Teuchos::SerialSpdDenseSolver<int, Real> solver;
  solver.setMatrix( rcp(&chol_factor, false) );
  solver.setVectors( rcp(&result, false), rcp(&rhs, false) );
  solver.factorWithEquilibration( true );
  int info(0);
  info = solver.factor();
  //solver.solveToRefinedSolution( true );
  info = solver.solve();
  info = solver.reciprocalConditionEstimate( rcond );
  return info;
}

int cholesky_solve(const RealMatrix& A, const RealMatrix& B,
                   RealMatrix& result, Real &rcond)
{
  Teuchos::LAPACK<int, Real> la;

  int m( A.numRows() );
  RealMatrix L;
  int info = cholesky( A, L, Teuchos::LOWER_TRI, true );

  if ( info != 0 ) return info;

  // Compute the reciprocal of the condition number of A = L*L**T or A = U**T*U 
  // from its cholesky decompostion computed by POTRF. 
  if ( rcond < 0 )
    {
      Real *work = new Real [3*m]; // workspace array 
      int *iwork = new int [m];  // workspace array
      la.POCON( Teuchos::EUploChar[Teuchos::LOWER_TRI], m, L.values(), 
		L.stride(), A.normOne(), &rcond, work, iwork, &info );
      delete [] work;
      delete [] iwork;
      if ( info < 0 ) 
	{
	  std::cout << "cholesky_solve() Incorrect arguments specified to ";
	  std::cout << "POCON()\n";
	  return info;
	}
    }

  info = solve_using_cholesky_factor( L, B, result, Teuchos::LOWER_TRI );

  return info;
};

void qr_factorization( const RealMatrix &A, RealMatrix &Qfactor, 
		       RealMatrix &Rfactor ){
  Teuchos::LAPACK<int, Real> la;
  int M( A.numRows() ), N( A.numCols() );
  int K( std::min(M,N) );
  RealMatrix qr_data( Teuchos::Copy, A, M, N );
  /*RealMatrix qr_data( M, std::max( M, N ), false );
  // Copy A into qr_data but if overdetermined allocate enough memory
  // to store Q when passed to dorgrq.
  for (int j=0; j<N; j++)
    for (int i=0; i<M; i++)
    qr_data(i,j) = A(i,j);*/
  
  //qr_data.shapeUninitialized( M, N );
  //qr_data.assign( A );

  //---------------------------------//
  // Get the optimal work array size //
  //---------------------------------//
  
  int lwork;     // Size of Teuchos::LAPACK work array
  Real *work;  // Teuchos::LAPACK work array
  int info;      // Teuchos::LAPACK output flag 
  lwork = -1;             // special code for workspace query
  work  = new Real [1]; // temporary work array
  //tau.sizeUninitialized( K ) ;// scalar factors of reflectors 
  RealVector tau( K, false );
  la.GEQRF( M, N, qr_data.values(), qr_data.stride(), tau.values(), 
	    work, lwork, &info );
  lwork = (int)work[0];  // optimal work array size returned by query
  delete [] work;
 
  //------------------------------//
  // Compute the QR factorization //
  //------------------------------//
  work  = new Real [lwork]; // Optimal work array
  la.GEQRF( M, N, qr_data.values(), qr_data.stride(), tau.values(), 
	    work, lwork, &info );
  delete [] work;

  //-------------------//
  // Form the R factor //
  //-------------------//

  Rfactor.shape( M, N ); // initialize to zero
  for (int j=0; j<N; j++)
    for (int i=0; i<=std::min(j,M-1); i++)
      Rfactor(i,j) = qr_data(i,j);

  //-------------------//
  // Form the Q factor //
  //-------------------//

  /*No economoy version.
    "economy size" R. If m>n, R has only n rows. 
    If m<=n, this is the same as R = qr(A).

  // Note ORGQR give the equivalent of Q from qr(econ) in matlab 
  // that is only the first N columns.
  // But I want all M columns so I use ORMQR to multiply the MxM identity by Q
  Qfactor.shape( M, M );; // initialize to zero
  for (int i=0; i<M; i++)
    Qfactor(i,i) = 1.;

  // query optimal work size
  lwork = -1;             // special code for workspace query
  work  = new Real [1]; // temporary work array
  int LDA = qr_data.stride();
  la.ORMQR( 'L', 'N', M, M, K, qr_data.values(), LDA, 
	    tau.values(), Qfactor.values(), Qfactor.stride(),
	    work, lwork, &info );
  lwork = (int)work[0];  // optimal work array size returned by query
  delete [] work;

  work  = new Real [lwork]; // Optimal work array
  la.ORMQR( 'L', 'N', M, M, K, qr_data.values(), LDA, 
	    tau.values(), Qfactor.values(), Qfactor.stride(),
	    work, lwork, &info );
  delete [] work; 
*/

  // NESTA needs economy version
  Qfactor.shape( M, N );; // initialize to zero
  for (int i=0; i<M; i++)
    Qfactor(i,i) = 1.;

  // query optimal work size
  lwork = -1;             // special code for workspace query
  work  = new Real [1]; // temporary work array
  la.ORGQR( M, N, tau.length(), qr_data.values(), qr_data.stride(), 
	    tau.values(), work, lwork, &info );
  lwork = (int)work[0];  // optimal work array size returned by query
  delete [] work;

  work  = new Real [lwork]; // Optimal work array
  la.ORGQR( M, N, tau.length(), qr_data.values(), qr_data.stride(), 
	    tau.values(), work, lwork, &info );
  delete [] work;

  Qfactor.shapeUninitialized( M, N );
  for (int j=0; j<N; j++)
    for (int i=0; i<M; i++)
      Qfactor(i,j) = qr_data(i,j);
  


}

/*void qr_solve( const RealMatrix &B, const RealMatrix &qr_data, 
	       const RealVector &tau, RealMatrix &result ){
  Teuchos::LAPACK<int, Real> la;
  int M( B.numRows() ), num_rhs( B.numCols() );
  int K = tau.length();
  RealMatrix QtB( Teuchos::Copy, B, M, num_rhs );
  
  // query optimal work size
  int info ;
  int lwork = -1;             // special code for workspace query
  Real *work  = new Real [1]; // temporary work array
  la.ORMQR( 'L', 'T', M, num_rhs, K, qr_data.values(), 
	    qr_data.stride(), tau.values(), QtB.values(), QtB.stride(),
	    work, lwork, &info );
  lwork = (int)work[0];  // optimal work array size returned by query
  delete [] work;
  work  = new Real [lwork]; // Optimal work array
  la.ORMQR( 'L', 'T', M, num_rhs, K, qr_data.values(), 
	    qr_data.stride(), tau.values(), QtB.values(), QtB.stride(),
	    work, lwork, &info );

  substitution_solve( qr_data, QtB, result, Teuchos::NO_TRANS, 
		      Teuchos::UPPER_TRI, Teuchos::NON_UNIT_DIAG );
  delete [] work;
		      }*/

void qr_solve( const RealMatrix &A, const RealMatrix &B, RealMatrix &result, 
	       Teuchos::ETransp trans )
{
  Teuchos::LAPACK<int, Real> la;

  RealMatrix A_copy( Teuchos::Copy, A, A.numRows(), A.numCols() );
  int M( A.numRows() ), N( A.numCols() ), num_rhs( B.numCols() );
  RealMatrix B_copy( Teuchos::Copy, B, B.numRows(), B.numCols() );;
  B_copy.reshape( std::max( M, N ), num_rhs );
  //---------------------------------//
  // Get the optimal work array size //
  //---------------------------------//
  
  int lwork;     // Size of Teuchos::LAPACK work array
  Real *work;  // Teuchos::LAPACK work array
  int info=0;      // Teuchos::LAPACK output flag 
  int lda = A_copy.stride();
  int ldb = B_copy.stride();

  lwork = -1;             // special code for workspace query
  work  = new Real [1]; // temporary work array
  la.GELS( Teuchos::ETranspChar[trans], M, N, num_rhs, A_copy.values(), 
	   lda, B_copy.values(), ldb, work, lwork, &info );
  lwork = (int)work[0];  // optimal work array size returned by query
  delete [] work;
  work  = new Real [lwork]; // Optimal work array

  //---------------------------------//
  // Solve Ax = b                    //
  //---------------------------------//  
  la.GELS( Teuchos::ETranspChar[trans], M, N, num_rhs, A_copy.values(), lda, 
	   B_copy.values(), ldb, work, lwork, &info );
  if ( info < 0 )
    {
      std::stringstream msg;
      msg << "qr_solve() dgels failed. ";
      msg << "The " << std::abs( info ) << "-th argument had an ";
      msg << "illegal value";
      throw( std::runtime_error( msg.str() ) );
    }
  if ( info > 0 )
    {
      std::stringstream msg;
      msg << "QR Solve dgels failed. ";
      msg << "The " << info << "-th diagonal element of the ";
      msg << "triangular factor of A is zero, so that A does not have";
      msg << "full rank; the least squares solution could not be computed.";
      throw( std::runtime_error( msg.str() ) );
    }
  delete [] work;
  result.reshape( N, num_rhs );
  for ( int j = 0; j < num_rhs; j++ )
    for ( int i = 0; i < N; i++ )
      result(i,j) = B_copy(i,j);
};


void svd_solve( const RealMatrix &A, const RealMatrix &B, RealMatrix &result_0,
		RealVector &result_1, int &rank, Real rcond )
{
  Teuchos::LAPACK<int, Real> la;

  //-----------------//
  // Allocate memory //
  //-----------------//

  int M( A.numRows() ),  N( A.numCols() ), num_rhs( B.numCols() );

  if (num_rhs<1)
    throw(std::runtime_error("B has zero columns"));
  
  RealMatrix A_copy( Teuchos::Copy, A, A.numRows(), A.numCols() );
  result_1.sizeUninitialized( std::min( M, N ) );

  //---------------------------------//
  // Get the optimal work array size //
  //---------------------------------//
  
  int lwork;     // Size of Teuchos::LAPACK work array
  Real *work;    // Teuchos::LAPACK work array
  int info;      // Teuchos::LAPACK output flag
  int lda = A_copy.stride();
  int ldb = std::max( std::max( B.stride(), lda ), N );

  result_0.shapeUninitialized( M, num_rhs );
  result_0.assign( B );
  result_0.reshape( ldb, num_rhs );

  lwork = -1;             // special code for workspace query
  work  = new Real [1]; // temporary work array
  la.GELSS( M, N, num_rhs, A_copy.values(), lda, result_0.values(), ldb, 
	    result_1.values(), rcond, &rank, work, lwork, &info );
  lwork = (int)work[0];  // optimal work array size returned by query

  delete [] work;
  work  = new Real [lwork]; // Optimal work array

  //---------------------------------//
  // Solve Ax = b                    //
  //---------------------------------//
  la.GELSS( M, N, num_rhs, A_copy.values(), lda, result_0.values(), ldb, 
	    result_1.values(), rcond, &rank, work, lwork, &info );

  result_0.reshape( N, num_rhs );
  delete [] work;
};

int cholesky_factorization_update_insert_column( RealMatrix &A, RealMatrix &U, 
						 RealMatrix &col, int iter,
						 Real delta )
{
  int info( 0 );
  Real col_norm = col.normFrobenius();
  if ( iter == 0 )
    {
      U(0,0) = std::sqrt( col_norm * col_norm + delta );
    }
  else
    {
      if ( iter >= U.numRows() )
	throw( std::runtime_error("cholesky_factorization_update_insert_column: iter out of bounds") );

      RealMatrix w;//( iter, 1, false );
      RealMatrix U_old( Teuchos::View, U, iter, iter, 0, 0 );
      // compute column k in gramian matrix A'A
      RealMatrix gramian_col( iter, 1, false );
      gramian_col.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
			    1.0, A, col, 0.0 );
      substitution_solve( U_old, gramian_col, w, Teuchos::TRANS,
			  Teuchos::UPPER_TRI );
      Real w_norm = w.normFrobenius();
      
      if ( col_norm * col_norm + delta - w_norm*w_norm <= 
	   std::numeric_limits<Real>::epsilon() )
	{
	  // New column is colinear. That is, it is in the span of the active
	  // set
	  info = 1;
	}
      else
	{
	  U(iter,iter) = std::sqrt(( col_norm*col_norm + delta )-w_norm*w_norm );
	  RealMatrix U_col( Teuchos::View, U, iter, 1, 0, iter );
	  U_col.assign( w );
	}
    }
  return info;
};

void givens_rotation( RealVector &x, RealVector &x_rot, RealMatrix &givens_matrix )
{
  givens_matrix.reshape( 2, 2 );
  x_rot.sizeUninitialized( x.length() );

  if ( x[1] == 0 )
    {
      //givens_matrix = I
      givens_matrix(0,0) = 1.0; givens_matrix(1,1) = 1.0;
      x_rot.assign( x );
    }
  else
    {
      Real x_norm = x.normFrobenius();

      givens_matrix(0,0) = x[0] / x_norm;
      givens_matrix(0,1) = x[1] / x_norm;
      givens_matrix(1,0) = -x[1] / x_norm;
      givens_matrix(1,1) = x[0] / x_norm;

      x_rot[0] = x_norm;
      x_rot[1] = 0.0;
    }
};

void cholesky_factorization_update_delete_column( RealMatrix &U, 
						  int col_index,
						  int N )
{
  if ( col_index != N - 1 )
    {
      // delete column but do not resize the matrix U
      delete_column( col_index, U, false );
    };
  
  RealVector x( 2, false );
  for ( int n = col_index; n < N-1; n++ )
    {
      RealMatrix givens_matrix;
      RealVector x_rot;
      x[0] = U(n,n); x[1] = U(n+1,n);
      givens_rotation( x, x_rot, givens_matrix );
      U(n,n) = x_rot[0]; U(n+1,n) = x_rot[1];
      if ( n < N - 2 )
	{
	  RealMatrix U_sub( Teuchos::View, U, 2, N - n - 1, n, n + 1 );
	  RealMatrix U_sub_rot( U_sub.numRows(), U_sub.numCols(), false );
	  U_sub_rot.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
			      1.0, givens_matrix, U_sub, 0.0 );
	  U_sub.assign( U_sub_rot );
	}
    }

  // Zero out last row and column of U
  for ( int m = 0; m < N; m++ ) U(m,N-1) = 0.0;
  for ( int n = 0; n < N; n++ ) U(N-1,n) = 0.0;
};

int conjugate_gradients_solve( const RealMatrix &A, const RealVector &b, RealVector &x, 
			       Real &relative_residual_norm,
			       int &iters_taken,
			       Real cg_tol, int max_iter,
			       int verbosity )
{
  int info( 0 );

  int M( A.numRows() ), N( A.numCols() );

  if ( x.length() != N )
    // Initialize to zero
    x.size( N );

  RealVector current_x(Teuchos::Copy, x.values(), x.length());


  Real b_norm = b.normFrobenius();
  RealVector residual(Teuchos::Copy, b.values(), b.length());
  residual.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		     -1.0, A, current_x, 1.0 );

  if ( b_norm < std::numeric_limits<Real>::epsilon() ) 
    b_norm = 1.0;

  relative_residual_norm = residual.normFrobenius() / b_norm ;

  if ( Teuchos::ScalarTraits<Real>::isnaninf( relative_residual_norm ) )
    {
      // At least one of the matrix inputs contains nan of inf
      if ( verbosity > 2 )
	{
	  std::stringstream msg;
	  msg << "conjugate_gradient_solve() Warning: at least one of the ";
	  msg << "matrix inputs contains nan and/or inf.\n";
	  std::cout << msg.str();
	}
      info = 3;
      return info;
    }
  
  if ( relative_residual_norm <= cg_tol )
    {
      info = 0;
      return info;
    }

  if ( verbosity > 2 )
    {
      std::cout << "CG iteration: " << 0 << ", residual: ";
      std::cout << relative_residual_norm << "\n";
    }

  RealVector p( residual );
  RealVector Ap( A.numRows() );
  
  iters_taken = 0;
  bool done( false );
  Real rtr_old( residual.dot( residual ) );

  while ( !done )
    {
      Ap.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		   1.0, A, p, 0.0 );

      Real ptAp = p.dot( Ap );

      Real alpha = rtr_old / ptAp;

      if ( Teuchos::ScalarTraits<Real>::isnaninf( alpha ) || ( ptAp <= 0 ) )
	{
	  if ( verbosity > 2 )
	    {
	      // The matrix is not positive definite
	      std::stringstream msg;
	      msg << "conjugate_gradient_solve() Warning: A is not postive ";
	      msg << "definite.\n";
	      std::cout << msg.str();
	    }
	  info = 2;
	  return info;
	}

      for ( int n = 0; n < N; n++ )
	current_x[n] += alpha * p[n];

      if ( (iters_taken+1)%50 == 0 )
	{
	  // Avoid numerical drift due to rounding error
	  residual.assign( b );
	  residual.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
			     -1.0, A, current_x, 1.0 );
	}
      else
	{
	  for ( int m = 0; m < M; m++ )
	    residual[m] -= alpha * Ap[m];
	}
      
      Real rtr = residual.dot( residual );

      for ( int m = 0; m < M; m++ )
	p[m] = residual[m] + rtr / rtr_old * p[m];

      rtr_old = rtr;

      Real current_relative_residual_norm = std::sqrt( rtr_old ) / b_norm;

      if ( current_relative_residual_norm < relative_residual_norm )
	{
	  relative_residual_norm = current_relative_residual_norm;
	  for ( int n = 0; n < N; n++ )
	    x[n] = current_x[n];
	}

      if ( verbosity > 2 )
	{
	  std::cout << "CG iteration: " << iters_taken + 1<< ", residual: ";
	  std::cout <<  current_relative_residual_norm << "\n";
	}

      if ( current_relative_residual_norm < cg_tol )
	{
	  if ( verbosity > 2 )
	    {
	      std::stringstream msg;
	      msg << "conjugate_gradient_solve() Exiting residual below ";
	      msg << "tolerance.\n";
	      std::cout << msg.str();
	    }
	  done = true;
	  info = 0;
	}

      iters_taken++;
      if ( iters_taken == max_iter )
	{
	  if ( verbosity > 2 )
	    {
	      std::stringstream msg;
	      msg << "conjugate_gradient_solve() Exiting maximum number of ";
	      msg << "iterations reached.\n";
	      std::cout << msg.str();
	    }
	  done = true;
	  info = 1;
	}
    }    
  return info;
};


void equality_constrained_least_squares_solve( const RealMatrix &A, 
					       const RealVector &b,
					       const RealMatrix &C, 
					       const RealVector &d,
					       RealVector &x, 
					       int verbosity )
{
  RealMatrix A_copy( Teuchos::Copy, A, A.numRows(), A.numCols() ),
    C_copy( Teuchos::Copy, C, C.numRows(), C.numCols() );
  RealVector b_copy( Teuchos::Copy, b.values(), b.length() ), 
    d_copy( Teuchos::Copy, d.values(), d.length() );

  int M( A_copy.numRows() ), N( A_copy.numCols() ), lda( A_copy.stride() ), 
    ldc( C_copy.stride() );

  x.sizeUninitialized( N );

  Teuchos::LAPACK<int, Real> la;
  double* work;    // LAPACK work array
  int info( 0 );   // LAPACK output flag
  int lwork; // size of LAPACK work array
  int num_cons( C_copy.numRows() ); // number of equality constraints

  // Get the optimal work array size
  lwork = -1; // special code for workspace query
  work  = new double [1]; // temporary work array
  la.GGLSE( M, N, num_cons, A_copy.values(), lda, 
	    C_copy.values(), ldc, b_copy.values(), d_copy.values(), x.values(), 
	    work, lwork, &info );
  lwork = (int)work[0]; // optimal work array size returned by query
  delete [] work;
  work  = new double [lwork]; // Optimal work array

  // Least squares computation using LAPACK's DGGLSE subroutine which uses
  // a GRQ factorization method for solving the eq-constrained LLS problem
  info = 0;
  la.GGLSE( M, N, num_cons, A_copy.values(), lda, C_copy.values(),
	    ldc, b_copy.values(), d_copy.values(), x.values(), 
	    work, lwork, &info );
  delete [] work;

  if ( info < 0 )
    {
      std::stringstream msg;
      msg << "equality_constrained_least_squares() dgglse failed. ";
      msg << "The " << std::abs( info ) << "-th argument had an ";
      msg << "illegal value";
      throw( std::runtime_error( msg.str() ) );
    }
  if ( info == 1 )
    {
      std::stringstream msg;
      msg << "the upper triangular factor R associated with C in the ";
      msg << "generalized RQ factorization of the pair (C, A) is ";
      msg << "singular, so that rank(C) < num_cons; the least squares ";
      msg << "solution could not be computed.";
      throw( std::runtime_error( msg.str() ) );
    }
  if ( info == 2 )
    {
      std::stringstream msg;
      msg << "the (N-P) by (N-P) part of the upper trapezoidal factor ";
      msg << "T associated with A in the generalized RQ factorization ";
      msg << "of the pair (C, A) is singular, so that\n";
      msg << "rank( (A) ) < N; the least squares solution could not\n";
      msg << "    ( (C) )\n";
      msg << "be computed.";
      throw( std::runtime_error( msg.str() ) );
    }
}

void cholesky_inverse( RealMatrix &U, RealMatrix &result, Teuchos::EUplo uplo  )
{
  Teuchos::LAPACK<int, Real> la;
  int N = U.numRows();
  // copy U into result. Lapack will compute result inplace
  result.shapeUninitialized( N, N );
  result.assign( U );
  
  // extract arguments needed for lapack call
  int lda( result.stride() );
  int info( 0 );
  la.POTRI( Teuchos::EUploChar[uplo], N, result.values(), lda, &info );

  if ( info < 0 )
    {
      std::stringstream msg;
      msg << "cholesky_inverse() dpotri failed. ";
      msg << "The " << std::abs( info ) << "-th argument had an ";
      msg << "illegal value";
      throw( std::runtime_error( msg.str() ) );
    }
  else if ( info > 0 )
    {
      std::stringstream msg;
      msg << "cholesky_inverse() dpotri failed. ";
      msg << "The (" << info << "," << info << ") element of the factor U or L is ";
      msg << "zero and the inverse could not be computed";
      throw( std::runtime_error( msg.str() ) );
    }

  // lapack returns the lower/upper triangle of the (symmetric) inverse of A
  // so we must use symmetry to fill in the entries above/below the diagonal
  if ( uplo == Teuchos::LOWER_TRI )
    {
      for ( int j = 1; j < N; j++ )
	for ( int i = 0; i < j; i++ )
	  result(i,j) = result(j,i);
    }
  else
    {
      for ( int j = 1; j < N; j++ )
	for ( int i = 0; i < j; i++ )
	  result(j,i) = result(i,j);
    }
};

void pivoted_qr_factorization( const RealMatrix &A, RealMatrix &Q, RealMatrix &R,
                               IntVector &p )
{
  Teuchos::LAPACK<int, Real> la;

  RealMatrix A_copy( Teuchos::Copy, A, A.numRows(), A.numCols() );
  int M( A.numRows() ), N( A.numCols() ), K( std::min( M, N ) );

  //Q.shapeUninitialized( M, K );
  Q.shape( M, K );
  R.shape( K, N ); // must be initialized to 0
  p.size( N ); // must be initialized to 0

  //---------------------------------//
  // Get the optimal work array size //
  //---------------------------------//
  
  int lwork;   // Size of Teuchos::LAPACK work array
  Real *work;  // Teuchos::LAPACK work array
  int info;    // Teuchos::LAPACK output flag 
  int lda = std::max( 1, A_copy.stride() );
  RealVector tau( K, true );

  lwork = -1;           // special code for workspace query
  work  = new Real [1]; // temporary work array
  
  DGEQP3_F77( &M, &N, A_copy.values(), &lda, p.values(), tau.values(), 
	      work, &lwork, &info );

  lwork = (int)work[0];  // optimal work array size returned by query
  delete [] work;
  work  = new Real [lwork]; // Optimal work array

  //---------------------------------//
  // Compute the QR factorization    //
  //---------------------------------//

  DGEQP3_F77( &M, &N, A_copy.values(), &lda, p.values(), tau.values(), 
	      work, &lwork, &info );

  if ( info < 0 )
    {
      std::stringstream msg;
      msg << "privoted_qr_factorization() dgeqp3 failed. ";
      msg << "The " << std::abs( info ) << "-th argument had an ";
      msg << "illegal value";
      throw( std::runtime_error( msg.str() ) );
    }

  delete [] work;

  //---------------------------------//
  // Form the Q and R matrices       //
  //---------------------------------//
  
  for ( int m = 0; m < K; m++ )
    {
      for ( int n = m; n < N; n++ )
	R(m,n) = A_copy(m,n);
    }

  lwork = -1;           // special code for workspace query
  work  = new Real [1]; // temporary work array

  la.ORGQR( M, K, K, A_copy.values(), lda, tau.values(), 
	   work, lwork, &info );
  
  lwork = (int)work[0];  // optimal work array size returned by query
  delete [] work;

  work  = new Real [lwork]; // Optimal work array

  la.ORGQR( M, K, K, A_copy.values(), lda, tau.values(), 
	   work, lwork, &info );

  for ( int n = 0; n < K; n++ )
    {
      for ( int m = 0; m < M; m++ )
	{
	  Q(m,n) = A_copy(m,n);
	}
    }

  // fortran returns indices 1,...,N
  // c++ requires 0,...,N-1
  for ( int n = 0; n < N; n++ )
    p[n]--;
 
  delete [] work;
}

void lu_inverse( const RealMatrix &L, const RealMatrix &U, const IntVector &p,
		 RealMatrix &LU_inv )
{
  int M = L.numRows(), N = U.numCols();
#ifdef DEBUG
  if ( M != N ) 
    {
      std::string msg="LU_inverse() The dimensions of L and U are inconsistent";
      throw( std::runtime_error( msg ) );
    }
#endif
  LU_inv.shape( M, M );
  RealMatrix I;
  eye( M, I );
  if ( p.length() != 0 )
    permute_matrix_columns( I, p );
  RealMatrix X;
  substitution_solve( L, I, X, Teuchos::NO_TRANS,
		      Teuchos::LOWER_TRI, Teuchos::NON_UNIT_DIAG ); 
  substitution_solve( U, X, LU_inv, Teuchos::NO_TRANS,
		      Teuchos::UPPER_TRI, Teuchos::NON_UNIT_DIAG );
};

void lu_solve( RealMatrix &A, 
	       const RealMatrix &B, 
	       RealMatrix &result,
	       bool copy,
	       Teuchos::ETransp trans ){
  Teuchos::LAPACK<int, Real> la;
  int M = A.numRows(), N = A.numCols();

  RealMatrix A_copy;
  if ( copy )
    {
      A_copy.shapeUninitialized( M, N );
      A_copy.assign( A );
    }

  IntVector ipiv( std::min( M, N ), false );
  int info;
  if ( copy )
    la.GETRF( M, N, A_copy.values(), A_copy.stride(), ipiv.values(), &info );
  else
    la.GETRF( M, N, A.values(), A.stride(), ipiv.values(), &info );

  if ( info < 0 )
    {
      std::stringstream msg;
      msg << "GETRF: The " << std::abs(info) <<  "ith argument had " <<
	"an illegal value";
      throw( std::runtime_error( msg.str() ) );
    }
  if ( info > 0 )
    {
      std::stringstream msg;
	msg << "U(" <<  info << "," << info << ") is exactly zero. " << 
	"The factorization has been completed, but the factor U is exactly " <<
	"singular, and division by zero will occur if it is used " <<
	"to solve a system of equations";
      throw( std::runtime_error( msg.str() ) );
    }
  result.shapeUninitialized( B.numRows(), B.numCols() );
  result.assign( B );

  if ( copy )
    la.GETRS( Teuchos::ETranspChar[trans], M, result.numCols(), A_copy.values(),
	      A_copy.stride(), ipiv.values(), result.values(), result.stride(), 
	      &info );
  else
    la.GETRS( Teuchos::ETranspChar[trans], M, result.numCols(), A.values(),
	      A.stride(), ipiv.values(), result.values(), result.stride(), 
	      &info );


  if ( info < 0 )
    {
      std::stringstream msg;
      msg << "GETRS: The " <<std::abs(info) <<  "ith argument had " <<
	"an illegal value";
      throw( std::runtime_error( msg.str() ) );
    }
};

void pivoted_lu_factorization( const RealMatrix &A, 
			       RealMatrix &L, 
			       RealMatrix &U,
			       IntVector &pivots )
{
  Teuchos::LAPACK<int, Real> la;

  RealMatrix A_copy( Teuchos::Copy, A, A.numRows(), A.numCols() );
  int M( A.numRows() ), N( A.numCols() ), K( std::min( M, N ) );  

  int info;
  IntVector la_pivots( K, false );
  la.GETRF( M, N, A_copy.values(), A_copy.stride(), la_pivots.values(), &info );

  if ( info < 0 )
    {
      std::stringstream msg;
      msg << "GETRF: The " << std::abs(info) <<  "ith argument had " <<
	"an illegal value";
      throw( std::runtime_error( msg.str() ) );
    }
  if ( info > 0 )
    {
      std::stringstream msg;
      msg << "U(" <<  info << "," << info << ") is exactly zero. " << 
	"The factorization has been completed, but the factor U is exactly " <<
	"singular, and division by zero will occur if it is used " <<
	"to solve a system of equations";
      throw( std::runtime_error( msg.str() ) );
    }

  L.shape( M, K );
  U.shape( K, N );

  // Fill U
  for ( int n = 0; n < N; n++ )
    {
      if ( n < K ) U(n,n) = A_copy(n,n);
      for ( int m = 0; m < std::min(n,K); m++ )
	U(m,n) = A_copy(m,n);
    }

  for ( int n = 0; n < K; n++ )
    {
      L(n,n) = 1.;
      for ( int m = n+1; m < M; m++ )
	L(m,n) = A_copy(m,n);
    }

  // Lapack does not support integers so must cast. 
  // TODO write implementation for integers
  RealVector dbl_pivots( M, false ); 
  for ( int n = 0; n < dbl_pivots.length(); n++ )
    dbl_pivots[n] = (double)n;

  // convert pivots to a format that supports slicing (e.g python) to swap rows
  RealVector row_pivots;
  pivot_matrix_rows( dbl_pivots, la_pivots, 1, true,  row_pivots );

  pivots.sizeUninitialized( K );
  for ( int n = 0; n < K; n++ )
    pivots[n] = (int)row_pivots[n];
  
}


void pivoted_cholesky_factorization( RealMatrix &A, RealMatrix &L, IntVector &p,
				     int &rank, double tol )
{
  Teuchos::LAPACK<int, Real> la;

  RealMatrix A_copy( Teuchos::Copy, A, A.numRows(), A.numCols() );
  int N( A.numRows() );

  L.shape( N, N );
  p.size( N ); // must be initialized to 0

  //------------------------------------//
  // Compute the cholesky factorization //
  //------------------------------------//
  
  Real *work;  // Teuchos::LAPACK work array
  int info;    // Teuchos::LAPACK output flag 
  int lda = std::max( 1, A_copy.stride() );
  rank = 0;

  work  = new Real [2*N]; // temporary work array
  
  //Teuchos::EUplo uplo
  //Teuchos::EUploChar[uplo]
  char uplo = 'L';
  DPSTRF_F77( &uplo, &N, A_copy.values(),
	      &lda, p.values(), &rank, &tol, work, &info );
  delete [] work;

  if ( info < 0 )
    {
      std::stringstream msg;
      msg << "privoted_cholesky_factorization() dpstrf failed. ";
      msg << "The " << std::abs( info ) << "-th argument had an ";
      msg << "illegal value";
      throw( std::runtime_error( msg.str() ) );
    }
  else if ( info > 0 )
    {
      std::stringstream msg;
      msg << "privoted_cholesky_factorization() dpstrf failed. ";
      msg << "The matrix A is either rank deficient with computed rank " << rank
	  << " , or is indefinite.  See Section 7 of "
	  << "LAPACK Working Note #161 for further information.\n";
      //std::cout << msg.str();
    }
  
  // Fill upper triangular part of cholesky factor
  for ( int m = 0; m < N; m++ )
    {
    for ( int n = 0; n < m+1; n++ )
      L(m,n) = A_copy(m,n);
    }

  // fortran returns indices 1,...,N
  // c++ requires 0,...,N-1
  for ( int n = 0; n < N; n++ )
    p[n]--;
}

Real cholesky_condition_number( RealMatrix &L )
{
  Teuchos::LAPACK<int, Real> la;
  int m( L.numRows() );
  Real *work = new Real [3*m]; // workspace array 
  int *iwork = new int [m];  // workspace array
  Real rcond = 0.;
  int info;
  la.POCON( Teuchos::EUploChar[Teuchos::LOWER_TRI], m, L.values(), 
	    L.stride(), L.normOne(), &rcond, work, iwork, &info );
  delete [] work;
  delete [] iwork;
  if ( info < 0 ) 
    {
      std::stringstream msg; 
      msg << "cholesky_condition_number() Incorrect arguments specified to "
	  << "POCON()\n";
      throw( std::runtime_error( msg.str() ) );
    }
  return  rcond;
}

void symmetric_eigenvalue_decomposition( const RealMatrix &A,
                                         RealVector &eigenvalues,
                                         RealMatrix &eigenvectors )
{
  Teuchos::LAPACK<int, Real> la;

  int N( A.numRows() );
  eigenvectors.shapeUninitialized( N, N );
  eigenvectors.assign( A );

  char jobz = 'V'; // compute eigenvectors

  char uplo = 'U'; // assume only upper triangular part of matrix is stored

  eigenvalues.sizeUninitialized( N );

  int info;        // Teuchos::LAPACK output flag
  RealVector work; // Teuchos::LAPACK work array;
  int lwork = -1;  // Size of Teuchos::LAPACK work array
  
  // Compute optimal size of work array
  work.sizeUninitialized( 1 ); // temporary work array
  la.SYEV( jobz, uplo, N, eigenvectors.values(), eigenvectors.stride(), 
	   eigenvalues.values(), work.values(), lwork, &info );

  lwork = (int)work[0];
  work.sizeUninitialized( lwork );
  
  la.SYEV( jobz, uplo, N, eigenvectors.values(), eigenvectors.stride(), 
	   eigenvalues.values(), work.values(), lwork, &info );

  if ( info > 0 )
    {
      std::stringstream msg;
      msg << "The algorithm failed to converge." << info
	  << " off-diagonal elements of an intermediate tridiagonal "
	  << "form did not converge to zero.";
      throw( std::runtime_error( msg.str() ) );
    }
  else if ( info < 0 )
    {
      std::stringstream msg;
      msg << " The " << std::abs( info ) << " argument had an illegal value.";
      throw( std::runtime_error( msg.str() ) );
    }
};

void complete_pivoted_lu_factorization( const RealMatrix &A, 
					RealMatrix &L_factor, 
					RealMatrix &U_factor,
					IntVector &row_pivots,
					IntVector &column_pivots,
					int max_iters )
{
  int num_rows = A.numRows(), num_cols = A.numCols();
  max_iters = std::min( max_iters, std::min( num_rows, num_cols ) );

  RealMatrix A_copy( Teuchos::Copy, A, A.numRows(), A.numCols() );
  row_pivots.sizeUninitialized( num_rows );
  for ( int i = 0; i < num_rows; i++ )
    row_pivots[i] = i;
  column_pivots.sizeUninitialized( num_cols );
  for ( int j = 0; j < num_cols; j++ )
    column_pivots[j] = j;

  RealVector max_col_values( num_cols, false );
  IntVector argmax_col_values( num_cols, false );
  RealVector temp_row( num_cols, false );
  RealVector temp_col( num_rows, false );
  for ( int k = 0; k < std::min( num_rows-1, num_cols ); k++ )
    {
      Real max_pivot = -1.;
      int pivot_column = k, pivot_row = k;
      max_col_values = -1.; 
      for ( int j = k; j < num_cols; j++ )
	{
	  for ( int i = k; i < num_rows; i++ )
	    {
	      Real pivot = std::abs( A_copy(i,j) );
	      if ( pivot > max_col_values[j-k] )
		{
		  max_col_values[j-k] = pivot;
		  argmax_col_values[j-k] = i-k;
		  if ( pivot > max_pivot )
		    {
		      max_pivot = pivot;
		      pivot_column = j-k;
		    }
		}
	    }
	}
      pivot_row = argmax_col_values[pivot_column]+k;
      pivot_column += k;

      // update pivot vectors
      int temp_index = row_pivots[k];
      row_pivots[k] = row_pivots[pivot_row];
      row_pivots[pivot_row] = temp_index;
      temp_index = column_pivots[k];
      column_pivots[k] = column_pivots[pivot_column];
      column_pivots[pivot_column] = temp_index;

      // swap rows
      for ( int j = 0; j < num_cols; j++ )
	{
	  Real temp_value = A_copy(k,j);
	  A_copy(k,j) = A_copy(pivot_row,j);
	  A_copy(pivot_row,j) = temp_value;
	}
      
      // swap columns
      for ( int i = 0; i < num_rows; i++ )
	{
	  Real temp_value = A_copy(i,k);
	  A_copy(i,k) = A_copy(i,pivot_column);
	  A_copy(i,pivot_column) = temp_value;
	}

      // check for singularity
      if ( std::abs( A_copy(k,k) ) < 1e-10 )
	{
	  std::cout << "pivot " << std::abs( A_copy(k,k) ) << " is to small. "
		    << "Stopping factorization.\n";
	  break;
	}

      // update LU factoization
      for ( int i = k+1; i < num_rows; i++ )
	A_copy(i,k) /= A_copy(k,k);

      RealMatrix sub_matrix( Teuchos::View, A_copy, num_rows-k-1, num_cols-k-1,
			     k+1, k+1 );
      RealMatrix col_vector( Teuchos::View, A_copy, num_rows-k-1, 1, k+1, k );
      RealMatrix row_vector( Teuchos::View, A_copy, 1, num_cols-k-1, k, k+1 );

      sub_matrix.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
			   -1.0, col_vector, row_vector, 1.0 );

      if ( k >= max_iters )
	break;
    }

  // final column swap if num_cols > num_rows
  if ( num_cols > num_rows && num_cols >= max_iters )
    {
      // determine pivot which is the maximum entry in A_copy
      Real max_entry = 0.;
      int pivot_column = num_rows-1;
      for ( int j = num_rows-1; j < num_cols; j++ )
	{
	  Real pivot = std::abs( A_copy(num_rows-1,j) );
	  if ( pivot >= max_entry )
	    {
	      max_entry = pivot;
	      pivot_column = j-(num_rows-1);
	    }
	}
      pivot_column += num_rows-1;
        

      // update pivot column vector
      int temp_index = column_pivots[num_rows-1];
      column_pivots[num_rows-1] = column_pivots[pivot_column];
      column_pivots[pivot_column] = temp_index;

      // swap column
      for ( int i = 0; i < num_rows; i++ )
	{
	  Real temp_value = A_copy(i,num_rows-1);
	  A_copy(i,num_rows-1) = A_copy(i,pivot_column);
	  A_copy(i,pivot_column) = temp_value;
	}
    }

  // Build L and U matrices
  int max_num_rows_cols = std::min( num_rows, num_cols );
  L_factor.shape( num_rows, max_num_rows_cols );
  U_factor.shape( max_num_rows_cols, num_cols );
  for ( int j = 0; j < max_num_rows_cols; j++ )
    {
      L_factor(j,j) = 1.;
      for ( int i = j+1; i < num_rows; i++ )
	L_factor(i,j) = A_copy(i,j);
    }
  
  for ( int j = 0; j < num_cols; j++ )
    {
      for ( int i = 0; i < std::min( num_rows, j+1 ); i++ )
	U_factor(i,j) = A_copy(i,j);
    }
};


/*
Does not seem to work
void complete_pivoted_lu_factorization_lapack( RealMatrix &A, 
					RealMatrix &L_factor, 
					RealMatrix &U_factor,
					IntVector &row_pivots,
					IntVector &column_pivots )
{
  RealMatrix A_copy( Teuchos::Copy, A, A.numRows(), A.numCols() );
  int M( A.numRows() ), N( A.numCols() );
  
  if ( M != N )
    throw( std::runtime_error("complete_pivoted_lu_factorization: matrix must be square" ) );

  row_pivots.sizeUninitialized( N );
  column_pivots.sizeUninitialized( N );
  row_pivots = -2;
  // dgetc2_ does not update pivots correctly seems to exit early.
  
  int lda = A_copy.stride(), info;
  disp( N );
  disp( lda );
  DGETC2_F77( &N, A_copy.values(), &lda, row_pivots.values(), 
	      column_pivots.values(), &info );
  
  disp( info );
  disp( row_pivots );
  if ( info > 0 )
    {
      std::stringstream msg;
      msg << "U(" <<  info-1 << "," << info-1 << ") is likely to produce " << 
	"overflow if we try to solve for x in Ax=b. So U is peturbed to avoid "
	"overflow.\n";
      std::cout << msg.str();
    }

  L_factor.shape( N, N );
  U_factor.shape( N, N );

  // Fill U
  for ( int n = 0; n < N; n++ )
    {
      for ( int m = 0; m < n; m++ )
	U_factor(m,n) = A_copy(m,n);
    }

  // Fill L
  for ( int n = 0; n < N; n++ )
    {
      L_factor(n,n) = 1.;
      for ( int m = n+1; m < N; m++ )
	L_factor(m,n) = A_copy(m,n);
    }

  // fortran returns indices 1,...,N
  // c++ requires 0,...,N-1
  for ( int n = 0; n < N; n++ )
    {
      row_pivots[n]--;
      column_pivots[n]--;
    }
}
*/

void pivot_matrix_rows( const RealMatrix &A, const IntVector &pivots, int dir, 
			bool fortran_indexing, RealMatrix &result )
{
  result.shapeUninitialized( A.numRows(), A.numCols() );
  result.assign( A );
  // Add 1 to pivots because we are using fotran indexing
  IntVector pivots_copy( pivots.length(), false );
  int shift = 1;
  if ( fortran_indexing ) shift = 0;
  for ( int i = 0; i < pivots.length(); i++ )
    pivots_copy[i] = pivots[i]+shift;
  // if dir is negative pivots are applied in reverse order
  int N = result.numCols(), LDA = result.stride(), K1 = 1, 
    K2 = pivots_copy.numRows(), INCX = dir;
  DLASWP_F77( &N, result.values(), &LDA, &K1, &K2, pivots_copy.values(), &INCX );
}


void truncated_pivoted_lu_factorization( const RealMatrix &A,
					 RealMatrix &L_factor, 
					 RealMatrix &U_factor,
					 IntVector &pivots,
					 int max_iters,
					 int num_initial_rows )
{

  Teuchos::BLAS<int, Real> blas;
  int num_rows = A.numRows(), num_cols = A.numCols();
  int min_num_rows_cols = std::min( num_rows, num_cols );
  max_iters = std::min( max_iters, num_rows );
  if ( A.numCols() < max_iters ){
    std::string msg = "truncated_pivoted_lu_factorization: ";
    msg += " A is inconsistent with max_iters. Try deceasing max_iters or ";
    msg += " increasing the number of columns of A";
    throw(std::runtime_error( msg ) );
  }

  // Use L to store both L and U during factoriation then copy out U in post
  // proccesing
  L_factor.shapeUninitialized( num_rows, min_num_rows_cols ); 
  L_factor.assign( A );
  pivots.sizeUninitialized( num_rows );
  for ( int i = 0; i < num_rows; i++ )
    pivots[i] = i;

  int iter = 0;
  for ( iter = 0; iter < min_num_rows_cols; iter++ )
    {
      // find best pivot
      int pivot =-1;
      if ( iter < num_initial_rows ){
	pivot = blas.IAMAX( num_initial_rows-iter, &(L_factor[iter][iter]), 1 )
	  -1+iter;	
      }else{
	pivot = blas.IAMAX( num_rows-iter, &(L_factor[iter][iter]), 1 )-1+iter;	
      }
      
      // update pivots vector
      int temp_index = pivots[iter];
      pivots[iter] = pivots[pivot];
      pivots[pivot] = temp_index;
      // apply pivots(swap rows) in L factorization
      for ( int j = 0; j < num_cols; j++ ){
	Real temp_value = L_factor(iter,j);
	L_factor(iter,j) = L_factor(pivot,j);
	L_factor(pivot,j) = temp_value;
      }

      // check for singularity
      if ( std::abs(L_factor(iter,iter))<std::numeric_limits<double>::epsilon()){
	std::cout << "pivot " << std::abs( L_factor(iter,iter) ) 
		  << " is to small. Stopping factorization.\n";
	break;
      }
      // update L_factor factoization
      for ( int i = iter+1; i < num_rows; i++ )
	L_factor(i,iter) /= L_factor(iter,iter);
      
      RealMatrix sub_matrix( Teuchos::View, L_factor, num_rows-iter-1, 
			     num_cols-iter-1,
			     iter+1, iter+1 );
      RealMatrix col_vector( Teuchos::View, L_factor, num_rows-iter-1, 1, 
			     iter+1, iter );
      RealMatrix row_vector( Teuchos::View, L_factor, 1, num_cols-iter-1, iter, 
			     iter+1 );

      sub_matrix.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
			   -1.0, col_vector, row_vector, 1.0 );
      if ( iter >= max_iters-1 )
	break;
    }

  // build L and U matrices

  // extract entries of U from L
  U_factor.shape( min_num_rows_cols, num_cols );
  for ( int j = 0; j < num_cols; j++ )
    {
      for ( int i = 0; i < std::min( min_num_rows_cols, j+1 ); i++ )
	U_factor(i,j) = L_factor(i,j);
    }

  // zero out U entries in L
  L_factor.reshape( iter+1, min_num_rows_cols );
  for ( int j = 0; j < min_num_rows_cols; j++ )
    {
      if ( j < iter+1 ) L_factor(j,j) = 1.;
      for ( int i = 0; i<std::min(j,iter+1); i++ )
	L_factor(i,j) = 0.;
    }
  // remove unused pivot entries
  pivots.resize( iter+1 );
}

}  // namespace util
}  // namespace Pecos
