#include "LinearAlgebra.hpp"

namespace Pecos {
int cholesky( RealMatrix &A, RealMatrix &result, Teuchos::EUplo uplo, 
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
      PCout << "cholesky() The matrix A is not positive definite\n";
      return info;
    }
  if ( info < 0 ) 
    {
      PCout << "cholesky() Incorrect arguments specified to POTRF()\n";
      return info;
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

int solve_using_cholesky_factor( RealMatrix &L, RealMatrix& B, 
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

int cholesky_solve( RealMatrix& A, RealMatrix& B, RealMatrix& result, Real &rcond )
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
	  PCout << "cholesky_solve() Incorrect arguments specified to ";
	  PCout << "POCON()\n";
	  return info;
	}
    }

  info = solve_using_cholesky_factor( L, B, result, Teuchos::LOWER_TRI );

  return info;
};

void qr_solve( RealMatrix &A, RealMatrix &B, RealMatrix &result, 
	       Teuchos::ETransp trans )
{
  Teuchos::LAPACK<int, Real> la;

  RealMatrix A_copy(Teuchos::Copy, A, A.numRows(), A.numCols());
  int M( A.numRows() ), N( A.numCols() ), num_rhs( B.numCols() );
  result.reshape( N, num_rhs );
  result.assign( B );

  //---------------------------------//
  // Get the optimal work array size //
  //---------------------------------//
  
  int lwork;     // Size of Teuchos::LAPACK work array
  Real *work;  // Teuchos::LAPACK work array
  int info;      // Teuchos::LAPACK output flag 
  int lda = A_copy.stride();
  int ldb = result.stride();

  lwork = -1;             // special code for workspace query
  work  = new Real [1]; // temporary work array
  la.GELS( Teuchos::ETranspChar[trans], M, N, num_rhs, A_copy.values(), 
	   lda, result.values(), ldb, work, lwork, &info );
  // Note la.GELS does not work because line 1358 in Teuchos_Teuchos::LAPACK.hpp 
  // last argument to DGELSS_F77 uses a & which should not be there
  lwork = (int)work[0];  // optimal work array size returned by query
  delete [] work;
  work  = new Real [lwork]; // Optimal work array

  //---------------------------------//
  // Solve Ax = b                    //
  //---------------------------------//

  la.GELS( Teuchos::ETranspChar[trans], M, N, num_rhs, A_copy.values(), lda, 
	   result.values(), ldb, work, lwork, &info );
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
};


void svd_solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0,
		RealVector &result_1, int &rank, Real rcond )
{
  Teuchos::LAPACK<int, Real> la;

  //-----------------//
  // Allocate memory //
  //-----------------//

  int M( A.numRows() ),  N( A.numCols() ), num_rhs( B.numCols() );
  RealMatrix A_copy(Teuchos::Copy, A, A.numRows(), A.numCols());
  result_1.sizeUninitialized( std::min( M, N ) );

  //---------------------------------//
  // Get the optimal work array size //
  //---------------------------------//
  
  int lwork;     // Size of Teuchos::LAPACK work array
  Real *work;  // Teuchos::LAPACK work array
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

void substitution_solve( RealMatrix &A, 
			 RealMatrix &B, 
			 RealMatrix &result,
			 Teuchos::ETransp trans,
			 Teuchos::EUplo uplo,
			 Teuchos::EDiag diag )
{
  int M( A.numRows() ), num_rhs( B.numCols() );

  Teuchos::LAPACK<int, Real> lapack;
  result.reshape( M, num_rhs );
  result.assign( B );

  int info;
  lapack.TRTRS( Teuchos::EUploChar[uplo], Teuchos::ETranspChar[trans], 
		Teuchos::EDiagChar[diag], 
		M, num_rhs, A.values(), A.stride(), 
		result.values(), result.stride(), &info );

  if ( info < 0 )
    {
      std::stringstream msg;
      msg << "substitution_solve() dtrtrs failed. ";
      msg << "The " << std::abs( info ) << "-th argument had an ";
      msg << "illegal value";
      throw( std::runtime_error( msg.str() ) );
    }
  if ( info > 0 )
    {
      std::stringstream msg;
      msg << "substitution_solve() dtrtrs failed. ";
      msg << "The " << info << "-th diagonal element of A is zero ";
      msg << "indicating that the matrix is singular and the solutions ";
      msg << "X have not been computed.";
      throw( std::runtime_error( msg.str() ) );
    }
   
}

int qr_factorization_update_insert_column( RealMatrix &Q, RealMatrix &R, 
					   RealMatrix &col, int iter )
{
  int info( 0 );
  int M( col.numRows() );
  Real col_norm = col.normFrobenius();

  if ( iter == 0 )
    {
      R(0,0) = col_norm;
      for ( int m = 0; m < M; m++ )
	Q(m,0) = col(m,0) / col_norm;
    }
  else
    {
      RealMatrix Q_old( Teuchos::View, Q, M, iter, 0, 0 );
      RealMatrix w( iter, 1, false ); 
      w.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
		  1.0, Q_old, col, 0.0 );
      Real w_norm = w.normFrobenius();

      if ( col_norm * col_norm - w_norm*w_norm <= 
	   std::numeric_limits<Real>::epsilon() )
	{
	  // New column is colinear. That is, it is in the span of the active
	  // set
	  info = 1;
	}
      else
	{
	  R(iter,iter) = std::sqrt( col_norm * col_norm - w_norm * w_norm );
	  RealMatrix R_col( Teuchos::View, R, iter, 1, 0, iter );
	  // must use assign below because operator= will not work
	  // because it will call deleteArrays when w is a copy 
	  // ( which it is here )
	  R_col.assign( w ); 
  
	  RealMatrix Qw( M, 1, false );
	  Qw.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		       1.0, Q_old, w, 0.0 );

	  for ( int m = 0; m < M; m ++ )
	    Q(m,iter) = ( col(m,0) - Qw(m,0) ) / R(iter,iter);
	}
    }
  return info;
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
      RealMatrix w( iter, 1, false );
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

int conjugate_gradients_solve( RealMatrix &A, RealVector &b, RealVector &x, 
			       Real &relative_residual_norm,
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
	  PCout << msg.str();
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
      PCout << "CG iteration: " << 0 << ", residual: ";
      PCout << relative_residual_norm << "\n";
    }

  RealVector p( residual );
  RealVector Ap( A.numRows() );
  
  int iter( 0 );
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
	      PCout << msg.str();
	    }
	  info = 2;
	  return info;
	}

      for ( int n = 0; n < N; n++ )
	current_x[n] += alpha * p[n];

      if ( (iter+1)%50 == 0 )
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
	  PCout << "CG iteration: " << iter + 1<< ", residual: ";
	  PCout <<  current_relative_residual_norm << "\n";
	}

      if ( current_relative_residual_norm < cg_tol )
	{
	  if ( verbosity > 2 )
	    {
	      std::stringstream msg;
	      msg << "conjugate_gradient_solve() Exiting residual below ";
	      msg << "tolerance.\n";
	      PCout << msg.str();
	    }
	  done = true;
	  info = 0;
	}

      iter++;
      if ( iter == max_iter )
	{
	  if ( verbosity > 2 )
	    {
	      std::stringstream msg;
	      msg << "conjugate_gradient_solve() Exiting maximum number of ";
	      msg << "iterations reached.\n";
	      PCout << msg.str();
	    }
	  done = true;
	  info = 1;
	}
    }    
  return info;
};


void equality_constrained_least_squares_solve( RealMatrix &A, 
					       RealVector &b,
					       RealMatrix &C, 
					       RealVector &d,
					       RealMatrix &x, 
					       int verbosity )
{
  RealMatrix A_copy(Teuchos::Copy, A, A.numRows(), A.numCols()), 
    C_copy(Teuchos::Copy, C, C.numRows(), C.numCols());
  RealVector b_copy(Teuchos::Copy, b.values(), b.length()), 
    d_copy(Teuchos::Copy, d.values(), d.length());

  int M( A_copy.numRows() ), N( A_copy.numCols() ), lda( A_copy.stride() ), 
    ldc( C_copy.stride() );

  x.shapeUninitialized( N, 1 );

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

void pivoted_qr_factorization( RealMatrix &A, RealMatrix &Q, RealMatrix &R,
			       IntVector &p )
{
  Teuchos::LAPACK<int, Real> la;

  RealMatrix A_copy(Teuchos::Copy, A, A.numRows(), A.numCols());
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
  //dgeqp3_( &M, &N, A_copy.values(), &lda, p.values(), tau.values(), 
  //work, &lwork, &info );

  lwork = (int)work[0];  // optimal work array size returned by query
  delete [] work;
  work  = new Real [lwork]; // Optimal work array

  //---------------------------------//
  // Compute the QR factorization    //
  //---------------------------------//

  DGEQP3_F77( &M, &N, A_copy.values(), &lda, p.values(), tau.values(), 
	      work, &lwork, &info );
  //dgeqp3_( &M, &N, A_copy.values(), &lda, p.values(), tau.values(), 
  //work, &lwork, &info );

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

  la.ORGQR ( M, K, K, A_copy.values(), lda, tau.values(),
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

void lu_inverse( RealMatrix &L, RealMatrix &U, IntVector &p, RealMatrix &LU_inv )
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


void truncated_pivoted_lu_factorization( RealMatrix &A,
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
      if ( std::abs( L_factor(iter,iter) ) < 1e-15 ){
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
      if ( iter >= max_iters )
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
  L_factor.reshape( iter, min_num_rows_cols );
  for ( int j = 0; j < min_num_rows_cols; j++ )
    {
      if ( j < iter ) L_factor(j,j) = 1.;
      for ( int i = 0; i<std::min(j,iter); i++ )
	L_factor(i,j) = 0.;
    }
  // remove unused pivot entries
  pivots.resize( iter );
}

void lu_solve( RealMatrix &A, 
	       RealMatrix &B, 
	       RealMatrix &result,
	       bool copy,
	       Teuchos::ETransp trans )
{
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


} // namespace Pecos
