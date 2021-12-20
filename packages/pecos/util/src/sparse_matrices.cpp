/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "sparse_matrices.hpp"

namespace Pecos {
namespace util {

void BlockDiagonalMatrix::get_row( int row_num, RealMatrix &row ) const {
  int num_total_rows = 0;
  int num_total_rows_prev = 0;
  int block_row = 0;
  int i=0;
  for (i=0; i<numBlocks_; i++ ){
    num_total_rows_prev = num_total_rows;
    num_total_rows += blocks_[i].numRows();
    if ( row_num < num_total_rows ){
      block_row = row_num - num_total_rows_prev;
      break;
    }
  }
  if ( ( row.numRows() !=1 ) || (row.numCols() != blocks_[i].numCols() ) ) 
    row.shapeUninitialized( 1, blocks_[i].numCols() );
  for ( int j=0; j<blocks_[i].numCols(); j++)
    row(0,j) = blocks_[i](block_row,j);
}

void BlockDiagonalMatrix::pre_multiply( const RealMatrix &matrix, 
					RealMatrix &result, 
					Teuchos::ETransp block_trans ) const {
  int block_num_cols = num_cols();
  int result_num_rows = num_rows();
  if ( block_trans==Teuchos::TRANS ){
     block_num_cols = num_rows();
     result_num_rows = num_cols();
  }

  if ( block_num_cols != matrix.numRows() ){
    std::string msg = "BlockDiagonalMatrix::pre_multiply() ";
    msg += "Matrices sizes are inconsistent\n";
    throw( std::runtime_error( msg ) );
  }

  result.shapeUninitialized( result_num_rows, matrix.numCols() );
  int sub_matrix_start_row = 0; 
  int sub_result_start_row = 0; 
  for (int i=0; i<numBlocks_; i++ ){
    int num_block_rows = blocks_[i].numRows();
    int num_block_cols = blocks_[i].numCols();
    if ( block_trans == Teuchos::TRANS ){
      num_block_rows = blocks_[i].numCols();
      num_block_cols = blocks_[i].numRows();
    }
    int num_submatrix_rows = num_block_cols;
    RealMatrix sub_matrix( Teuchos::View, matrix, 
			   num_submatrix_rows, matrix.numCols(), 
			   sub_matrix_start_row, 0 );
    RealMatrix sub_result( Teuchos::View, result, 
			   num_block_rows, matrix.numCols(), 
			   sub_result_start_row, 0 );
    sub_result.multiply( block_trans, Teuchos::NO_TRANS, 
			 1.0, blocks_[i], sub_matrix, 0.0 );
    sub_matrix_start_row += num_submatrix_rows;
    sub_result_start_row += num_block_rows;
  }
}


void BlockDiagonalMatrix::post_multiply( const RealMatrix &matrix, 
					 RealMatrix &result, 
					 Teuchos::ETransp block_trans,
					 int subblock_num_rows ) const { 
  int block_num_rows = num_rows();
  int result_num_cols = num_cols();
  if ( block_trans==Teuchos::TRANS ){
     block_num_rows = num_cols();
     result_num_cols = num_rows();
  }

  if ( subblock_num_rows < 0 )
      subblock_num_rows = block_num_rows;
    
  if ( subblock_num_rows != matrix.numCols() ){
    std::string msg = "BlockDiagonalMatrix::post_multiply() ";
    msg += "Matrices sizes are inconsistent\n";
    throw( std::runtime_error( msg ) );
  }

  if ( subblock_num_rows > block_num_rows ){
    std::string msg = "BlockDiagonalMatrix::post_multiply() ";
    msg += "The number of subset rows (subblock_num_rows) is to large\n";
    throw( std::runtime_error( msg ) );
  }


  result.shape( matrix.numRows(), result_num_cols );
  int sub_matrix_start_col = 0; 
  int sub_result_start_col = 0; 
  
  bool done = false;
  int num_total_rows = 0;
  int num_total_rows_prev = 0;
  for (int i=0; i<numBlocks_; i++ ){
    int num_block_rows = blocks_[i].numRows();
    int num_block_cols = blocks_[i].numCols();
    if ( block_trans == Teuchos::TRANS ){
      num_block_rows = blocks_[i].numCols();
      num_block_cols = blocks_[i].numRows();
    }

    num_total_rows_prev = num_total_rows;
    num_total_rows += num_block_rows;
    int num_sub_block_rows = num_block_rows;
    int num_sub_block_cols = num_block_cols;
    if ( num_total_rows > subblock_num_rows ){
      num_sub_block_rows = subblock_num_rows - num_total_rows_prev;
      done = true;
    }
    int num_no_trans_sub_block_rows =  num_sub_block_rows;
    int num_no_trans_sub_block_cols =  num_sub_block_cols;
    if  ( block_trans == Teuchos::TRANS ){
      num_no_trans_sub_block_rows = num_sub_block_cols;
      num_no_trans_sub_block_cols = num_sub_block_rows;
    }
    RealMatrix sub_block( Teuchos::View, blocks_[i], num_no_trans_sub_block_rows,
 			  num_no_trans_sub_block_cols, 0, 0 );

    int num_submatrix_cols = num_sub_block_rows;
    RealMatrix sub_matrix( Teuchos::View, matrix, 
			   matrix.numRows(), num_submatrix_cols, 
			   0, sub_matrix_start_col );
    RealMatrix sub_result( Teuchos::View, result, 
			   matrix.numRows(), num_sub_block_cols,
			   0, sub_result_start_col );
    sub_result.multiply( Teuchos::NO_TRANS, block_trans, 
			 //1.0, sub_matrix, blocks_[i], 0.0 );
			 1.0, sub_matrix, sub_block, 0.0 );

    sub_matrix_start_col += num_submatrix_cols;
    sub_result_start_col += num_sub_block_cols;

    if ( done ) break;
  }
}

void BlockDiagonalMatrix::post_multiply_block( int block_num, 
					       const RealMatrix &matrix,
					       RealMatrix & result ) const
  {
    if ( block_num >= num_blocks() ){
    std::string msg = "BlockDiagonalMatrix::post_multiply_block() ";
    msg += "block num exceeds the number of blocks\n";
    throw( std::runtime_error( msg ) );
  }

    int block_num_rows = blocks_[block_num].numRows();
    int block_num_cols = blocks_[block_num].numCols();

    int matrix_num_rows = matrix.numRows();
    int matrix_num_cols = matrix.numCols();

    if ( block_num_rows != matrix_num_cols ){
    std::string msg = "BlockDiagonalMatrix::post_multiply_block() ";
    msg += "Matrices sizes are inconsistent\n";
    throw( std::runtime_error( msg ) );
    }

    result.shapeUninitialized( matrix_num_rows, block_num_cols );
    result.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		     1.0, matrix, blocks_[block_num], 0.0 );
  };

}  // namespace util
}  // namespace Pecos
