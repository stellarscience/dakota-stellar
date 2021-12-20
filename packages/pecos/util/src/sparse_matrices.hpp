/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_SPARSE_MATRICES_HPP
#define PECOS_UTIL_SPARSE_MATRICES_HPP

#include "linear_algebra.hpp"

namespace Pecos {
namespace util {

class BlockDiagonalMatrix
{
 private:
  int numBlocks_;
  std::vector<RealMatrix> blocks_;

 public:

  /// Default constructor
  BlockDiagonalMatrix() : numBlocks_(0) {};

  /// Deconstructor
  ~BlockDiagonalMatrix(){};

  /// Resize the matrix to have j number of diagonal blocks
  void resize( int num_blocks ){
    blocks_.resize( num_blocks );
    numBlocks_ = num_blocks;
  };
  
  /// Return the total number of blocks in the block diagonal matrix
  int num_blocks() const{
    return blocks_.size();
  }
  
  /// Return the total number of rows in the block diagonal matrix
  int num_rows() const{
    int num_rows = 0;
    for (int i=0; i<numBlocks_; i++ )
      num_rows += blocks_[i].numRows();
    return num_rows;
  }

  /// Return the total number of columns in the block diagonal matrix
  int num_cols() const{
    int num_cols = 0;
    for (int i=0; i<numBlocks_; i++ )
      num_cols += blocks_[i].numCols();
    return num_cols;
  }

  /**
   * \brief Set a block of the block diagonal matrix. Only use a view if 
   * the block being stored is a view
   */
  void set_matrix( int block_num, RealMatrix &block, 
		   Teuchos::DataAccess CV=Teuchos::View ){
    if ( block_num >= numBlocks_ ){
      std::string msg = "BlockDiagonalMatrix::set_matrix() ";
      msg += "Index out of range\n";
      throw( std::runtime_error( msg ) );
    }
    if ( CV == Teuchos::View ){
      RealMatrix empty_matrix;
      // Must set equal to empty matrix before copy otherwise
      // if first value of array points to the same data for both
      // the current block and the new block then the current block
      // will be returned and will not be update
      blocks_[block_num] = empty_matrix;
      blocks_[block_num] = block;
    }
    else{
      blocks_[block_num].shapeUninitialized( block.numRows(), block.numCols() );
      blocks_[block_num].assign( block );
    }
  }

  /**
   * \brief return C = B*A or C = B'*A, where B is the block diagonal matrix
   * and A is a dense matrix
   */
  void pre_multiply( const RealMatrix &matrix, RealMatrix &result, 
		     Teuchos::ETransp block_trans = Teuchos::NO_TRANS ) const;

  /**
   * \brief return C = A*B or C = A*B', where B is the block diagonal matrix
   * and A is a dense matrix
   */
  void post_multiply( const RealMatrix &matrix, RealMatrix &result, 
		      Teuchos::ETransp block_trans = Teuchos::NO_TRANS,
		      int subblock_num_rows = -1 ) const;

  /** 
   * \brief Get the k-th row of the j-th diagonal block
   */
  void get_row( int row_num, RealMatrix &row ) const;

  /// Print the blocks of the block diagonal matrix
  void print( std::ostream& os ) const{
    os << std::endl;
    os << "Block diagonal matrix" << std::endl;
    os << "Blocks: " << num_blocks() << std::endl;
    os << "---------------------" << std::endl;
    for (int i=0; i<numBlocks_; i++ ){
      os << "Block: " << i << "\n";
      blocks_[i].print( os );
    }
    os << "---------------------" << std::endl;
  };

  /// Convert the block diagonal matrix to a dense matrix
  void get_dense_matrix( RealMatrix &result) const{
    int sub_result_row_start = 0;
    int sub_result_col_start = 0;
    result.shape( num_rows(), num_cols() );
    for (int i=0; i<numBlocks_; i++ ){
      for (int k=0; k<blocks_[i].numCols(); k++ ){
	for( int j=0; j<blocks_[i].numRows(); j++ )
	  result(sub_result_row_start+j,sub_result_col_start+k) = 
	    blocks_[i](j,k);
      }
      sub_result_row_start += blocks_[i].numRows();
      sub_result_col_start += blocks_[i].numCols();
    }
  };

  /// Return the number of rows in the ith diagonal block 
  int num_rows( int block_num ) const{
    return blocks_[block_num].numRows();
  }

  /// Return the number of cols in the ith diagonal block 
  int num_cols( int block_num ) const{
    return blocks_[block_num].numCols();
  }

  

  /**
   * \brief multiply a matrix A by the nth  block Bn on the diagonal
   */
  void post_multiply_block( int block_num, const RealMatrix &matrix,
			    RealMatrix & result ) const;
};

}  // namespace util
}  // namespace Pecos

#endif  // include guard
