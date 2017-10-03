#include "MathTools.hpp"

namespace Pecos {

void ind2sub( const IntVector &sizes, int index, int num_elems,
	      IntVector &result )
{

  int denom = num_elems;
  int num_sets = sizes.length();
  resize( result, num_sets );
  for ( int i = num_sets-1; i >= 0; i-- )
    {
      denom /= sizes[i];
      result[i] = index / denom;
      index = index % denom;
    }
};

int sub2ind( const IntVector &sizes, const IntVector &multi_index )
{
  int num_dims = sizes.length();
  int index = 0, shift = 1;
  for ( int i = 0; i < num_dims; i++ )
    {
      index += shift * multi_index[i];
      shift *= sizes[i];
    }
  return index;
};

void linspace( RealVector &result, Real a, Real b, int n )
{
  result.resize( n );
  Real h = ( b - a ) / (Real)(n-1);
  for ( int i = 0; i < n ; i++ )
    {
      result[i] = a + (Real)i * h;
    }
};

Real round( Real r ) 
{
  return (r > 0.0) ? std::floor(r + 0.5) : std::ceil(r - 0.5);
};

bool abs_compare( Real a, Real b )
{
  return ( std::fabs( a ) < std::fabs( b ) );
};

int nearest_integer( Real x )
{
  int value;
  
  if ( x < 0.0 )
    {
      value = - ( int ) ( std::fabs ( x ) + 0.5 );
    }
  else
    {
      value =   ( int ) ( std::fabs ( x ) + 0.5 );
    }

  return value;
};

Real factorial(int n)
{
  Real temp = 1.0;
  for ( int i = 1; i <= n; i++ )
    {
      temp *= (Real)(i);
    }
  return temp;
};

int nchoosek(int n, int k)
/* Determine the number of possible combinations of k elements that can be
   chosen from n elements. */
{
  Real value = 1;
  for ( int i = 0; i < n-k; i++ )
    {
      value *= (Real)(n-i) / (Real)(n-k-i);
    }
  return (int)round( value );
  //return (int)(Factorial(n)/(Factorial(n-k)*Factorial(k)));
};

int get_total_degree_from_num_samples( int num_dims,
				       int num_samples ){
  int degree = 0;
  while ( true ){
    int num_terms = nchoosek(degree+num_dims,num_dims);
    if ( num_terms >= num_samples ){
      break;
    }
    degree++;
  }
  return degree;
}

Real rescale( Real x, Real a, Real b, int dir )
{
  if ( dir == 0 )
    {
      return a + ( b-a ) * x; // [0,1] -> [a,b]
    }
  else
    {
      return ( x - a ) / ( b - a); // [a,b] -> [0,1]
    }
};

void rescale( RealMatrix &x, const RealVector &domain, int dir )
{
  for ( int i = 0; i < x.numCols(); i++ )
    {
      for ( int d = 0; d < x.numRows(); d++ )
	{
	  x(d,i) = rescale( x(d,i), domain[2*d], domain[2*d+1], dir );
	}
    }
};

void hypercube_map( RealMatrix &x, const RealVector &domain_in,  
		    const RealVector &domain_out,
		    RealMatrix &result )
{
  result.shapeUninitialized( x.numRows(), x.numCols() );
  result.assign( x );
  rescale( result, domain_in, 1 );
  rescale( result, domain_out, 0 );
};

void mesh_grid( const IntVector &num_pts_1d, 
		const RealVector &domain, 
		RealMatrix &result )
{
#ifdef DEBUG
  if ( num_pts_1d.length() % 2 != 0 )
    {
      
      throw(std::invalid_argument("Ensure num_pts_1d.length() %2 == 0."));
    }
#endif

  int num_dims( domain.length() / 2 );

  std::vector< RealVector > pointSets1D;
  pointSets1D.resize( num_dims );  
  for ( int d = 0; d < num_dims; d++ )
    {
      linspace( pointSets1D[d], domain[2*d], domain[2*d+1], num_pts_1d[d] );
    }
  cartesian_product( pointSets1D, result, 1 );
};

void get_multi_dimensional_polynomial_subspace_indices( IntMatrix &B, 
							int* elems, 
							int len_elems, 
							int* pos, int len_pos, 
							int choices_made, 
							int first_pos, 
							int order, int &col )
{
  int start,temp;
  // Have we selected the number of required elements?
  if ( choices_made >= len_pos )
    {
      start = 0;
      temp = 0;
      for ( int j = 0; j < len_pos; j++ )
	{
	  if ( pos[j]-start == 0 )
	    {
	      B(j,col) = 0;
	    }
	  else
	    {
	      B(j,col) = pos[j] - start;
	      temp = temp + ( pos[j] - start );
	    }
	  start = pos[j] + 1;
	}
      B(len_pos,col) = order - temp;
      col++;
      return;
    }
  
  // Are there enough remaining elements to be selected?
  if ( ( len_elems - first_pos ) < ( len_pos - choices_made ) )
    {
      return;
    }

  // Try to select new elements in the same position or
  // to the right of the last selected one.
  for ( int j = first_pos; j < len_elems; j++ )
    {
      pos[choices_made] = j;
      get_multi_dimensional_polynomial_subspace_indices( B, elems, len_elems, 
							 pos, len_pos, 
							 choices_made+1, j+1, 
							 order, col );
      //change to ii if items can be selected repeatedly: Do Not for gPC
    }
  return;
};

void get_multi_dimensional_polynomial_indices( int num_dims, int degree, 
					       IntMatrix &result )
{
  // m: degree of polynomial m=0,...,degree
  // i: index of basis function i=0,...,P-1
  // j: index of random dimension j=0,...,num_dims-1
  
  int len_elems, len_pos, choices_made,first_pos;
  int *elems,*pos;
  int P = nchoosek( num_dims + degree - 1, num_dims - 1 );
  if ( ( result.numRows() != num_dims ) || ( result.numCols() != P ) )
    {
      //std::string msg = "get_multi_dimensional_polynomial_indices() ";
      //msg += "The size of result must be set before entry to function.";
      //throw( std::runtime_error( msg ) );
      result.shape( num_dims, P );
    }

  int col = 0;
  len_pos = num_dims - 1;
  pos = new int [len_pos];
  
  // Generate all multimonomial indices of degree degree
  len_elems = ( num_dims - 1 ) + degree;
  elems = new int [len_elems];

  for ( int i = 0; i < len_elems; i++ )
    {
      elems[i] = i;
      if ( i < len_pos )
	{
	  pos[i] = 0;
	}
    }
  choices_made = 0;
  first_pos = 0;
  get_multi_dimensional_polynomial_subspace_indices( result, elems, len_elems, 
						     pos, len_pos, choices_made,
						     first_pos, degree, col );
  delete[] elems;
  delete[] pos;
};

void set_hypercube_domain( RealVector &result, int num_dims, Real a, Real b )
{
  result.resize( 2*num_dims );
  for ( int d = 0; d < num_dims; d++ )
    {
      result[2*d] = a;
      result[2*d+1] = b;
    }
};

void lp_error( RealMatrix &referenceValues, RealMatrix &approximateValues,
	      std::vector<lp_norm> error_norms, RealMatrix &error, 
	      IntVector &activeColumns, bool normalise )
{
  Teuchos::BLAS<int, Real> blas;

  if ( referenceValues.numRows() != approximateValues.numRows() )
    throw( std::runtime_error( "lp_error() Matrix sizes do not match." ) );

  if ( activeColumns.length() == 0 )
    range( activeColumns, 0, approximateValues.numCols(), 1 );
 
  RealMatrix diff( referenceValues );
  diff -= approximateValues;

  int M( diff.numRows() );

  reshape( error, activeColumns.length(), (int)error_norms.size() );

  for ( int j = 0; j < (int)error_norms.size(); j++ )
    {
      switch ( error_norms[j] )
	{
	case mt_linf_norm:
	  {
	    // Infinity norm
	    for ( int i = 0; i < activeColumns.length(); i++ )
	      {
		// Access the active columns
		int colNumber( activeColumns[i] );
		Real* diffCol = diff[colNumber];
		// Find the index of the element of the col with the 
		// maximum magnitude. 
		int amax = blas.IAMAX( M, diffCol, 1 ) - 1;
		error(i,j) = std::abs(diffCol[amax]);
		if ( normalise )
		  {
		    // Normalise by the half the range of the data in
		    // col
		    Real minima, maxima;
		    Real *refCol = referenceValues[colNumber];
		    minima = refCol[argmin( M, refCol )];
		    maxima = refCol[argmax( M, refCol )];
		    error(i,j) /= ( 0.5 * ( maxima - minima ) );
		  }
	      }
	    break;
	  }
	case mt_l1_norm:
	  {
	    // L1 norm
	    for ( int i = 0; i < activeColumns.length(); i++ )
	      {
		// Access the active column
		int colNumber( activeColumns[i] );
		Real* diffCol = diff[colNumber];
		//Sum the absolute values of the entries of col. 
		error(i,j) = blas.ASUM( M, diffCol, 1 ) / (Real)M;
	      }
	    break;
	  }
	case mt_l2_norm:
	  {
	    // L2 norm
	    for ( int i = 0; i < activeColumns.length(); i++ )
	      {
		// Access the active column
		int colNumber( activeColumns[i] );
		Real *diffCol = diff[colNumber];
		//Sum the absolute values of the entries of col. 
		error(i,j) = blas.NRM2( M, diffCol, 1 ) / std::sqrt( (Real)M );
		if ( normalise )
		  {
		    // Normalise by the standard deviation of the data in
		    // col
		    Real *refCol = referenceValues[colNumber];
		    Real stddev = std::sqrt( variance( M, refCol ) );
		    error(i,j) /= stddev;
		  }
	      }
	    break;
	  }
	}
    }
};

void latin_hypercube_design( int num_pts, int num_dims, RealMatrix &result, int seed )
{
  IntMatrix permutations;
  result.shapeUninitialized( num_dims, num_pts );
  get_permutations( permutations, num_pts, num_dims, seed );
  for ( int j = 0; j < num_pts; j++ )
    {
      for ( int d = 0; d < num_dims; d++ )
	{
	  result(d,j) = ( ( ( Real ) permutations(j,d) ) + 0.5 ) / 
	    ( ( Real ) num_pts );
	}
    }
};

Real mean( int n, Real *x )
{
  return sum( n, x ) / (Real)n;
};

Real variance( int n, Real *x, int dof )
{
  Real mu = mean( n, x );

  Real variance( 0.0 );
  for ( int i = 0; i < n; i++ )
    {
      variance += ( x[i] - mu ) * ( x[i] - mu );
    }
  return variance / (Real)( n - dof );
}

void get_permutations( IntMatrix &permutations, 
		       int M , int N, unsigned int seed )
{
  std::srand( seed );
  
  permutations.reshape( M, N );
  IntMatrix numbers;
  for ( int j = 0; j < N; j++ )
    {
      std::vector<int> random( M );
      for ( int i = 0; i < M; i++ ) 
	{
	  random[i] = i;
	}

      std::random_shuffle( random.begin(), random.end() );

      for ( int i = 0; i < M; i++ )
	{
	  permutations(i,j) = random[i];
	}
    }
};

void compute_next_combination( int num_dims, int level, IntVector &index, 
			       bool &extend, int &h, int &t )
{
   if ( !extend )
     {
       t = level;
       h = 0;
       index[0] = level;
       for ( int d = 1; d < num_dims; d++ ) 
	 index[d] = 0;
     }
   else
     {
       if ( 1 < t )
	 h = 0;

       t = index[h];
       index[h] = 0;
       index[0] = t - 1;
       index[h+1] = index[h+1] + 1;
       h += 1;
     }
   extend = ( index[num_dims-1] != level );
}

void compute_combinations( int num_dims, int level, IntMatrix &indices )
{
  int num_indices;
  if ( level > 0 )
    {
      num_indices = nchoosek( num_dims + level, num_dims ) - 
	nchoosek( num_dims + level-1, num_dims );
      indices.shapeUninitialized( num_indices, num_dims );
      bool extend = false; 
      int h = 0, t = 0; 
      IntVector index( num_dims ); // important this is initialized to zero
      int i = 0;
      while ( true )
	{
	  compute_next_combination( num_dims, level, index, 
				    extend, h, t );
	  for ( int d = 0; d < num_dims; d++ ) 
	    indices( i, d ) = index[d];
	  i++;

	  if ( !extend ) break;
	}
    }
  else
    {
      indices.shape( 1, num_dims );
    }
}

void compute_hyperbolic_level_subdim_indices( int num_dims, int level, 
					      int num_active_dims, 
					      Real p,
					      IntMatrix &indices )
{

  Real eps = 100 * std::numeric_limits<double>::epsilon();

  // assumes p <= 1
  //int max_num_indices =  nchoosek( num_dims + level, num_dims ) - 
  //  nchoosek( num_dims + level-1, num_dims );
  //indices.shapeUninitialized( num_active_dims, max_num_indices );
  indices.shapeUninitialized( num_active_dims, 1000 );
  int l = num_active_dims;
  int num_indices = 0;
  while ( true )
    {
      IntMatrix tmp;
      compute_combinations( num_active_dims, l, tmp );
      IntMatrix level_data( tmp, Teuchos::TRANS );
      for ( int i = 0; i < level_data.numCols(); i++ )
	{
	  IntVector index( Teuchos::View, level_data[i], num_active_dims );
	  if ( num_nonzeros( index ) == num_active_dims )
	    {
	      Real p_norm = pnorm( index, p );
	      if ( ( p_norm > level-1 + eps ) &&( p_norm < level + eps ) )
		{
		  if ( num_indices >= indices.numCols() )
		    indices.reshape( indices.numRows(), num_indices+1000 );
		  IntVector col( Teuchos::View, indices[num_indices],
				 num_active_dims );
		  col.assign( index );
		  num_indices++;
		}
	    }
	}
      l++;
      if ( l > level ) break;
    }
  // Remove unwanted memory
  indices.reshape( num_active_dims, num_indices );
  // Return row major matrix
  IntMatrix transpose( indices, Teuchos::TRANS );
  indices = transpose;
}

void compute_hyperbolic_level_indices( int num_dims, int level, Real p, 
				       IntMatrix &indices )
{
  if ( level == 0 )
    indices.reshape( 1, num_dims );
  else
    {
      indices.shapeUninitialized( num_dims, num_dims );
      for ( int  i = 0; i < num_dims; i++ )
	{
	  IntVector index( num_dims );
	  index[i] = level;
	  for ( int  d = 0; d < num_dims; d++ )
	    indices(d,i) = index[d];
	}

      for ( int d = 2; d < std::min( level+1, num_dims+1 ); d++ )
	{
	  IntMatrix level_comb;
	  compute_hyperbolic_level_subdim_indices( num_dims, level,
						   d, p, level_comb );
	  IntMatrix level_indices( level_comb, Teuchos::TRANS );
	  
	  if ( level_comb.numRows() == 0  )  break;

	  IntMatrix tmp;
	  compute_combinations( num_dims, d, tmp );
	  IntMatrix dim_indices( tmp.numCols(), tmp.numRows() ), 
	    dims_comb( tmp, Teuchos::TRANS );
	  int num_dim_indices = 0;
	  for ( int i = 0; i < dims_comb.numCols(); i++ )
	    {
	      IntVector index( Teuchos::View, dims_comb[i], num_dims );
	      if ( num_nonzeros( index ) == d )
		{
		  IntVector col( Teuchos::View, dim_indices[num_dim_indices],
				 num_dims );
		  col.assign( index );
		  num_dim_indices++;
		}
	    }
	  // Chop off unused memory;
	  dim_indices.reshape( num_dims, num_dim_indices );
	  
	  IntMatrix new_indices( num_dims, 
				 dim_indices.numCols()*level_indices.numCols() );
	  int num_new_indices = 0;
	  for ( int i = 0; i < dim_indices.numCols(); i++ )
	    {
	      IntVector dim_index( Teuchos::View, dim_indices[i], num_dims );
	      IntVector I;
	      nonzero( dim_index, I );
	      for ( int j = 0; j < level_indices.numCols(); j++ )
		{
		  IntVector index( Teuchos::View, new_indices[num_new_indices], 
				   indices.numRows() );
		  for ( int k = 0; k < I.length(); k++ )
		    {
		      index[I[k]] = level_indices(k,j);
		    }
		  num_new_indices++;
		}
	    }

	  column_append( new_indices, indices );
	}
      // Return row major matrix
      IntMatrix transpose( indices, Teuchos::TRANS );
      indices = transpose;
    } 
}

void compute_hyperbolic_indices( int num_dims, int level, Real p, 
				 IntMatrix &indices )
{
  compute_hyperbolic_level_indices( num_dims, 0, p, indices );
  for ( int l = 1; l < level+1; l++ )
    {
      IntMatrix level_indices;
      compute_hyperbolic_level_indices( num_dims, l, p, level_indices );
      row_append( level_indices, indices );
    }
}

void compute_anisotropic_hyperbolic_indices( int num_dims, int level, Real p, 
					     RealVector &weights,
					     IntMatrix &indices )
{
  compute_hyperbolic_level_indices( num_dims, 0, p, indices );
  int num_indices = indices.numRows();
  for ( int l = 1; l < level+1; l++ )
    {
      IntMatrix level_indices;
      compute_hyperbolic_level_indices( num_dims, l, p, level_indices );
      if ( num_indices + level_indices.numRows() >= indices.numRows() )
	indices.reshape( num_indices + level_indices.numRows(), num_dims );
      for ( int k = 0; k < level_indices.numRows(); k++ )
	{
	  Real weighted_level = 0.;
	  for ( int d = 0; d < num_dims; d++ )
	    weighted_level += std::pow((Real)level_indices(k,d),p)*weights[d];
	    //weighted_level += (Real)level_indices(k,d) * weights[d];
	  weighted_level = std::pow( weighted_level, 1.  / (Real)p );
	  if ( weighted_level <= (double)level )
	    {
	      for( int d = 0; d < num_dims; d++ )
		indices(num_indices,d) = level_indices(k,d);
	      num_indices++;
	    }	    
	}
    } 
  indices.reshape( num_indices, num_dims );
}

void get_column_norms( RealMatrix &A, RealVector &result )
{
  int M = A.numRows(), N = A.numCols();
  result.sizeUninitialized( N );
  for ( int i = 0; i < N; i++ )
    {
      RealVector col( Teuchos::View, A[i], M );
      result[i] = col.normFrobenius();
    }
}

} // namespace Pecos
