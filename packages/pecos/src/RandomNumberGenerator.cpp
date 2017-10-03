#include "RandomNumberGenerator.hpp"

namespace Pecos {

RandomNumberGenerator::RandomNumberGenerator()
{};

RandomNumberGenerator::~RandomNumberGenerator()
{};

void RandomNumberGenerator::uniform( IntMatrix &A, int M, int N, 
				     unsigned int seed )
{ 
  const int range_min = 0; 
  const int range_max = 1; 
  typedef boost::uniform_int<> NumberDistribution; 
  typedef boost::mt19937 RNG; 
  typedef boost::variate_generator<RNG&,NumberDistribution> Generator; 
 
  NumberDistribution distribution( range_min, range_max ); 
  RNG generator; 
  Generator numberGenerator( generator, distribution ); 
  generator.seed( seed ); // seed with the current time 
  
  A.reshape( M, N );
  for ( int i = 0; i < M; i++ )
    {
      for ( int j = 0; j < N; j++ )
	{
	  A(i,j) = numberGenerator();
	}
    }
};

void RandomNumberGenerator::uniform( RealMatrix &A, int M, int N, 
				     unsigned int seed )
{ 
  const Real range_min = 0.0; 
  const Real range_max = 1.0; 
  typedef boost::uniform_real<> NumberDistribution; 
  typedef boost::mt19937 RNG; 
  typedef boost::variate_generator<RNG&,NumberDistribution> Generator; 
 
  NumberDistribution distribution( range_min, range_max ); 
  RNG generator; 
  Generator numberGenerator( generator, distribution ); 
  generator.seed( seed ); // seed with the current time 
  
  A.reshape( M, N );
  for ( int i = 0; i < M; i++ )
    {
      for ( int j = 0; j < N; j++ )
	{
	  A(i,j) = numberGenerator();
	}
    }
};

void RandomNumberGenerator::gaussian( RealMatrix &A, int M, int N, 
				      Real mean, Real variance, 
				      unsigned int seed )
{ 
  const Real mu = mean; 
  const Real sigma2 = variance; 
  typedef boost::normal_distribution<Real> NumberDistribution; 
  typedef boost::mt19937 RNG; 
  typedef boost::variate_generator<RNG&,NumberDistribution> Generator; 
 
  NumberDistribution distribution( mu, sigma2 ); 
  RNG generator; 
  generator.seed( seed ); // seed with the current time 
  Generator numberGenerator( generator, distribution ); 

  A.reshape( M, N );
  for ( int i = 0; i < M; i++ )
    {
      for ( int j = 0; j < N; j++ )
	{
	  A(i,j) = numberGenerator();
	}
    }
};

void RandomNumberGenerator::beta( RealMatrix &A, int M, int N, 
				  Real alpha, Real beta, 
				  unsigned int seed )
{ 
  const Real range_min = 0; 
  const Real range_max = 1; 
  typedef boost::uniform_real<> NumberDistribution; 
  typedef boost::mt19937 RNG; 
  typedef boost::variate_generator<RNG&,NumberDistribution> Generator; 
 
  NumberDistribution distribution( range_min, range_max ); 
  RNG generator; 
  Generator numberGenerator( generator, distribution ); 
  generator.seed( seed ); // seed with the current time 
  
  boost::math::beta_distribution<> betaDistribution( alpha, beta );

  A.reshape( M, N );
  for ( int i = 0; i < M; i++ )
    {
      for ( int j = 0; j < N; j++ )
	{
	  // Generate Random Number on [0,1)
	  Real u = numberGenerator();
	  // Evaluate Beta CDF at u;
	  A(i,j) = boost::math::quantile( betaDistribution, u );
	}
    }
};

void RandomNumberGenerator::permutation( IntMatrix &permutations, 
					 int M, int N,
					 unsigned int seed )
{
  boost::mt19937 generator( seed );
  boost::uniform_int<> dist( 0, M-1 );
  permutations.shapeUninitialized( M, N );
  for ( int j = 0; j < N; j++ )
    {
      for ( int i = 0; i < M; i++ )
	permutations(i,j) = i;

      for ( int i = 0; i < M; i++ )
	{
	  int index = dist( generator );
	  int tmp = permutations(i,j);
	  permutations(i,j) = permutations(index,j);
	  permutations(index,j) = tmp;
	}
    }
};

}
