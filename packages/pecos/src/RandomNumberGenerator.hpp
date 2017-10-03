/**
 * \file RandomNumberGenerator.hpp
 * \brief Easy to use wrapper to the Boost psuedo random number generator
 * \author John D. Jakeman
 */


#ifndef RAND_NUMBER_GENERATOR_HPP
#define RAND_NUMBER_GENERATOR_HPP

#include <boost/random.hpp> 
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions/beta.hpp>

#include "MathTools.hpp"

namespace Pecos {

/**
 * \class RandomNumberGenerator
 * \brief Easy to use wrapper to the Boost psuedo random number generator 
 */
class RandomNumberGenerator 
{
 public:
  /**
   * Default constructor
   */
  RandomNumberGenerator();

  /**
   * Deconstructor
   */
  ~RandomNumberGenerator();
  
  void uniform( IntMatrix& A, int M, int N, unsigned int seed );

  void uniform( RealMatrix& A, int M, int N, unsigned int seed );

  void gaussian( RealMatrix& A, int M, int N, Real mean, Real variance,
		 unsigned int seed );

  void beta( RealMatrix& A, int M, int N, Real alpha, Real beta, 
	     unsigned int seed );

  void permutation( IntMatrix &result, int M, int N, unsigned int seed );
};

}
#endif
