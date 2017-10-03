/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       DakotaPsuade
//- Description: Implementation code for the DakotaPsuade class
//- Owner:       Brian M. Adams, Sandia National Laboratories

#include "DakotaPsuade.H"
#include <cstddef>
#include <algorithm>
#include <vector>
#include <cstdlib>

static const char rcsId[]="@(#) $Id$";

// a lightweight definition of rand[0,N-1] that will work with random_shuffle
struct rand_for_random_shuffle_std {
  explicit rand_for_random_shuffle_std()
  { }

  ptrdiff_t operator()(ptrdiff_t upperBound)
  { return std::rand() % upperBound; }
};

/// Base class for the DakotaPsuade-provided RNG, with default implementation
class DakotaPsuadeBaseRNG
{
public:

  DakotaPsuadeBaseRNG() 
  { /* empty constructor */ }

  // TODO: probably don't want to do this due to random seed manipulation
  // in DAKOTA; TBD
  /// Constructor that accepts a seed
  DakotaPsuadeBaseRNG(int seed) 
  { std::srand(seed); }

  virtual ~DakotaPsuadeBaseRNG() 
  { }

  /// get random number in [0,1)
  virtual double uniform_0_1() 
  { return (double)std::rand()/(double)RAND_MAX; }

  /// return an RNG suitable for used with std::random_shuffle
  rand_for_random_shuffle_std& stl_rng() 
  { return stlRNG; };

protected:

  // copy/assign is NOT allowed
  DakotaPsuadeBaseRNG(const DakotaPsuadeBaseRNG&);
  DakotaPsuadeBaseRNG& operator=(const DakotaPsuadeBaseRNG&);

private:

  /// RNG compatible with STL
  rand_for_random_shuffle_std stlRNG;
};

#ifdef HAVE_BOOST

#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

// a lightweight definition of rand[0,N-1] that will work with random_shuffle
// Example taken from
// http://www.boost.org/doc/libs/1_40_0/libs/random/random_test.cpp
// (couldn't get boost::random_number_generator to work)
template <typename EngineT>
struct rand_for_random_shuffle {
  explicit rand_for_random_shuffle(EngineT &engine): m_engine(engine)
  { }

  template <typename IntT>
  IntT operator()(IntT upperBound) {
    if (upperBound == 1)
      return 0;
    typedef boost::uniform_int<IntT> distribution_type;
    typedef boost::variate_generator<EngineT &, distribution_type> generator_type;
    return generator_type(m_engine, distribution_type(0, upperBound - 1))();
  }

  EngineT &m_engine;
};

/// Specialized class for the DakotaPsuade-provided RNG, using Boost
class DakotaPsuadeBoostRNG: public DakotaPsuadeBaseRNG
{
public:

  DakotaPsuadeBoostRNG():
    rngSeed(41u), rnumGenerator(rngSeed), uniRealDist(0, 1), 
    uniRealVar(rnumGenerator, uniRealDist), stlRNG(rnumGenerator)
  { /* empty constructor */ }


  /// Constructor that accepts a seed
  DakotaPsuadeBoostRNG(int seed) :
    rngSeed(seed), rnumGenerator(rngSeed), uniRealDist(0, 1), 
    uniRealVar(rnumGenerator, uniRealDist), stlRNG(rnumGenerator)
  { /* empty constructor */ }

  /// get random number in [0,1)
  virtual double uniform_0_1()
  { return uniRealVar(); }

  /// return an RNG suitable for std::random_shuffle
  rand_for_random_shuffle<boost::mt19937>& stl_rng()
  { return stlRNG; }

private:

  /// seed for shared generator
  unsigned int rngSeed;
  /// shared Mersenne Twister generator
  boost::mt19937 rnumGenerator;
  
  /// distribution of U[0,1) random numbers
  boost::uniform_real<> uniRealDist;
  /// stream of U[0,1) random numbers
  boost::variate_generator<boost::mt19937&,boost::uniform_real<> > uniRealVar;

  /// RNG compatible with STL
  rand_for_random_shuffle<boost::mt19937> stlRNG;
};
#endif  // HAVE_BOOST


//
//- Implementation of DakotaPsuade
//

/** Constructor using default seed */
DakotaPsuade::DakotaPsuade():
#ifdef HAVE_BOOST
  rngPtr( new DakotaPsuadeBoostRNG() )
#else
  rngPtr( new DakotaPsuadeBaseRNG() )
#endif
{ /* empty constructor */ }

/** Constructor using DAKOTA-specified seed */
DakotaPsuade::DakotaPsuade(int seed):
#ifdef HAVE_BOOST
  rngPtr( new DakotaPsuadeBoostRNG(seed) )
#else
  rngPtr( new DakotaPsuadeBaseRNG(seed) )
#endif
{ /* empty constructor */ }

DakotaPsuade::~DakotaPsuade()
{ delete rngPtr; }

double DakotaPsuade::PSUADE_drand()
{ return rngPtr->uniform_0_1(); }

/** emulation of PSUADE's integer vector shuffler generateRandomIvector
    presumes permute has been allocated
    populates with [0:num_inputs-1] and permutes */
void DakotaPsuade::generateRandomIvector(int num_inputs, int *permute)
{
  // TODO: make more efficient by using original data in permute instead of 
  // copying
  std::vector<int> p;
  for (int i=0; i<num_inputs; i++) p.push_back(i);
  std::random_shuffle(p.begin(), p.end(), rngPtr->stl_rng());
  for (int i=0; i<num_inputs; i++) permute[i] = p[i];
}
