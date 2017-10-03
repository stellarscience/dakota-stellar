/*  ______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */
 
//- Class:       BoostRNG_Monostate
//- Description: Wrapper for various implementations of Random Number Generators
//- Owner:       Laura Swiler, Dave Gay, and Bill Bohnhoff 
//- Checked by:
//- Version: $Id$

#ifndef BOOST_RNG_MONOSTATE_HPP
#define BOOST_RNG_MONOSTATE_HPP

#include "pecos_data_types.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>


// BoostRNG_Monostate helps replace the default LHS rnumlhs1/rnumlhs2
// RNGs. In standalone LHS, these functions are provided by LHS (see 
// Lhs_drv.f90), but BoostRNG_Monostate provides them when LHS is compiled as a
// library for use by Dakota. This permits user-selectable RNG, either LHS
// default or mt19937 provided by Boost.
//
// Summary of random number generator replacement call chain:
//
// * LHS calls to RNG call RNUMLHS1, RNUMLHS2, which return double
//
// * BoostRNG_Monostate defines global symbols RNUMLHS1, RNUMLHS2 which
//   forward calls to BoostRNG_Monostate::randomNum(), randomNum2().
//   These replace those provided by LHS, so symbol resolution at link
//   time remains important.
//
// * randomNum/randomNum2 default to point to
//   BoostRNG_Monostate::mt19937 (Boost mt19937), but can be overridden by 
//   LHSDriver to point to defaultrnum1, defaultrnum2, provided by LHS
//


namespace Pecos {

/// RNG function pointer type taking no arguments and returning Real
typedef Real (*Rfunc)();

/** The monostate BoostRNG class provides static pointers randomNum
    and randomNum2 that at runtime will either point to
    BoostRNG_Monostate::mt19937, or be overridden to point to
    defaultnum1 / defaultrnum2.  At construct time they point to the
    BoostRNG_Monostate generators. */
class BoostRNG_Monostate
{
private:

  //
  //- Heading: Data
  //

  static unsigned int rngSeed;
  static boost::mt19937 rnumGenerator;
  static boost::uniform_real<> uniDist;
  static boost::variate_generator<boost::mt19937&,boost::uniform_real<> > uniMT;

public:

  //
  //- Heading: Member functions
  //

  /// set randomSeed
  static void seed(unsigned int rng_seed);
  /// return randomSeed
  static unsigned int seed();

  static Real mt19937();

  //
  //- Heading: Data
  //

  /// Cached function pointers
  static Rfunc randomNum;
  static Rfunc randomNum2;
};


inline void BoostRNG_Monostate::seed(unsigned int rng_seed)
{ rngSeed = rng_seed; rnumGenerator.seed(rng_seed); }

inline unsigned int BoostRNG_Monostate::seed()
{ return rngSeed; }

inline Real BoostRNG_Monostate::mt19937()
{ return uniMT(); }

unsigned int BoostRNG_Monostate::rngSeed(41u); // 41 used in the Boost examples

boost::mt19937 BoostRNG_Monostate::rnumGenerator( BoostRNG_Monostate::seed() );

boost::uniform_real<> BoostRNG_Monostate::uniDist(0, 1);

boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
  BoostRNG_Monostate::uniMT(BoostRNG_Monostate::rnumGenerator,
                            BoostRNG_Monostate::uniDist);
 
Real (*BoostRNG_Monostate::randomNum)()  = BoostRNG_Monostate::mt19937;
Real (*BoostRNG_Monostate::randomNum2)() = BoostRNG_Monostate::mt19937;
} // namespace Pecos


// This section defines Fortran 90 functions rnumlhs1 and rnumlhs2,
// which forward to the BoostRNG_Monostate randomNum and randomNum2
// functions (which may further delegate to Boost or back to LHS
// defaultrnum1/defaultrnum2).  LHS will call Fortran-mangled versions of
// rnumlhs1 and rnumlhs2 throughout its code.

#ifdef HAVE_LHS
#include "LHS.h"
#define rnum1 LHS_GLOBAL(rnumlhs1,RNUMLHS1)
#define rnum2 LHS_GLOBAL(rnumlhs2,RNUMLHS2)
#endif  // HAVE_LHS

extern "C" Pecos::Real rnum1(void), rnum2(void);

// inline leads to DIFFs in the Regression suite
/* inline */ Pecos::Real rnum1(void)
{
  //PCout << "running Boost MT" << "\n";
  return Pecos::BoostRNG_Monostate::randomNum();
}

/* inline */ Pecos::Real rnum2(void)
{
  // clone of rnum1
  //PCout << "running Boost MT" << "\n";
  return Pecos::BoostRNG_Monostate::randomNum2();
}

#endif  // BOOST_RNG_MONOSTATE_HPP

