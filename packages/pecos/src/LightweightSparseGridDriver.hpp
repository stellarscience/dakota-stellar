/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 LightweightSparseGridDriver
//- Description: Lightweight sparse grid implementation for adapting index sets
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef LIGHTWEIGHT_SPARSE_GRID_DRIVER_HPP
#define LIGHTWEIGHT_SPARSE_GRID_DRIVER_HPP

#include "SparseGridDriver.hpp"

namespace Pecos {


/// Lighweight sparse grid driver class that generates index sets for
/// N-dimensional Smolyak sparse grids.

/** This class is used for basis adaptation using sparse grid index sets and
    omits many of the details needed for computing expectation integrals. */

class LightweightSparseGridDriver: public SparseGridDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  LightweightSparseGridDriver();
  /// constructor
  LightweightSparseGridDriver(unsigned short ssg_level,
			      const RealVector& dim_pref = RealVector(),
			      short growth_rate = MODERATE_RESTRICTED_GROWTH,
			      short refine_control = NO_CONTROL);
  /// destructor
  ~LightweightSparseGridDriver();

  //
  //- Heading: Virtual function redefinitionss
  //

  void initialize_sets();
  void push_trial_set(const UShortArray& set);
  void pop_trial_set();
  const UShortArray& trial_set() const;
  void print_smolyak_multi_index() const;

  //
  //- Heading: Member functions
  //

  /// initialize a lightweight instance that only generates index sets
  /** this lightweight sparse grid mode only generates index sets
      without points, weights, or collocation bookkeeping. */
  void initialize_grid(size_t num_v, unsigned short ssg_level);

  /// prune sets from oldMultiIndex and redefine activeMultiIndex
  void prune_sets(const SizetSet& save_tp);

  /// return smolyakMultiIndex
  const UShort2DArray& smolyak_multi_index() const;

private:

  //
  //- Heading: Data
  //

  /// numSmolyakIndices-by-numVars array for identifying the index to use
  /// within the polynomialBasis for a particular variable
  /** The index sets correspond to j (0-based) for use as indices, which
      are offset from the i indices (1-based) normally used in the Smolyak
      expressions.  The indices correspond to levels, one for each
      anisotropic tensor-product grid within a Smolyak recursion. */
  UShort2DArray smolyakMultiIndex;
};


inline LightweightSparseGridDriver::LightweightSparseGridDriver():
  SparseGridDriver()
{ }


inline LightweightSparseGridDriver::
LightweightSparseGridDriver(unsigned short ssg_level,
			    const RealVector& dim_pref,
			    short growth_rate, short refine_control):
  SparseGridDriver(ssg_level, dim_pref, growth_rate, refine_control)
{ }


inline LightweightSparseGridDriver::~LightweightSparseGridDriver()
{ }


inline void LightweightSparseGridDriver::push_trial_set(const UShortArray& set)
{ smolyakMultiIndex.push_back(set); }


inline void LightweightSparseGridDriver::pop_trial_set()
{ smolyakMultiIndex.pop_back(); }


inline const UShortArray& LightweightSparseGridDriver::trial_set() const
{ return smolyakMultiIndex.back(); }


inline void LightweightSparseGridDriver::print_smolyak_multi_index() const
{
  size_t i, sm_mi_len = smolyakMultiIndex.size();
  for (i=0; i<sm_mi_len; ++i) {
    PCout << "Smolyak index set " << i << ':';
    print_index_set(PCout, smolyakMultiIndex[i]);
  }
}


inline const UShort2DArray& LightweightSparseGridDriver::
smolyak_multi_index() const
{ return smolyakMultiIndex; }

} // namespace Pecos

#endif
