/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       DakotaPsuade
//- Description: This class provides necessary interface for DAKOTA with PSUADE
//- Owner:       Brian M. Adams, Sandia National Laboratories

#ifndef DAKOTA_PSUADE_H
#define DAKOTA_PSUADE_H

#define PSUADE_UNDEFINED -DBL_MAX

class DakotaPsuadeBaseRNG; // forward declaration

/// glue class providing RNG services to PSUADE
class DakotaPsuade
{

public:

  DakotaPsuade();

  /// constructor that accepts a random seed
  DakotaPsuade(int seed);

  ~DakotaPsuade();

  // for safety for now, only allowing element-wise access to sampleMatrix
  /// return the i,j entry of the sample matrix
  double sample_matrix(int i, int j);

protected:

  /// disallow copy/assignment due to pointer
  DakotaPsuade(const DakotaPsuade&);
  DakotaPsuade& operator=(const DakotaPsuade&);

  //-- Member functions

  /// could contain initializer for RNG -- not used (as named in MOAT*)
  void PSUADE_randInit() { /* no-op */ }

  /// generate a random number in [0,1) (half open) (as named in MOAT*) 
  double PSUADE_drand();

  /// create a vector with integers [0:num_inputs-1] and permute
  /// (as named in MOAT*)
  void generateRandomIvector(int num_inputs, int *permute);

  //-- Data

  /// pointer to a random number generator engine (pImpl pattern)  
  DakotaPsuadeBaseRNG *rngPtr;

  /// matrix of samples to evaluate (name as given in client MOAT* files)
  double **sampleMatrix_;

};

inline double DakotaPsuade::sample_matrix(int i, int j)
{ return sampleMatrix_[i][j]; }

#endif
