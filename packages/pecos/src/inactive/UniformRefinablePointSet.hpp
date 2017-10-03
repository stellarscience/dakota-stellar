/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        UniformRefinablePointSet
//- Description:  Class implementing a refinable point set using uniformly
//	              spaced points on [a,b]
//-               
//- Owner:        Christopher Miller: University of Maryland at College Park
//- Contact:      cmiller@math.umd.edu

#ifndef UNIFORMREFINABLEPOINTSET_HPP
#define UNIFORMREFINABLEPOINTSET_HPP

#include "RefinablePointSet.hpp"

namespace Pecos {

/// Class implementing RefinablePointSet for uniformly spaced points on [a,b]. 

class UniformRefinablePointSet : public RefinablePointSet
{
public:

  //
  //- Heading: Constructor and Destructor
  //

  /// Standard constructor
  UniformRefinablePointSet(const Real left_end_ = 0, 
                           const Real right_end_ = 1,
                           const unsigned int starting_level = 1,
			   size_t map_starting_capacity = 100);
    
  /// Destructor
  virtual ~UniformRefinablePointSet();

  //
  //- Heading: Virtual function redefinitions
  //

  /// Implements RefinablePointSet::refine().
  virtual void refine(const BoolDeque& points);
  /// Implements RefinablePointSet::get_left_neighbor().
  virtual const Real get_left_neighbor(const unsigned int idx) const;
  /// Implements RefinablePointSet::get_right_neighbor().
  virtual const Real get_right_neighbor(const unsigned int idx) const;
	
protected:
  
  //
  //- Heading: Protected member functions
  //
  
  /// Implements RefinablePointSet::get_interp_point().
  virtual const Real get_interp_point(const IntArray& point) const;  

  //
  //- Heading: Data
  //
  
  /// Left endpoint of interpolation.
  Real left_end;
  /// Right endpoint of interpolation.
  Real right_end;
  /// (left_end + right_end)/2
  Real midpoint;
  

private:

}; //End class definition

} //End namespace Pecos

#endif
