/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        RefinablePointSet
//- Description:  Virtual base class for grids used by HierarchicalBasis
//-               
//- Owner:        Christopher Miller: University of Maryland at College Park
//- Contact:      cmiller@math.umd.edu

#ifndef REFINABLEPOINTSET_HPP
#define REFINABLEPOINTSET_HPP

#include "pecos_data_types.hpp"

namespace Pecos {

/// Virtual base class for grids used by HierarchicalBasis.
/** Conceptually each grid point is specified by a [level,index] pair that 
    describes the location of the point.  The most important data member is the 
    std::vector<IntArray> int_index_to_level_index_map.  This vector asscociates 
    an integer index with the [level,index] pair that describes the location of 
    the associated point.  Subclasses need to implement 
    get_interp_point(IntArray) to associate an actual value with a [level,index]
    pair.  Subclasses also need to implement refine(BoolDeque) which packs the
    children of the points indicated in the BoolDeque into the map. Subclasses
    need to implement get_left_neighbor() and get_right_neighbor() which
    return the location of the hierarchical neighbors of a given point.

    For additional details see:
    For additional information see:
    X. Ma and N. Zabras, "An adaptive hierarchical sparse grid 
    collocation algorithm for the solution of stochastic 
    differential equations", J. Comput. Phys. 228:3084--3113,
    2009.
*/

class RefinablePointSet
{
public:

  //
  //- Heading: Constructor and Destructor
  //

  /// Standard constructor
  RefinablePointSet(size_t map_starting_capacity = 100);
    
  /// Destructor
  virtual ~RefinablePointSet();

  //
  //- Heading: Function declarations
  //

  /// Refines all of the points in the highest level.
  /** Calls refine(). */
  virtual void refine_all();

  /// Returns the point associated with the integer index.
  /**  The ordering of points is in the hierarchical sense
   i.e. ordering is lexiographic with lower levels first then left to right
   within a given level. */
  virtual const Real get_interp_point(unsigned int point) const;

  /// Returns the entire set of interpolation points ordered from left to right.
  virtual const RealArray& get_interp_points();

  /// Returns the points at the current highest level ordered from left to right
  virtual const RealArray& get_highest_level_points();

  /// Returns the current level.
  virtual const unsigned int get_current_level() const;

  /// Returns the [level/index] pair associated with a given index.
  virtual const IntArray& get_level_index_pair(unsigned int i) const;

  /// Returns the size of the current level.
  virtual const unsigned int get_current_level_size() const;

  /// Returns the total number of interpolation points
  virtual const unsigned int get_num_interp_points() const;

  /// Prints the hierarchical grid.
  friend std::ostream& 
  operator<<(std::ostream& ostr, const RefinablePointSet& pointSet);

  //
  //- Heading: Pure virtual function definitions
  //

  /// Refines the indicated points in the highest level.
  /** Refine handles the expansion of the int_index_to_level_index_map map
      each new point gets packed into the map.  Note that refine does NOT
      update the vector interpPts or highestLevelPoints.*/
  virtual void refine(const BoolDeque& points) = 0;

  /// Returns the nearest point in a lower level to the left of the given point 
  virtual const Real get_left_neighbor(const unsigned int idx) const = 0;
  
  /// Returns the nearest point in a lower level to the right of the given point
  virtual const Real get_right_neighbor(const unsigned int idx) const = 0;  

protected:

  //
  //- Heading: Protected member functions
  //

  /// Updates the vector of interpPts.
  /** Method resizes the interpPts vector, and loops through  
      int_index_to_level_index_map packing the associated points into the
      vector.  The vector is then sorted.

      Called by get_interp_points().*/
  virtual void update_interpPts();

  /// Updates highestLevelPoints
  virtual void update_highestLevelPoints();

  /// Returns the point associated with a given [level,index] pair.
  /** Note that
     this method does not care if the point is actually in the grid or not
     this function simply performs the translation (i.e. it could
     be static).  Information on what points are actually in the grid is
     provided by the int_index_to_level_index_map data member which
     associates an int in the range [0,num_interp_points) with its 
     corresponding [level,index] pair. */
  virtual const Real get_interp_point(const IntArray& point) const = 0;

  
  //
  //- Heading: Data
  //

  /// The current level of the grid.
  unsigned int current_level;
  /// The number of points at the current level.
  unsigned int current_level_size;
  /// The total number of interpolation points.
  unsigned int num_interp_points;
  /// Flag for checking if interpPts is current. 
  bool interpPointsValid;
  /// Flag for checking if highestLevelPoints is current.
  bool highestLevelPointsValid;
  /// Vector containing the interpolation points.
  RealArray interpPts;
  /// Vector containing the points at the current level.
  RealArray highestLevelPoints;

  /// Map associates integers with [level,index] pairs.  
  /** It is up to subclasses
   to be able to convert from a [level, index] to the actual point location.  
   The construction of this map should be handled by the refine() function.*/
  std::vector<IntArray> int_index_to_level_index_map;
  
  

private:

}; //End class definition

} //End namespace Pecos

#endif
