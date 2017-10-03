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

#include "UniformRefinablePointSet.hpp"

namespace Pecos{

UniformRefinablePointSet::
UniformRefinablePointSet(
  const Real left_end_,
  const Real right_end_,
  const unsigned int starting_level,
  size_t map_starting_capacity):
    RefinablePointSet(map_starting_capacity),
    left_end(left_end_),
    right_end(right_end_),
    midpoint(left_end + (right_end - left_end)/2.0)
{
  if( left_end >= right_end ){
    PCerr << "Error: left_end_ should be <= right_end_." 
	  << " Got left = " << left_end_ << " right = " << right_end_ 
	  << ". Throwing an exception." << std::endl;
    throw ( std::invalid_argument("") );
  }

  if( starting_level < 1 ){
    PCerr << "Error: starting_level shoud be >= 1." 
          << " Got starting_level = " << starting_level 
          << ". Throwing an exception." << std::endl;
    throw ( std::invalid_argument("") );
  }
  interpPts[0] = midpoint;
  highestLevelPoints[0] = midpoint;
  IntArray firstIndex(2);
  firstIndex[0] = 1;
  firstIndex[1] = 0;
  int_index_to_level_index_map[0] = firstIndex;
  while ( current_level != starting_level ) {
    refine_all();
  }
}

UniformRefinablePointSet::~UniformRefinablePointSet()
{
}

void UniformRefinablePointSet::refine(const BoolDeque& points)
{
  //check for improperly sized BoolDeque
  if ( points.size() != current_level_size ) {
    PCerr << "Error: Refinement selector vector 'points' needs to be the "
          << "same size as current_level_size. Got points.size() == " 
          << points.size() << " and current_level_size == "
          << current_level_size << "." << std::endl;
    throw ( std::invalid_argument("") );
  }

  const unsigned int currentGridSize = num_interp_points; 
  unsigned int masterIndex;
  IntArray levelIndex;
  IntArray childIndex1(2);
  IntArray childIndex2(2);
  //Each marked parent spawns two children (unless on boundary, handled later)
  unsigned int current_level_size_new = 2*
    (int) std::count(points.begin(),points.end(),true);
  int_index_to_level_index_map.resize(int_index_to_level_index_map.size() + 
				      current_level_size_new);
  for ( unsigned int idx = 0; idx < current_level_size; idx++ ) {
    if ( points[idx] ) {
      masterIndex = currentGridSize - current_level_size + idx;
      levelIndex = int_index_to_level_index_map[masterIndex];
      
      /* Need to handle special case where the boundary points at lvl
	 2 each only have one child */
      if ( current_level + 1 == 3 ) {
        childIndex1[0] = current_level + 1;
        childIndex1[1] = levelIndex[1];
        int_index_to_level_index_map[num_interp_points++] = childIndex1;
        current_level_size_new -= 1;
	int_index_to_level_index_map.resize(int_index_to_level_index_map.size() - 1);
      } else { //general case
        childIndex1[0] = current_level + 1;
        childIndex1[1] = 2*levelIndex[1];
        childIndex2[0] = current_level + 1;
        childIndex2[1] = 2*levelIndex[1] + 1;
        int_index_to_level_index_map[num_interp_points++] = childIndex1;
        int_index_to_level_index_map[num_interp_points++] = childIndex2;
      }
    }
  }
  interpPointsValid = false;
  highestLevelPointsValid = false;
  ++current_level;
  current_level_size = current_level_size_new;
}

  
const Real UniformRefinablePointSet::
get_interp_point(const IntArray& point) const
{
  //This seems bad but it's not so rough.  I handle
  //the level 1 and 2 cases manually.  The first
  //point in level 3 is half way between left_end
  //and the middle.  The first point in level 4
  //is a quarter way, level 5 an eighth and so on.
  //This is the offset of the first point from the
  //left end.  Subsequent points on the same level
  //are spaced in increments 2*offset.

  const unsigned int level = point[0];
  const unsigned int index = point[1];

  if( level == 1 ) return (left_end + right_end) / 2.0;
  else if( level == 2 ){
    if( index == 0 ) return left_end;
    else return right_end;
  } else {
    Real offset = ( (left_end + right_end) / 2.0 - left_end) 
      / std::pow(2.0,(double)level-2);
    return left_end + (1 + 2*index)*offset;
  }
}
 
const Real UniformRefinablePointSet::get_left_neighbor(const unsigned int idx) const
{
  const IntArray& levelIndex = get_level_index_pair(idx);
  const unsigned int level = levelIndex[0];
  const unsigned int index = levelIndex[1];
  if ( level == 1 ) {
    return left_end;
  }
  else if ( level == 2 ) {
    if ( index == 0 ) return left_end;
    else return midpoint;
  } else {
    Real offset = (midpoint - left_end)/std::pow(2.0,(double)level-2);
    return get_interp_point( levelIndex ) - offset;
  } 
}

const Real UniformRefinablePointSet::get_right_neighbor(const unsigned int idx) const{
    const IntArray& levelIndex = get_level_index_pair(idx);
  const unsigned int level = levelIndex[0];
  const unsigned int index = levelIndex[1];
  if ( level == 1 ) {
    return right_end;
  }
  else if ( level == 2 ) {
    if ( index == 0 ) return midpoint;
    else return right_end;
  } else {
    Real offset = (midpoint - left_end)/std::pow(2.0,(double)level-2);
    return get_interp_point( levelIndex ) + offset;
  }
}
  
} //End namespace Pecos
