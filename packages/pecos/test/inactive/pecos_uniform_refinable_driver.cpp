/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_uniform_refinable_driver.cpp
    \brief A test program for UniformRefinablePointSet class. */

#include "UniformRefinablePointSet.hpp"

using namespace Pecos;
int main(int argc, char** argv)
{

  RefinablePointSet *a, *b, *c;
  // Test constructor
  a = new UniformRefinablePointSet();
  b = new UniformRefinablePointSet(-1.2,2.4,3);
  bool caughtException1 = false;
  bool caughtException2 = false;

  // constructors should throw exceptions if left >= right...
  try{
    c = new UniformRefinablePointSet(-1.2,-1.2,1);
  }catch ( std::invalid_argument& ex ){
    caughtException1 = true;
  }

  // or if the desired level <= 1.
  try{
    c = new UniformRefinablePointSet(-1.2, 1.2, 0);
  }catch ( std::invalid_argument& ex ){
    caughtException2 = true;
  }

  if ( !caughtException1 ) {
    std::cout << "UniformRefinablePointSet(-1.2,-1.2,1) should have thrown"
              << " an exception." << std::endl;
    return EXIT_FAILURE;
  }

  if ( !caughtException2 ) {
    std::cout << "UniformRefinablePointSet(-1.2,1.2,0) should have thrown" 
              << " an exception." << std::endl;
    return EXIT_FAILURE;
  }

  caughtException1 = false;
  caughtException2 = false;

  if ( b->get_current_level() != 3 ) {
    std::cout << "Constructor failed to set b to a third level grid." 
              << std::endl;
    return EXIT_FAILURE;
  }

  // Default constructor should give a grid with one point = 0.5.
  RealArray a_interp_points = a->get_interp_points();
  RealArray a_highest_interp_points = a->get_highest_level_points();

  if ( ( a_interp_points != a_highest_interp_points ) ||
       ( a_interp_points[0] != 0.5 ) ||
       ( a_interp_points.size() != 1 ) ||
       ( a->get_current_level() != 1 ) ) {
    std::cout << "Constructor failed to set a to a first level grid on [0,1]"
              << std::endl;
    std::cout << "a_interp_points != a_highest_interp_points: "
              << (a_interp_points != a_highest_interp_points) << std::endl;
    std::cout << "a_interp_points[0] =  "
              << a_interp_points[0] << std::endl;
    std::cout << "a_interp_points.size() =  "
              << a_interp_points.size() << std::endl;
    std::cout << "a->get_current_level() =  "
              <<  a->get_current_level()  << std::endl;
    return EXIT_FAILURE;
  }

  delete a, b;

  // Let's try a different interval.  Also has a midpoint of 0.5.
  a = new UniformRefinablePointSet(-1,2);
  a_interp_points = a->get_interp_points();
  a_highest_interp_points = a->get_highest_level_points();

  if ( ( a_interp_points != a_highest_interp_points ) ||
       ( a_interp_points[0] != 0.5 ) ||
       ( a_interp_points.size() != 1 ) ||
       ( a->get_current_level() != 1 ) ) {
    std::cout << "Constructor failed to set a to a first level grid on [-1,2]"
              << std::endl;
    return EXIT_FAILURE;
  }

  // Now lets do some refinements.
  a->refine_all();

  RealArray correct_interp_pts(3);
  correct_interp_pts[0] = -1;
  correct_interp_pts[1] = .5;
  correct_interp_pts[2] = 2;

  RealArray correct_highest_interp_pts(2);
  correct_highest_interp_pts[0] = -1;
  correct_highest_interp_pts[1] = 2;
  // Highest level should contain just the endpoints of [-1,2]
  a_highest_interp_points = a->get_highest_level_points();
  // and the whole thing should be {-1,.5,2}
  a_interp_points = a->get_interp_points();

  if ( ( a_interp_points != correct_interp_pts ) ||
       ( a_highest_interp_points != correct_highest_interp_pts ) ||
       ( a->get_current_level() != 2 ) ) {
    std::cout << "a not refined to level 2 correctly." << std::endl;
    bool check1 = ( a_interp_points != correct_interp_pts );
    bool check2 = ( a_highest_interp_points != correct_highest_interp_pts );
    bool check3 = ( a->get_current_level() != 2 );
    return EXIT_FAILURE;
  }

  a->refine_all();

  correct_interp_pts.resize(5);
  correct_interp_pts[0] = -1;
  correct_interp_pts[1] = -.25;
  correct_interp_pts[2] = .5;
  correct_interp_pts[3] = 1.25;
  correct_interp_pts[4] = 2;

  correct_highest_interp_pts[0] = -.25;
  correct_highest_interp_pts[1] = 1.25;
  // Highest level should be {-.25, 1.25}
  a_highest_interp_points = a->get_highest_level_points();
  // and the whole thing should be {-1,-.25,.5,1.25,2}
  a_interp_points = a->get_interp_points();

  if ( ( a_interp_points != correct_interp_pts ) ||
       ( a_highest_interp_points != correct_highest_interp_pts ) ||
       ( a->get_current_level() != 3 ) ) {
    std::cout << "a not refined to level 3 correctly." << std::endl;
    return EXIT_FAILURE;
  }

  // One more uniform refinement and then a bitset refinement.
  a->refine_all();
  BoolDeque selective_refinement(4,false);
  selective_refinement[0] = true;
  selective_refinement[2] = true;
  a->refine(selective_refinement);

  correct_interp_pts.resize(13);
  correct_interp_pts[0] = -1;
  correct_interp_pts[1] = -.8125;
  correct_interp_pts[2] = -.625;
  correct_interp_pts[3] = -.4375;
  correct_interp_pts[4] = -.25;
  correct_interp_pts[5] = .125;
  correct_interp_pts[6] = .5;
  correct_interp_pts[7] = .6875;
  correct_interp_pts[8] = .875;
  correct_interp_pts[9] = 1.0625;
  correct_interp_pts[10] = 1.25;
  correct_interp_pts[11] = 1.625;
  correct_interp_pts[12] = 2;

  correct_highest_interp_pts.resize(4);
  correct_highest_interp_pts[0] = -.8125;
  correct_highest_interp_pts[1] = -.4375;
  correct_highest_interp_pts[2] = .6875;
  correct_highest_interp_pts[3] = 1.0625;

  a_highest_interp_points = a->get_highest_level_points();
  a_interp_points = a->get_interp_points();

  if ( ( a_interp_points != correct_interp_pts ) ||
       ( a_highest_interp_points != correct_highest_interp_pts ) ||
       ( a->get_current_level() != 5 ) ) {
    std::cout << "a not refined to level 5 correctly." << std::endl;
    return EXIT_FAILURE;
  }
  //Try to refine with an improperly sized booldeque.
  selective_refinement.resize(a_highest_interp_points.size() + 1);
  try {
    a->refine(selective_refinement);
  } catch ( std::invalid_argument& ex ){
    caughtException1 = true;
  }

  if ( !caughtException1 ) {
    std::cout << "a->refine() should have thrown an exception when called with"
              << " an improperly sized BoolDeque." << std::endl;
    return EXIT_FAILURE;
  }

  // test get_interp_point
  if ( ( a->get_interp_point(0) != .5 ) ||
       ( a->get_interp_point(1) != -1 ) ||
       ( a->get_interp_point(2) != 2 ) ||
       ( a->get_interp_point(6) != .125 ) ||
       ( a->get_interp_point(12) != 1.0625 ) ) {
    std::cout << "a->get_interp_points isn't working" << std::endl;
    return EXIT_FAILURE;
  }

  // test get_level_index_pair
  const IntArray& levelIndex = a->get_level_index_pair(11);
  if ( ( levelIndex[0] != 5 ) || ( levelIndex[1] != 4 ) ) {
    std::cout << "a->get_level_index_pair isn't working" << std::endl;
    return EXIT_FAILURE;
   }

  // test get_left_neighbor get_right_neighbor
  if ( ( a->get_left_neighbor(0) != -1 ) || 
       ( a->get_right_neighbor(0) != 2 ) ) {
    std::cout << "get_left_neighbor, get_right_neighbor not working on "
              << "root node." << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->get_left_neighbor(1) != -1 ) ||
       ( a->get_right_neighbor(1) != .5 ) ) {
    std::cout << "get_left_neighbor, get_right_neighbor not working on "
              << "left boundary node." << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->get_left_neighbor(2) != .5 ) ||
       ( a->get_right_neighbor(2) != 2 ) ) {
    std::cout << "get_left_neighbor, get_right_neighbor not working on "
              << "right boundary node." << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->get_left_neighbor(4) != .5 ) ||
       ( a->get_right_neighbor(4) != 2 ) ) {
    std::cout << "get_left_neighbor, get_right_neighbor not working on "
              << "node 4." << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->get_left_neighbor(10) != -.625 ) ||
       ( a->get_right_neighbor(10) != -.25 ) ) {
    std::cout << "get_left_neighbor, get_right_neighbor not working on "
              << "node 10." << std::endl;
    return EXIT_FAILURE;
  }

  //grid printing test
  std::cout << *a << std::endl;

  delete a;
    
  return EXIT_SUCCESS;

}
