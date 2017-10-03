/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_hierarchical_driver.cpp
    \brief A test program for HierarchicalLinearBasis class. */


#include "HierarchPWInterpPolynomial.hpp"
#include "UniformRefinablePointSet.hpp"
using namespace Pecos;
int main(int argc, char** argv)
{
  
  RefinablePointSet *pointSet = new UniformRefinablePointSet(-1,1);
  HierarchPWInterpPolynomial *a = new HierarchPWInterpPolynomial(*pointSet);

  if ( ( a->type1_value(0.1, 0) != 1 ) ||
       ( a->type1_value(-1.1,0) != 0 ) ){
    std::cout << "test failure:  basis[0] did not evaluate correctly."
              << " Got f(0.1) = " << a->type1_value(0.1, 0)
              << " and f(-1.1) = " << a->type1_value(-1.1, 0) 
              << std::endl;
    return EXIT_FAILURE;
  }

  pointSet->refine_all();
  if ( ( a->type1_value(-1, 1) != 1.0 ) ||
       ( a->type1_value(-.5,1) != 0.5 ) ||
       ( a->type1_value(2.0,1) != 0.0 ) ||
       ( a->type1_value(-1.1,1) != 0.0 ) ||
       ( a->type1_value(0.1,1) != 0.0 ) ){
    std::cout << "test failure:  basis[1] did not evaluate correctly."
              << " Got f(-1) = " << a->type1_value(-1, 1)
              << " and f(-.5) = " << a->type1_value(-.5, 1) 
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type1_value(1, 2) != 1.0 ) ||
       ( a->type1_value(.5,2) != 0.5 ) ||
       ( a->type1_value(0.0,2) != 0.0 ) ||
       ( a->type1_value(-.1,2) != 0.0 ) ||
       ( a->type1_value(1.1,2) != 0.0 ) ){
    std::cout << "test failure:  basis[2] did not evaluate correctly."
              << " Got f(1) = " << a->type1_value(1, 2)
              << " and f(.5) = " << a->type1_value(.5, 2) 
              << std::endl;
    return EXIT_FAILURE;
  }

  pointSet->refine_all();
  BoolDeque bitset(2,false);
  bitset[1] = true;
  pointSet->refine(bitset);

  if ( ( a->type1_value(.25, 5) != 1.0 ) ||
       ( a->type1_value(.5,5) != 0.0 ) ||
       ( a->type1_value(0.0,5) != 0.0 ) ||
       ( a->type1_value(.75,5) != 0.0 ) ||
       ( a->type1_value(.125,5) != 0.5 ) ||
       ( a->type1_value(.25 + .125,5) != .5 ) ||
       ( a->type1_value(-.1,5) != 0 ) ||
       ( a->type1_value(.51,5) != 0 ) ){
    std::cout << "test failure:  basis[5] did not evaluate correctly."
              << std::endl;
    return EXIT_FAILURE;
  } 

  // test gradient

  if ( ( a->type1_gradient(-.5,0) != 0 ) ||
       ( a->type1_gradient(.5,0) != 0 ) ||
       ( a->type1_gradient(1.5,0) != 0 ) ) {
    std::cout << "test failure: basis[0] did not evaluate gradient."
	      << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type1_gradient(-1, 1) != 0.0 ) ||
       ( a->type1_gradient(0.0,5) != 0.0 ) ||
       ( a->type1_gradient(-.5,1) != -1.0 ) ||
       ( a->type1_gradient(-1.1,1) != 0.0 ) ||
       ( a->type1_gradient(0.1,1) != 0.0 ) ){
    std::cout << "test failure:  basis[1] did not evaluate gradient."
	      << " Got f'(-1) == " << a->type1_gradient(-1, 1)
	      << ", f'(-.5) == " << a->type1_value(-.5,1)
              << ", f'(1.1) == " << a->type1_value(1.1,1)
              << ", f'(0.1) == " << a->type1_value(0.1,1)
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type1_gradient(1, 2) != 0.0 ) ||
       ( a->type1_value(0.0,5) != 0.0 ) ||
       ( a->type1_gradient(.5,2) != 1.0 ) ||
       ( a->type1_gradient(-.1,2) != 0.0 ) ||
       ( a->type1_gradient(1.1,2) != 0.0 ) ){
    std::cout << "test failure:  basis[2] did not evaluate gradient." 
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type1_gradient(.25, 5) != 0.0 ) ||
       ( a->type1_gradient(.5,5) != 0.0 ) ||
       ( a->type1_gradient(0.0,5) != 0.0 ) ||
       ( a->type1_gradient(.75,5) != 0.0 ) ||
       ( a->type1_gradient(.125,5) != 4 ) ||
       ( a->type1_gradient(.25 + .125,5) != -4 ) ||
       ( a->type1_gradient(-.1,5) != 0 ) ||
       ( a->type1_gradient(.51,5) != 0 ) ){
    std::cout << "test failure:  basis[5] did not evaluate gradient."
              << std::endl;
    return EXIT_FAILURE;
  }

  //test upcasting to a nodal basis.
  PiecewiseInterpPolynomial b = HierarchPWInterpPolynomial(*pointSet);
  /*  Will add this in once issue with boundary points is settled.
  if ( ( b.type1_value(-1, 0) != 1.0 ) ||
       ( b.type1_value(-.75,0) != 0.5 ) ||
       ( b.type1_value(-0.5,0) != 0.0 ) ||
       ( b.type1_value(-.5,2) != 0.0 ) ||
       ( b.type1_value(-.25,2) != .5 ) ||
       ( b.type1_value(0.0,2) != 1 ) ||
       ( b.type1_value(.125,2) != 0.5 ) ) {
    std::cout << "test failure:  upcasting to nodal basis unsuccessful. VALUE"
              << std::endl;
    return EXIT_FAILURE;
  }
  */

  delete a;

  // test cubic Hermite evaluations

  a = new HierarchPWInterpPolynomial(*pointSet, PIECEWISE_CUBIC_INTERP);
  
  if ( ( a->type1_value(0.1, 0) != 1 ) ||
       ( a->type1_value(-1.1,0) != 0 ) ){
    std::cout << "test failure:  basis[0] did not evaluate correctly."
              << " Got f(0.1) = " << a->type1_value(0.1, 0)
              << " and f(-1.1) = " << a->type1_value(-1.1, 0) 
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type1_value(-1, 1) != 1.0 ) ||
       ( a->type1_value(-.5,1) != 0.5 ) ||
       ( a->type1_value(0.0,1) != 0.0 ) ||
       ( a->type1_value(-1.1,1) != 0.0 ) ||
       ( a->type1_value(0.1,1) != 0.0 ) ){
    std::cout << "test failure:  basis[1] did not evaluate correctly."
              << " Got f(-1) = " << a->type1_value(-1, 1)
              << " and f(-.5) = " << a->type1_value(-.5, 1) 
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type1_value(1, 2) != 1.0 ) ||
       ( a->type1_value(0.0,2) != 0.0 ) ||
       ( a->type1_value(.5,2) != 0.5 ) ||
       ( a->type1_value(-.1,2) != 0.0 ) ||
       ( a->type1_value(1.1,2) != 0.0 ) ){
    std::cout << "test failure:  basis[2] did not evaluate correctly."
              << " Got f(1) = " << a->type1_value(1, 2)
              << " and f(.5) = " << a->type1_value(.5, 2) 
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type1_value(.25, 5) != 1.0 ) ||
       ( a->type1_value(.5,5) != 0.0 ) ||
       ( a->type1_value(0.0,5) != 0.0 ) ||
       ( a->type1_value(.75,5) != 0.0 ) ||
       ( a->type1_value(.125,5) != 0.5 ) ||
       ( a->type1_value(.25 + .125,5) != .5 ) ||
       ( a->type1_value(-.1,5) != 0 ) ||
       ( a->type1_value(.51,5) != 0 ) ){
    std::cout << "test failure:  basis[5] did not evaluate correctly."
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type2_value(0.1, 0) != 0 ) ||
       ( a->type2_value(-1.1,0) != 0 ) ){
    std::cout << "test failure:  basis[0] did not evaluate type 2 correctly."
              << " Got f(0.1) = " << a->type2_value(0.1, 0)
              << " and f(-1.1) = " << a->type2_value(-1.1, 0) 
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type2_value(-1, 1) != 0.0 ) ||
       ( a->type2_value(0.0,1) != 0.0 ) ||
       ( a->type2_value(-.5,1) != 1/8.0 ) ||
       ( a->type2_value(-1.1,1) != 0.0 ) ||
       ( a->type2_value(0.1,1) != 0.0 ) ){
    std::cout << "test failure:  basis[1] did not evaluate type 2 correctly."
              << " Got f(-1) = " << a->type2_value(-1, 1)
              << " and f(-.5) = " << a->type2_value(-.5, 1) 
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type2_value(1, 2) != 0.0 ) ||
       ( a->type2_value(0.0,2) != 0.0 ) ||
       ( a->type2_value(.5,2) != -1/8.0 ) ||
       ( a->type2_value(-.1,2) != 0.0 ) ||
       ( a->type2_value(1.1,2) != 0.0 ) ){
    std::cout << "test failure:  basis[2] did not evaluate type 2 correctly."
              << " Got f(1) = " << a->type2_value(1, 2)
              << " and f(.5) = " << a->type2_value(.5, 2) 
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type2_value(.25, 5) != 0.0 ) ||
       ( a->type2_value(.5,5) != 0.0 ) ||
       ( a->type2_value(0.0,5) != 0.0 ) ||
       ( a->type2_value(.75,5) != 0.0 ) ||
       ( a->type2_value(.125,5) != -1/32.0 ) ||
       ( a->type2_value(.25 + .125,5) != 1/32.0 ) ||
       ( a->type2_value(-.1,5) != 0 ) ||
       ( a->type2_value(.51,5) != 0 ) ){
    std::cout << "test failure:  basis[5] did not evaluate type 2 correctly."
              << " Got f(.25) == " << a->type2_value(.25, 5)
	      << ", got f(.5) == " << a->type2_value(.5, 5)
	      << ", got f(.75) == " << a->type2_value(.25, 5)
	      << ", got f(.125) == " << a->type2_value(.125, 5)
	      << ", got f(.25+.125) == " << a->type2_value(.25+.125, 5)
	      << ", got f(-.1) == " << a->type2_value(-.1, 5)
              << ", got f(.51) == " << a->type2_value(.51, 5)
              << std::endl;
    return EXIT_FAILURE;
  }

  // test cubic hermite gradient

  if ( ( a->type1_gradient(-.5,0) != 0 ) ||
       ( a->type1_gradient(.5,0) != 0 ) ||
       ( a->type1_gradient(1.5,0) != 0 ) ) {
    std::cout << "test failure: cubic_basis[0] did not evaluate type 1 gradient."
	      << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type1_gradient(-1, 1) != 0.0 ) ||
       ( a->type1_gradient(0.0,1) != 0.0 ) ||
       ( a->type1_gradient(-.5,1) != -1.5 ) ||
       ( a->type1_gradient(-1.1,1) != 0.0 ) ||
       ( a->type1_gradient(0.1,1) != 0.0 ) ){
    std::cout << "test failure:  cubic_basis[1] did not evaluate type 1 gradient."
	      << " Got f'(-1) == " << a->type1_gradient(-1, 1)
	      << ", f'(-.5) == " << a->type1_gradient(-.5,1)
              << ", f'(1.1) == " << a->type1_gradient(1.1,1)
              << ", f'(0.1) == " << a->type1_gradient(0.1,1)
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type1_gradient(1, 2) != 0.0 ) ||
       ( a->type1_gradient(0.0,2) != 0.0 ) ||
       ( a->type1_gradient(.5,2) != 1.5 ) ||
       ( a->type1_gradient(-.1,2) != 0.0 ) ||
       ( a->type1_gradient(1.1,2) != 0.0 ) ){
    std::cout << "test failure:  cubic_basis[2] did not evaluate type 1 gradient." 
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type1_gradient(.25, 5) != 0.0 ) ||
       ( a->type1_gradient(.5,5) != 0.0 ) ||
       ( a->type1_gradient(0.0,5) != 0.0 ) ||
       ( a->type1_gradient(.75,5) != 0.0 ) ||
       ( a->type1_gradient(.125,5) != 1.5*4 ) ||
       ( a->type1_gradient(.25 + .125,5) != -1.5*4 ) ||
       ( a->type1_gradient(-.1,5) != 0 ) ||
       ( a->type1_gradient(.51,5) != 0 ) ){
    std::cout << "test failure:  cubic_basis[5] did not evaluate type 1 gradient."
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type2_gradient(-.5,0) != 0 ) ||
       ( a->type2_gradient(.5,0) != 0 ) ||
       ( a->type2_gradient(1.5,0) != 0 ) ) {
    std::cout << "test failure: cubic_basis[0] did not evaluate type 2 gradient."
	      << " Got f(-.5) = " << a->type2_gradient(-.5,0)
              << " Got f(.5) = " << a->type2_gradient(.5,0)
              << " Got f(1.5) = " << a->type2_gradient(1.5,0)
	      << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type2_gradient(-1, 1) != 1.0 ) ||
       ( a->type2_gradient(0.0,1) != 0.0 ) ||
       ( a->type2_gradient(-.5,1) != -.25 ) ||
       ( a->type2_gradient(-1.1,1) != 0.0 ) ||
       ( a->type2_gradient(0.1,1) != 0.0 ) ){
    std::cout << "test failure:  cubic_basis[1] did not evaluate type 2 gradient."
	      << " Got f'(-1) == " << a->type2_gradient(-1, 1)
	      << ", f'(-.5) == " << a->type2_gradient(-.5,1)
              << ", f'(1.1) == " << a->type2_gradient(1.1,1)
              << ", f'(0.1) == " << a->type2_gradient(0.1,1)
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type2_gradient(0, 2) != 0.0 ) ||
       ( a->type2_gradient(0.0,2) != 0.0 ) ||
       ( a->type2_gradient(.5,2) != -.25 ) ||
       ( a->type2_gradient(-.1,2) != 0.0 ) ||
       ( a->type2_gradient(1,2) != 1.0 ) ){
    std::cout << "test failure:  cubic_basis[2] did not evaluate type 2 gradient."
	      << " Got f'(0) == " << a->type2_gradient(0, 1)
	      << ", f'(.5) == " << a->type2_gradient(.5,1)
              << ", f'(-.1) == " << a->type2_gradient(-.1,1)
              << ", f'(1) == " << a->type2_gradient(1,1)
              << std::endl;
    return EXIT_FAILURE;
  }

  if ( ( a->type2_gradient(.25, 5) != 1.0 ) ||
       ( a->type2_gradient(.5,5) != 0.0 ) ||
       ( a->type2_gradient(0.0,5) != 0.0 ) ||
       ( a->type2_gradient(.75,5) != 0.0 ) ||
       ( a->type2_gradient(.125,5) != -.25 ) ||
       ( a->type2_gradient(.25 + .125,5) != -.25 ) ||
       ( a->type2_gradient(-.1,5) != 0 ) ||
       ( a->type2_gradient(.51,5) != 0 ) ){
    std::cout << "test failure:  cubic_basis[5] did not evaluate type 2 gradient."
              << std::endl;
    return EXIT_FAILURE;
  }


  delete a;
    

  return EXIT_SUCCESS;
  
}
