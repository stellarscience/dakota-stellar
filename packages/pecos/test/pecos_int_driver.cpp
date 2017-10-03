/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_int_driver.cpp
    \brief A driver program for PECOS */

#include "CombinedSparseGridDriver.hpp"
#include "TensorProductDriver.hpp"
#include "CubatureDriver.hpp"
//#include "LocalRefinableDriver.hpp"


/// A driver program for PECOS.

/** Generates sparse, tensor, and cubature grids for numerical integration. */

int main(int argc, char* argv[])
{
  std::cout << "Instantiating basis:\n";
  size_t i, num_vars = 4;
  std::vector<Pecos::BasisPolynomial> poly_basis(num_vars);
  for (i=0; i<num_vars; ++i)
    poly_basis[i] = Pecos::BasisPolynomial(Pecos::HERMITE_ORTHOG);
  Pecos::RealMatrix variable_sets;

  // Smolyak sparse grid
  std::cout << "Instantiating CombinedSparseGridDriver:\n";
  unsigned short level = 3; //Pecos::RealVector dimension_pref; // isotropic
  Pecos::CombinedSparseGridDriver csg_driver(level);//, dimension_pref);
  csg_driver.initialize_grid(poly_basis);
  csg_driver.compute_grid(variable_sets);
  std::cout << "Sparse grid points:\n";
  Pecos::write_data(std::cout, variable_sets, false, true, true);
  std::cout << "Sparse grid weights:\n";
  Pecos::write_data(std::cout, csg_driver.type1_weight_sets());

  // Tensor-product quadrature
  std::cout << "Instantiating TensorProductDriver:\n";
  Pecos::UShortArray quad_order(num_vars, 3);
  Pecos::TensorProductDriver tp_driver(quad_order);
  tp_driver.initialize_grid(poly_basis);
  tp_driver.compute_grid(variable_sets);
  std::cout << "Tensor grid points:\n";
  Pecos::write_data(std::cout, variable_sets, false, true, true);
  std::cout << "Tensor grid weights:\n";
  Pecos::write_data(std::cout, tp_driver.type1_weight_sets());

  // Cubature
  std::cout << "Instantiating CubatureDriver:\n";
  unsigned short integrand_order = 5;
  Pecos::CubatureDriver c_driver(integrand_order);
  c_driver.initialize_grid(poly_basis);
  c_driver.compute_grid(variable_sets);
  std::cout << "Cubature points:\n";
  Pecos::write_data(std::cout, variable_sets, false, true, true);
  std::cout << "Cubature weights:\n";
  Pecos::write_data(std::cout, c_driver.type1_weight_sets());

  /*
  // Local refinable grid
  std::cout << "Instantiating LocalRefinableDriver" << std::endl;
  Pecos::LocalRefinableDriver l_driver;
  unsigned short refinable_level = 4;
  Pecos::Real2DArray bounds(2);
  bounds[0] = Pecos::RealArray(2,-1);
  bounds[1] = Pecos::RealArray(2,1);

  // Check correct initialization
  l_driver.initialize_grid(bounds[0],bounds[1],refinable_level,
			   Pecos::PIECEWISE_LINEAR_INTERP);
  l_driver.compute_grid(variable_sets);
  assert( l_driver.get_current_level() == 4 );
  assert( l_driver.grid_size() == 29 );
  assert( l_driver.highest_level_size() == 16 );

  // Attempt local refinement
  Pecos::BoolDeque refinement_selector(16,false);
  refinement_selector[5] = true;
  refinement_selector[10] = true;
  refinement_selector[15] = true;
  l_driver.refine_locally(refinement_selector);
  assert( l_driver.get_current_level() == 5 );
  assert( l_driver.grid_size() == 40 );
  assert( l_driver.highest_level_size() == 11 );
  Pecos::const_ColPtIterator highestLevelPoints = 
    l_driver.get_highest_level_start();
  Pecos::RealVector pointCheck(2);
  pointCheck[0] = -.375;
  pointCheck[1] = 0.0;
  assert( *highestLevelPoints == pointCheck );
  pointCheck[0] = -.125;
  pointCheck[1] = 0.0;
  assert( *(highestLevelPoints+3) == pointCheck );
  pointCheck[0] = 0.5;
  pointCheck[1] = 0.5;
  assert( *(highestLevelPoints+8) == pointCheck );

  std::ostringstream output_check;
  output_check << l_driver;
  assert( output_check.good() );
  Pecos::RealMatrix GridPoints_check(2,29);
  std::cout << "Quadrature points:\n";
  Pecos::write_data(std::cout, variable_sets, false, true, true);
  
  std::cout << "Quadrature weights:\n";
  assert(l_driver.type1_weight_sets()[0] == 4);
  assert(l_driver.type1_weight_sets()[5] == .25);
  Pecos::write_data(std::cout, l_driver.type1_weight_sets());
  */

  return 0;
}
