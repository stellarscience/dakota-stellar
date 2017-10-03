/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_hierarch_approx_driver.cpp
    \brief A test program for HierarchInterpPolyApproximation class. */

#include "LRDHierarchInterpPolyApproximation.hpp"
#include "LocalRefinableDriver.hpp"

using namespace Pecos;

//Forward declarations of the test functions.
Real test_function(const RealVector& x);
RealVector test_function_grad(const RealVector& x);

Real test_function_2(const RealVector& x);
void test_function_2_grad(const RealVector& x, RealVector& grad);


int main(int argc, char** argv)
{

  /*Declare a new LRDHierarchInterpPolyApproximation object
    of basis_type 0 in one variable.  This approximation
    will consist of piecewise linear functions since
    use_derivs is false.
  */
  LRDHierarchInterpPolyApproximation *a = 
    new LRDHierarchInterpPolyApproximation(0,1,false);
  
  /*Construct an integration driver for the approximation.
    The domain of the approximation is [0,1] and the approximation is 
    intialized to level 1.
  */
  LocalRefinableDriver l_driver;
  l_driver.initialize_grid(RealArray(1,0),RealArray(1,1),1);
  const std::vector<CollocationPoint>& col_pts = 
    l_driver.get_collocation_points();
  a->integration_driver_rep(&l_driver);

  /*Perform a few simple tests to make sure everything is working. First we 
    will push a few dummy values into the approximation, compute the 
    coefficients and test that the correct values and gradients 
    are computed by the approximation.
  */
  RealVector x(1);
  Real fn_val;
  RealVector fn_grad;
  RealSymMatrix fn_hess;
  fn_val = 2.718;
  x[0] = col_pts[0].get_point()[0];

  SurrogateData points;
  points.push_back( SurrogateDataVars(x),
		    SurrogateDataResp(fn_val,fn_grad,fn_hess,1) );

  a->surrogate_data(points);

  a->compute_coefficients();

  std::cout << "The value is... " << a->value(x) << std::endl;
  std::cout << "The gradient is... " << a->gradient_basis_variables(x)[0] << std::endl;

  assert(a->value(x) == 2.718);
  assert(a->gradient_basis_variables(x)[0] == 0);

  l_driver.refine_globally();

  x[0] = col_pts[1].get_point()[0];
  fn_val = 0.0;
  points.push_back(SurrogateDataVars(x),
		   SurrogateDataResp(fn_val,fn_grad,fn_hess,1) );

  x[0] = col_pts[2].get_point()[0];
  fn_val = 5.6;
  points.push_back(SurrogateDataVars(x),
		   SurrogateDataResp(fn_val,fn_grad,fn_hess,1) );

  a->surrogate_data(points);
  a->compute_coefficients();
  
  x[0] = col_pts[1].get_point()[0];

  std::cout << "The value is... " << a->value(x) << std::endl;
  std::cout << "The gradient is... " << a->gradient_basis_variables(x)[0] << std::endl;

  assert(a->value(x) == 0);
  assert(a->gradient_basis_variables(x)[0] == 0);

  x[0] = col_pts[2].get_point()[0];

  std::cout << "The value is... " << a->value(x) << std::endl;
  std::cout << "The gradient is... " << a->gradient_basis_variables(x)[0] << std::endl;

  assert(a->value(x) == 5.6);
  assert(a->gradient_basis_variables(x)[0] == 0);

  /*These next two test check that the linera interpolation is correct. At
    this point the 'unknown' function is known to be equal to 0 at 0, 2.718
    at .5 and 5.6 at 1.  Since the approximation is piecewise linear, at .25
    the interpolant should be equal to (0 + 2.718)/2 = 1.359.  The gradient
    should be (2.718 - 0)/(.5 - 0) = 5.436.  A similar check is done for
    x = .75.
  */
  x[0] = .25;

  std::cout << "The value is... " << a->value(x) << std::endl;
  std::cout << "The gradient is... " << a->gradient_basis_variables(x)[0] << std::endl;

  assert(std::abs(a->value(x) - 1.359) < 1e-10);
  assert(std::abs(a->gradient_basis_variables(x)[0] - 5.436) < 1e-10);

  x[0] = .75;

  std::cout << "The value is... " << a->value(x) << std::endl;
  std::cout << "The gradient is... " << a->gradient_basis_variables(x)[0] << std::endl;

  assert(std::abs(a->value(x) - 4.159) < 1e-10);
  assert(std::abs(a->gradient_basis_variables(x)[0] - 5.764) < 1e-10);

  /*Now we do a local refinement of the grid.  Only the collocation point at
    zero is refined.  The new grid has points {0, .25, .5, 1} a few tests are
    checked at this level.  The interested (or skeptical) reader can check that
    things are working as intended.
  */
  BoolDeque refinementSelector(2,false);
  refinementSelector[0] = true;
  l_driver.refine_locally(refinementSelector);
  x[0] = (col_pts[3].get_point())[0];
  fn_val = -1;
  
  points.push_back( SurrogateDataVars(x),
		    SurrogateDataResp(fn_val,fn_grad,fn_hess,1) );

  a->surrogate_data(points);
  a->increment_coefficients();
  
  x[0] = .25;
  std::cout << "The value is... " << a->value(x) << std::endl;
  std::cout << "The gradient is... " << a->gradient_basis_variables(x)[0] << std::endl;
  
  /*2D example.  The new approximation domain is [-1,1]^2 and the grid is
    initialized to level 8.
  */
  points.clear_data();
  l_driver.initialize_grid(RealArray(2,-1),RealArray(2,1),8);

  RealMatrix col_pts_mat;
  l_driver.compute_grid(col_pts_mat);
  
  /*Loop over the collocation points and evaluate the test function at each of
    the collocation points.  Again we're testing the linear interpolant so the
    gradients don't matter.
  */
  for (unsigned int idx = 0; idx < l_driver.grid_size() ; ++idx) {
    RealVector x(Teuchos::Copy,col_pts_mat[idx],2);
    Real fn_val = test_function(x);
    points.push_back( SurrogateDataVars(x),
		      SurrogateDataResp(fn_val,fn_grad,fn_hess,1) );
  }
  std::cout << "Grid size is: " << l_driver.grid_size() << std::endl;
  a->surrogate_data(points);
  
  //This is a rather expensive operation.
  a->compute_coefficients();

  //Test that the interpolant is exact at the grid points.
  for (unsigned int idx = 0; idx < l_driver.grid_size(); ++idx ) {
    assert( std::abs( a->value( points.continuous_variables(idx) ) -  
		      test_function(points.continuous_variables(idx)) < 1e-9));
  }

  /*Test the interpolant at a dummy point.  This grid is fine enough that the
    error should be small.
  */
  x.size(2);
  x[0] = .69324;
  x[1] = .84529;
  std::cout << "Value =  " << test_function(x) << std::endl;
  std::cout << "Approximate = " << a->value(x) << std::endl;
  std::cout << "Error = " << std::abs(a->value(x) - test_function(x)) 
	    << std::endl;

  //Check the mean computation also.  This should evaluate to near zero.
  std::cout << "Mean = " << a->mean() << std::endl;


  /*2D example with gradients. For this test we're going to use a function
    with gradients.  First construct a new approximation to use gradients and
    initialize a new grid.  The approximation domain is [-1,1]^2 and the grid
    is initialized to level 9.
  */ 
  delete a;
  a = new LRDHierarchInterpPolyApproximation(0,2,true);
  a->integration_driver_rep(&l_driver);
  points.clear_data();
  l_driver.initialize_grid(RealArray(2,-1),RealArray(2,1),9,PIECEWISE_CUBIC_INTERP,true);

  l_driver.compute_grid(col_pts_mat);
  
  //Evaluate the function and the gradient at the collocation points.
  fn_grad.size(2);
  for (unsigned int idx = 0; idx < l_driver.grid_size() ; ++idx) {
    RealVector x(Teuchos::Copy,col_pts_mat[idx],2);
    Real fn_val = test_function_2(x);
    test_function_2_grad(x,fn_grad);
    points.push_back( SurrogateDataVars(x),
		      SurrogateDataResp(fn_val,fn_grad,fn_hess,3) );
  }
  
  a->surrogate_data(points);
  std::cout << "Grid size is: " << l_driver.grid_size() << std::endl;
  a->compute_coefficients();
  
  //Check that the interpolation is exact at the grid points.
  for (unsigned int idx = 0; idx < l_driver.grid_size(); ++idx ) {
    assert(std::abs( a->value( points.continuous_variables(idx) ) -  
		   test_function_2(points.continuous_variables(idx)) < 1e-9) );
  }

  //Check the interpolation at a non-grid point.
  x.size(2);
  x[0] = .69324;
  x[1] = .84529;
  RealVector grad2(2);
  std::cout << "Value =  " << test_function_2(x) << std::endl;
  std::cout << "Approximate = " << a->value(x) << std::endl;
  std::cout << "Error = " 
	    << std::abs(a->value(x) - test_function_2(x))/test_function_2(x) 
	    << std::endl;
  assert( std::abs(test_function_2(x) - a->value(x) ) / test_function_2(x) < 
	  1e-2);

  test_function_2_grad(x,grad2);
  std::cout << "Gradient = [" << grad2[0] 
	    << "," << grad2[1] << "]" <<std::endl;
  std::cout << "Approx_gradient = [" << a->gradient_basis_variables(x)[0] 
	    << "," << a->gradient_basis_variables(x)[1] << "]" <<std::endl;
  
  SizetArray dvv(1);
  dvv[0] = 1;
  double approx_df_dx1 = a->gradient_basis_variables(x)[0];
  double approx_df_dx2 = a->gradient_basis_variables(x)[1];
  std::cout << "All variables gradient dx1 = " << a->gradient_basis_variables(x,dvv)[0] << std::endl;
  assert( approx_df_dx1 = a->gradient_basis_variables(x,dvv)[0] ); 
  dvv[0] = 2;
  std::cout << "All variables gradient dx2 = " << a->gradient_basis_variables(x,dvv)[0] << std::endl;
  assert( approx_df_dx2 = a->gradient_basis_variables(x,dvv)[0] );
  assert( std::abs(a->mean()) < 1e-10 );
  std::cout << "Mean = " << a->mean() << std::endl;
  

  return EXIT_SUCCESS;

}

Real test_function(const RealVector& x) {

  unsigned int dimension = x.length();
  Real return_val = 1;
  for ( unsigned int idx = 0; idx < dimension ; ++idx ) {
    if ( x[idx] == 0 ) return 0;
    else return_val *= std::abs(x[idx])*std::sin(1/x[idx]);
  }
  return return_val;
}

RealVector test_function_grad(const RealVector& x){
  
  unsigned int dimension = x.length();
  RealVector grad(dimension,1);
  for ( unsigned int idx = 0; idx < dimension; ++idx ){
    for ( unsigned int idx2 = 0; idx2 < dimension; ++idx2) {
      if ( idx != idx2 ) {
	if ( x[idx2] != 0 ) {
	  grad[idx] *= std::abs(x[idx2])*std::sin(1/x[idx2]);
	} else  grad[idx] *= 0;  
      }
    }
  }
  for ( unsigned int idx = 0; idx< dimension; ++idx ){
    if ( x[idx] < 0 ) {
      grad[idx] *= std::abs(x[idx])*cos(1/x[idx])*(-x[idx]*x[idx]) - sin(1/x[idx]);
    } else if (x[idx] > 0) {
      grad[idx] *= std::abs(x[idx])*cos(1/x[idx])*(-x[idx]*x[idx]) + sin(1/x[idx]);
    } else {
      grad[idx] *= 0;
    }
  }
  return grad;
}

Real test_function_2(const RealVector& x) {

  unsigned int dimension = x.length();
  Real return_val = 1;
  for ( unsigned int idx = 0; idx < dimension ; ++idx ) {
    if ( x[idx] == 0 ) return 0;
    else return_val *= std::sin(x[idx]);
  }
  return return_val;
}

void test_function_2_grad(const RealVector& x, RealVector& grad){
  
  unsigned int dimension = x.length();
  for (unsigned int idx = 0; idx < dimension; ++idx) grad[idx] = 1;
  for (unsigned int idx = 0; idx < dimension; ++idx) {
    for ( unsigned int idx2 = 0; idx2 < dimension; ++idx2 ) {
      if ( idx != idx2) grad[idx] = grad[idx] * std::sin(x[idx2]);
      else grad[idx] = grad[idx] * ( std::cos(x[idx]) );
    }
  }
}
