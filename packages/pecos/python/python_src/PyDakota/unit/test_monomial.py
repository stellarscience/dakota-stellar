#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

import unittest, numpy
from PyDakota.approximation import *
from PyDakota.regression import *
from PyDakota.math_tools import *
from PyDakota.options_list import *
class TestMonomialApproximation(unittest.TestCase):
    def setUp(self):
        pass

    def test_pyfunction(self):
        """
        Make a python function a PyDakota function.
        """
        num_vars = 2; num_samples = 10
        additive_quadratic_function = \
          lambda x: numpy.sum(x**2)*numpy.arange(1,3) + \
          numpy.sum(x)*numpy.arange(2,4) + numpy.arange(1,3)
        function = PyFunction(additive_quadratic_function)
        samples = numpy.random.uniform(-1,1,(num_vars,num_samples))
        values = function.value(samples)

        # Check function.value produces the same result as the function
        # it is wrapping
        true_values = numpy.empty((num_samples,2),float)
        for i in xrange(num_samples):
            true_values[i,:] = additive_quadratic_function(samples[:,i])
        assert numpy.allclose(true_values, values)

        # Check that function.value can be called with a 1D array
        assert numpy.allclose(
            function.value(samples[:,0]), true_values[0,:])

        # check that function.value can be used with a function that
        # returns a scalar
        additive_quadratic_scalar_valued_function = \
          lambda x: numpy.sum(x**2) + numpy.sum(x)*2 +1.
        function = PyFunction(additive_quadratic_scalar_valued_function)
        values = function.value(samples)
        true_values = numpy.empty((num_samples,1),float)
        for i in xrange(num_samples):
            true_values[i,0] = additive_quadratic_scalar_valued_function(
                samples[:,i])
        assert numpy.allclose(true_values, values)

    def test_define_homogeneous_ranges(self):
        """Generate a hypercube with the same bounds for each dimension"""
        num_vars = 3
        ranges = define_homogeneous_ranges(num_vars, 0., 1.);
        true_ranges = numpy.ones((2*num_vars),float)
        true_ranges[::2]=0.
        assert numpy.allclose(ranges,true_ranges)

        ranges = define_homogeneous_ranges(num_vars, -3., 2.);
        true_ranges = numpy.ones((2*num_vars),float)*2.
        true_ranges[::2]=-3.
        assert numpy.allclose(ranges,true_ranges)

    def test_compute_hyperbolic_indices(self):
        """Generate total-degree polynomial indices"""
        num_vars = 2; degree = 3
        basis_indices = compute_hyperbolic_indices(num_vars, degree, 1.)
        num_indices = cardinality_of_total_degree_polynomial(
            num_vars,degree)
        assert basis_indices.shape==(num_vars,num_indices)
        true_indices = numpy.array(
            [[0,0],[1,0],[0,1],[2,0],[1,1],
             [0,2],[3,0],[2,1],[1,2],[3,0]]).T
        for i in xrange(num_indices):
            found = False
            for j in xrange(num_indices):
                if numpy.allclose(basis_indices[:,j],true_indices[:,j]):
                    found = True
                    break
            assert found

    def test_set_basis_indices(self):
        """
        Check exception is thrown in basis_indices is set before
        variable transformation
        """
        # Define the function variables
        num_vars = 2; degree = 3
        variables = BoundedVariables()
        ranges = define_homogeneous_ranges(num_vars, 0., 1.);
        variables.set_ranges(ranges)
        
        # Define the variable transformation used to covert data in
        # the user define space into the space native to the approximation
        var_transform = AffineVariableTransformation()
        var_transform.set_variables(variables)
        
        approx = Monomial()
        basis_indices=compute_hyperbolic_indices(num_vars, degree, 1.)
        self.assertRaises(
            RuntimeError,approx.set_basis_indices,basis_indices)
        approx.set_variable_transformation(var_transform)
        approx.set_basis_indices(basis_indices)
    
    def initialize_homogeneous_polynomial_01(self, num_vars, degree,
                                             use_tensor_product_indices,
                                             poly_type):

        # Define the function variables
        variables = BoundedVariables()
        ranges = define_homogeneous_ranges(num_vars, 0., 1.);
        variables.set_ranges(ranges)
        
        # Define the variable transformation used to covert data in
        # the user define space into the space native to the approximation
        var_transform = AffineVariableTransformation()
        var_transform.set_variables(variables)

        from PyDakota.univariate_polynomials import LEGENDRE_ORTHOG
        basis_types = numpy.asarray([LEGENDRE_ORTHOG]*num_vars,dtype=numpy.int32)
        opts = {'poly_type':poly_type,'basis_types':basis_types}
        approx = polynomial_approximation_factory(var_transform, opts)

        if use_tensor_product_indices:
            degrees = numpy.array([degree]*num_vars,dtype=numpy.int32)
            basis_indices = tensor_product_indices(degrees)
        else:
            basis_indices=compute_hyperbolic_indices(num_vars, degree, 1.)
        approx.set_basis_indices(basis_indices)

        return approx

    def test_generate_monomial_basis_matrix(self):
        """
        Generate basis matrix envoking variable transformation.

        Samples x are generated in the user-domain [0,1]^d. Canonical
        domain of monomial is [-1,1]^d. So in 2D basis_matrix will
        be (2*x[0,:]-1)**j*(2*x[1,:]-1)**k for a polynomial index
        [j,k].
        """
        num_vars = 2; degree = 3; use_tensor_product_indices=True
        approx = self.initialize_homogeneous_polynomial_01(
            num_vars, degree, use_tensor_product_indices, MONOMIAL)
        num_samples = 10

        index_sets_1d = [numpy.arange(degree+1)]*num_vars
        tensor_product_indices = cartesian_product(index_sets_1d, 1)
        assert numpy.allclose(
            tensor_product_indices,approx.get_basis_indices())

        samples = numpy.random.uniform(0,1,(num_vars,num_samples))
        basis_matrix = approx.generate_basis_matrix(samples)

        for i in range(num_samples):
            basis_vals_1d = [(2*samples[0,i]-1.)**numpy.arange(degree+1)]
            for d in range(1,num_vars):
                basis_vals_1d+=[(2*samples[d,i]-1)**numpy.arange(degree+1)]
            assert numpy.allclose(
                basis_matrix[i,:],outer_product(basis_vals_1d,1))
        

    def test_train_polynomial_using_fine_grained_regression_tools(self):
        num_vars = 2; degree = 3; use_tensor_product_indices=False
        approx = self.initialize_homogeneous_polynomial_01(
            num_vars, degree, use_tensor_product_indices, poly_type=PCE)
        num_training_samples = 2*approx.num_terms()

        # Define the function to approximate
        additive_quadratic_function = \
          lambda x: numpy.sum(x**2)*numpy.arange(1,3) + \
          numpy.sum(x)*numpy.arange(2,4) + numpy.arange(1,3)
        function = PyFunction(additive_quadratic_function)

        # Check that if value is called before coefficients is set an
        # exception is thrown
        samples = numpy.random.uniform(-1,1,(num_vars,10))
        self.assertRaises(RuntimeError,approx.value,samples)

        # Construct the approximation
        training_samples = numpy.random.uniform(0,1,(num_vars,num_training_samples))
        training_function_vals = function.value(training_samples)
        basis_matrix = approx.generate_basis_matrix(training_samples)
        opts_dict = {'regression_type':SVD_LEAST_SQ_REGRESSION}
        linsys_opts = OptionsList(opts_dict)
        solver = regression_solver_factory(linsys_opts);
        solver.solve(
            basis_matrix, training_function_vals, linsys_opts);
        coeffs = solver.get_final_solutions();
        approx.set_coefficients(coeffs)

        # Check that approximation is an interpolant
        approx_vals = approx.value(training_samples)
        assert numpy.allclose(approx_vals,training_function_vals)

        # Check that approximation is exact everywhere
        num_samples = 100
        samples = numpy.random.uniform(0,1,(num_vars,num_samples))
        approx_vals = approx.value(samples)
        function_vals = function.value(samples)
        assert numpy.allclose(approx_vals,function_vals)

    def test_train_polynomial_using_regression_builder(self):
        num_vars = 2; degree = 3; use_tensor_product_indices=False
        approx = self.initialize_homogeneous_polynomial_01(
            num_vars, degree, use_tensor_product_indices, poly_type=PCE)
        num_training_samples = 2*approx.num_terms()

        # Define the function to approximate
        additive_quadratic_function = \
          lambda x: numpy.sum(x**2)*numpy.arange(1,3) + \
          numpy.sum(x)*numpy.arange(2,4) + numpy.arange(1,3)
        function = PyFunction(additive_quadratic_function)


        opts_dict = {'num_samples':num_training_samples,
                     'regression_type':SVD_LEAST_SQ_REGRESSION,
                     'sample_type':'probabilistic_MC'}
        regression_opts = OptionsList(opts_dict)
        builder = RegressionBuilder()
        builder.set_target_function(function)
        result = OptionsList()
        builder.build(regression_opts, approx, result)
        training_samples = result['training_samples']
        training_function_vals = result['training_values']

        # Check that approximation is an interpolant
        approx_vals = approx.value(training_samples)
        assert numpy.allclose(approx_vals,training_function_vals)

        # Check that approximation is exact everywhere
        num_samples = 100
        samples = numpy.random.uniform(-1,1,(num_vars,num_samples))
        approx_vals = approx.value(samples)
        function_vals = function.value(samples)
        assert numpy.allclose(approx_vals,function_vals)

if __name__ == '__main__':
    unittest.main()

    suite = unittest.TestSuite()
    suite.addTest( TestMonomialApproximation( "test_generate_monomial_basis_matrix" ) ) 
    runner = unittest.TextTestRunner()
    runner.run( suite )
