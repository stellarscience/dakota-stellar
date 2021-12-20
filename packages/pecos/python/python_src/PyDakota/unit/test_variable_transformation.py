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
class TestMonomialApproximation(unittest.TestCase):
    def setUp(self):
        pass

    def test_affine_variable_transformation(self):
        """
        Apply affine transformation of bounded variables.
        E.g. from [a,b] to [-1,1].
        """
        # Define the function variables
        num_vars = 2; num_samples = 10
        variables = BoundedVariables()
        ranges = define_homogeneous_ranges(num_vars, 0., 1.);
        variables.set_ranges(ranges)

        var_transform = AffineVariableTransformation()
        var_transform.set_variables(variables)

        samples = numpy.random.uniform(0.,1.,(num_vars,num_samples))
        transformed_samples = var_transform.map_samples_from_user_space(
            samples)

        # Check all samples are in the canonical domain [-1,1]
        assert numpy.all(transformed_samples.min(axis=1)>=-1.)
        assert numpy.all(transformed_samples.max(axis=1)<=1.)

        # Check that the canonical samples can be transformed
        # back into user domain
        user_space_samples = var_transform.map_samples_to_user_space(
            transformed_samples)
        assert numpy.allclose(user_space_samples,samples)

    def test_askey_probabilistic_transformation(self):
        """Transform probabilistic variables into their canonical form. 
        For example Beta(1,2)[0,1] to Beta(1,2)[-1,1]. 
        """
        assert False, 'Test not yet implemented'

    def test_nataff_probabilistic_transformation(self):
        """Transform correlated probabilistic variables into uncorrelated
        random variables. 
        """
        assert False, 'Test not yet implemented'


    def test_rosenblatt_probabilistic_transformation(self):
        """Transform probabilistic variables into uniform random 
        variables. 
        """
        assert False, 'Test not yet implemented'
    

if __name__ == '__main__':
    unittest.main()
