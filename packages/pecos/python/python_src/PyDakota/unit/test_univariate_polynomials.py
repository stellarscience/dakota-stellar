#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

import unittest, numpy, os
from PyDakota.univariate_polynomials import *

def test_orthogonal_polynomial(samples, poly, true_values, true_derivatives,
                               true_norms=None):
    """
    Test the following functions of a polynomial basis :
        type1_value
        type1_gradient
        norm_squared (if true_norms is not None and poly is an
                      OrthogonalPolynomial)

    Parameters
    ----------
    true_values : np.ndarray (num_terms x num_samples)
       Exact values of the polynomial for each degree

    true_derivaties : np.ndarray (num_terms x num_samples)
       Exact values of the 1st derivative of the polynomial for each degree

    true_norms : np.ndarray (num_terms-1)
        l2_norm of basis from degree 1 and up
    """
    assert samples.ndim==1
    nsamples = samples.shape[0]
    assert true_values.shape[1]==nsamples
    nterms = true_values.shape[0]
    assert true_derivatives.shape[1]==nsamples

    for i in range(nterms):
      for j in range(nsamples):
        assert numpy.allclose(poly.type1_value(samples[j],i), true_values[i,j])

    degree = nterms-1
    assert numpy.allclose(poly.values(samples, degree), true_values.T)

    gradient_nterms=true_derivatives.shape[0]
    for i in range(gradient_nterms):
      for j in range(nsamples):
        assert numpy.allclose(
          poly.type1_gradient(samples[j],i), true_derivatives[i,j])

    if true_norms is not None:
        for i in xrange( true_norms.shape[0] ):
            assert numpy.allclose( true_norms[i], poly.norm_squared( i+1 ) )




def test_orthogonality(poly, degree, eps=2*numpy.finfo( numpy.float ).eps):
    x = poly.collocation_points(degree+1)
    w = poly.type1_collocation_weights(degree+1)
    V = poly.values(x,degree)
    nterms = V.shape[1]
    for i in xrange(nterms):
        for j in xrange(nterms):
            integral = numpy.sum(V[:,i]*V[:,j]*w)
            if i!=j:
                assert numpy.absolute(integral)<eps,(i,j,numpy.absolute(integral))
            else:
                l2_norm = poly.norm_squared(i)
                msg = (j,numpy.absolute(l2_norm-integral))
                assert numpy.absolute(l2_norm-integral)<eps, msg

class TestUnivariatePolynomials(unittest.TestCase):
    def setUp( self ):
        pass

    def test_legendre_polynomial(self):
        """Test Legendre which are orthogonal to uniform random variables"""
        poly_type=LEGENDRE_ORTHOG
        poly = BasisPolynomial(poly_type)
        x = numpy.linspace(-1., 1., 100)
        true_values = numpy.array( [0.*x+1., x, ( 3. * x**2 - 1. ) / 2.,
                                    ( 5. * x**3 - 3.*x ) / 2.,
                                    ( 35. * x**4 - 30 * x**2 + 3. ) / 8.,
                                    x * ( 63. * x**4 - 70. * x**2 + 15) / 8.])

        true_derivatives = numpy.array( [0.*x, 0.*x+1., 3.*x] )
        true_norms = numpy.array([1./3.,1./5.,1./7.,1./9.])
        test_orthogonal_polynomial(
            x, poly, true_values, true_derivatives, true_norms )

        degree = 100
        test_orthogonality(poly, degree, eps=1e-14)

    def test_jacobi_polynomials( self ):
        """ Jacobi( 0.5, 0.5 ) which are orthogonal to Beta(1.5,1.5)
        random variables"""
        poly_type=JACOBI_ORTHOG
        poly = BasisPolynomial(poly_type)
        poly.alpha_stat(1.5); poly.beta_stat(1.5)
        x = numpy.linspace(-1., 1., 100)
        true_values = numpy.array(
            [0.*x+1.,
            3. * x / 2, 5. * ( 4. * x**2 - 1. ) / 8.,
            35. * x * (2. * x**2 - 1. ) / 16.,
            63. * ( 16. * x**4 - 12 * x**2 + 1. )/128.,
            231. * x * ( 16. * x**4 - 16.*x**2+3)/256.])
        true_derivatives = numpy.array([0.*x,0.*x+3./2.,5.*x])

        # answers found using wolfram alpha
        #N[int[(JacobiP[5,1/2,1/2,2*x-1])^2*(x)^(1/2)*(1-x)^(1/2)/Gamma[1.5]^2*Gamma[3],{x,0,1}],16]
        true_norms = numpy.array([0.5625,0.390625, 0.299072265625,
                                  0.242248535156245, 0.2035560607910156])
        test_orthogonal_polynomial(
            x, poly, true_values, true_derivatives, true_norms )

        degree = 98# probem with degree > 98
        test_orthogonality(poly, degree, eps=1e-14)

if __name__ == '__main__':
    unittest.main()
