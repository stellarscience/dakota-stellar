#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

import unittest, numpy
from PyDakota.math_tools import *
class TestMathTools(unittest.TestCase):
    def setUp(self):
        pass

    def test_range( self ):
        assert numpy.allclose(numpy.arange(0, 3),range_double(0,3,1))
        assert numpy.allclose(numpy.arange(1, 10, 3), range_double(1, 10, 3))
        assert numpy.allclose(numpy.arange(1, 9, 3), range_double(1, 9, 3))

    def test_cartesian_product(self):
        """
        """
        # test when num elems = 1
        s1 = numpy.arange( 0, 3 )
        s2 = numpy.arange( 3, 5 )

        sets = numpy.array( [[0,3], [1,3], [2,3], [0,4],
                                [1,4], [2,4]], numpy.int )
        output_sets = cartesian_product( [s1,s2], 1 )
        assert numpy.array_equal( output_sets.T, sets )

        # test when num elems > 1
        s1 = numpy.arange( 0, 6 )
        s2 = numpy.arange( 6, 10 )

        sets = numpy.array( [[ 0, 1, 6, 7], [ 2, 3, 6, 7],
                                [ 4, 5, 6, 7], [ 0, 1, 8, 9],
                                [ 2, 3, 8, 9], [ 4, 5, 8, 9]], numpy.int )
        output_sets = cartesian_product( [s1,s2], 2 )
        assert numpy.array_equal( output_sets.T, sets )

    def test_outer_product( self ):
        s1 = numpy.arange( 0, 3 )
        s2 = numpy.arange( 3, 5 )

        test_vals = numpy.array( [0.,3.,6.,0.,4.,8.])
        output = outer_product( [s1,s2] )
        assert numpy.allclose( test_vals, output )

        output = outer_product( [s1] )
        assert numpy.allclose( output, s1 )

    def test_binary_search( self ):
        a = numpy.arange( 5., dtype = numpy.double )
        assert binary_search_double( 2.2, a ) == 2
        assert binary_search_double( 1.9, a ) == 1
        assert binary_search_double( 2., a ) == 2
        assert binary_search_double( 0., a ) == 0
        assert binary_search_double( 4., a ) == 4
        assert binary_search_double( 5., a ) == 4
        assert binary_search_double( -1., a ) == 0

        a = numpy.arange( 6., dtype = numpy.double )
        assert binary_search_double( 2.2, a ) == 2
        assert binary_search_double( 1.9, a ) == 1
        assert  binary_search_double( 3.1, a ) == 3
        assert binary_search_double( 0., a ) == 0
        assert binary_search_double( 4., a ) == 4
        assert binary_search_double( 5., a ) == 5
        assert binary_search_double( 5.5, a ) == 5
        assert binary_search_double( -1., a ) == 0

    def test_polynomial_space_indices( self ):

        num_dims = 2
        max_order = 2
        output_indices = \
          get_multi_dimensional_polynomial_indices( num_dims, max_order )
        indices = numpy.array( [[0,2], [1,1], [2,0]],
                                   numpy.int )
        assert numpy.array_equal( output_indices.T, indices )

        num_dims = 3
        max_order = 2
        output_indices = \
          get_multi_dimensional_polynomial_indices( num_dims, max_order )
        indices = numpy.array(
            [[0,0,2], [0,1,1], [0,2,0], [1,0,1], [1,1,0], [2,0,0]],
            numpy.int )
        assert numpy.array_equal( output_indices.T, indices )

    def test_multi_dimensional_matrix_indexing( self ):

        sizes = numpy.array( [2,3] )
        num_indices = numpy.prod( sizes )
        for i in range( sizes[0] ):
            for j in range( sizes[1] ):
                assert sub2ind( sizes, [i,j] ) == sizes[0]*j + i
                assert numpy.array_equal(
                    ind2sub( sizes,  sizes[0]*j + i, num_indices), [i,j] )

        sizes = numpy.array( [2,3,2] )
        num_indices = numpy.prod( sizes )
        for i in range( sizes[0] ):
            for j in range( sizes[1] ):
                for k in range( sizes[2] ):
                    assert sub2ind( sizes, [i,j,k] ) == \
                      sizes[0]*sizes[1]*k + sizes[0]*j + i
                    assert numpy.array_equal(
                        ind2sub( sizes, sizes[0]*sizes[1]*k + sizes[0]*j+i,
                                num_indices), [i,j,k] )

from scipy.linalg import lu as scipy_lu
class TestLinearAlgebra(unittest.TestCase):
    def test_svd_solve( self ):
        M = 3; N = 20
        A = numpy.random.normal( 0., 1., ( M, N ) )
        b = numpy.random.normal( 0., 1., ( M ) )

        out = numpy.linalg.lstsq( A , b )
        x_true = out[0]
        x = svd_solve( A, b.reshape(b.shape[0],1))[0].squeeze()
        assert numpy.allclose( x, x_true )

    def test_cholesky_solve( self ):
        M = 10; N = 5
        A = numpy.random.normal( 0., 1., ( M, N ) )
        b = numpy.random.normal( 0., 1., ( M ) )

        AtA = numpy.dot( A.T, A )
        assert numpy.allclose( numpy.linalg.cholesky( AtA ), cholesky( AtA ) )

        out = numpy.linalg.lstsq( A , b )
        x_true = out[0]
        x = cholesky_solve( AtA, numpy.dot( A.T, b ) ).squeeze()
        assert numpy.allclose( x, x_true )

        x = solve_using_cholesky_factor(
            numpy.linalg.cholesky( AtA ),
            numpy.dot( A.T, b ).reshape(A.shape[1],1), 1 ).squeeze()
        assert numpy.allclose( x, x_true )

    def test_qr_solve(self):
        M = 10; N = 5
        A = numpy.random.normal( 0., 1., ( M, N ) )
        b = numpy.random.normal( 0., 1., ( M ) )

        out = numpy.linalg.lstsq( A , b )
        x_true = out[0]

        x = qr_solve(
            numpy.dot( A.T, A ),
            numpy.dot( A.T, b ).reshape(A.shape[1],1), 1 ).squeeze()
        assert numpy.allclose( x, x_true )

        A =  numpy.array([[1,2,3],[4,5,6]]).T
        x_true = numpy.arange( 2 ).reshape( 2, 1 )
        b = numpy.dot( A, x_true )
        x = qr_solve(A, b.reshape(b.shape[0],1),0 ).squeeze()
        assert numpy.allclose( x, x_true.squeeze() )

    def test_pivoted_qr_factorization(self):
        rng = numpy.random.RandomState( 0 )
        A = rng.randn( 9, 6 )
        q, r, p = pivoted_qr_factorization( A )
        P = numpy.zeros( ( p.shape[0], p.shape[0] ) )
        P[numpy.arange(p.shape[0],dtype=int),p] = 1.
        numpy.allclose( A, numpy.dot( q, r ) )
        d = abs(numpy.diag(r))
        assert all( d[1:] <= d[:-1] )
        assert numpy.allclose( A[:, p], numpy.dot( q, r ) )

        A = rng.randn( 6, 9 )
        q, r, p = pivoted_qr_factorization( A )
        numpy.allclose( A, numpy.dot( q, r ) )
        d = abs(numpy.diag(r))
        assert all( d[1:] <= d[:-1] )
        assert numpy.allclose( A[:, p], numpy.dot( q, r ) )

        A = rng.randn( 6, 6 )
        q, r, p = pivoted_qr_factorization( A )
        numpy.allclose( A, numpy.dot( q, r ) )
        d = abs(numpy.diag(r))
        assert all( d[1:] <= d[:-1] )
        assert numpy.allclose( A[:, p], numpy.dot( q, r ) )

        A = rng.randn( 6, 6 )
        P, L, U = scipy_lu( A )
        p = numpy.asarray( numpy.nonzero( P )[1],  numpy.int32 )
        LU_inv = lu_inverse( L, U, p )
        assert numpy.allclose( LU_inv, numpy.linalg.inv( A ) )

        A = rng.randn( 3, 3 )
        A = numpy.dot( A.T, A )
        P, L, U = scipy_lu( A )
        if numpy.allclose( P,  numpy.eye( 3 ) ):
            p = numpy.empty( (0), numpy.int32 )
        else:
            p = numpy.asarray( numpy.nonzero( P )[1],  numpy.int32 )
        LU_inv = lu_inverse( L, U, p )
        assert numpy.allclose( LU_inv, numpy.linalg.inv( A ) )

    def test_complete_pivoted_lu_factorization(self):
        A = numpy.arange( 1., 13. ).reshape( 3, 4 )
        L, U, rp, cp = complete_pivoted_lu_factorization( A )
        permuted_A = A[:,cp]
        permuted_A = permuted_A[rp,:]
        assert numpy.allclose( numpy.dot( L, U ), permuted_A )

        A = numpy.arange( 1., 13. ).reshape( 4, 3 )
        L, U, rp, cp = complete_pivoted_lu_factorization( A )
        permuted_A = A[:,cp]
        permuted_A = permuted_A[rp,:]
        assert numpy.allclose( numpy.dot( L, U ), permuted_A )

        A = numpy.arange( 1., 13. ).reshape( 3, 4 ).T
        L, U, rp, cp = complete_pivoted_lu_factorization( A )
        permuted_A = A[:,cp]
        permuted_A = permuted_A[rp,:]
        assert numpy.allclose( numpy.dot( L, U ), permuted_A )

        A = numpy.random.RandomState(1).normal( 0., 1., ( 100, 20 ) )
        L, U, rp, cp = complete_pivoted_lu_factorization( A )
        permuted_A = A[:,cp]
        permuted_A = permuted_A[rp,:]
        assert numpy.allclose( numpy.dot( L, U ), permuted_A )

        A = numpy.random.RandomState(1).normal( 0., 1., ( 20, 100 ) )
        L, U, rp, cp = complete_pivoted_lu_factorization( A )
        permuted_A = A[:,cp]
        permuted_A = permuted_A[rp,:]
        assert numpy.allclose( numpy.dot( L, U ), permuted_A )

    def test_symmetric_eigenvalue_decomposition(self):
        A = numpy.arange( 1., 13. ).reshape( 3, 4 )
        A = numpy.dot( A.T, A )
        eigvals, eigvecs = symmetric_eigenvalue_decomposition( A )
        assert numpy.allclose(
            A, numpy.dot(numpy.dot(eigvecs,numpy.diag(eigvals)),eigvecs.T) )

        A = numpy.random.RandomState(1).normal( 0., 1., ( 20, 100 ) )
        A = numpy.dot( A.T, A )
        eigvals, eigvecs = symmetric_eigenvalue_decomposition( A )
        assert numpy.allclose( A, numpy.dot(numpy.dot(eigvecs,numpy.diag(eigvals)),eigvecs.T) )

    def test_pivoted_lu_factorization(self):
        # test pivoted lu factorization
        A = numpy.array([1.8,2.88,2.05,-0.89,5.25,-2.95,-0.95,-3.8,1.58,
                             -2.69,-2.90,-1.04,-1.11,-0.66,-0.59,0.8]).reshape(4,4)
        P, L, U = scipy_lu(A)
        l,u,p = pivoted_lu_factorization(A)
        assert numpy.allclose( L, l )
        assert numpy.allclose( U, u )
        PIV = numpy.where(P==1)[1]
        assert numpy.allclose( PIV, p )

    def test_truncated_pivoted_lu_factorization(self):

        # test truncated_pivoted lu factorization
        A = numpy.random.RandomState(1).normal( 0, 1, (4,3) )
        P, L, U = scipy_lu(A)
        num_pivots = 3
        l,u,p = truncated_pivoted_lu_factorization( A, num_pivots )
        assert numpy.allclose( A[p,:], numpy.dot( l, u ) )
        assert numpy.allclose( L[:num_pivots,:], l )
        assert numpy.allclose( U[:num_pivots,:], u[:num_pivots,:] )
        assert numpy.allclose( p, [2,3,1] )

        # test truncated_pivoted lu factorization which enforces first
        # n rows to be chosen
        tmp = A[p[0],:].copy()
        A[p[0],:] = A[0,:].copy()
        A[0,:] = tmp
        num_pivots = 3
        num_initial_rows = 1
        l,u,p = truncated_pivoted_lu_factorization( A, num_pivots,
                                                        num_initial_rows )
        assert numpy.allclose( A[p,:], numpy.dot( l, u ) )
        assert numpy.allclose( p, [0,3,1] )

        # Modify the above test to first factorize 4,3 A then factorize
        # B = [A; C] where C is 2*3 and if B was factorized without enforcing
        # A then the factors would be different. Then check that first 
        # 4 rows of LU factors of B are the same as when A was factored.

    def test_pivot_matrix_rows( self ):
        p = numpy.array( [1,1,2,3],dtype=numpy.int32 )
        A = numpy.arange(p.shape[0])
        A_p = A.copy()
        for i in xrange( p.shape[0] ):
            tmp = A_p[i]
            A_p[i] = A_p[p[i]]
            A_p[p[i]] = tmp

        # todo template pivot_matrix rows so I do not have to specify dtype
        A = numpy.arange(p.shape[0],dtype=float)
        A = A.reshape(A.shape[0],1)
        assert numpy.allclose( A_p, pivot_matrix_rows(A,p,1,False).squeeze() )

if __name__ == '__main__':
    unittest.main()
