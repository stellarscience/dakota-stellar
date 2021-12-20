#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

import unittest, numpy, os
from PyDakota.regression import *
from PyDakota.options_list import OptionsList
class TestLASSO(unittest.TestCase):
    def setUp( self ):
        self.diabetes_matrix, self.diabetes_rhs = self.load_diabetes_data()
        self.pce_matrix, self.pce_rhs, self.pce_exact_coef = self.load_pce_data()

    def load_diabetes_data( self ):
        base_dir = os.path.join(os.path.dirname(__file__), 'data')
        matrix = numpy.loadtxt(os.path.join(base_dir, 'diabetes_data.csv.gz'))
        rhs = numpy.loadtxt(os.path.join(base_dir, 'diabetes_target.csv.gz'))
        return matrix, rhs

    def load_pce_data( self ):
        base_dir = os.path.join(os.path.dirname(__file__), 'data')
        matrix = numpy.loadtxt(os.path.join(base_dir, 'pce_data.csv.gz'))
        rhs = numpy.loadtxt(os.path.join(base_dir, 'pce_target.csv.gz'))
        exact_coef = numpy.loadtxt(os.path.join(base_dir, 'pce_coef.csv.gz'))
        return matrix, rhs, exact_coef
    
    def test_lar_last_step( self ):
        """ 
        Test that the last step of least angle regression returns the
        least squares solution
        """
        solver = LARSolver()
        matrix = self.diabetes_matrix
        rhs = self.diabetes_rhs
        for store_history in [False,True]:
            regression_opts = {'regression_type':LEAST_ANGLE_REGRESSION,
                               'verbosity':0,'store-history':store_history}
            regression_opts = OptionsList(regression_opts)
            solver.solve(matrix, rhs, regression_opts)
            coef_lasso = solver.get_solutions_for_all_regularization_params(0) 
            coef_lstsq = numpy.linalg.lstsq(matrix, rhs)[0]
            assert numpy.allclose(coef_lasso[:,-1].squeeze(), coef_lstsq.squeeze())

    def test_store_history_consistent_with_no_history(self):
        solver = LARSolver()
        matrix = self.diabetes_matrix
        rhs = self.diabetes_rhs

        coef_lasso = []; residuals=[]
        for store_history in [False,True]:
            regression_opts = {'regression_type':LEAST_ANGLE_REGRESSION,
                               'verbosity':0,'store-history':store_history}
            solver.solve(matrix, rhs, regression_opts)
            coef_lasso.append(solver.get_solutions_for_all_regularization_params(0))
            residuals.append(solver.get_residuals_for_all_regularization_params(0))
        assert numpy.allclose(coef_lasso[0][:,0],coef_lasso[1][:,-1],atol=1e-15)
        assert numpy.allclose(residuals[0][0],residuals[1][-1],atol=1e-12)
        

    def test_lar_factory( self ):
        """ 
        Test that the regression factory returns a lar solver and that the 
        solver works correctly. That is the last step of least angle 
        regression returns the least squares solution  
        """
        factory_opts = OptionsList({'regression_type':LEAST_ANGLE_REGRESSION})
        print type(factory_opts)
        solver = regression_solver_factory(factory_opts)
        matrix = self.diabetes_matrix
        rhs = self.diabetes_rhs
        for store_history in [False,True]:
            regression_opts = {'regression_type':LEAST_ANGLE_REGRESSION,
                               'verbosity':0,'store-history':store_history}
            regression_opts = OptionsList(regression_opts)
            solver.solve(matrix, rhs, regression_opts)
            coef_lasso = solver.get_solutions_for_all_regularization_params(0) 
            coef_lstsq = numpy.linalg.lstsq(matrix, rhs)[0]
            assert numpy.allclose(coef_lasso[:,-1].squeeze(), coef_lstsq.squeeze())

    def test_lar_memory_management(self):
        """ 
        LAR internally allocates memory for solutions and QR factorization in 
        blocks so that if residual tolerance is met before maximum number of
        solutions is reached min (num_rows, num_cols) then memory has not been 
        wasted. Memory is allocated in chunks of size 100. Test algorithm works
        when more than one chunk is needed.
        """
        num_rows = 200
        num_cols = 200
        sparsity = 100
        memory_chunk_size = 20
        solver = LARSolver()
        matrix = numpy.random.normal(0.,1.,(num_rows,num_cols))
        x = numpy.zeros((num_cols),float)
        I = numpy.random.permutation(num_cols)[:sparsity]
        x[I] = numpy.random.normal(0.,1.,(sparsity))
        rhs = numpy.dot(matrix, x)
        regression_opts = {'regression_type':LEAST_ANGLE_REGRESSION,
                           'verbosity':0,'memory-chunk-size':memory_chunk_size}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_lasso = solver.get_final_solutions()
        assert numpy.allclose(coef_lasso[:,0].squeeze(), x.squeeze())
        

    def test_lasso_last_step( self ):
        """
        Test that the last step of the lasso variant of 
        least angle regression returns the least squares solution
        """
        solver = LARSolver()
        matrix = self.diabetes_matrix # use unnormalized data
        rhs = self.diabetes_rhs
        regression_opts = {'regression_type':LASSO_REGRESSION,'verbosity':0,
                           'normalize-inputs':False}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_lasso = solver.get_solutions_for_all_regularization_params(0) 
        coef_lstsq = numpy.linalg.lstsq(matrix, rhs)[0]
        #print coef_lasso[:,-1].squeeze(), coef_lstsq.squeeze()
        assert numpy.allclose(coef_lasso[:,-1].squeeze(), coef_lstsq.squeeze())

    def test_lar_path( self ):
        """ 
        Test max covariance is tied and decreasing
        """
        solver = LARSolver()
        matrix = self.diabetes_matrix.copy()
        rhs = self.diabetes_rhs
        regression_opts = {'regression_type':LASSO_REGRESSION}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_lasso = solver.get_solutions_for_all_regularization_params(0) 
        max_covariance_prev = numpy.finfo(float).max
        for i in xrange( coef_lasso.shape[1] ):
            coef = coef_lasso[:,i]
            residual = rhs - numpy.dot(self.diabetes_matrix, coef)
            covariance = numpy.dot(self.diabetes_matrix.T, residual)
            
            max_covariance = numpy.max(numpy.absolute(covariance))
            assert max_covariance < max_covariance_prev
            max_covariance_prev = max_covariance
            eps = 1e-3
            num_non_zeros = len(covariance[max_covariance-eps<
                                           numpy.absolute(covariance)])
            if i < self.diabetes_matrix.shape[1]-1:
                assert num_non_zeros == i+2
            else:
                # no more than max_pred variables can go into the active set
                assert num_non_zeros == self.diabetes_matrix.shape[1]

    def test_lasso_path( self ):
        """ 
        Test max covariance is tied and decreasing
        """

        solver = LARSolver()
        matrix = self.diabetes_matrix
        rhs = self.diabetes_rhs
        regression_opts = {'regression_type':LASSO_REGRESSION}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_lasso = solver.get_solutions_for_all_regularization_params(0) 
        max_covariance_prev = numpy.finfo(float).max
        for i in xrange( coef_lasso.shape[1] ):
            coef = coef_lasso[:,i]
            residual = rhs - numpy.dot(self.diabetes_matrix, coef)
            covariance = numpy.dot(self.diabetes_matrix.T, residual)
            
            max_covariance = numpy.max(numpy.absolute(covariance))
            assert max_covariance < max_covariance_prev
            max_covariance_prev = max_covariance
            eps = 1e-3
            num_non_zeros = len(covariance[max_covariance-eps<
                                           numpy.absolute(covariance)])
            if i < self.diabetes_matrix.shape[1]-1:
                assert num_non_zeros == i+2
            else:
                # no more than max_pred variables can go into the active set
                assert num_non_zeros == self.diabetes_matrix.shape[1]

    def test_non_negative_lasso_path( self ):
        """ 
        Test when enforcing non-negative coefficients the 
        max covariance is tied and decreasing
        """
        solver = LARSolver()
        assert False, 'test of non-negative lasso not implemented'
            
    def test_lasso_early_exit_tol( self ):
        """
        Test that the algorithm terminates correctly when tolerance is set
        """
        tol = 3.40e3
        solver = LARSolver()
        matrix = self.diabetes_matrix
        rhs = self.diabetes_rhs
        regression_opts = {'regression_type':LASSO_REGRESSION,
                           'residual-tolerance':tol}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_lasso = solver.get_solutions_for_all_regularization_params(0) 
        coef = coef_lasso[:,-1]
        residual = rhs - numpy.dot(self.diabetes_matrix, coef)
        assert numpy.linalg.norm( residual ) < tol

    def test_lasso_early_exit_num_non_zeros( self ):
        """
        Test that the algorithm terminates correctly when the number of nonzeros
        is set
        """
        solver = LARSolver()
        matrix = self.diabetes_matrix
        rhs = self.diabetes_rhs
        max_nnz = 9
        regression_opts = {'regression_type':LASSO_REGRESSION,
                           'verbosity':0,'max-num-non-zeros':max_nnz}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_lasso = solver.get_solutions_for_all_regularization_params(0) 
        assert numpy.count_nonzero( coef_lasso[:,-1] )==max_nnz

    def test_lasso_pce_exact_recovery( self ):
        base_dir = os.path.join(os.path.dirname(__file__), 'data')
        matrix = numpy.loadtxt(os.path.join(base_dir, 'pce_data.csv.gz'))
        rhs = numpy.loadtxt(os.path.join(base_dir, 'pce_target.csv.gz'))
        exact_coef = numpy.loadtxt(os.path.join(base_dir, 'pce_coef.csv.gz'))
        solver = LARSolver()
        regression_opts = {'regression_type':LASSO_REGRESSION,
                           'verbosity':0}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_lasso = solver.get_solutions_for_all_regularization_params(0) 
        assert numpy.allclose( coef_lasso[:,-1], exact_coef )

        
class TestOMP(unittest.TestCase):
    def setUp( self ):
        self.diabetes_matrix, self.diabetes_rhs = self.load_diabetes_data()
        self.pce_matrix, self.pce_rhs, self.pce_exact_coef = self.load_pce_data()

    def load_diabetes_data( self ):
        base_dir = os.path.join(os.path.dirname(__file__), 'data')
        matrix = numpy.loadtxt(os.path.join(base_dir, 'diabetes_data.csv.gz'))
        rhs = numpy.loadtxt(os.path.join(base_dir, 'diabetes_target.csv.gz'))
        return matrix, rhs

    def load_pce_data( self ):
        base_dir = os.path.join(os.path.dirname(__file__), 'data')
        matrix = numpy.loadtxt(os.path.join(base_dir, 'pce_data.csv.gz'))
        rhs = numpy.loadtxt(os.path.join(base_dir, 'pce_target.csv.gz'))
        exact_coef = numpy.loadtxt(os.path.join(base_dir, 'pce_coef.csv.gz'))
        return matrix, rhs, exact_coef
    
    def test_omp_last_step( self ):
        """ 
        Test that the last step of orthogonal matching pursuit returns the
        least squares solution
        """
        solver = OMPSolver()
        #solver.set_verbosity( 3 )
        matrix = self.diabetes_matrix
        rhs = self.diabetes_rhs
        regression_opts = {'regression_type':ORTHOG_MATCH_PURSUIT,'verbosity':0}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_omp = solver.get_solutions_for_all_regularization_params(0) 
        coef_lstsq = numpy.linalg.lstsq(matrix, rhs)[0]
        #print coef_omp[:,-1].squeeze(), coef_lstsq.squeeze()
        assert numpy.allclose(coef_omp[:,-1].squeeze(), coef_lstsq.squeeze())

    def test_omp_memory_management( self ):
        """ 
        OMP internally allocates memory for solutions and QR factorization in 
        blocks so that if residual tolerance is met before maximum number of
        solutions is reached min (num_rows, num_cols) then memory has not been 
        wasted. Memory is allocated in chunks of size 100. Test algorithm works
        when more than one chunk is needed.
        """
        num_rows = 200
        num_cols = 200
        sparsity = 100
        memory_chunk_size = 100
        solver = OMPSolver()
        matrix = numpy.random.normal(0.,1.,(num_rows,num_cols))
        x = numpy.zeros((num_cols),float)
        I = numpy.random.permutation(num_cols)[:sparsity]
        x[I] = numpy.random.normal(0.,1.,(sparsity))
        rhs = numpy.dot(matrix, x)
        regression_opts = {'regression_type':ORTHOG_MATCH_PURSUIT,'verbosity':0,
                           'memory-chunk-size':memory_chunk_size}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_omp = solver.get_final_solutions()
        assert numpy.allclose(coef_omp[:,0].squeeze(), x.squeeze())
        
    def test_omp_path( self ):
        """ 
        Test residual is always decreasing
        """
        solver = OMPSolver()
        #solver.set_verbosity( 3 )
        matrix = self.diabetes_matrix.copy()
        rhs = self.diabetes_rhs
        regression_opts = {'regression_type':ORTHOG_MATCH_PURSUIT,'verbosity':0}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_omp = solver.get_solutions_for_all_regularization_params(0) 
        coef_lstsq = numpy.linalg.lstsq(matrix, rhs)[0]
        residual_norm_prev = numpy.finfo(float).max
        for i in xrange( coef_omp.shape[1] ):
            coef = coef_omp[:,i]
            residual = rhs - numpy.dot(self.diabetes_matrix, coef)
            residual_norm = numpy.linalg.norm( residual )
            assert residual_norm < residual_norm_prev

    def test_omp_early_exit_tol( self ):
        """
        Test that the algorithm terminates correctly when tolerance is set
        """
        tol = 3.40e3
        solver = OMPSolver()
        matrix = self.diabetes_matrix
        rhs = self.diabetes_rhs
        regression_opts = {'regression_type':ORTHOG_MATCH_PURSUIT,'verbosity':0,
                           'residual-tolerance':tol}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_omp = solver.get_solutions_for_all_regularization_params(0) 
        coef = coef_omp[:,-1]
        residual = rhs - numpy.dot(self.diabetes_matrix, coef)
        assert numpy.linalg.norm( residual ) < tol

    def test_omp_early_exit_num_non_zeros( self ):
        """
        Test that the algorithm terminates correctly when the number of nonzeros
        is set
        """
        solver = OMPSolver()
        matrix = self.diabetes_matrix
        rhs = self.diabetes_rhs
        max_nnz = 9
        # setting max_iters will set max_nnz as no columns of A can be removed
        # once added
        regression_opts = {'regression_type':ORTHOG_MATCH_PURSUIT,'verbosity':0,
                           'max-iters':max_nnz}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_omp = solver.get_solutions_for_all_regularization_params(0) 
        assert numpy.count_nonzero( coef_omp[:,-1] )==max_nnz

    def test_omp_pce_exact_recovery( self ):
        base_dir = os.path.join(os.path.dirname(__file__), 'data')
        matrix = numpy.loadtxt(os.path.join(base_dir, 'pce_data.csv.gz'))
        rhs = numpy.loadtxt(os.path.join(base_dir, 'pce_target.csv.gz'))
        exact_coef = numpy.loadtxt(os.path.join(base_dir, 'pce_coef.csv.gz'))
        solver = OMPSolver()
        regression_opts = {'regression_type':ORTHOG_MATCH_PURSUIT,'verbosity':0}
        regression_opts = OptionsList(regression_opts)
        solver.solve(matrix, rhs, regression_opts)
        coef_omp = solver.get_solutions_for_all_regularization_params(0) 
        assert numpy.allclose( coef_omp[:,-1], exact_coef )

def single_test_suite():
    suite = unittest.TestSuite()
    suite.addTest( TestLASSO( "test_lar_last_step" ) ) 
   
    return suite      

if __name__ == '__main__':
    unittest.main()
    #runner = unittest.TextTestRunner()
    #single = single_test_suite()
    #runner.run( single )
