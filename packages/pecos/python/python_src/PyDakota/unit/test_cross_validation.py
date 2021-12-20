#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

import numpy, unittest
from scipy import interpolate
from PyDakota.regression import *
from PyDakota.options_list import OptionsList
def get_linear_system(num_pts, num_eq_per_pt, num_cols, num_rhs,
                      num_nonzeros, noise_std):
    """
    Build a random linear system used to test solvers of AX=B.
    Returns
    -------
    A : matrix (num_pts*num_eq_per_pt x num_cols)
        The matrix A consiting of entries drawn randomly from N(0,1)
    B : matrix (num_pts*num_eq_per_pt x num_rhs) 
        The data generated from A*X + eps, eps~N(0,noise_std**2)
    X : matrix (num_cols x num_rhs)
        The solutions used to generate the data B. X will be sparse with
        num_nonzeros number of nonzero entries
        
    """
    A = numpy.random.normal(0.,1.,(num_pts*num_eq_per_pt,num_cols))
    X = numpy.zeros((num_cols, num_rhs),float)
    for i in xrange(num_rhs):
        non_zero_indices = \
          numpy.random.permutation(num_cols)[:num_nonzeros]
        X[non_zero_indices,i]=1.
    B = numpy.dot(A,X)
    if noise_std>0:
        B += numpy.random.normal(0.,noise_std,B.shape)
    return A, B, X

def split_sequence(seq, num_chunks):
    """
    Split a sequence into chunks or roughly equal size. Ensures each chunk 
    is of a minimum size and adds remainding indices to first chunks.

    Parameters
    ----------
    seq : list of ndarray
        The sequence to be split

    num_chunks : integer
        The number of chunks desired 

    Returns
    -------
    newseq : list of each chunks
    """
    assert len(seq)>=num_chunks
    newseq = []
    n = len(seq) / num_chunks    # min items per subsequence
    r = len(seq) % num_chunks    # remaindered items
    b,e = 0, n + min(1, r)  # first split
    for i in range(num_chunks):
        newseq.append(seq[b:e])
        r = max(0, r-1)  # use up remainders
        b,e = e, e + n + min(1, r)  # min(1,r) is always 0 or 1

    return newseq

def create_random_fault_data(num_pts, num_failures, num_secondary_failures=0):
    """
    Create random faults data.

    Parameters
    ----------
    num_pts : integer
        The total number of points with data
        
    num_failures : integer
        The number of primary equation failures. There is only one primary 
        equation per point

    num_secondary_failures : integer
        The number of secondary equation failures. There may be more than
        one secondary equation per point, but if one fails we assume all fail

    Returns
    -------
    fault_data : vector (num_pts x 1)
        Vector indicating no failures (0), primary equation failure (1)
        and secondary equation failure (2). On exit total number of non-zeros
        will be <= num_failures + num_secondary_failures. Some secondary 
        failures may occur with a primary failure. In our code we assume 
        primary failures imply secondary failures so when both occur we 
        return the value (1)
    """
    assert num_failures<=num_pts
    assert num_secondary_failures<=num_pts
    
    fault_data = numpy.zeros((num_pts),dtype=numpy.int32)
    # secondary failures
    J = numpy.random.permutation(num_pts)[:num_secondary_failures]
    fault_data[J]=2
    # primary failures. This may overwrite secondary failures as we assume
    # primary failures imply secondary failures
    I = numpy.random.permutation(num_pts)[:num_failures]
    fault_data[I]=1
    
    assert numpy.count_nonzero(fault_data)<=num_failures+num_secondary_failures
    return fault_data

def get_active_row_indices(indices, fault_data, num_eq_per_pt, num_pts):
    """
    For a matrix A with N primary equations and N*(num_eq_per_pt-1)
    secondary equations. Return the sub indices of the rows of A 
    indexed by indices that correspond to non-faulty equations. 
    We assume primary equation failure forces seconday equation failure.

    Parameters
    ----------
    indices : vector (num_indices)
        Row indices of the rows of A to consider. num_indices <=N

    fault_data : vector (A.shape[0])
        Vector indicating no failures (0), primary equation failure (1)
        and secondary equation failure (2).

    num_eq_per_pt : integer
        The number of equations per point. num_eq_per_pt*num_pts = A.shape[0]

    num_pts : integer
        The number of points used to construct A. 
        num_eq_per_pt*num_pts = A.shape[0]

    Returns
    -------
    I : vector (num_faultless_eq)
       The indices of A assocaited with the subset of the rows of A 
       that correspond to non-faulty equations.    
    """
    num_indices = indices.shape[0]
    if fault_data is not None:
        # recall we assume primary eq failure (1) implies secondary eq
        # failure (2)
        num_indices_with_secondary_eq =\
          numpy.where(fault_data[indices]==0)[0].shape[0]
    else:
        num_indices_with_secondary_eq = num_indices
    num_active_matrix_rows=\
      num_indices+(num_eq_per_pt-1)*num_indices_with_secondary_eq
    I = numpy.zeros( (num_active_matrix_rows), numpy.int )
    I[:num_indices] = indices
    K = numpy.arange(num_eq_per_pt-1)
    cntr = 0
    for i in xrange(indices.shape[0]):
        if ( (fault_data is not None) and not (fault_data[indices[i]]) ):
            start_secondary_eq = num_indices+cntr*(num_eq_per_pt-1)
            end_secondary_eq = num_indices+(cntr+1)*(num_eq_per_pt-1)
            I[start_secondary_eq:end_secondary_eq]=\
              num_pts+indices[i]*(num_eq_per_pt-1)+K
            cntr += 1
    return I

def partition_data(A, rhs, cv_iterator, i, num_pts, num_eq_per_pt, fault_data):
    """
    For a matrix A with N primary equations and N*(num_eq_per_pt-1)
    secondary equations. Return the subset of rows of A and b
    indexed by indices that correspond to non-faulty equations. 
    We assume primary equation failure forces seconday equation failure.

    Parameters
    ----------
    A : matrix (num_eq_per_pt*num_pts x num_cols)
        The A matrix of the linear system Ax=rhs

    rhs : vector (num_eq_per_pt*num_pts)
        The rhs of the linear system Ax=rhs

    cv_iterator : CrossValidationIterator 
        A cross validation iterator used to partition the data

    i : integer
        The id of the fold left out for validation

    fault_data : vector (A.shape[0])
        Vector indicating no failures (0), primary equation failure (1)
        and secondary equation failure (2).

    num_eq_per_pt : integer
        The number of equations per point. num_eq_per_pt*num_pts = A.shape[0]

    num_pts : integer
        The number of points used to construct A. 
        num_eq_per_pt*num_pts = A.shape[0]

    Returns
    -------
    A_sub : matrix (num_faultless_eq x num_cols)
       The subset of the rows of A that correspond to non-faulty equations.
    
    b_sub : vector (num_faultless_eq)
       The subset of the rows of b that correspond to non-faulty equations.

    validation_indices : vector (num_faultless_primary_eq)
       The indices of the subset of the rows of A that correspond to 
       non-faulty primary equations.
    """
    training_indices, validation_indices = cv_iterator.get_fold_indices( i )
    I = get_active_row_indices(
        training_indices, fault_data, num_eq_per_pt, num_pts)
    A = A[I,:]
    b = rhs[I]
    b = b.reshape( b.shape[0] )
    return A, b, validation_indices


class TestLinearSystemCrossValidation(unittest.TestCase):
    """
    Testing of cross validation for solving least squares problem Ax=b. 
    """
    def fold_partitioning_test(self, num_pts, num_folds, fault_data=None):
        """
        Helper function to test partitioning of a vector of indices into
        training and validation sets without the presence of faults

        Tests:
          * Intersection of training and validation indices is empty

          * Union training and validation indices contains all the original
            faultless indices
        """
        cv_iterator = CrossValidationIterator()
        cv_iterator.set_seed( -1 ) # always use deterministic partioning
        cv_iterator.set_num_points( num_pts )
        cv_iterator.set_num_folds( num_folds )
        if fault_data is not None:
            cv_iterator.set_fault_data(fault_data)
        cv_iterator.create_partitions()
        indices = numpy.arange(num_pts)
        if fault_data is not None:
            indices = indices[numpy.where(fault_data==0)[0]]
        v_indices = split_sequence(indices,num_folds)
        for it in xrange( cv_iterator.num_folds() ):
            training_indices, validation_indices = \
              cv_iterator.get_fold_indices( it )
            t_indices = numpy.setdiff1d(indices, v_indices[it])
            assert numpy.allclose( v_indices[it], validation_indices )
            assert numpy.allclose( t_indices, training_indices )
            assert numpy.allclose(
                numpy.unique(numpy.hstack((training_indices,
                                           validation_indices))),
                indices )

    def test_fold_partitioning( self ):
        """ 
        Test the creation of training and validation indices with and without 
        the presence of faults.
        """
        for num_pts in [15,16,17]:
            for num_folds in [4,5,6]:
                # Test partitioning when faults ARE NOT present
                self.fold_partitioning_test(num_pts, num_folds)
                for num_failures in xrange(1,4):
                    # Test partitioning when faults ARE present
                    I = numpy.random.permutation(num_pts)[:num_failures]
                    fault_data = numpy.zeros((num_pts),dtype=numpy.int32)
                    fault_data[I]=1
                    self.fold_partitioning_test(num_pts, num_folds, fault_data)

    def test_fold_partitioning_error_catching( self ):
        """
        Tests runtime errors are thrown when assumed conditions of partioning
        are violated.
        """
        # Test error is thrown when number of folds exceeds number of points
        num_pts = 10; num_folds = 11
        self.assertRaises(RuntimeError, self.fold_partitioning_test,
                          num_pts, num_folds)

        # Test error is thrown when number of folds exceeds number of points
        # after faults are removed.
        num_pts = 10; num_folds = 8; num_failures = 3
        fault_data = create_random_fault_data(num_pts, num_failures)
        self.assertRaises(RuntimeError, self.fold_partitioning_test,
                          num_pts, num_folds, fault_data)

    def fold_data_extraction_test(self, num_folds, num_pts, num_eq_per_pt,
                                      seed, fault_data=None):
        """
        Helper function for testing extracting training vand validation data
        for solving AX=B accounting for the presence of faults.

        Tests:
         * union of training indices and validation indices of rows of A
           correspond to the non-faulty indices of the primary equations of A
           Also test that their intersection is empty.
        
         * Submartices of A corresponding to validation and training indices 
           are extract correctly

         * Submartices of B corresponding to validation and training indices 
           are extract correctly
        """

        cv_iterator = LinearSystemCrossValidationIterator()
        # The following functions MUST be called in the order given
        cv_iterator.set_seed( seed ) # negative means indices are not randomized
        cv_iterator.set_num_points( num_pts )
        cv_iterator.set_num_folds( num_folds )
        cv_iterator.set_num_equations_per_point( num_eq_per_pt )
        if fault_data is not None:
            cv_iterator.set_fault_data(fault_data)
        cv_iterator.create_partitions()
        indices = numpy.arange(num_pts, dtype=numpy.int32)
        if fault_data is not None:
            indices = indices[numpy.where(fault_data!=1)[0]]
        A = numpy.random.normal( 0., 1., ( num_pts*num_eq_per_pt, 20 ) )
        b = numpy.arange( num_pts*num_eq_per_pt, dtype=float )
        for iter in xrange( cv_iterator.num_folds() ):
            training_indices, validation_indices = \
                cv_iterator.get_fold_indices( iter )
            assert numpy.allclose(
                numpy.unique(numpy.hstack((training_indices,validation_indices))), indices )

            I = get_active_row_indices(
                training_indices, fault_data, num_eq_per_pt, num_pts)
            A_train = cv_iterator.extract_matrix( A, training_indices )
            assert numpy.allclose( A_train, A[I,:] )
            b_train = cv_iterator.extract_values( b, training_indices )
            assert numpy.allclose( b_train, b[I] )
            J = get_active_row_indices(
                validation_indices, fault_data, num_eq_per_pt, num_pts)
            A_valid = cv_iterator.extract_matrix( A, validation_indices )
            assert numpy.allclose( A_valid, A[J,:] )
            b_valid = cv_iterator.extract_values( b, validation_indices )
            assert numpy.allclose( b_valid, b[J] )

    def test_fold_data_extraction(self):
        """
        Test the extraction of training vand validation data
        for solving AX=B accounting with and without the presence of faults.
        """
        num_pts = 10
        for num_folds in [4,6,10]:
            for num_eq_per_pt in [1,3]:
                #self.fold_data_extraction_test(
                #    num_folds, num_pts, num_eq_per_pt, seed=1)
                for num_failures in [1,2,3]:
                    for num_secondary_failures in [0,1,2,10]:
                        fault_data = create_random_fault_data(
                            num_pts, num_failures, num_secondary_failures)
                        if num_folds<=num_pts-num_failures:
                            self.fold_data_extraction_test(
                                num_folds, num_pts, num_eq_per_pt, seed=1,
                                fault_data=fault_data)
                    
    def linear_system_cross_validation_test(self, vand, rhs, num_folds,
                                            num_eq_per_pt, solver,
                                            regression_opts,
                                            residual_tols = None,
                                            fault_data=None ):
        """
        Helper function for teseting K fold cross validation for linear systems

        Tests:
         * The validation residuals on each fold are correct
        
         * The cross validation scores are correct. This implicilty test
           the interpolation of cross validation scores on each fold from
           an arbitray set of residual tolerances onto a common (unique)
           set of tolerances.
        """
        #rhs = rhs.squeeze()
        #assert rhs.ndim==1
        assert vand.shape[0]%num_eq_per_pt==0
        num_build_pts = vand.shape[0]/num_eq_per_pt
        #################################
        # Use in built cross validation #
        #################################
        cv_iterator = LinearSystemCrossValidationIterator()
        cv_iterator.set_linear_system_solver( solver )
        cv_opts = {'num-points':num_build_pts,'num-folds':num_folds}
        cv_opts = OptionsList(cv_opts)
        cv_opts.set("regression-opts", regression_opts)
        cv_iterator.run( vand, rhs.squeeze(), cv_opts )
        cv_scores = cv_iterator.get_scores()
        cv_fold_diffs = cv_iterator.get_fold_validation_residuals()
        assert len(cv_fold_diffs)==num_build_pts
        cv_fold_tols = cv_iterator.get_fold_tolerances()
        cv_unique_tols = cv_iterator.get_unique_tolerances()
        cv_fold_scores = cv_iterator.get_fold_scores()

        #################################
        # Use brute force cross validation and use as reference `truth'
        #################################
        # cv_iterator randomizes points incase there is a pattern, e.g
        # the points are from a tensor grid, to do test we must know what
        # random permutation of the points was used

        for rhs_num in range(rhs.shape[1]):        
            fold_tols = []
            fold_diffs = []
            for i in xrange( num_folds ):
                # partition the data into a training and validation set
                A_train, b_train, valid_ind = partition_data(
                    vand, rhs[:,rhs_num], cv_iterator, i, num_build_pts, num_eq_per_pt,
                    fault_data)
                # Compute the solution on the training data
                solver.solve(A_train, b_train, regression_opts)
                x = solver.get_solutions_for_all_regularization_params(0)
                # Test that the validation residual is stored correctly for
                # each linear model.
                fold_diffs.append((numpy.tile( rhs[:,rhs_num].reshape(rhs.shape[0],1),
                                      (1,x.shape[1]) )-numpy.dot( vand, x ) )[valid_ind])
                if rhs_num==rhs.shape[1]-1:
                    #cv_fold_diffs is only stored for last RHS
                    assert numpy.allclose(fold_diffs[i], cv_fold_diffs[i] )
                fold_tols.append(solver.get_residuals_for_all_regularization_params(0))

            # Determine the unique tolerances at which to compute cross
            # validation errors
            max_num_path_steps = 0
            unique_tols = numpy.empty((0), numpy.double)
            for i in xrange(num_folds):
                fold_tols[i]=fold_tols[i][::-1]
                if rhs_num==rhs.shape[1]-1:
                    # cv_fold_tols is only stored for last RHS
                    assert numpy.allclose(fold_tols[i],cv_fold_tols[i])

                unique_tols = numpy.concatenate((unique_tols,fold_tols[i]))
                num_path_steps = fold_diffs[i].shape[1]
                max_num_path_steps = max( max_num_path_steps,
                                           num_path_steps )
            unique_tols = numpy.unique( unique_tols )

            # There are often thousands of unique parameter values so take
            # a thinned subset of these. The size of the subset is controlled by
            # max_num_path_steps
            num_unique_res = unique_tols.shape[0]
            stride = num_unique_res / max_num_path_steps
            unique_tols = unique_tols[::stride]
            assert numpy.allclose(cv_unique_tols[rhs_num],unique_tols)

            unique_errors = numpy.empty( ( unique_tols.shape[0],
                                           num_folds ), numpy.double )
            for i in xrange( num_folds ):
                # for the current fold compute the cross validation error
                fold_tols_i = fold_tols[i]
                errors_i = numpy.sum(fold_diffs[i]**2, axis = 0)

                if unique_tols.shape[0] > 1:
                    # fold_tols is reversed internally by c++ code
                    # so we must reverse errors_i to be consistent
                    errors_i = errors_i[::-1]
                    if rhs_num==rhs.shape[1]-1:
                        # cv_errors only stored for last RHS
                        assert numpy.allclose(errors_i,cv_fold_scores[i])

                    # Enforce constant interpolation when interpolation is
                    # outside the range of fold_tols_i
                    if ( fold_tols_i[0] > unique_tols[0] ):
                        fold_tols_i = numpy.r_[ unique_tols[0], fold_tols_i ]
                        errors_i = numpy.r_[ errors_i[0], errors_i ]
                    if  ( fold_tols_i[-1] < unique_tols[-1] ):
                        fold_tols_i = numpy.r_[ fold_tols_i, unique_tols[-1] ]
                        errors_i = numpy.r_[ errors_i, errors_i[-1] ]

                    # Interpolate the cross validation errors onto a unique
                    # set of tolerances
                    poly = interpolate.interp1d( fold_tols_i, errors_i )
                    unique_errors[:,i] = poly( unique_tols )
                else:
                    unique_errors[:,i] = errors_i

            # Test that the cross validation scores computed by the `truth'
            # and the internal linear solver cross validation module are the same
            scores = numpy.sum( unique_errors, axis = 1 ) / num_build_pts
            assert numpy.allclose( scores, cv_scores[rhs_num] )


    def test_linear_system_cross_validation(self):
        """
        Test cross validation applied to least squares regression using
        the various supported methods: QR, LU and SVD regression.

        This test does not use the fast LSQCrossValidationIterator
        which solves least squares problem Ax=b 
        using x = inv(A'A)*A'*b computing inv(A'A) using cholesky factorization
        """
        num_pts = 11; num_eq_per_pt = 1 
        num_nonzeros = 3; num_cols = 10; num_rhs=2; noise_std=0.1
        num_folds = num_pts
        regression_types = [SVD_LEAST_SQ_REGRESSION, QR_LEAST_SQ_REGRESSION,
                            ORTHOG_MATCH_PURSUIT, LEAST_ANGLE_REGRESSION,
                            LASSO_REGRESSION]
        #LU_LEAST_SQ_REGRESSION] # not working
        for regression_type in regression_types:
            opts_dict = {'regression_type':regression_type,'verbosity':0}
            regression_opts = OptionsList(opts_dict)
            solver = regression_solver_factory(regression_opts)

            A, rhs, x_truth = get_linear_system(
                num_pts, num_eq_per_pt, num_cols, num_rhs,
                num_nonzeros, noise_std)
            solver.solve(A, rhs, regression_opts)
            x = solver.get_solutions_for_all_regularization_params(0)
            # todo make my own allclose function that calls squeeze before
            # comparison, like i do here
            #assert numpy.allclose(x.squeeze(),x_truth)
            opts = OptionsList()
            self.linear_system_cross_validation_test(
                A, rhs, num_folds, num_eq_per_pt, solver, regression_opts)

    def xtest_equality_constrained_least_squares_cross_validation(self):
        """
        Test cross validation applied to equality constrained regression
        """
        num_pts = 11; num_eq_per_pt = 1 
        num_nonzeros = 3; num_cols = 10; num_rhs=1; noise_std=0.1
        num_folds = num_pts; num_primary_eq = 9
        regression_type = EQ_CONS_LEAST_SQ_REGRESSION
        opts_dict = {'regression_type':regression_type,
                     'num-primary-equations':num_primary_eq}
        regression_opts = OptionsList(opts_dict)
        solver = regression_solver_factory(regression_opts)

        A, rhs, x_truth = get_linear_system(
            num_pts, num_eq_per_pt, num_cols, num_rhs,
            num_nonzeros, noise_std)
        solver.solve(A, rhs, regression_opts)
        x = solver.get_solutions_for_all_regularization_params(0)
        # todo make my own allclose function that calls squeeze before
        # comparison, like i do here
        #assert numpy.allclose(x.squeeze(),x_truth)
        opts = OptionsList()
        self.linear_system_cross_validation_test(
            A, rhs, num_folds, num_eq_per_pt, solver, regression_opts)


    def test_equality_constrained_least_squares_cross_validation_error_catching(self):
        """
        Test that error is thrown when cross validation applied to 
        equality constrained regression induces an underdetermined system
        """
        num_pts = 10; num_eq_per_pt = 1 
        num_nonzeros = 3; num_cols = 10; num_rhs=1; noise_std=0.001
        num_folds = num_pts; num_primary_eq = 10
        regression_type = EQ_CONS_LEAST_SQ_REGRESSION
        opts_dict = {'regression_type':regression_type,
                     'num-primary-equations':num_primary_eq}
        regression_opts = OptionsList(opts_dict)
        solver = regression_solver_factory(regression_opts)

        A, rhs, x_truth = get_linear_system(
            num_pts, num_eq_per_pt, num_cols, num_rhs,
            num_nonzeros, noise_std)
        solver.solve(A, rhs, regression_opts)
        x = solver.get_solutions_for_all_regularization_params(0)
        # todo make my own allclose function that calls squeeze before
        # comparison, like i do here
        assert numpy.allclose(x,x_truth,atol=noise_std*10)
        self.assertRaises(RuntimeError,self.linear_system_cross_validation_test,
                          A, rhs, num_folds, num_eq_per_pt, solver,
                          regression_opts)

class TestLSQCrossValidation(unittest.TestCase):
    """
    Tests of fast cross validation for solving least squares problem Ax=b 
    using x = inv(A'A)*A'*b computing inv(A'A) using cholesky factorization.
    
    When using leave one out cross validation the validation residual on
    each point are:
    
    H = numpy.dot(A,numpy.dot(numpy.linalg(numpy.dot(A.T,A)),A.T))
    validation_rediduals = residuals/(1.-numpy.diag(H)[:,numpy.newaxis])

    When using leave k out cross validation the validation residuals on the ith
    fold are:

    validation_indices = cv_iterator.get_fold_validation_indices(i)
    A_valid = cv_iterator.extract_matrix(A, validation_indices)
    residuals_valid=cv_iterator.extract_values(
    residuals, validation_indices)
    I = numpy.eye(validation_indices.shape[0])
    H_valid = I - numpy.dot(A_valid,numpy.dot(AtA_inv,A_valid.T))
    H_valid_inv = numpy.linalg.inv(H_valid)
    validation_residuals_i = numpy.dot(H_valid_inv, residuals_valid)

    """
    
    def lsq_cross_validation_test(self, A, rhs, num_pts, num_folds, seed=0):
        """
        Helper function for testing LSQCrossValidationIterator.

        Tests:
         * Test validation residuals on each fold are correct

         * Test validation residuals on each fold are consistent with
           those are produced by LinearSystemCrossValidationIterator using 
           SVD regression

         * Test cross validation scores are consistent with
           those are produced by LinearSystemCrossValidationIterator using 
           SVD regression

        """
        lsq_cv_iterator = LSQCrossValidationIterator()
        cv_opts = {'num-points':num_pts, 'num-folds':num_folds, 'seed':seed}
        cv_opts = OptionsList(cv_opts)
        regression_opts = {'regression_type':SVD_LEAST_SQ_REGRESSION}
        regression_opts = OptionsList(regression_opts)
        cv_opts.set("regression-opts", regression_opts)
        lsq_cv_iterator.run(A, rhs, cv_opts)
        fold_diffs = lsq_cv_iterator.get_fold_validation_residuals()
        x = lsq_cv_iterator.get_coefficients()

        solver = regression_solver_factory(regression_opts)
        cv_iterator = LinearSystemCrossValidationIterator()
        cv_iterator.set_linear_system_solver( solver )
        cv_iterator.run( A, rhs, cv_opts )
        brute_force_fold_diffs = cv_iterator.get_fold_validation_residuals()
        assert len(fold_diffs)==num_folds
        for i in xrange(len(fold_diffs)):
            # Get fold validation residuals using analytical formula
            # which is implemented by lsq_cv_iterator. I.e.
            AtA_inv = numpy.linalg.inv(numpy.dot(A.T,A))
            residuals = rhs - numpy.dot(A,x)
            if num_pts==num_folds:
                H = numpy.dot(A,numpy.dot(AtA_inv,A.T))
                cv_diffs = \
                  residuals[i:i+1]/(1.-numpy.diag(H)[i:i+1,numpy.newaxis])
            else:
                validation_indices = cv_iterator.get_fold_validation_indices(i)
                A_valid = cv_iterator.extract_matrix(A, validation_indices)
                residuals_valid=cv_iterator.extract_values(
                    residuals, validation_indices)
                I = numpy.eye(validation_indices.shape[0])
                H_valid = I - numpy.dot(A_valid,numpy.dot(AtA_inv,A_valid.T))
                H_valid_inv = numpy.linalg.inv(H_valid)
                cv_diffs = numpy.dot(H_valid_inv, residuals_valid)

            # Test the validation residuals on this fold are correct
            assert numpy.allclose(cv_diffs, fold_diffs[i])
            
            # Test validation residuals on this fold are consistent with
            # those are produced by LinearSystemCrossValidationIterator
            # cv_iterator.get_fold_validation_residuals() only contains diffs
            # from last rhs column considered
            assert numpy.allclose(
                fold_diffs[i][0,-1],brute_force_fold_diffs[i][0])
        # Test cross validation scores are consistent with
        # those are produced by LinearSystemCrossValidationIterator
        assert numpy.allclose(
            cv_iterator.get_scores(),lsq_cv_iterator.get_scores())
    
    def test_leave_one_out_cross_validation(self):
        """
        Test LSQCrossValidationIterator leave one our cross validation
        with mutiple RHS
        """
        num_pts = 11; num_eq_per_pt = 1
        num_nonzeros = 3; num_cols = 10; num_rhs = 2; noise_std=0.1; 
        num_folds = num_pts; 

        A, rhs, x_truth = get_linear_system(
            num_pts, num_eq_per_pt, num_cols, num_rhs, num_nonzeros, noise_std)
        self.lsq_cross_validation_test(A, rhs, num_pts, num_folds, seed=-1)

    def test_leave_K_out_cross_validation(self):
        """
        Test LSQCrossValidationIterator K fold cross validation
        with mutiple RHS
        """
        num_pts = 20; num_eq_per_pt = 1
        num_nonzeros = 3; num_cols = 10; num_rhs = 2; noise_std=0.1; 
        num_folds = 4; 

        A, rhs, x_truth = get_linear_system(
            num_pts, num_eq_per_pt, num_cols, num_rhs, num_nonzeros, noise_std)
        self.lsq_cross_validation_test(A, rhs, num_pts, num_folds, seed=-1)

class TestCrossValidatedSolver(unittest.TestCase):
    """
    Test that when cross-validation is requested for a linear system solver
    a CrossValidedSolver object is created and the object finds the best
    cross validated regularization parameters and returns the solution built
    on the entier data set corresponding to the best parameters.
    """

    # TODO: this test will break with older NumPy, e.g., RHEL7's NumPy 1.7.1
    # Should work with NumPy 1.8 or newer
    def cross_validated_solver_wrappers_of_lscv_iterator_helper(
            self, regression_type, regression_opts, A, rhs):
        
        # Get solutions using cross validated solver
        cv_solver = CrossValidatedSolver()
        cv_solver.set_linear_system_solver(regression_type)
        cv_solver.solve(A, rhs, regression_opts)
        solutions = cv_solver.get_final_solutions()
        scores = cv_solver.get_best_scores()

        # make sure the best scores returned by cross validated solver
        # are the same as the best scores computed by its internal
        # cross validation iterator
        cv_iterator = cv_solver.get_cross_validation_iterator()
        cv_scores = cv_iterator.get_scores()
        best_score_indices = cv_iterator.get_best_score_indices()
        num_rhs = rhs.shape[1]
        for i in xrange(num_rhs):
            assert numpy.allclose(cv_scores[i].min(),scores[i])
            assert numpy.argmin(cv_scores[i])==best_score_indices[i]

        # For each rhs extract residual tolerances of
        # best cross validated solution and run a new solver instance
        # (of the same type wrapped by the cross validated solver)
        # and check it produces the same solution when run
        # on the entire data set
        for i in xrange(num_rhs):

            # run new instance of correct linear system solver
            local_regression_opts = {'regression_type':regression_type,
                                '   verbosity':0}
            local_regression_opts = OptionsList(local_regression_opts)
            solver = regression_solver_factory(local_regression_opts)
            solver.solve(A, rhs, local_regression_opts)

            # check residuals of entire data set solutions are the same
            cv_solutions = cv_solver.get_final_solutions()
            cv_residual_norms = cv_solver.get_final_residuals()
            residual_norms = numpy.linalg.norm(
                rhs-numpy.dot(A,solutions),axis=0)
            assert numpy.allclose(residual_norms, cv_residual_norms)

            # Check the residuals used to compute final solution
            # (adjusted residuals) account
            # difference in size between the training sets
            # used in cross validation and the entire data set
            cv_iterator = cast_linear_cv_iterator(
                cv_iterator, regression_type)
            cv_adjusted_residuals = \
              cv_iterator.get_adjusted_best_residual_tolerances()
            cv_best_residual_tols=cv_iterator.get_best_residual_tolerances()
            cv_opts = regression_opts.get("cv-opts")
            num_folds = cv_opts.get("num-folds")
            assert numpy.allclose(
                cv_adjusted_residuals,
                cv_best_residual_tols*num_folds/(num_folds-1.0))

            # check that the residual tolerances of the best cv solution
            # are set correctly in the linear system solver
            # contained in the cross validated solver instance
            derived_class = cast_linear_system_solver(
                cv_iterator.get_linear_system_solver(),regression_type)
            assert numpy.allclose(derived_class.get_residual_tolerances(),
                                  cv_adjusted_residuals)

            # Check that the l2 norm of the residual of the last solution
            # computed using the entire data set is the largest
            # residual norm that is smaller than the adjusted residuals
            residuals = \
              solver.get_residuals_for_all_regularization_params(i)
            I = numpy.where(residuals>=cv_adjusted_residuals[i])[0]
            I = I[-1]+1
            assert cv_residual_norms[i]>=residuals[I]
            print  cv_residual_norms[i],residuals[I],residuals[I-1],I
            solver_solutions = \
              solver.get_solutions_for_all_regularization_params(i)[:,I]

            # check that the final solutions obtained by the cross
            # validated solver and the new solver instance are the same
            assert numpy.allclose(
                solver_solutions, cv_solver.get_final_solutions()[:,i])

    
    def test_cross_validated_solver_wrappers_of_lscv_iterator(self):
        """
        Make sure the CrossValidatedSolver returns the solution associated
        with the best cross validation error found using linear system
        solvers that envoke LinearSystemCrossValidationIterator
        """
        num_pts = 20; num_eq_per_pt = 1
        num_nonzeros = 3; num_cols = 10; num_rhs = 2; noise_std=0.1; 
        num_folds = 4; 
        
        A, rhs, x_truth = get_linear_system(
            num_pts, num_eq_per_pt, num_cols, num_rhs, num_nonzeros, noise_std)

        cv_opts = {'num-points':num_pts,'num-folds':num_folds}
        cv_opts = OptionsList(cv_opts)
        regression_opts = {'verbosity':0}
        regression_opts = OptionsList(regression_opts)
        regression_opts.set("cv-opts", cv_opts)
        regression_types = [ORTHOG_MATCH_PURSUIT, LEAST_ANGLE_REGRESSION,
                            LASSO_REGRESSION]
        for regression_type in regression_types:
            self.cross_validated_solver_wrappers_of_lscv_iterator_helper(
                regression_type, regression_opts, A, rhs)

    def test_cross_validated_solver_wrappers_of_lsqcv_iterator(self):
        #        SVD_LEAST_SQ_REGRESSION, QR_LEAST_SQ_REGRESSION,
        assert False, "test not implemented"
def single_test_suite():
    suite = unittest.TestSuite()
    suite.addTest(TestCrossValidatedSolver( "test_cross_validated_solver_wrappers_of_lscv_iterator" ) )
    return suite      

            
if __name__ == '__main__':
    numpy.random.seed(2)
    runner = unittest.TextTestRunner()
    #single = single_test_suite()
    #runner.run( single )

    unittest.main()
