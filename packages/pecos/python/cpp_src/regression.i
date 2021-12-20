/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/* regression.i */

%define %regression_docstring
"
PyDakota.regression is the python interface to the Dakota utilities for
solving  regression problems.

This package provides tools for
  * Least squares regression, using SVD, QR, CHOLESKY, and LU decompositions
  * Orthgonal matching pursuit (OMP)
  * Least angle regession (LAR)
"
%enddef

%module(package="PyDakota",
        directors="1",
        implicitconv="1",
        autodoc="1",
        docstring=%regression_docstring) regression

%feature("director") LinearSystemSolver;
%{
  #define SWIG_FILE_WITH_INIT
  //unlike math_tools we must include following file here
  #include "numpy_include.hpp"

// Local includes
#include "linear_solvers.hpp"
#include "OMPSolver.hpp"
#include "LARSolver.hpp"
#include "LSQSolver.hpp"
#include "EqConstrainedLSQSolver.hpp"
#include "CrossValidatedSolver.hpp"
#include "CrossValidationIterator.hpp"
#include "LinearSystemCrossValidationIterator.hpp"
#include "LSQCrossValidationIterator.hpp"
#include "CrossValidatedSolver.hpp"

  // Following needed to allow OptionsList to be used and function such as
  // pyDictToNewOptionsList to be found
  #include "python_helpers.hpp"


using namespace Pecos;
%}

// We utilize very little of numpy.i, so some of the fragments
// do not get automatically instantiated.  This forces the issue.
%include "numpy.i"
// Following is needed otherwise will get segfault when wrapped functions in */
// module is run
%init %{
  import_array();
%}

// Standard exception handling
%include "exception.i"
%include "std_except.i"
%exception
{
  try{
    // use these print statements to debug
    //printf("Entering function : $name\n");
    $action
      if (PyErr_Occurred()) SWIG_fail;
    //else    printf("Exiting function : $name\n");
  }
  SWIG_CATCH_STDEXCEPT
    catch (...) {
    SWIG_exception(SWIG_UnknownError, "unknown C++ exception");
    throw(std::runtime_error(""));
  }
}

// %ignore and rename must be included before the function decleration, i.e.
// before the %include
%ignore *::operator[];
%ignore *::operator=;
%ignore *::print;

// include Teuchos enums, such as TRANS, NO_TRANS
%include "Teuchos_BLAS_types.hpp"

%rename(extract_values_cpp) extract_values;
// note following is applied to that class and all derived classes
%rename(solve_cpp) Pecos::util::LinearSystemSolver::solve;
%rename(run_cpp) Pecos::util::LinearSystemCrossValidationIterator::run;

// importing math_tools.i avoids need to %include
// Teuchos_SerialDenseVector?Matrix.i and data_Structures.i
%import "math_tools.i"
%import "OptionsList.i"

%shared_ptr(Pecos::util::LinearSystemSolver)
%shared_ptr(Pecos::util::SparseSolver)
%shared_ptr(Pecos::util::OMPSolver)
%shared_ptr(Pecos::util::LARSolver)
%shared_ptr(Pecos::util::LSQSolver)
%shared_ptr(Pecos::util::EqConstrainedLSQSolver)
%shared_ptr(Pecos::util::CrossValidatedSolver)

%shared_ptr(Pecos::util::CrossValidatedSolver)
%shared_ptr(Pecos::util::CrossValidationIterator)
%shared_ptr(Pecos::util::LinearSystemCrossValidationIteratorBase)
%shared_ptr(Pecos::util::LinearSystemCrossValidationIterator)
%shared_ptr(Pecos::util::LSQCrossValidationIterator)

// %%include of a base class needs to be called before %include for derived
// classes. LinearSolver base also needed to get enums. perhaps move enums to a type definitions file
%include "LinearSystemSolver.hpp" 
%include "OMPSolver.hpp"
%include "LARSolver.hpp"
%include "LSQSolver.hpp"
%include "EqConstrainedLSQSolver.hpp"

%include "CrossValidationIterator.hpp"
%include "LinearSystemCrossValidationIterator.hpp"
%include "LSQCrossValidationIterator.hpp"

%include "CrossValidatedSolver.hpp"
%include "linear_solvers.hpp"

%extend Pecos::util::CrossValidationIterator {
  %pythoncode
    %{
    def extract_values( self, values, indices ):
        if ( values.ndim==1):
            values  = values.reshape(values.shape[0],1)
        return self.extract_values_cpp( values, indices ).squeeze()
    %}
}

%extend Pecos::util::LinearSystemSolver {
  %pythoncode
    %{
     def solve( self, A, rhs, opts ):
        if (rhs.ndim==1):
            rhs = rhs.reshape(rhs.shape[0],1)
        return self.solve_cpp( A, rhs, opts )
    %}
}

%extend Pecos::util::LinearSystemCrossValidationIterator {
  %pythoncode
    %{
     def run( self, A, rhs, opts ):
        if (rhs.ndim==1):
            rhs = rhs.reshape(rhs.shape[0],1)
        return self.run_cpp( A, rhs, opts )
    %}
}

%pythoncode
%{import numpy
def cast_linear_cv_iterator(cv_iterator, regression_type):
    if regression_type in [ORTHOG_MATCH_PURSUIT, LEAST_ANGLE_REGRESSION,
                          LASSO_REGRESSION]:
        return cast_to_linear_system_cross_validation_iterator(cv_iterator)
    else:
        return cast_to_least_squares_cross_validation_iterator(cv_iterator)

def cast_linear_system_solver(solver, regression_type):
    if regression_type==ORTHOG_MATCH_PURSUIT:
        return cast_linear_system_solver_to_ompsolver(solver)
    elif regression_type in [LEAST_ANGLE_REGRESSION,LASSO_REGRESSION]:
        return cast_linear_system_solver_to_larssolver(solver)
    elif regression_type == EQ_CONS_LEAST_SQ_REGRESSION:
        return cast_linear_system_solver_to_equalityconstrainedlsqsolver(
            solver)
    elif regression_type in [SVD_LEAST_SQ_REGRESSION,QR_LEAST_SQ_REGRESSION,
                             LU_LEAST_SQ_REGRESSION]:
        return cast_linear_system_solver_to_lsqsolver(solver)
    else:
        raise Exception, 'incorrect regression_type specified'

%}
