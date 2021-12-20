/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/* math_tools.i */

%define %math_tools_docstring
"
PyDakota.math_tools is the python interface to the Dakota miscelaneous math
tools and utilities

The purpose of this package is to provide a number of utilities often
needed for mainpulating operatring on matrices and vectors
"
%enddef

%module(package="PyDakota",
        directors=1,
        implicitconv="1",
        autodoc="1",
        docstring=%math_tools_docstring) math_tools
%{
  #define SWIG_FILE_WITH_INIT
  
  // #define NO_IMPORT_ARRAY causes missing symbol error when loading
  // _math_tools module when used with %fragment("NumPy_Fragments");
  // Folloing not needed when using %fragment("NumPy_Fragments");
  //#include "numpy_include.hpp"
  
// Local includes
#include "math_tools.hpp"
#include "linear_algebra.hpp"
#include "teuchos_data_types.hpp"

using namespace Pecos;
%}

// Global swig features
%feature("autodoc", "1");
%feature("compactdefaultargs");

// %ignore and rename must be included before the function decleration, i.e.
// before the %include
%ignore *::operator[];
%ignore *::operator=;
%ignore *::print;

// We utilize very little of numpy.i, so some of the fragments
// do not get automatically instantiated.  This forces the issue.
%include "numpy.i"
 //Following needed to include obj_to_array_fortran_allow_conversion
%fragment("NumPy_Fragments");
// Following is needed otherwise will get segfault when wrapped functions in
// module are run
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

// include typemaps for scalar outputs etc
%include <typemaps.i>
%apply( int &OUTPUT ){ int &rank };
%apply( double &INOUT ){ double &rcond };

// include Teuchos enums, such as TRANS, NO_TRANS
%include "Teuchos_BLAS_types.hpp"

// include typemaps to convert to from NumPy Arrays to
// Teuchos::SerialDenseVector/Matrix
%include "Teuchos_SerialDenseVector.i"
%include "Teuchos_SerialDenseMatrix.i"
%include "data_structures.i"

%ignore qr_solve( const RealMatrix &, const RealMatrix &, RealMatrix & );
%ignore svd_solve( RealMatrix &A, RealMatrix &B, RealMatrix &result_0, RealVector &result_1, int &rank );
%rename(tensor_product_indices_cpp) tensor_product_indices;
%rename(svd_solve_cpp) svd_solve;
%rename (cholesky_solve_cpp) cholesky_solve;
%rename (cholesky_cpp) cholesky;
%rename (solve_using_cholesky_factor_cpp) solve_using_cholesky_factor;
%ignore equality_constrained_least_squares_solve( RealMatrix &A, RealVector &b, RealMatrix &C, RealVector &d, RealMatrix &result);
%rename(complete_pivoted_lu_factorization_cpp) complete_pivoted_lu_factorization;
%rename (truncated_pivoted_lu_factorization_cpp) truncated_pivoted_lu_factorization;
%rename (ind2sub_cpp) ind2sub;
%rename (sub2ind_cpp) sub2ind;


// typemaps are applied to special instances of SerialDenseVecor/Matrix and
// std::vectors of these classes in the Teuchos_SerialDenseVector/Matrix.i
// and data_strtuctures.i files. argout typemaps applied to
// any non-const argument passed by reference named argout, result, result_0,
// result_1, for IntVector/Matrix, std::vector<IntVector/Matrix> similarly for
// RealVector/Matrix, std::vector<RealVector/Matrix>. Additional applys can be
// applied here between ----or in the .i files
// ----
// %apply() go here
// ----

// Must specify here to ensure that functions involving the function with
// parameters renamed using typedef can be wrapped.
// If no match is found using the above rules SWIG applies a typedef
// reduction to the type and repeats the typemap search for the reduced type.
// The original %apply statements (in this file or in data_structures.i,
// Teuchos_SerialDenseVector/Matrix.i) should not use these typemaps. If they
// do the functions that are explictly written with these typemaps will be
//wrapped correctly, but other functions (such as templated functions,
// which cannot use typemaps), will not be wrapped correctly
typedef double Real;
typedef Teuchos::SerialDenseVector<int,double> RealVector;
typedef Teuchos::SerialDenseVector<int,int> IntVector;
typedef Teuchos::SerialDenseVector<int,Complex> ComplexVector;
typedef Teuchos::SerialDenseMatrix<int,double> RealMatrix;
typedef Teuchos::SerialDenseMatrix<int,int> IntMatrix;
typedef Teuchos::SerialDenseMatrix<int,Complex> ComplexMatrix;
// ----
// add extra typedefs here
// ----

%include "math_tools.hpp"
%include "linear_algebra.hpp"

namespace Pecos {
namespace util {
%template(cartesian_product_int) cartesian_product<int,int>;
%template(cartesian_product_double) cartesian_product<int,double>;
%template(outer_product_int) outer_product<int,int>;
%template(outer_product_double) outer_product<int,double>;
%template(binary_search_double) binary_search<int,double>;
%template(binary_search_int) binary_search<int,int>;
%template(range_double) range<int,double>;
%template(range_int) range<int,double>;
}
}


%pythoncode %{
import numpy
def remove_common_columns( A, B ):
    """
    Remove columns from A that are also in B
    """
    D = numpy.hstack( ( B, A ) )
    order = numpy.lexsort( D )
    C = D.copy()[:,order]
    I = numpy.hstack((numpy.arange(B.shape[1])+1,-numpy.arange(A.shape[1]) ))
    I = I[order]
    diff = numpy.diff( C, axis = 1 )
    ui = numpy.ones( C.shape[1], 'bool' )

    ui[1:] = ( numpy.absolute( diff ) >
               numpy.finfo( numpy.double ).eps ).any( axis = 0 )
    ui = ( I[ui] <= 0 )
    return A[:,ui]

def unique_matrix_rows(A):
    """
    Remove duplicate columns from A
    """
    return numpy.vstack(set(tuple(row) for row in A))

def unique_matrix_cols(A):
    return unique_matrix_cols(A.T).T

def cartesian_product(input_sets,elem_size=1):
    """Wrapper of cpp function that converts to accepted types."""
    if input_sets[0].dtype==float or input_sets[0].dtype==numpy.float64:
        return cartesian_product_double(input_sets,elem_size)
    elif (input_sets[0].dtype==numpy.int64 or input_sets[0].dtype==numpy.int32 or input_sets[0].dtype==int):
        for i in xrange(len(input_sets)):
            if input_sets[i].dtype!=numpy.int32:
                input_sets[i] = numpy.asarray(input_sets[i],dtype=numpy.int32)
        return cartesian_product_int(input_sets,elem_size)
    else:
        raise Exception, 'element type not supported'

def outer_product(input_sets,elem_size=1):
    """Wrapper of cpp function that converts to accepted types."""
    if input_sets[0].dtype==float or input_sets[0].dtype==numpy.float64:
        return outer_product_double(input_sets)
    elif (input_sets[0].dtype==numpy.int64 or input_sets[0].dtype==numpy.int32 or input_sets[0].dtype==int):
        for i in xrange(len(input_sets)):
            if input_sets[i].dtype!=numpy.int32:
                input_sets[i] = numpy.asarray(input_sets[i],dtype=numpy.int32)
        return outer_product_int(input_sets)
    else:
        raise Exception, 'element type not supported'

def tensor_product_indices(degrees):
    """Wrapper of cpp function that converts to accepted integer type."""
    assert degrees.dtype==numpy.int32 or degrees.dtype==numpy.int64
    if degrees.dtype!=numpy.int32:
        degrees = numpy.asarray(degrees,dtype=numpy.int32)
    return tensor_product_indices_cpp(degrees)

def cholesky_solve( A, b, rcond=-1. ):
    return cholesky_solve_cpp( A, b.reshape(b.shape[0],1), rcond )[1];

def cholesky(A):
    return cholesky_cpp( A, 1, False )[1];

def solve_using_cholesky_factor(L, b, uplo):
    return solve_using_cholesky_factor_cpp(L, b, uplo)[1];

def svd_solve( A, b, rcond=-1. ):
    return svd_solve_cpp( A, b.reshape(b.shape[0],1), rcond );

def complete_pivoted_lu_factorization(A, max_iters=None):
      if max_iters is None:
          max_iters = min( A.shape[0], A.shape[1] )
      return complete_pivoted_lu_factorization_cpp( A, max_iters )

def truncated_pivoted_lu_factorization(A, max_iters=None, num_initial_rows=None):
      if max_iters is None: max_iters = min( A.shape[0], A.shape[1] )
      if num_initial_rows is None: num_initial_rows = 0
      return truncated_pivoted_lu_factorization_cpp(
          A, max_iters, num_initial_rows)

def ind2sub(sizes, index, num_elems):
    # TODO: discuss with wfspotz how to flexibly support various integer arrays
    assert sizes.dtype==numpy.int32 or sizes.dtype==numpy.int64
    sizes = numpy.asarray(sizes, dtype=numpy.int32)
    if isinstance(index,(int,numpy.integer)):
        index = numpy.int(index)
    else:
        raise Exception('index is not an integer')
    if isinstance(num_elems,(int,numpy.integer)):
        num_elems = numpy.int(num_elems)
    else:
        raise Exception('num_elems is not an integer')
    return ind2sub_cpp(sizes,index,num_elems)

def sub2ind(sizes, multi_index):
    # TODO: discuss with wfspotz how to flexibly support various integer arrays
    assert sizes.dtype==numpy.int32 or sizes.dtype==numpy.int64
    sizes = numpy.asarray(sizes, dtype=numpy.int32)
    return sub2ind_cpp(sizes, multi_index)

   %}
