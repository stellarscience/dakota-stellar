/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/* fundamentals.i */
%{
// Required for interfacing with NumPy
// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
#include <vector> 
// Namespace flattening
using std::string;
%}

// Namespace flattening
using std::string;

// Standard exception handling
%include "exception.i"

// Global swig features
%feature("autodoc", "1");
%feature("compactdefaultargs");

// C++ STL support.  If the wrapped class uses standard template
// library containers, the following %include wraps the containers
// and makes certain conversions seamless, such as between std::string
// and python strings.
%include "stl.i"
%include "std_except.i"

 // General exception handling
%feature("director:except")
{
  if ($error != NULL) {
    throw Swig::DirectorMethodException();
  }
}
%exception
{
  try{
   $action
   if (PyErr_Occurred()) SWIG_fail;
 }
 catch (Swig::DirectorException &e){
   SWIG_fail;
 }
 SWIG_CATCH_STDEXCEPT
 catch (...) {
   SWIG_exception(SWIG_RuntimeError, "unknown C++ exception");
 }
}

// General ignore directives
%ignore *::operator[];
%ignore *::operator=;
%ignore *::print;

// We utilize very little of numpy.i, so some of the fragments
// do not get automatically instantiated.  This forces the issue.
%include "numpy.i"
%fragment("NumPy_Fragments");
%init %{
  import_array();
%}

%include "std_complex.i"
%include "std_vector.i"

%include <std_shared_ptr.i>
%shared_ptr(std::basic_ostream)
%shared_ptr(std::ostream)

%include "Teuchos_SerialDenseVector.i"
%include "Teuchos_SerialDenseMatrix.i"
%include "data_structures.i"

// Must specify here to ensure that functions involving the function with
// parameters renamed using typedef can be wrapped.
// If no match is found using the above rules SWIG applies a typedef
// reduction to the type and repeats the typemap search for the reduced type
typedef double Real;
typedef Teuchos::SerialDenseVector<int,double> RealVector;
typedef Teuchos::SerialDenseVector<int,int> IntVector;
typedef Teuchos::SerialDenseVector<int,Complex> ComplexVector;
typedef Teuchos::SerialDenseMatrix<int,double> RealMatrix;
typedef Teuchos::SerialDenseMatrix<int,int> IntMatrix;
typedef Teuchos::SerialDenseMatrix<int,Complex> ComplexMatrix;

%apply IntVector &argout { IntVector &result }
%apply IntVector &argout { IntVector &result_0 }
%apply IntVector &argout { IntVector &result_1 }
%apply IntMatrix &argout { IntMatrix &result }
%apply IntMatrix &argout { IntMatrix &result_0 }
%apply IntMatrix &argout { IntMatrix &result_1 }

%apply RealVector &argout { RealVector &result }
%apply RealVector &argout { RealVector &result_0 }
%apply RealVector &argout { RealVector &result_1 }
%apply RealMatrix &argout { RealMatrix &result }
%apply RealMatrix &argout { RealMatrix &result_0 }
%apply RealMatrix &argout { RealMatrix &result_1 }

%apply std::vector<RealVector> &argout {std::vector<RealVector> &result}
%apply std::vector<RealVector> &argout {std::vector<RealVector> &result_0}
%apply std::vector<RealVector> &argout {std::vector<RealVector> &result_1}
%apply std::vector<RealMatrix>  &argout {std::vector<RealMatrix> &result}
%apply std::vector<RealMatrix>  &argout {std::vector<RealMatrix> &result_0}
%apply std::vector<RealMatrix>  &argout {std::vector<RealMatrix> &result_1}

%apply std::vector<IntVector> &argout {std::vector<IntVector> &result}
%apply std::vector<IntVector> &argout {std::vector<IntVector> &result_0}
%apply std::vector<IntVector> &argout {std::vector<IntVector> &result_1}
%apply std::vector<IntMatrix>  &argout {std::vector<IntMatrix> &result}
%apply std::vector<IntMatrix>  &argout {std::vector<IntMatrix> &result_0}
%apply std::vector<IntMatrix>  &argout {std::vector<IntMatrix> &result_1}

