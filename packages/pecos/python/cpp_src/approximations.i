/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/* approximations.i */

%define %approximation_docstring
"
PyDakota.approximations is the python interface to the Dakota approximation
classes and training methods.
"
%enddef

%module(package      = "PyDakota",
        directors    = "1",
        autodoc      = "1",
        implicitconv = "1",
        docstring    = %approximation_docstring) approximation

%feature("director") Function;

// The following code is inserted directly into the wrapper file
%{
#include "numpy_include.hpp"//gets rid of deprecated warnings

// Approximation includes
#include "Function.hpp"
#include "CppFunction.hpp"
#include "Approximation.hpp"
#include "PolyApproximation.hpp"
#include "Monomial.hpp"
#include "PolynomialChaosExpansion.hpp"
#include "Variables.hpp"
#include "BoundedVariables.hpp"
#include "VariableTransformation.hpp"
#include "AffineVariableTransformation.hpp"
#include "polynomial_approximation_drivers.hpp"
#include "RegressionBuilder.hpp"
#include "PCEFactory.hpp"
  //  using std::string;
using namespace Pecos;

#include "typedefs_for_python_wrapper.hpp"
#include <vector> 

  // Following needed to allow OptionsList to be used and function such as
  // pyDictToNewOptionsList to be found
  #include "python_helpers.hpp"
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

// include typemaps to convert to from NumPy Arrays to
// Teuchos::SerialDenseVector/Matrix
/* %include "Teuchos_SerialDenseVector.i" */
/* %include "Teuchos_SerialDenseMatrix.i" */
/* %include "data_structures.i" */

%include "std_vector.i"
%template() std::vector<short>;
%include "typedefs_for_python_wrapper.hpp"

%include <std_shared_ptr.i>
%shared_ptr(std::basic_ostream)
%shared_ptr(std::ostream)

%shared_ptr(Pecos::surrogates::Function)
%shared_ptr(Pecos::surrogates::CppFunction)
%shared_ptr(Pecos::surrogates::Approximation)
%shared_ptr(Pecos::surrogates::PolyApproximation)
%shared_ptr(Pecos::surrogates::Monomial)
%shared_ptr(Pecos::surrogates::PolynomialChaosExpansion)
%shared_ptr(Pecos::surrogates::Variables)
%shared_ptr(Pecos::surrogates::BoundedVariables)
%shared_ptr(Pecos::surrogates::VariableTransformation)
%shared_ptr(Pecos::surrogates::AffineVariableTransformation)

// importing math_tools.i avoids need to %include
// Teuchos_SerialDenseVector?Matrix.i and data_Structures.i
%import "math_tools.i"

%import "OptionsList.i"

%include "Function.hpp"
%include "CppFunction.hpp"
%include "Variables.hpp"
%include "BoundedVariables.hpp"
%include "VariableTransformation.hpp"
%include "AffineVariableTransformation.hpp"
// If classes uses other classes then the classes it use must
// be declared first. E.g. Approximation.hpp calls VariableTransformation.hpp
// Also if derived class does not have redefinition of virtual function in
// base class The base class must be %include so that python can see the
 // virtual function
%include "Approximation.hpp"
%include "PolyApproximation.hpp"
%include "Monomial.hpp"
%include "PolynomialChaosExpansion.hpp"
%include "SurrogateBuilder.hpp"
%include "RegressionBuilder.hpp"

%include "polynomial_approximation_drivers.hpp"
%include "PCEFactory.hpp"

%extend Pecos::surrogates::Function{
%pythoncode %{
def __reduce__(self):
    return (type(self), (None, ))
%}
}

%pythoncode %{
import numpy
class PyFunction(Function):
    def __init__(self,target_function):
        """
        Parameters
        ----------
        target_function : callable function
            Calls to target funcation are assumed to follow
            vals = target_function(sample). Where vals
            is a 1D array and sample is a 1D array or scalar.
        """
        Function.__init__(self)
        self.target_function = target_function

    def value(self,samples):
        """
        Evaluate the function at a set of samples.

        The number of QoI of the vectored_valued function is determined
        by probing the target_function with the first sample in the set
        of samples.

        Parameters
        ----------
        samples : (num_vars x num_samples)
            The coordinates of the samples.
            If samples is a 1D array it will be assumed that
            samples is only one sample and be converted to
            a matrix (num_vars x 1)

        Returns
        -------
        values : (num_samples x num_qoi) matrix
            The vector-valued function value at the samples
        """
	if samples.ndim==1:
            samples = samples.reshape(samples.shape[0],1)
        num_samples = samples.shape[1]
        values = numpy.empty((num_samples),float)
        values_0 = self.target_function(samples[:,0])
        if numpy.isscalar(values_0):
            values_0 = numpy.array([values_0])
        assert values_0.ndim==1
        num_qoi = values_0.shape[0]
        values = numpy.empty((num_samples,num_qoi),float)
        values[0,:] = values_0
        for i in xrange(1,samples.shape[1]):
            values_i = self.target_function(samples[:,i])
            if numpy.isscalar(values_i):
                values_i = numpy.array([values_i])
            values[i,:] = values_i
        return values
 %}
