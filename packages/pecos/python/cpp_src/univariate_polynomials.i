/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/* univariate_polynomials.i */

%define %univariate_polynomials_docstring
"
PyDakota.univariate_polynomials is the python interface to the Dakota
univariate polynomial classes.
"
%enddef

%module(package="PyDakota",
        directors=1,
        autodoc=1,
        implicitconv = "1",
        doc_string=%univariate_polynomials_docstring) univariate_polynomials

%{
#include "numpy_include.hpp"//gets rid of deprecated warnings

  //#include "pecos_data_types.hpp"
#include <vector>
#include "BasisPolynomial.hpp"

  // needed because of confusion induced by surrogates and Pecos namespaces
  // when defining Real
  using namespace Pecos;
%}

%include "numpy.i"
%fragment("NumPy_Fragments");
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

// importing math_tools.i avoids need to %include
// Teuchos_SerialDenseVector?Matrix.i and data_Structures.i
%import "math_tools.i"


// custom std vector typemap needed for returning values of 
// Pecos::ApproximationBasis
%define %std_vector_typemaps(SCALAR_TYPE, TYPECODE)
%typemap(out) std::vector<SCALAR_TYPE> const &
{
  npy_intp dims[1] = {static_cast<npy_intp>($1->size())};
  $result = PyArray_SimpleNew( 1, dims, TYPECODE );
  if (!$result) SWIG_fail;
  SCALAR_TYPE *array = (SCALAR_TYPE *)PyArray_DATA((PyArrayObject*)$result );
  for (std::vector<SCALAR_TYPE>::iterator it=$1->begin() ; it!=$1->end(); ++it){
    *array++ = *it;
  }
}
%enddef
%std_vector_typemaps( double            , NPY_DOUBLE   )

%include "pecos_global_defs.hpp"
%include "BasisPolynomial.hpp"

typedef double Real;
typedef std::vector<double>                RealArray;

%pythoncode %{
    import numpy as np
%}

%extend Pecos::BasisPolynomial{  
%pythoncode
%{
    def values(self,samples,degree):
        """
        Evaluate the polynomial basis at a set of samples for all degrees
        d=0,....,degree.

        TODO: currently this wraps self.type1_value. In future create
        a new C++ function
        RealMatrix& values(const RealVector& samples, degree)
        and use swig to wrap this. rename values in .i file and use
        this function to deal with cases when samples is a scalar and a
        np.ndarray

        Parameters
        ----------
        samples : np.ndarray (nsamples) or double
            The samples at which to evaluate the polynomial.

        degree : integer
            The maximum degree.

        Returns
        -------
        values : np.ndarray (nsamples x nterms)
            The basis value, for each degree, at each sample.
        """
        if np.isscalar(samples):
            samples = np.asarray([samples])
        assert samples.ndim==1
        nsamples = samples.shape[0]; nterms=degree+1
        values = np.empty((nsamples, nterms),dtype=float)
        for i in range(nsamples):
            for j in range(nterms):
                values[i,j] = self.type1_value(samples[i], j)
        return values
%}
}
