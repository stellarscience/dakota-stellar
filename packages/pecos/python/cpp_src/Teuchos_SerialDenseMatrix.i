/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

// Teuchos_SerialDenseMatrix.i is a SWIG interface file that provides SWIG
// directives to handle Teuchos SerialDenseMatrix types.  These class is not
// wrapped, but instead typemaps are defined so that the python user
// can use NumPy arrays or Python sequences instead.
%{
#include "Teuchos_SerialDenseMatrix.hpp"
using Teuchos::SerialDenseMatrix;
%}
%import  "Teuchos_SerialDenseMatrix.hpp"

////////////////////////////////////////////////////////////////////////
// The philosophy is that wherever Teuchos SerialDenseMatrix classes are used in
// C++, NumPy arrays will be used in python.  Thus we need the NumPy
// SWIG directives.
%include "numpy.i"

////////////////////////////////////////////////////////////////////////
// Define a macro that takes a C++ data ordinal and scalar type
// (ORDINAL_TYPE,SCALAR_TYPE) and a
// corresponding NumPy typecode (TYPECODE) and define all of the
// typemaps needed to handle that (ORDINAL_TYPE,SCALAR_TYPE) array.
%define %teuchos_sdm_typemaps(ORDINAL_TYPE, SCALAR_TYPE, TYPECODE)

// If an SerialDenseMatrix argument has a template parameter argument that is
// a const TYPE, then we know that the argument is input only.
// Therefore we allow any type of sequence to be converted to a
// PyArrayObject and then extract the resulting data pointer to
// construct the SerialDenseMatrix.  If the conversion creates a new
// PyArrayObject, then we have to be sure to decrement its reference
// count once the SerialDenseMatrix has been used.

//////////////////////////////////////
// Teuchos::SerialDenseMatrix<const ORDINAL_TYPE,const SCALAR_TYPE> //
//////////////////////////////////////
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (Teuchos::SerialDenseMatrix<const ORDINAL_TYPE,const SCALAR_TYPE>)
{
  $1 = is_array($input) || PySequence_Check($input);
}

%typemap(in) Teuchos::SerialDenseMatrix<const ORDINAL_TYPE,const SCALAR_TYPE>
(int is_new = 0,
 PyArrayObject * npArray = NULL)
{
  // ENFORCE that $input is contiguous and has fortran (column-major ordering)
  // If it is not the matrix is copied into fortran format, else a view is taken
  // and the reference count is incremented
  npArray = obj_to_array_fortran_allow_conversion($input, TYPECODE, &is_new);
  if ( npArray == NULL ) SWIG_fail;

  // Now we need to check that the NumPy array that we have is 2D.
  if ( PyArray_NDIM( npArray ) != 2 ) {
    PyErr_SetString(PyExc_ValueError, "Array data must be two dimensional");
    SWIG_fail;
  }

  // Get array shape
  ORDINAL_TYPE m = (ORDINAL_TYPE)PyArray_DIM(npArray,0),
    n = (ORDINAL_TYPE)PyArray_DIM(npArray,1),
    stride = (ORDINAL_TYPE)( PyArray_STRIDE(npArray,1) /
				 PyArray_ITEMSIZE(npArray) );

  $1 = Teuchos::SerialDenseMatrix<const ORDINAL_TYPE,const SCALAR_TYPE>(
		Teuchos::View,
		(SCALAR_TYPE*)PyArray_DATA(npArray),
		 stride, m, n);
}

%typemap(freearg) Teuchos::SerialDenseMatrix<const ORDINAL_TYPE,const SCALAR_TYPE>
{
  if (is_new$argnum) Py_DECREF(npArray$argnum);
}

%typemap(out) Teuchos::SerialDenseMatrix<const ORDINAL_TYPE,const SCALAR_TYPE>
{
  ORDINAL_TYPE m = $1->numRows(), n = $1->numCols();
  npy_intp dims[2] = { m, n };
  $result = PyArray_EMPTY( 2, dims, TYPECODE,  NPY_FORTRANORDER );
  if (!$result) SWIG_fail;
  SCALAR_TYPE *buffer = (SCALAR_TYPE *)PyArray_DATA( $result );
  for ( int j = 0; j < n; j++ ){
    for ( int i = 0; i < m; i++ ){
      buffer[j*m+i] = $1->operator()(i,j);
    }
  }
}

//////////////////////////////////////////////
// Teuchos::SerialDenseMatrix<const ORDINAL_TYPE,const SCALAR_TYPE> const & //
//////////////////////////////////////////////
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (Teuchos::SerialDenseMatrix<const ORDINAL_TYPE,const SCALAR_TYPE> const &)
{
  $1 = is_array($input) || PySequence_Check($input);
}

%typemap(in) Teuchos::SerialDenseMatrix<const ORDINAL_TYPE,const SCALAR_TYPE> const &
(int is_new = 0,
 PyArrayObject * npArray = NULL,
 Teuchos::SerialDenseMatrix<const ORDINAL_TYPE,const SCALAR_TYPE> temp)
{
  // ENFORCE that $input is contiguous and has fortran (column-major ordering)
  // If it is not the matrix is copied into fortran format, else a view is taken
  // and the reference count is incremented
  npArray = obj_to_array_fortran_allow_conversion($input, TYPECODE, &is_new);
  if ( npArray == NULL ) SWIG_fail;

  // Now we need to check that the NumPy array that we have is 2D.
  if ( PyArray_NDIM( npArray ) != 2 ) {
    PyErr_SetString(PyExc_ValueError, "Array data must be two dimensional");
    SWIG_fail;
  }

  // Get array shape
  ORDINAL_TYPE m = (ORDINAL_TYPE)PyArray_DIM(npArray,0),
    n = (ORDINAL_TYPE)PyArray_DIM(npArray,1),
    stride = (ORDINAL_TYPE)( PyArray_STRIDE(npArray,1) /
				 PyArray_ITEMSIZE(npArray) );

  temp = Teuchos::SerialDenseMatrix<const ORDINAL_TYPE, const SCALAR_TYPE>(
		  Teuchos::View,
		  (SCALAR_TYPE*)PyArray_DATA(npArray),
		  stride, m, n);

  // The use of & is important
  $1 = &temp;
}

%typemap(freearg) Teuchos::SerialDenseMatrix<const ORDINAL_TYPE,const SCALAR_TYPE> const &
{
  if (is_new$argnum) Py_DECREF(npArray$argnum);
}

// If an SerialDenseMatrix argument has template parameter argument that is a
// non-const TYPE, then the default behavior is to assume that the
// array is input/output.  Therefore the input python argument must be
// a NumPy array.

////////////////////////////////
// Teuchos::SerialDenseMatrix<ORDINAL_TYPE, SCALAR_TYPE> //
////////////////////////////////
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (Teuchos::SerialDenseMatrix<ORDINAL_TYPE, SCALAR_TYPE>)
{
  $1 = is_array($input);
}

%typemap(in) Teuchos::SerialDenseMatrix<ORDINAL_TYPE, SCALAR_TYPE>
{
  // ENFORCE that $input is contiguous and has fortran (column-major ordering)
  // If it is not the matrix is copied into fortran format, else a view is taken
  // and the reference count is incremented
  PyObject *npArray = PyArray_FROM_OTF($input, TYPECODE, NPY_ARRAY_F_CONTIGUOUS);
  if ( npArray == NULL ) SWIG_fail;

  // Now we need to check that the NumPy array that we have is 2D.
  if ( PyArray_NDIM( (PyArrayObject*)npArray ) != 2 ) {
    PyErr_SetString(PyExc_ValueError, "Array data must be two dimensional");
    SWIG_fail;
  }

  // Get array shape
  ORDINAL_TYPE m = (ORDINAL_TYPE)PyArray_DIM((PyArrayObject*)npArray,0),
    n = (ORDINAL_TYPE)PyArray_DIM((PyArrayObject*)npArray,1),
    stride = (ORDINAL_TYPE)( PyArray_STRIDE((PyArrayObject*)npArray,1) /
				 PyArray_ITEMSIZE((PyArrayObject*)npArray) );

  $1 = Teuchos::SerialDenseMatrix<ORDINAL_TYPE,SCALAR_TYPE>(
		Teuchos::View,
		(SCALAR_TYPE*)PyArray_DATA((PyArrayObject*)npArray),
		stride, m, n);
}

%typemap(out) Teuchos::SerialDenseMatrix<ORDINAL_TYPE, SCALAR_TYPE>
{
  ORDINAL_TYPE m = $1->numRows(), n = $1->numCols();
  npy_intp dims[2] = { m, n };
  $result = PyArray_EMPTY( 2, dims, TYPECODE,  NPY_FORTRANORDER );
  if (!$result) SWIG_fail;
  SCALAR_TYPE *buffer = (SCALAR_TYPE *)PyArray_DATA( $result );
  for ( int j = 0; j < n; j++ ){
    for ( int i = 0; i < m; i++ ){
      buffer[j*m+i] = $1->operator()(i,j);
    }
  }
}

////////////////////////////////////////
// Teuchos::SerialDenseMatrix<ORDINAL_TYPE, SCALAR_TYPE> const & //
////////////////////////////////////////
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (Teuchos::SerialDenseMatrix<ORDINAL_TYPE, SCALAR_TYPE> const &)
{
  $1 = is_array($input);
}

%typemap(in) Teuchos::SerialDenseMatrix<ORDINAL_TYPE, SCALAR_TYPE> const &
(Teuchos::SerialDenseMatrix<ORDINAL_TYPE, SCALAR_TYPE> temp)
{
  // ENFORCE that $input is contiguous and has fortran (column-major ordering)
  // If it is not the matrix is copied into fortran format, else a view is taken
  // and the reference count is incremented
  PyObject *npArray = PyArray_FROM_OTF($input, TYPECODE, NPY_ARRAY_F_CONTIGUOUS);
  if ( npArray == NULL ) SWIG_fail;

  // Now we need to check that the NumPy array that we have is 2D.
  if ( PyArray_NDIM( (PyArrayObject*)npArray ) != 2 ) {
    PyErr_SetString(PyExc_ValueError, "Array data must be two dimensional");
    SWIG_fail;
  }

  // Get array shape
  ORDINAL_TYPE m = (ORDINAL_TYPE)PyArray_DIM((PyArrayObject*)npArray,0),
    n = (ORDINAL_TYPE)PyArray_DIM((PyArrayObject*)npArray,1),
    stride = (ORDINAL_TYPE)( PyArray_STRIDE((PyArrayObject*)npArray,1) /
				 PyArray_ITEMSIZE((PyArrayObject*)npArray) );

  temp = Teuchos::SerialDenseMatrix<ORDINAL_TYPE,SCALAR_TYPE>(
		Teuchos::View,
		(SCALAR_TYPE*)PyArray_DATA((PyArrayObject*)npArray),
		stride, m, n);

  // The use of & is important
  $1 = &temp;
}

%typemap(out) Teuchos::SerialDenseMatrix<ORDINAL_TYPE, SCALAR_TYPE> const &
{
  ORDINAL_TYPE m = $1->numRows(), n = $1->numCols();
  npy_intp dims[2] = { m, n };
  $result = PyArray_EMPTY( 2, dims, TYPECODE,  NPY_FORTRANORDER );
  if (!$result) SWIG_fail;
  SCALAR_TYPE *buffer = (SCALAR_TYPE *)PyArray_DATA( $result );
  for ( int j = 0; j < n; j++ ){
    for ( int i = 0; i < m; i++ ){
      buffer[j*m+i] = $1->operator()(i,j);
    }
  }
}

// If an SerialDenseMatrix is an output argument

////////////////////////////////
// Teuchos::SerialDenseMatrix<ORDINAL_TYPE, SCALAR_TYPE> & result//
////////////////////////////////

// Specify how to return result to python
%typemap(argout) Teuchos::SerialDenseMatrix<ORDINAL_TYPE,SCALAR_TYPE> &argout
{
  ORDINAL_TYPE m = $1->numRows(), n = $1->numCols();
  npy_intp dims[2] = { m, n };
  PyObject *npArray = PyArray_EMPTY( 2, dims, TYPECODE, 1 );
  if (!npArray) SWIG_fail;
  SCALAR_TYPE *buffer = (SCALAR_TYPE *)PyArray_DATA( (PyArrayObject*)npArray );
  for ( int j = 0; j < n; j++ ){
    for ( int i = 0; i < m; i++ ){
      buffer[j*m+i] = $1->operator()(i,j);
    }
  }
  $result = SWIG_Python_AppendOutput($result,npArray);
}

// Remove result from the python function call.
%typemap(in,numinputs=0) Teuchos::SerialDenseMatrix<ORDINAL_TYPE,SCALAR_TYPE> &argout
  ( Teuchos::SerialDenseMatrix<ORDINAL_TYPE,SCALAR_TYPE> TSDM )
{
  // Must declare the variable that needs to be passed to the c++ function
  TSDM = Teuchos::SerialDenseMatrix<ORDINAL_TYPE,SCALAR_TYPE>();
  $1 = &TSDM;
};

%typemap(freearg) Teuchos::SerialDenseMatrix<ORDINAL_TYPE,SCALAR_TYPE> &argout
{
  // Do nothing. Needed so that default free arg is not used
  // Check with Bill Spotz
};

//Allow serial dense matrix be passed back to a python class derived from a c++ director class
%typemap(directorin) Teuchos::SerialDenseMatrix<int,double> &%{
  npy_intp dims$argnum[2] = { $1_name.numRows(),$1_name.numCols() };
  double * data$argnum = $1_name.values();
  PyArray_Descr * dtype$argnum = PyArray_DescrFromType( NPY_DOUBLE );
  PyArrayObject* py_array$argnum =
    (PyArrayObject*) PyArray_NewFromDescr( &PyArray_Type, dtype$argnum, 2,
					   dims$argnum, NULL, (void*)data$argnum,
					   NPY_ARRAY_F_CONTIGUOUS, NULL );
  // Cast to (PyObject*) is necessary
  $input = (PyObject*)py_array$argnum;
  %}

%typemap(directorargout) Teuchos::SerialDenseMatrix<int,double> &argout%{

  // ENFORCE that $result is contiguous and has fortran (column-major ordering)
  // If it is not the matrix is copied into fortran format, else a view is taken
  // and the reference count is incremented
  PyObject *npArray$argnum=PyArray_FROM_OTF($result,TYPECODE,NPY_ARRAY_F_CONTIGUOUS);
  if ( npArray$argnum == NULL )
    throw( std::runtime_error( "Conversion to fortran array failed" ) );

  // Now we need to check that the NumPy array that we have is 2D.
  if ( PyArray_NDIM( (PyArrayObject*)npArray$argnum ) != 2 ) {
    throw( std::runtime_error( " Values out must be 2 dimensional" ) );
    //PyErr_SetString(PyExc_ValueError, "Array data must be two dimensional");
    //SWIG_fail;
  }

  // Get array shape
  ORDINAL_TYPE m$argnum = (ORDINAL_TYPE)PyArray_DIM((PyArrayObject*)npArray$argnum,0),
    n$argnum = (ORDINAL_TYPE)PyArray_DIM((PyArrayObject*)npArray$argnum,1),
    stride$argnum = (ORDINAL_TYPE)( PyArray_STRIDE((PyArrayObject*)npArray$argnum,1) /
				    PyArray_ITEMSIZE((PyArrayObject*)npArray$argnum) );

  Teuchos::SerialDenseMatrix<ORDINAL_TYPE,SCALAR_TYPE> temp$argnum(
		Teuchos::View,
		(SCALAR_TYPE*)PyArray_DATA((PyArrayObject*)npArray$argnum),
		stride$argnum, m$argnum, n$argnum);

  $1.shapeUninitialized(m$argnum,n$argnum);
  $1.assign(temp$argnum);
%}

%enddef

////////////////////////////////////////////////////////////////////////
// Call the %teuchos_sdm_typemaps() macro for specific data types
// that are supported by NumPy
%teuchos_sdm_typemaps(int, signed char       , NPY_BYTE     )
%teuchos_sdm_typemaps(int, unsigned char     , NPY_UBYTE    )
%teuchos_sdm_typemaps(int, short             , NPY_SHORT    )
%teuchos_sdm_typemaps(int, unsigned short    , NPY_USHORT   )
%teuchos_sdm_typemaps(int, int               , NPY_INT      )
%teuchos_sdm_typemaps(int, unsigned int      , NPY_UINT     )
%teuchos_sdm_typemaps(int, long              , NPY_LONG     )
%teuchos_sdm_typemaps(int, unsigned long     , NPY_ULONG    )
%teuchos_sdm_typemaps(int, long long         , NPY_LONGLONG )
%teuchos_sdm_typemaps(int, unsigned long long, NPY_ULONGLONG)
%teuchos_sdm_typemaps(int, float             , NPY_FLOAT    )
%teuchos_sdm_typemaps(int, double            , NPY_DOUBLE   )


// ---------------------------------------------------------------------//
// have typemaps be applied to classes named result, result_0, result_1 //
// ---------------------------------------------------------------------//
%apply Teuchos::SerialDenseMatrix<int,int> &argout { Teuchos::SerialDenseMatrix<int,int> &result }
%apply Teuchos::SerialDenseMatrix<int,int> &argout { Teuchos::SerialDenseMatrix<int,int> &result_0 }
%apply Teuchos::SerialDenseMatrix<int,int> &argout { Teuchos::SerialDenseMatrix<int,int> &result_1 }
%apply Teuchos::SerialDenseMatrix<int,double> &argout { Teuchos::SerialDenseMatrix<int,double> &result }
%apply Teuchos::SerialDenseMatrix<int,double> &argout { Teuchos::SerialDenseMatrix<int,double> &result_0 }
%apply Teuchos::SerialDenseMatrix<int,double> &argout { Teuchos::SerialDenseMatrix<int,double> &result_1 }
