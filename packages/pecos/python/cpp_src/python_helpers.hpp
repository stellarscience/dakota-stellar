/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PYTHON_HELPERS_HPP
#define PYTHON_HELPERS_HPP

// NumPy include
// Python.h must be included before any headers to avoid 
// warning: "_POSIX_C_SOURCE" redefined. 
// Python.h is included in numpy_include.hpp
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
#include "Teuchos_SerialDenseVector.hpp"

#include "OptionsList.hpp"
#include <memory>

namespace Pecos {

/** \brief Copy a numpy ndarray into a Teuchos SerialDenseVector (SDV).
 * T is the scalar type of the SDV and S is the scalar type of the numpy array
 * Having these two template parameters instead of just one T allows the
 # conversion of both NPY_INT and NPY_LONG into A SDV<int,int>.
 */
template< typename T, typename S >
void copyNumPyToTeuchosVector(PyObject * pyArray,
			      Teuchos::SerialDenseVector<int,T > & tvec);

template< typename T >
PyObject * copyTeuchosVectorToNumPy(Teuchos::SerialDenseVector< int,T > &tvec);


template< typename T, typename S >
void copyNumPyToTeuchosMatrix(PyObject * pyArray,
			      Teuchos::SerialDenseMatrix<int,T > & tmat);

template< typename T >
PyObject * copyTeuchosMatrixToNumPy(Teuchos::SerialDenseMatrix< int,T > &tmat);

bool setPythonParameter(util::OptionsList & opts_list,
			const std::string      & name,
			PyObject               * value);

bool updateOptionsListWithPyDict(PyObject    * dict,
				   util::OptionsList & opts_list);

util::OptionsList * pyDictToNewOptionsList(PyObject* dict);

PyObject * getPythonParameter(const util::OptionsList & plist,
			      const std::string & name);


bool updatePyDictWithOptionsList(PyObject          * dict,
				 const util::OptionsList & opts_list);

PyObject * optionsListToNewPyDict(const util::OptionsList & opts_list);

template< typename TYPE >
int NumPy_TypeCode();

void pyListToNewStdVector(PyObject *pylist, std::vector< util::OptionsList > &list);

PyObject * copyStdVectorToPyList(std::vector<util::OptionsList> &list);

bool PyListOfOptionsList_check(PyObject* pylist);

  // Functions to help debugging
void print_pyobject_string_rep(PyObject * value);
  
void print_PyType(PyObject *value);

} //namespace Pecos
#endif // PYTHON_HELPERS_HPP
