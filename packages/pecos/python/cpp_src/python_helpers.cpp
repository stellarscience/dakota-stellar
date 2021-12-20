/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "python_helpers.hpp"
#include "swigpyrun.h"

using Pecos::util::OptionsList;

namespace Pecos {

void print_pyobject_string_rep(PyObject * value){
  PyObject* objectsRepresentation = PyObject_Repr(value);
  const char* s = PyString_AsString(objectsRepresentation);
  std::printf(s);
  std::printf("\n");
}

  
void print_PyType(PyObject *value){
  std::printf("The object: ");
  print_pyobject_string_rep(value);
  if (PyBool_Check(value))
    std::printf("Is int\n");
  else if (PyInt_Check(value))
    std::printf("Is int\n");
  else if (PyFloat_Check(value))
    std::printf("Is float\n");
  else if (PyString_Check(value))
    std::printf("Is string\n");
  else if (value == Py_None)
    std::printf("Is None\n");
  else if (PyDict_Check(value))
    std::printf("Is dict\n");
  else if (PyListOfOptionsList_check(value))
    std::printf("Is list of options\n");
  else if (PyList_Check(value))
    std::printf("Is list\n");
  else if (PyArray_Check(value))
    std::printf("Is array\n");
  else if (PyTuple_Check(value))
    std::printf("Is tuple\n");
  else if (PySequence_Check(value))
    std::printf("Is sequence\n");
  else
    std::printf("value unknown\n");
}

bool PyListOfOptionsList_check(PyObject* pylist){
  void * argp;
  static swig_type_info * swig_TPL_ptr =
    SWIG_TypeQuery("std::shared_ptr< OptionsList >*");
  Py_ssize_t len = PyList_Size(pylist);
  for(Py_ssize_t i=0; i<len; ++i){
    PyObject *value = PyList_GetItem(pylist, i);
    if (!SWIG_CheckState(SWIG_Python_ConvertPtr(value,
						&argp,
						swig_TPL_ptr,
						0)) &&
	(!PyDict_Check(value)))
      return false;
  }
  return true;
}

template<>
int NumPy_TypeCode<int>(){
  return NPY_INT;
}

template<>
int NumPy_TypeCode<long>(){
  return NPY_LONG;
}


template<>
int NumPy_TypeCode<double>(){
  return NPY_DOUBLE;
}

template< typename T, typename S >
void copyNumPyToTeuchosVector(PyObject * pyArray,
			      Teuchos::SerialDenseVector<int,T > & tvec){
  int length = PyArray_DIM((PyArrayObject*) pyArray, 0);
  tvec.sizeUninitialized(length);
  S * data = (S*) PyArray_DATA((PyArrayObject*) pyArray);
  for (int i=0; i<length; ++i){
    // static cast need because we are casting npy_long and npy_int to int
    tvec[i] = static_cast<T>(*(data++));
  }
}

template< typename T, typename S >
void copyNumPyToTeuchosMatrix(PyObject * pyArray,
			      Teuchos::SerialDenseMatrix<int,T > & tmat){
  int num_rows = PyArray_DIM((PyArrayObject*) pyArray, 0),
    num_cols = PyArray_DIM((PyArrayObject*) pyArray, 1);
  tmat.shapeUninitialized(num_rows,num_cols);
  S * data = (S*) PyArray_DATA((PyArrayObject*) pyArray);
  for (int j=0; j<num_cols; ++j){
    for (int i=0; i<num_rows; ++i){
      // static cast need because we are casting npy_long and npy_int to int
      // unlike all other matrix based functions here we 
      // assume c ordering of numpy array
      tmat(i,j) = static_cast<T>(data[i*num_cols+j]);
    }
  }
}

template< typename T >
PyObject * copyTeuchosVectorToNumPy(Teuchos::SerialDenseVector< int,T > &tvec)
{
  int typecode = NumPy_TypeCode< T >();
  npy_intp dims[] = { tvec.length() };
  PyObject * pyArray = PyArray_SimpleNew(1, dims, typecode);
  T * data = (T*) PyArray_DATA((PyArrayObject*) pyArray);
  for (int i=0; i<tvec.length(); ++i)
    *(data++) = tvec[i];
  return pyArray;
}

template< typename T >
PyObject * copyTeuchosMatrixToNumPy(Teuchos::SerialDenseMatrix< int,T > &tmat)
{
  int m = tmat.numRows(), n = tmat.numCols();
  int typecode = NumPy_TypeCode< T >();
  npy_intp dims[] = { m,n };
  PyObject * pyArray = PyArray_EMPTY( 2, dims, typecode,  NPY_FORTRANORDER);
  T * data = (T*) PyArray_DATA((PyArrayObject*) pyArray);
  for ( int j = 0; j < n; j++ ){
    for ( int i = 0; i < m; i++ ){
      data[j*m+i] = tmat(i,j);
    }
  }
  return pyArray;
}

bool setPythonParameter(OptionsList & opts_list,
			const std::string      & name,
			PyObject               * value)
{ 
  static swig_type_info * swig_TPL_ptr =
    SWIG_TypeQuery("std::shared_ptr< OptionsList >*");
  void * argp;
  int newmem = 0;

  // Boolean values
  if (PyBool_Check(value))
  {
    if (value == Py_True) opts_list.set(name,true );
    else                  opts_list.set(name,false);
  }
  
  // Integer values
  else if (PyInt_Check(value))
  {
    opts_list.set(name, (int)PyInt_AsLong(value));
  }

  // Floating point values
  else if (PyFloat_Check(value))
  {
    opts_list.set(name, PyFloat_AsDouble(value));
  }

  // String values
  else if (PyString_Check(value))
  {
    opts_list.set(name, std::string(PyString_AsString(value)));
  }

  // None object not allowed: this is a python type not usable by
  // Trilinos solver packages, so we reserve it for the
  // getPythonParameter() function to indicate that the requested
  // parameter does not exist in the given OptionsList.
  // For logic reasons, this check must come before the check for
  // OptionsList
  else if (value == Py_None)
  {
    return false;
  }

  // Dictionary values.  This must come before the check for Python
  // sequences, because the Python ditionary is a sequence.
  else if (PyDict_Check(value))
  {  
    // Convert the python dictionary to a OptionsList
    OptionsList * sublist = pyDictToNewOptionsList(value);

    // Store the OptionsList
    opts_list.set(name,*sublist);
    delete sublist;
  }

  // OptionsList values
  else if (SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(value,
                                                        &argp,
                                                        swig_TPL_ptr,
                                                        0,
                                                        &newmem)))
  {
    if (newmem & SWIG_CAST_NEW_MEMORY)
    { 
      std::shared_ptr< OptionsList > tempshared =
	*reinterpret_cast< std::shared_ptr< OptionsList > * >(argp);
      delete reinterpret_cast< std::shared_ptr< OptionsList > * >(argp);
      opts_list.set(name, *(tempshared.get()));
    }
    else
    {
      std::shared_ptr< OptionsList > * smartarg =
	reinterpret_cast< std::shared_ptr< OptionsList > * >(argp);
      if (smartarg) opts_list.set(name, *(smartarg->get()));
    }
  }
  
  else if (PyList_Check(value) && PyListOfOptionsList_check(value)) {
    // only allow lists of optionslist here
    // list of numerical entries are handled below and converted to
    // Teuchos matrices and vectors
    std::vector<OptionsList> list;
    pyListToNewStdVector(value, list);
    opts_list.set(name, list);
   }
  
  else if (PyArray_Check(value) || PySequence_Check(value)){
    PyObject * pyArray =
      PyArray_CheckFromAny(value,
                           NULL,
                           1,
                           2,
                           NPY_ARRAY_DEFAULT | NPY_ARRAY_NOTSWAPPED,
                           NULL);
    if (!pyArray) return false;
    if ( PyArray_NDIM( (PyArrayObject*)pyArray)==1 ) {
      if (PyArray_TYPE((PyArrayObject*) pyArray) == NumPy_TypeCode<int>()){
        Teuchos::SerialDenseVector<int,int> tvec;
        copyNumPyToTeuchosVector<int,int>(pyArray, tvec);
        opts_list.set(name, tvec);
      }
      else if (PyArray_TYPE((PyArrayObject*) pyArray) == NumPy_TypeCode<long>()){
        Teuchos::SerialDenseVector<int,int> tvec;
        copyNumPyToTeuchosVector<int,long>(pyArray, tvec);
        opts_list.set(name, tvec);
      }
      else if (PyArray_TYPE((PyArrayObject*) pyArray) == NumPy_TypeCode<double>()){
        Teuchos::SerialDenseVector< int,double > tvec;
        copyNumPyToTeuchosVector<double,double>(pyArray, tvec);
        opts_list.set(name, tvec);
      }else{
        // Unsupported data type
        if (pyArray != value) Py_DECREF(pyArray);
        return false;
      }
    }
    if ( PyArray_NDIM( (PyArrayObject*)pyArray)==2 ) {
      if (PyArray_TYPE((PyArrayObject*) pyArray) == NumPy_TypeCode<int>()){
        Teuchos::SerialDenseMatrix<int,int> tmat;
        copyNumPyToTeuchosMatrix<int,int>(pyArray, tmat);
        opts_list.set(name, tmat);
      }
      else if (PyArray_TYPE((PyArrayObject*) pyArray) == NumPy_TypeCode<long>()){
        Teuchos::SerialDenseMatrix<int,int> tmat;
        copyNumPyToTeuchosMatrix<int,long>(pyArray, tmat);
        opts_list.set(name, tmat);
      }
      else if (PyArray_TYPE((PyArrayObject*) pyArray) == NumPy_TypeCode<double>()){
        Teuchos::SerialDenseMatrix< int,double > tmat;
        copyNumPyToTeuchosMatrix<double,double>(pyArray, tmat);
        opts_list.set(name, tmat);
      }else{
        // Unsupported data type
        if (pyArray != value) Py_DECREF(pyArray);
        return false;
      }
    }
  }
  // All other value types are unsupported
  else
  {
    return false;
  }

  // Successful type conversion
  return true;
}    // setPythonParameter


bool updateOptionsListWithPyDict(PyObject    * dict,
				   OptionsList & opts_list)
{
  PyObject  * key    = NULL;
  PyObject  * value  = NULL;
  Py_ssize_t  pos    = 0;
  bool        result = true;
  std::string name;

  // The dict pointer must point to a dictionary
  if (!PyDict_Check(dict))
  {
    PyErr_SetString(PyExc_TypeError, "Expected a dictionary");
    goto fail;
  }

  // Iterate over all items in the python dictionary and ensure they
  // are synchronized with the OptionsList
  while (PyDict_Next(dict, &pos, &key, &value))
  {

    // If the key is not a string, we can't synchronize
    if (!PyString_Check(key))
    {
      PyErr_SetString(PyExc_TypeError, "Encountered non-string key in dictionary");
      goto fail;
    }

    name = std::string(PyString_AsString(key));
    if (!setPythonParameter(opts_list, name, value))
    {
      PyErr_SetString(PyExc_TypeError, "value could not be set");
      goto fail;
    }
  }
  return result;
  fail:
  return false;
}    // updateOptionsListWithPyDict


OptionsList *
pyDictToNewOptionsList(PyObject* dict)
{
  OptionsList * opts_list = 0;
  // The dict pointer must point to a dictionary
  if (!PyDict_Check(dict))
  {
    PyErr_SetString(PyExc_ValueError, "Expected a dictionary");
    goto fail;
  }

  // Create a new OptionsList and synchronize it with the python
  // dictionary
  opts_list = new OptionsList();
  if (!updateOptionsListWithPyDict(dict,*opts_list))
  {
      delete opts_list;
      goto fail;
  }
  return opts_list;
  fail:
  return NULL;
}    // pyDictToNewOptionsList


// **************************************************************** //

bool updatePyDictWithOptionsList(PyObject                     * dict,
				   const OptionsList & opts_list)
{
  PyObject   * value   = NULL;
  PyObject   * param   = NULL;
  bool         result  = true;
  const char * nameStr = NULL;
  std::string  name;
  std::vector<std::string> keys;

  // The dict pointer must point to a dictionary
  if (!PyDict_Check(dict))
  {
    PyErr_SetString(PyExc_TypeError, "Expected a dictionary");
    goto fail;
  }

  // Iterate over all entries in OptionsList and ensure they are
  // mirrored in the python dictionary
  opts_list.get_keys(keys);
  for (size_t i=0; i<keys.size(); ++i)
  {
    name    = keys[i];
    nameStr = name.c_str();
    param   = getPythonParameter(opts_list,nameStr);
    value   = PyDict_GetItemString(dict,nameStr);

    if (param == NULL)
    {
      PyErr_Format(PyExc_RuntimeError, "Unexpected enumeration encountered");
      goto fail;
    }
    else
      {
      // If param is a sublist, mirror with a dictionary by calling
      // this routine recursively
	if (opts_list.get_entry(name)->type()==typeid(OptionsList))
      {
	if (value == NULL) value = PyDict_New();
	else if (!PyDict_Check(value))
	{
	  Py_DECREF(value);
	  value = PyDict_New();
	}
	result = result and 
	  updatePyDictWithOptionsList(value,opts_list.get<OptionsList>(name));
	PyDict_SetItemString(dict,nameStr,value);
      }

      // Else synchronize the dictionary value to the parameter
      else PyDict_SetItemString(dict,nameStr,param);
    }
    Py_XDECREF(param);
  }
  return result;
  fail:
  return false;
}    // updatePyDictWithOptionsList

// **************************************************************** //

PyObject * optionsListToNewPyDict(const OptionsList & opts_list)
{
  
  // Create a new dictionary and synchronize it with the OptionsList
  PyObject * dict = PyDict_New();
  if (!updatePyDictWithOptionsList(dict,opts_list))
    {
      Py_XDECREF(dict);
      goto fail;
    }
  return dict;
 fail:
  return NULL;
}    // optionsListToNewPyDict



PyObject * getPythonParameter(const OptionsList & opts_list,
			      const std::string            & name)
{
  // used if returning optionslist wrapper and not python dict
  // which is the case now
  //static swig_type_info * swig_TPL_ptr =
  //      SWIG_TypeQuery("std::shared_ptr< OptionsList >*");
  
  // If parameter does not exist, return None
  if (!opts_list.exists(name)) return Py_BuildValue("");

  // Get the parameter entry.  I now deal with the Teuchos::ParameterEntry
  // objects so that I can query the OptionsList without setting
  // the "used" flag to true.
  const boost::any* entry = opts_list.get_entry(name);

  // Boolean parameter values
  if (entry->type()==typeid(bool))
  {
    bool value = boost::any_cast<bool>(*entry);
    return PyBool_FromLong((long)value);
  }
  // Integer parameter values
  else if (entry->type()==typeid(int))
  {
    int value = boost::any_cast< int >(*entry);
    return PyInt_FromLong((long)value);
  }
  // Double parameter values
  else if (entry->type()==typeid(double))
  {
    double value = boost::any_cast< double >(*entry);
    return PyFloat_FromDouble(value);
  }
  // String parameter values
  else if (entry->type()==typeid(std::string))
  {
    std::string value = boost::any_cast< std::string >(*entry);
    return PyString_FromString(value.c_str());
  }
  // Char * parameter values
  else if (entry->type()==typeid(char *))
  {
    char * value = boost::any_cast< char * >(*entry);
    return PyString_FromString(value);
  }
  // OptionsList values
  else if (entry->type()==typeid(OptionsList))
  {
    // OptionsList value = boost::any_cast< OptionsList >(*entry);
    // std::shared_ptr< OptionsList > * valuercp =
    //   new std::shared_ptr< OptionsList >(&value);
    // return SWIG_NewPointerObj((void*)valuercp, swig_TPL_ptr, SWIG_POINTER_NOSHADOW);
    // // Bill uses this but I get a segfault unless I used above. Why?
    // //return SWIG_NewPointerObj((void*)valuercp, swig_TPL_ptr, SWIG_POINTER_OWN);
    return optionsListToNewPyDict(boost::any_cast< OptionsList >(*entry));
  }
  else if (entry->type()==typeid(Teuchos::SerialDenseVector<int,int>))
  {
    try
      {
        Teuchos::SerialDenseVector< int,int > tvec =
          boost::any_cast< Teuchos::SerialDenseVector< int,int > >(*entry);
        return copyTeuchosVectorToNumPy(tvec);
      }
    catch(boost::bad_any_cast &e) {return NULL;}
  }
  else if (entry->type()==typeid(Teuchos::SerialDenseVector<int,double>))
  {
    try
      {
	Teuchos::SerialDenseVector< int,double > tvec =
	  boost::any_cast< Teuchos::SerialDenseVector< int,double > >(*entry);
	return copyTeuchosVectorToNumPy(tvec);
      }
    catch(boost::bad_any_cast &e) {return NULL;}
  }
  else if (entry->type()==typeid(Teuchos::SerialDenseMatrix<int,int>))
  {
    try
      {
        Teuchos::SerialDenseMatrix< int,int > tmat =
          boost::any_cast< Teuchos::SerialDenseMatrix< int,int > >(*entry);
        return copyTeuchosMatrixToNumPy(tmat);
      }
    catch(boost::bad_any_cast &e) {return NULL;}
  }
  else if (entry->type()==typeid(Teuchos::SerialDenseMatrix<int,double>))
  {
    try
      {
	Teuchos::SerialDenseMatrix< int,double > tmat =
	  boost::any_cast< Teuchos::SerialDenseMatrix< int,double > >(*entry);
	return copyTeuchosMatrixToNumPy(tmat);
      }
    catch(boost::bad_any_cast &e) {return NULL;}
  }
  else if (entry->type()==typeid(std::vector< OptionsList >))
  {
    try
      {
	std::vector< OptionsList > list =
	  boost::any_cast< std::vector< OptionsList > >(*entry);
	return copyStdVectorToPyList(list);
      }
    catch(boost::bad_any_cast &e) {return NULL;}
  }
  // All  other types are unsupported
  return NULL;
}

void pyListToNewStdVector(PyObject *pylist,
			  std::vector< OptionsList > &list){

  void * argp;
  int newmem = 0;
  static swig_type_info * swig_TPL_ptr =
    SWIG_TypeQuery("std::shared_ptr< OptionsList >*");
  //std::printf("pyListToNewStdVector\n");
  //print_pyobject_string_rep(pylist);
  if (PyList_Check(pylist)) {
    Py_ssize_t len = PyList_Size(pylist);
    list.reserve(len);
    for(Py_ssize_t i=0; i<len; ++i){
      PyObject *value = PyList_GetItem(pylist, i);
      //print_PyType(value);
      if (SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(value,
						       &argp,
						       swig_TPL_ptr,
						       0,
						       &newmem))){
	if (newmem & SWIG_CAST_NEW_MEMORY){ 
	  std::shared_ptr< OptionsList > tempshared =
	    *reinterpret_cast< std::shared_ptr< OptionsList > * >(argp);
	  delete reinterpret_cast< std::shared_ptr< OptionsList > * >(argp);
	  list.push_back(*(tempshared.get()));
	}else{
	  std::shared_ptr< OptionsList > * smartarg =
	    reinterpret_cast< std::shared_ptr< OptionsList > * >(argp);
	  if (smartarg) list.push_back(*(smartarg->get()));
	}
      }else if (PyDict_Check(value)){
	
	OptionsList * sublist = pyDictToNewOptionsList(value);
	list.push_back(*sublist);
	//delete sublist;
      }else
	throw(std::runtime_error("Passed PyObject pointer was not a Python dictionary or a swig wrapped instance of OptionsList"));
    }
  }else
    throw(std::runtime_error("Passed PyObject was not a list"));
}

PyObject * copyStdVectorToPyList(std::vector<OptionsList> &list){
  PyObject* pylist = PyList_New(list.size());
  if (!pylist)
    goto fail;
  for (size_t i=0; i<list.size(); ++i){
    PyObject *item =
      optionsListToNewPyDict(list[i]);
    if (!item) {
     Py_XDECREF(pylist);
     goto fail;
    }
    PyList_SET_ITEM(pylist, i, item);
  }
  return pylist;
 fail:
  return NULL;
}
  
}//  namespace Pecos
