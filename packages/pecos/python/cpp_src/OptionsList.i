/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

%module(implicitconv="1", autodoc="1",package="PyDakota") options_list
%{
#define SWIG_FILE_WITH_INIT
// Required for interfacing with NumPy
#include "numpy_include.hpp"//gets rid of deprecated warnings

#include "OptionsList.hpp"
#include "python_helpers.hpp"

using Pecos::util::OptionsList;
%}


%include <std_shared_ptr.i>
%shared_ptr(Pecos::util::OptionsList)


%feature("docstring") Pecos::util::OptionsList
"The ``OptionsList`` class is used for communicating arbitrary-type
options to functions."

// General ignore directives
%ignore *::operator=; // ignore any member operator=
%ignore *::print;     // ignore any member functions name print
%ignore operator<<;   // ignore any operator<< (including member functions)
%ignore operator==;   // ignore any operator!= (including member functions)
%ignore operator!=;   // ignore any operator!= (including member functions)

// We utilize very little of numpy.i, so some of the fragments
// do not get automatically instantiated.  This forces the issue.
%include "numpy.i"
%fragment("NumPy_Fragments");
%init %{
  import_array();
%}

%include "std_string.i"

/////////////////////////////////////////////////////////////////////
// Override typemaps for OptionsLists to allow PyDicts as input //
/////////////////////////////////////////////////////////////////////
%define %pydict_overrides(CONST, CLASS...)
// Input a plain reference
%typemap(in) CONST CLASS &
(void *argp=0, int res=0, bool cleanup=false, std::shared_ptr< CONST CLASS > tempshared)
{
  if (PyDict_Check($input))
  {
    $1 = Pecos::pyDictToNewOptionsList($input);
    if (!$1) SWIG_fail;
    cleanup = true;
  }
  else
  {
    int newmem = 0;
    res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(std::shared_ptr< CLASS > *),
				%convertptr_flags, &newmem);
    if (!SWIG_IsOK(res))
    {
      %argument_fail(res, "$type", $symname, $argnum);
    }
    if (!argp)
    {
      %argument_nullref("$type", $symname, $argnum);
    }
    if (newmem & SWIG_CAST_NEW_MEMORY)
    {
      tempshared = *%reinterpret_cast(argp, std::shared_ptr< CONST CLASS > *);
      delete %reinterpret_cast(argp, std::shared_ptr< CONST CLASS > *);
      $1 = %const_cast(tempshared.get(), $1_ltype);
    }
    else
    {
      $1 = %const_cast(%reinterpret_cast(argp, std::shared_ptr< CONST CLASS > *)->get(), $1_ltype);
    }
  }
}
// Perform type checking
%typemap(typecheck,precedence=SWIG_TYPECHECK_POINTER,noblock=1)
  CONST CLASS,
  CONST CLASS &,
  CONST CLASS *,
  CONST CLASS *&,
  std::shared_ptr< CONST CLASS >,
  std::shared_ptr< CONST CLASS > &,
  std::shared_ptr< CONST CLASS > *,
  std::shared_ptr< CONST CLASS > *&
{
  // Accept PyDicts or CLASS instances
  $1 = PyDict_Check($input);
  if (!$1)
  {
    int res = SWIG_ConvertPtr($input, 0, $descriptor(std::shared_ptr< CLASS > *), 0);
    $1 = SWIG_CheckState(res);
  }
}

// Cleanup
%typemap(freearg) CONST CLASS &
{
  if (cleanup$argnum && $1) delete $1;
}

// Input a std::shared_ptr by reference
%typemap(in) std::shared_ptr< CONST CLASS > &
(void *argp, int res = 0, $*1_ltype tempshared)
{
  if (PyDict_Check($input))
  {
    tempshared = std::shared_ptr<Pecos::util::OptionsList>(Pecos::pyDictToNewOptionsList($input));
    if (!tempshared) SWIG_fail;
    $1 = &tempshared;
  }
  else
  {
    int newmem = 0;
    res = SWIG_ConvertPtrAndOwn($input, &argp,
				$descriptor(std::shared_ptr< CLASS > *),
				%convertptr_flags, &newmem);
    if (!SWIG_IsOK(res))
    {
      %argument_fail(res, "$type", $symname, $argnum);
    }
    if (newmem & SWIG_CAST_NEW_MEMORY)
    {
      if (argp) tempshared = *%reinterpret_cast(argp, $ltype);
      delete %reinterpret_cast(argp, $ltype);
      $1 = &tempshared;
    }
    else
    {
      $1 = (argp) ? %reinterpret_cast(argp, $ltype) : &tempshared;
    }
  }
}
%enddef

%pydict_overrides(SWIGEMPTYHACK, Pecos::util::OptionsList)
%pydict_overrides(const        , Pecos::util::OptionsList)

%rename(_print) Pecos::util::OptionsList::print() const;

%include "OptionsList.hpp"
 //%include "test_options_list.hpp"

%extend Pecos::util::OptionsList
{
  /******************************************************************/
  // Dictionary constructor
  OptionsList(PyObject * dict)
  {
    OptionsList * opts_list;
    if (!PyDict_Check(dict)) {
      // setting PyErr_SetString and returning NULL does not seem
      // to raise an exception. so throw hard error
      throw(std::runtime_error("Argument must be a Python dictionary"));
    }
    opts_list =
      Pecos::pyDictToNewOptionsList(dict);
    if (opts_list == NULL) goto fail;
    return opts_list;
  fail:
    return NULL;
  }
  
  /******************************************************************/
  // Set method: accept only python objects as values
  void set(const std::string &name, PyObject *value)
  {
    if (!Pecos::setPythonParameter(*self,name,value))
    {
      PyErr_SetString(PyExc_TypeError, "OptionsList value type not supported");
    }
  }

  /******************************************************************/
  // Get method: return entries as python objects
  PyObject * get(const std::string &name, PyObject * default_value=NULL) const
  {
    PyObject * value = Pecos::getPythonParameter(*self, name);
    // Type not supported
    if (value == NULL)
    {
      PyErr_SetString(PyExc_TypeError, "OptionsList value type not supported");
      goto fail;
    }
    // Name not found
    else if (value == Py_None)
    {
      if (default_value == NULL)
      {
  	PyErr_Format(PyExc_KeyError, "'%s'", name.c_str());
  	goto fail;
      }
      Py_DECREF(value);
      Py_INCREF(default_value);
      return default_value;
    }
    // Type supported and name found
    else
      return value;
  fail:
    Py_XDECREF(value);
    return NULL;
  }

  /******************************************************************/
  // Length operator
  int __len__() const
  {
    return (*self).size();
  }

  /******************************************************************/
  // Equals operators
  
  PyObject * __eq__(const OptionsList & opts_list) const
  {
    if ((*self)==opts_list) return Py_True;
    else return Py_False;
  }

  PyObject * __eq__(PyObject * obj) const
  {
    PyObject * dict   = Pecos::optionsListToNewPyDict(*self);
    PyObject * result = 0;
    if (dict == NULL) goto fail;
    result = PyObject_RichCompare(dict,obj,Py_EQ);
    Py_DECREF(dict);
    return result;
  fail:
    return NULL;
  }

  /******************************************************************/
  // GetItem operator
  PyObject * __getitem__(const std::string & name) const
  {
    PyObject * value = Pecos::getPythonParameter(*self, name);
    // Type not supported
    if (value == NULL)
    {
      PyErr_SetString(PyExc_TypeError, "OptionsList value type not supported");
      goto fail;
    }
    // Name not found
    else if (value == Py_None)
    {
      PyErr_Format(PyExc_KeyError, "'%s'", name.c_str());
      goto fail;
    }
    // Type supported and name found
    else
      return value;
  fail:
    Py_XDECREF(value);
    return NULL;
  }

  /******************************************************************/
  // Contains operator
  int __contains__(const std::string & name) const
  {
    return (int)self->exists(name);
  }

  /******************************************************************/
  // SetItem operator
  void __setitem__(const std::string & name, PyObject * value)
  {
    if (!Pecos::setPythonParameter(*self,name,value))
    {
      PyErr_SetString(PyExc_TypeError, "OptionsList value type not supported");
    }
  }

  /******************************************************************/
  // String conversion method
  PyObject * __str__() const
  {
    PyObject * dict = Pecos::optionsListToNewPyDict(*self);
    PyObject * str  = 0;
    if (dict == NULL) goto fail;
    str = PyObject_Str(dict);
    Py_DECREF(dict);
    return str;
  fail:
    return NULL;
  }

  /******************************************************************/
  // String representation method
  PyObject * __repr__() const
  {
    std::string reprStr;
    PyObject * dict    = Pecos::optionsListToNewPyDict(*self);
    PyObject * dictStr = 0;
    PyObject * result = 0;
    if (dict == NULL) goto fail;
    dictStr = PyObject_Str(dict);
    reprStr = std::string("OptionsList(") +
              std::string(PyString_AsString(dictStr)) +
              std::string(")");
    result = PyString_FromString(reprStr.c_str());
    Py_DECREF(dict   );
    Py_DECREF(dictStr);
    return result;
  fail:
    Py_XDECREF(dict);
    return NULL;
  }
}
