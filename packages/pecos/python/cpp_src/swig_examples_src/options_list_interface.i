/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

%module(implicitconv="1", autodoc="1",package="PyDakota.swig_examples") options_list_interface
%{
#include "numpy_include.hpp"//gets rid of deprecated warnings
#include "options_list_interface.hpp"
#include "python_helpers.hpp"
%}
%include "numpy.i"
%fragment("NumPy_Fragments");
%init %{
  import_array();
%}

%import "OptionsList.i"
%include "options_list_interface.hpp"
