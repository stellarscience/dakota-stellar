/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

%module dot

%{
    #define SWIG_FILE_WITH_INIT
    #include "dot.hpp"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* vec1), (int len2, double* vec2)}


%include "dot.hpp"
%rename (dot) my_dot;

%inline %{
    double my_dot(int len1, double* vec1, int len2, double* vec2) {
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError, "Arrays of lengths (%d,%d) given", len1, len2);
        return 0.0;
    }
    return dot(len1, vec1, vec2);
}
%}
