/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef TEST_OPTIONS_LIST_HPP
#define TEST_OPTIONS_LIST_HPP

#include "OptionsList.hpp"

#define TEST_FUNC_PROTOS(TYPE, SNAME)				  \
								  \
  Pecos::util::OptionsList SNAME ## set_entry(Pecos::util::OptionsList &opts, \
		    const std::string &name, const TYPE & item);  \


TEST_FUNC_PROTOS(int, int)
TEST_FUNC_PROTOS(double, double)
TEST_FUNC_PROTOS(std::string, string)
TEST_FUNC_PROTOS(Pecos::util::OptionsList, optionslist)

#endif // TEST_OPTIONS_LIST_HPP

// based upon Vector.h Vector.cxx testVector.py at
// https://github.com/numpy/numpy/tree/master/tools/swig/test
