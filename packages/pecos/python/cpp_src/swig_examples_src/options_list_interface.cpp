/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "options_list_interface.hpp"

using Pecos::util::OptionsList;

#define TEST_FUNCS(TYPE, SNAME)					\
								\
  OptionsList SNAME ## set_entry(OptionsList &opts,	\
		    const std::string &name, const TYPE &item){	\
    opts.set(name,item);					\
    return opts;						\
  }								

TEST_FUNCS(int, int)
TEST_FUNCS(double, double)
TEST_FUNCS(std::string, string)
TEST_FUNCS(OptionsList, optionslist)
