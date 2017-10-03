/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef BOOST_TEST_DIST_HPP
#define BOOST_TEST_DIST_HPP

/// Class for testing Boost statistical functions


class boost_test_dist
{
public:

  /// default constructor
  boost_test_dist() { }
  /// destructor
  ~boost_test_dist() { }

  /// output quantities for comparison between GSL and Boost
  void print_comparison();
};

#endif // BOOST_TEST_DIST_HPP

