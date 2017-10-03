/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef BOOST_TEST_RNG_HPP
#define BOOST_TEST_RNG_HPP

/// Base class for testing Boost random number generators


class boost_test_rng
{
public:
 
  /// default constructor
  boost_test_rng() { }
 
  /// destructor
  ~boost_test_rng() { }

  //  void experiment(boost::mt19937& generator);

};

#endif  // BOOST_TEST_RNG_HPP

