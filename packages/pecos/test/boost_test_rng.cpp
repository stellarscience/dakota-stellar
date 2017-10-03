/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "boost_test_rng.hpp"
#include "pecos_data_types.hpp"
#include "pecos_global_defs.hpp"
#include <boost/random/linear_congruential.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <algorithm>
#include <iomanip>


// This is a reproducible simulation experiment.  See main().
void experiment(boost::mt19937& generator)
{
  using namespace boost;
  // Define a uniform random number distribution of integer values between
  // 1 and 6 inclusive.
  //typedef variate_generator<mt19937&, uniform_int<> > gen_type;
  variate_generator<mt19937&, uniform_int<> >
    die_gen(generator, uniform_int<>(1, 6));

  // If you want to use an STL iterator interface, use iterator_adaptors.hpp.
  generator_iterator<variate_generator<mt19937&, uniform_int<> > >
    die(&die_gen);
  for(int i = 0; i < 10; ++i)
    PCout << *die++ << " ";
  PCout << '\n';
}


int main(int argc, char* argv[])
{
  using namespace boost;
  //boost_test_rng bt;

  PCout.setf(std::ios::fixed);
  PCout.precision(16); 
  PCout.setf(std::ios::scientific);
  // PCout.setf(std::ios::showpoint);

  // Define a RNG and initialize it with a reproducible seed.
  // (The seed is unsigned, otherwise the wrong overload may be selected
  // when using mt19937 as the base_generator_type.)
  mt19937 generator(41u);
  lagged_fibonacci607 generator2(41u);

  PCout << "10 samples of a uniform distribution in [0..1):\n";

  // Define a uniform random number distribution which produces "double"
  // values between 0 and 1 (0 inclusive, 1 exclusive).
  uniform_real<> uni_dist(0,1);
  variate_generator<mt19937&, uniform_real<> > uni(generator, uni_dist);
  
  //  new mt19937() mers_twister;
  //new lagged_fibonacci607() lag_fib;
  
  // You can now retrieve random numbers from that distribution by means
  // of a STL Generator interface, i.e. calling the generator like a C-function
  // with no arguments.
  for(int i = 0; i < 10; ++i){
    PCout << "Uniform:  " << uni() << '\n';
    PCout << "Mersenne Twister:  " << generator() << '\n';
    PCout << "Lagged Fibonacci:  " << generator2() << '\n';
    //PCout << "Uniform:  " << uni() << '\n';
  }

  generator.seed(static_cast<unsigned int>(std::time(0)));

  PCout << "\nexperiment: roll a die 10 times:\n";

  // You can save a generator's state by copy construction/assignment.
  mt19937 saved_generator = generator;

  // When calling other functions which take a generator or distribution
  // as a parameter, make sure to always call by reference (or pointer).
  experiment(generator);

  PCout << "redo the experiment to verify it:\n";
  experiment(saved_generator);

  // After that, both generators are equivalent
  //assert(generator == saved_generator);

  // as a degenerate case, you can set min = max for uniform_int
  uniform_int<> degen_dist(4,4);
  variate_generator<mt19937&, uniform_int<> > deg(generator, degen_dist);
  PCout << deg() << " " << deg() << " " << deg() << std::endl;

  return 0;
}
