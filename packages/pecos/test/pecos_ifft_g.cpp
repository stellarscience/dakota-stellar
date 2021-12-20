/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_ifft_g.cpp
    \brief A driver program for PECOS */

#include "DataTransformation.hpp"
#include "pecos_data_types.hpp"

//#define TIMING

#ifdef TIMING
#include <ctime>
#endif // TIMING

using namespace Pecos;


/// A driver program for PECOS.

/** Generates samples from an analytic PSD for a MATLAB example
    problem using the Grigoriu iFFT algorithm. */

int main(int argc, char* argv[])
{
  // Instantiate/initialize the data transformation instance which manages
  // the ProbabilityTransformation and BasisFunction instances.
  DataTransformation ifft_transform("inverse_fourier_grigoriu");

  // Constants for this problem
  Real   vbar  = 10000.; // cut-off frequency (rad/s)
  Real   T     = 5.;     // stop time (sec)
  size_t nseed = 314;    // random number seed
  size_t ns    = 50;     // number of samples

  // Demonstrate 1st-order Markov internally-defined PSD
  ifft_transform.initialize(T, vbar, nseed);
  ifft_transform.power_spectral_density("first_order_markov", 1000.);

  // compute and return ns samples all at once
  const RealMatrix& samples = ifft_transform.compute_samples(ns);
#ifndef TIMING
  PCout << "Sample matrix:\n" << samples << std::endl;
#endif // TIMING

  // compute and return two more, one at a time
  for (size_t i=ns; i<ns+2; i++) {
    const RealVector& sample = ifft_transform.compute_sample();
#ifndef TIMING
    PCout << "Sample " << i+1 << ":\n" << sample << std::endl;
#endif // TIMING
  }

  // compute all at once and return one at a time (?)
  // requires Matrix(i,:) -> Vector conversion
  //ifft_transform.compute_samples(ns);
  //for (size_t i=0; i<ns; i++) {
  //  const RealVector& sample = ifft_transform.sample(i);
  //  PCout << "Sample " << i+1 << ":\n" << sample_i << std::endl;
  //}

#ifdef TIMING
  PCout << "Processor time: " << (double)clock()/(double)CLOCKS_PER_SEC
	<< " seconds." << std::endl;
#endif // TIMING

  return 0;
}
