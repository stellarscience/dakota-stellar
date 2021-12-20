/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_lhs_driver.cpp
    \brief A driver program for PECOS */

#include "LHSDriver.hpp"


/// A driver program for PECOS.

/** Generates an LHS sample set from a DistributionParams specification. */

int main(int argc, char* argv[])
{
  // data set corresponds to short column test problem

  // Instantiate/initialize random variables
  std::cout << "Instantiating distribution parameters:\n";
  std::vector<Pecos::RandomVariable> ran_vars(3);
  ran_vars[0] = Pecos::RandomVariable(Pecos::NORMAL);
  ran_vars[0].push_parameter(Pecos::N_MEAN,    500.);
  ran_vars[0].push_parameter(Pecos::N_STD_DEV, 100.);
  ran_vars[1] = Pecos::RandomVariable(Pecos::NORMAL);
  ran_vars[1].push_parameter(Pecos::N_MEAN,   2000.);
  ran_vars[1].push_parameter(Pecos::N_STD_DEV, 400.);
  ran_vars[2] = Pecos::RandomVariable(Pecos::LOGNORMAL);
  ran_vars[2].push_parameter(Pecos::LN_MEAN,    5.);
  ran_vars[2].push_parameter(Pecos::LN_STD_DEV, 0.5);

  // Instantiate/initialize correlations
  Pecos::RealSymMatrix uv_corr(3); // inits to 0 -> specify nonzero lwr triangle
  uv_corr(0,0) = uv_corr(1,1) = uv_corr(2,2) = 1.; uv_corr(1,0) = 0.5;
  
  // Instantiate/initialize LHSDriver
  std::cout << "Instantiating LHSDriver:\n";
  Pecos::LHSDriver lhs_driver("lhs"); // default sample_ranks_mode, reports
  lhs_driver.seed(1234567);

  // Compute and output samples
  int num_samples = 100;
  Pecos::RealMatrix samples_array, sample_ranks;
  lhs_driver.generate_samples(ran_vars, uv_corr, num_samples,
			      samples_array, sample_ranks);
  //std::cout << "Samples:\n"; // << samples_array << '\n';
  //Pecos::write_data(std::cout, samples_array, false, true, true);

  return 0;
}
