/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_lhs_driver.cpp
    \brief A driver program for PECOS */

#include "DistributionParams.hpp"
#include "LHSDriver.hpp"


/// A driver program for PECOS.

/** Generates an LHS sample set from a DistributionParams specification. */

int main(int argc, char* argv[])
{
  // data set corresponds to short column test problem
  Pecos::RealVector num_means(2); num_means[0] = 500.; num_means[1] = 2000.;
  Pecos::RealVector num_std_devs(2);
  num_std_devs[0] = 100.; num_std_devs[1] = 400.;
  Pecos::RealVector lnuv_means(1);    lnuv_means[0]    = 5.;
  Pecos::RealVector lnuv_std_devs(1); lnuv_std_devs[0] = 0.5;
  Pecos::RealSymMatrix uv_corr(3); // inits to 0 -> specify nonzero lwr triangle
  uv_corr(0,0) = uv_corr(1,1) = uv_corr(2,2) = 1.; uv_corr(1,0) = 0.5;

  // Instantiate/initialize DistributionParams
  std::cout << "Instantiating distribution parameters:\n";
  Pecos::AleatoryDistParams dp; // default ctor
  dp.normal_means(num_means);     dp.normal_std_deviations(num_std_devs);
  dp.lognormal_means(lnuv_means); dp.lognormal_std_deviations(lnuv_std_devs);
  dp.uncertain_correlations(uv_corr);

  // Instantiate/initialize LHSDriver
  std::cout << "Instantiating LHSDriver:\n";
  Pecos::LHSDriver lhs_driver("lhs"); // default sample_ranks_mode, reports
  lhs_driver.seed(1234567);

  // Compute and output samples
  int num_samples = 100;
  Pecos::RealMatrix samples_array;
  lhs_driver.generate_samples(dp, num_samples, samples_array);
  std::cout << "Samples:\n"; // << samples_array << '\n';
  Pecos::write_data(std::cout, samples_array, false, true, true);

  return 0;
}
