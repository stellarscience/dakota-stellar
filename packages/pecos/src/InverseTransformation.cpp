/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "InverseTransformation.hpp"
#include <cmath>

static const char rcsId[]="@(#) $Id: InverseTransformation.cpp 4768 2007-12-17 17:49:32Z mseldre $";

//#define DEBUG


namespace Pecos {

void InverseTransformation::
initialize(const Real& total_t, const Real& w_bar, size_t seed)
{
  bool err_flag = false;
  if (total_t < 0.) {
    PCerr << "Error: total time must be non-negative." << std::endl;
    err_flag = true;
  }
  if (w_bar <= 0.) {
    PCerr << "Error: cut-off frequency must be positive." << std::endl;
    err_flag = true;
  }
  if (err_flag)
    abort_handler(-1);

  totalTime = total_t;
  omegaBar  = w_bar;

  // omegaBar and totalTime specify deltaTime and m
  deltaTime  = 2.*PI/omegaBar;  // rad/sec -> sec
  size_t m   = 1 + (size_t)std::floor(totalTime/deltaTime);
  deltaOmega = omegaBar/(m-1);

  timeSequence.sizeUninitialized(m);
  omegaSequence.sizeUninitialized(m);
  for (size_t i=0; i<m; i++) {
    timeSequence[i]  = i*deltaTime;  // from 0 to final time <= totalTime
    omegaSequence[i] = i*deltaOmega; // from 0 to omegaBar
  }

  lhsSampler.seed(seed);
}


void InverseTransformation::
power_spectral_density(const String& psd_name, const Real& param)
{
  size_t i, m = omegaSequence.length();
  bool err_flag = false;
  if (!m) {
    PCerr << "Error: initialize() must be called prior to "
	  << "power_spectral_density()." << std::endl;
    err_flag = true;
  }
  if (err_flag)
    abort_handler(-1);

  psdSequence.sizeUninitialized(m);
  //if (psd_name == "white_noise")            // TO DO
  //else if (psd_name == "rectangular_pulse") // TO DO
  if (psd_name == "band_limited_white_noise")
    // One-sided PSD for (unit-variance) band-limited white noise:
    // param is upper bound wc (which differs from omegaBar)
    //        { 1/wc   w \in [0,wc]
    // g(w) = { 
    //        { 0      otherwise
    for (i=0; i<m; i++)
      psdSequence[i] = (omegaSequence[i] <= param) ? 1./param : 0.;
  else if (psd_name == "increasing_linear" || psd_name == "decreasing_linear") {
    // One-sided linear PSD with unit variance (s^2 = area) defining either the
    // lower right or lower left half portion of band_limited_white_noise:
    // param is upper bound wc (which differs from omegaBar)
    Real intercept, slope;
    if (psd_name == "increasing_linear") {
      // For positive slope (lower right triangle - emphasize high freq):
      // g(w) = 2/wc^2 w for w in [0,wc], 0 otherwise
      intercept = 0.;
      slope     = 2./param/param;
    }
    else {
      // For negative slope (lower left triangle - emphasize low freq):
      // g(w) = 2/wc - 2/wc^2 w for w in [0,wc], 0 otherwise
      intercept = 2./param;
      slope     = -intercept/param;
    }
    for (i=0; i<m; i++)
      psdSequence[i] = (omegaSequence[i] <= param) ?
	intercept + slope * omegaSequence[i] : 0.;
  }
  else if (psd_name == "first_order_markov") {
    // One-sided PSD for (unit-variance) 1st-order Markov process:
    // param is lambda
    //            2 * lam
    // g(w) = ----------------
    //        pi (w^2 + lam^2)
    Real p_2_over_pi = 2.*param/PI, param_sq = param*param;
    for (i=0; i<m; i++)
      psdSequence[i] = p_2_over_pi/(std::pow(omegaSequence[i], 2) + param_sq);
  }
  else if (psd_name == "second_order_markov") {
    // One-sided PSD for (unit-variance) 2st-order Markov process:
    // param is lambda
    //            4 * lam^3
    // g(w) = -------------------
    //        pi (w^2 + lam^2)^2
    Real param_sq = param*param, p_sq_4_over_pi = 4.*param_sq/PI;
    for (i=0; i<m; i++)
      psdSequence[i]
	= p_sq_4_over_pi/std::pow(std::pow(omegaSequence[i], 2) + param_sq, 2);
  }
}


/** Incoming psd if a set of (x,y) pairs that define a PSD curve (w, g(w)).  
    The incoming discretization cannot be assumed to match that selected
    for psdSequence; therefore, linear interpolation is performed in any
    combination of log or linear scales for the two axes. */
void InverseTransformation::power_spectral_density(const RealRealPairArray& psd)
{
  // Interpolation: psd[i].second -> psdSequence

  // interpolation: log-log , linear-log, log-linear, linear-linear
}

} // namespace Pecos
