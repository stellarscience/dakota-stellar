/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: pecos_stat_util.cpp 4768 2007-12-17 17:49:32Z mseldre $";

namespace Pecos {

/*
#ifndef HAVE_BOOST
#ifdef HAVE_GSL
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

//  Solve is performed in scaled space (for the standard beta distribution).
Real cdf_beta_Pinv(const Real& cdf, const Real& alpha, const Real& beta)
{
  // F(x) = Phi(z) = cdf
  // F(x) - cdf = 0
  // x^{i+1} = x - (F(x) - cdf)/f(x)

  // Initial guess: model as uniform
  // (linear CDF accumulation between [0,1] bounds)
  Real scaled_x = cdf;

  // evaluate residual F(x) - cdf
  Real res = gsl_cdf_beta_P(scaled_x, alpha, beta) - cdf;

  size_t newton_iters = 0, max_iters = 20;
  bool converged = false, terminate = false;
  Real convergence_tol = 1.e-4; // hardwired for now
  while (!terminate && !converged) {

    // evaluate residual derivative f(x)
    Real dres_dx = gsl_ran_beta_pdf(scaled_x, alpha, beta);

    // compute Newton step
    Real delta_scaled_x;
    if (std::fabs(dres_dx) > DBL_MIN) {
      delta_scaled_x = -res/dres_dx; // full Newton step
      if (std::fabs(delta_scaled_x) < convergence_tol)
	converged = true; // but go ahead and take the step, if beneficial
    }
    else
      terminate = true;

    // Simple backtracking line search globalization
    bool reduction = false;
    size_t backtrack_iters = 0;
    while (!reduction && !terminate) {
      Real scaled_x_step = scaled_x + delta_scaled_x;
      // scaled_x must lie in [0,1]
      if (scaled_x_step < 0.) scaled_x_step = 0.;
      if (scaled_x_step > 1.) scaled_x_step = 1.;

      // evaluate residual at scaled_x_step

      Real res_step = gsl_cdf_beta_P(scaled_x_step, alpha, beta) - cdf;

      // perform backtracking line search to enforce decrease in res
      if ( std::fabs(res_step) < std::fabs(res) ) { // accept step
	reduction = true;
	scaled_x = scaled_x_step;
	res      = res_step;
	//PCout << "residual = " << res << " delta = " << delta_scaled_x
	//      << " scaled_x = " << scaled_x << '\n';
      }
      else if (converged)
	terminate = true; // kick out of inner while
      else {
	//PCout << "Backtracking\n";
	delta_scaled_x /= 2.; // backtrack
	if (backtrack_iters++ >= max_iters) // backtrack iter must complete
	  terminate = true;
      }
    }
    if (++newton_iters >= max_iters) // Newton iteration has completed
      terminate = true;
  }

  return scaled_x;
}

#else

static Real erf_inverse(const Real& p)
{
  // Adapted from ltqnorm.m see URL below for more info

  // The algorithm uses a minimax approximation by rational functions and
  // the result has a relative error less than 1.15e-9.  A last refinement
  // by Halley's rational method is applied to achieve full machine precision.

  // Author:      Peter J. Acklam
  // Time-stamp:  2000-07-19 16:44:07
  // E-mail:      pjacklam@online.no
  // URL:         http://home.online.no/~pjacklam/notes/invnorm/

  Real a[7] = { 0.,         -3.969683028665376e+01,  2.209460984245205e+02,
    -2.759285104469687e+02,  1.383577518672690e+02, -3.066479806614716e+01,
     2.506628277459239e+00 };
  Real b[6] = { 0.,         -5.447609879822406e+01,  1.615858368580409e+02,
    -1.556989798598866e+02,  6.680131188771972e+01, -1.328068155288572e+01 };
  Real c[7] = { 0.,         -7.784894002430293e-03, -3.223964580411365e-01,
    -2.400758277161838e+00, -2.549732539343734e+00,  4.374664141464968e+00,
     2.938163982698783e+00 };
  Real d[5] = { 0.,          7.784695709041462e-03,  3.224671290700398e-01,
     2.445134137142996e+00,  3.754408661907416e+00 };

  // Modify p to since this is the error function
  // 1/std::sqrt(2.)*ltqnorm((p+1)/2)=erf_inverse(p)
  // redefine p here
  Real p_new = 0.5*(p + 1.);
  // scale at the end

  // Define break-points.
  Real plow  = 0.02425;
  Real phigh = 1. - plow;

  // Rational approximation for central region
  Real r, q, z;
  if (p_new > plow && p_new < phigh) { // -plow <= p_new <= phigh 
    q = p_new - 0.5;
    r = q*q;
    z = (((((a[1]*r+a[2])*r+a[3])*r+a[4])*r+a[5])*r+a[6])*q/ 
        (((((b[1]*r+b[2])*r+b[3])*r+b[4])*r+b[5])*r+1);
  }
  // Rational approximation for lower region
  else if (p_new < plow && p_new > 0.) { // plow > p_new > 0.
    q = std::sqrt(-2*std::log(p_new));
    z = (((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6])/
        ((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1);  
  }
  else if (p_new > phigh && p_new < 1.) { // phigh < p_new < 1.
    q = std::sqrt(-2*(std::log(1-p_new)));
    z = -(((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6])/
         ((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1); 
  }
  else if (p_new == 0.)
    z = -1.*std::pow(10., 150.);
  else if (p_new == 1.)
    z = std::pow(10., 150.); 
  else if (p_new < 0. || p_new > 1.) {
    PCerr << "Error: probability greater than 1 or less than 0 in "
	  << "erf_inverse()." << std::endl;
    abort_handler(-1);
  }
  // use erf instead of erfc
  Real e = 0.5*(1. - erf(-z/std::sqrt(2.))) - p_new;
  Real u = e*std::sqrt(2.*PI)*std::exp(z*z/2.);
  z = z - u/(1. + z*u/2.);
  // scale since this is the erf inverse not gaussian inverse
  // see above
  z = 1/std::sqrt(2.)*z;
  return z;
}

#endif // HAVE_GSL
#endif // HAVE_BOOST
*/

} // namespace Pecos
