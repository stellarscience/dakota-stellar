/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        KrawtchoukOrthogPolynomial
//- Description:  Implementation code for KrawtchoukOrthogPolynomial class
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#include "KrawtchoukOrthogPolynomial.hpp"
#include "pecos_stat_util.hpp"


namespace Pecos {

Real KrawtchoukOrthogPolynomial::type1_value(Real x, unsigned short order)
{
  Real t1_val, nt = (Real)numTrials;

  switch (order) {
  case 0:
    t1_val = 1.;
    break;

  case 1:
    t1_val = 1. - x/(probPerTrial*nt);
    break;

  case 2: {
    Real omN = 1. - nt, ppt2 = probPerTrial*probPerTrial;
    t1_val = 1. + x*(1.-2.*probPerTrial*omN - x)/(ppt2*nt*omN);
    break;
  }
  default: {
    // Support higher order polynomials using the 3 point recursion formula:
    Real omN = 1.-nt, ppt2 = probPerTrial*probPerTrial,
      Kc_nm1 = 1. - x/(probPerTrial*nt),                         //1
      Kc_n   = 1. + x*(1.-2.*probPerTrial*omN - x)/(ppt2*nt*omN);//2
    Real om1, A, C;
    for (size_t i=3; i<=order; i++) {
      om1 = (Real)i-1.;
      A = probPerTrial*(nt-om1), C = om1*(1.-probPerTrial);
      t1_val = ((A+C-x)*Kc_n - C*Kc_nm1)/A; // Kc_nplus1
      if (i < order)
	{ Kc_nm1 = Kc_n;  Kc_n = t1_val; }
    }
    break;
  }
  }

  return t1_val;
}

} // namespace Pecos
