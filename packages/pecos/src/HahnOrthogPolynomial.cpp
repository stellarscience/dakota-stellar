/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HahnOrthogPolynomial
//- Description:  Implementation code for HahnOrthogPolynomial class
//-               
//- Owner:        Russell Hooper, Sandia National Laboratories

#include "HahnOrthogPolynomial.hpp"
#include "pecos_stat_util.hpp"


namespace Pecos {

Real HahnOrthogPolynomial::type1_value(Real x, unsigned short order)
{
  Real t1_val, nd = (Real)numDrawn, tp = (Real)totalPop, sp = (Real)selectPop;

  switch (order) {
  case 0:
    t1_val = 1.;
    break;

  case 1:
    t1_val = 1. + (2.+tp+sp)/(-nd*(tp+1.))*x;
    break;

  case 2:
    t1_val = 1. - 2.*(3.+tp+sp)*x/(nd*(tp+1.))
           + (3.+tp+sp)*(4.+tp+sp)/((tp+1.)*(tp+2.)*nd*(nd-1.))*x*(x-1.);
    break;

  default: {
    // Support higher order polynomials using the 3 point recursion formula:
    Real Ha_nm1 = 1. + (2.+tp+sp)/(-nd*(tp+1.))*x,                        //1
         Ha_n   = 1. - 2.*(3.+tp+sp)*x/(nd*(tp+1.))
           + (3.+tp+sp)*(4.+tp+sp)/((tp+1.)*(tp+2.)*nd*(nd-1.))*x*(x-1.); //2
    Real om1, A, C;
    for (size_t i=3; i<=order; i++) {
      om1 = (Real)i-1.;
      A = (om1+tp+sp+1.)*(om1+tp+1.)*(nd-om1)
	/ ((2.*om1+tp+sp+1.)*(2.*om1+tp+sp+2.));
      C = om1*(om1+tp+sp+nd+1.)*(om1+sp)/((2.*om1+tp+sp)*(2.*om1+tp+sp+1.));
      t1_val = ((A+C-x)*Ha_n - C*Ha_nm1)/A; // Ha_nplus1
      if (i < order)
	{ Ha_nm1 = Ha_n;  Ha_n = t1_val; }
    }
    break;
  }
  }

  return t1_val;
}

} // namespace Pecos
