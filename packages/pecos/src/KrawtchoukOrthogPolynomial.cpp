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
  Real t1_val;
  Real omN = 1.0 - betaPoly;
  Real om1 = Real(order) - 1.0;

  switch (order) {
    case 0:
      t1_val = 1.;
      break;

    case 1:
      t1_val = (alphaPoly*Real(betaPoly) - x)/(alphaPoly*Real(betaPoly));
      break;

    case 2:
      t1_val = (alphaPoly*alphaPoly*betaPoly*omN + (1.0-2.0*alphaPoly*omN)*x - x*x)/(alphaPoly*alphaPoly*betaPoly*omN);
      break;

    default: {
      // Support higher order polynomials using the 3 point recursion formula:
      Real fm2 = type1_value(x, order-2);
      Real fm1 = type1_value(x, order-1);
      t1_val = ((alphaPoly*(betaPoly-om1)+om1*(1.0-alphaPoly)-x)*fm1 - om1*(1.0-alphaPoly)*fm2)/(alphaPoly*(betaPoly-om1));
      break;
    }
  }

  return t1_val;
}


} // namespace Pecos
