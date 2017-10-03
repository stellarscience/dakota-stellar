/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HermiteInterpPolynomial
//- Description:  Implementation code for HermiteInterpPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "HermiteInterpPolynomial.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.hpp"
#endif
//#define INTERPOLATION_TEST


namespace Pecos {

/** Pre-compute denominator products that are only a function of the
    interpolation points. */
void HermiteInterpPolynomial::precompute_data()
{
  // set up divided difference tables for type1/type2 matching conditions

  int i, num_pts = interpPts.size(), n2 = 2*num_pts, n2m1 = n2 - 1;
  RealArray values(num_pts, 0.), derivs(num_pts, 0.);

  xValDiffTab.resize(n2);        xGradDiffTab.resize(n2m1);      // type1/2
  yT1ValDiffTab.resize(num_pts); yT1GradDiffTab.resize(num_pts); // type 1
  yT2ValDiffTab.resize(num_pts); yT2GradDiffTab.resize(num_pts); // type 2
  for (i=0; i<num_pts; ++i) {
    // type 1
    RealArray& yt1v_tab_i = yT1ValDiffTab[i];  yt1v_tab_i.resize(n2);
    RealArray& yt1g_tab_i = yT1GradDiffTab[i]; yt1g_tab_i.resize(n2m1);
    values[i] = 1.; // i-th type1 matching condition
    webbur::hermite_interpolant(num_pts, &interpPts[0], &values[0], &derivs[0],
				&xValDiffTab[0],  &yt1v_tab_i[0],
				&xGradDiffTab[0], &yt1g_tab_i[0]);
    values[i] = 0.; // restore
    // type 2
    RealArray& yt2v_tab_i = yT2ValDiffTab[i];  yt2v_tab_i.resize(n2);
    RealArray& yt2g_tab_i = yT2GradDiffTab[i]; yt2g_tab_i.resize(n2m1);
    derivs[i] = 1.; // i-th type1 matching condition
    webbur::hermite_interpolant(num_pts, &interpPts[0], &values[0], &derivs[0],
				&xValDiffTab[0],  &yt2v_tab_i[0],
				&xGradDiffTab[0], &yt2g_tab_i[0]);
    derivs[i] = 0.; // restore
#ifdef INTERPOLATION_TEST
    PCout << "xv =\n" << xValDiffTab << "xg =\n" << xGradDiffTab
	  << "yt1v[" << i+1 << "] =\n" << yt1v_tab_i
	  << "yt1g[" << i+1 << "] =\n" << yt1g_tab_i
	  << "yt2v[" << i+1 << "] =\n" << yt2v_tab_i
	  << "yt2g[" << i+1 << "] =\n" << yt2g_tab_i;
#endif // INTERPOLATION_TEST
  }

#ifdef INTERPOLATION_TEST
  int wp = WRITE_PRECISION+7;
  for (i=0; i<num_pts; ++i) {
    Real pt_i = interpPts[i];
    PCout << "Interp pt " << std::setw(3) << i+1 << " = "
	  << std::setw(wp) << pt_i << ":\n";
    for (int j=0; j<num_pts; ++j)
      PCout << "  interpolant " << j+1
	    << " t1 val = "  << std::setw(wp) << type1_value(pt_i, j)
	    << " t1 grad = " << std::setw(wp) << type1_gradient(pt_i, j)
	    << " t2 val = "  << std::setw(wp) << type2_value(pt_i, j)
	    << " t2 grad = " << std::setw(wp) << type2_gradient(pt_i, j) <<'\n';
  }
#endif // INTERPOLATION_TEST
}


/** Compute value of the Hermite type 1 polynomial corresponding to
    interpolation point i. */
Real HermiteInterpPolynomial::type1_value(Real x, unsigned short i)
{
  int n2 = 2*interpPts.size(); Real x_copy = x, t1_val, t1_grad;
  webbur::hermite_interpolant_value(n2, &xValDiffTab[0], &yT1ValDiffTab[i][0],
    &xGradDiffTab[0], &yT1GradDiffTab[i][0], 1, &x_copy, &t1_val, &t1_grad);
  return t1_val;
}


/** Compute value of the Hermite type 2 polynomial corresponding to
    interpolation point i. */
Real HermiteInterpPolynomial::type2_value(Real x, unsigned short i)
{
  int n2 = 2*interpPts.size(); Real x_copy = x, t2_val, t2_grad;
  webbur::hermite_interpolant_value(n2, &xValDiffTab[0], &yT2ValDiffTab[i][0],
    &xGradDiffTab[0], &yT2GradDiffTab[i][0], 1, &x_copy, &t2_val, &t2_grad);
  return t2_val;
}


/** Compute derivative with respect to x of the Hermite type 1
    polynomial corresponding to interpolation point i. */
Real HermiteInterpPolynomial::type1_gradient(Real x, unsigned short i)
{ 
  int n2 = 2*interpPts.size(); Real x_copy = x, t1_val, t1_grad;
  webbur::hermite_interpolant_value(n2, &xValDiffTab[0], &yT1ValDiffTab[i][0],
    &xGradDiffTab[0], &yT1GradDiffTab[i][0], 1, &x_copy, &t1_val, &t1_grad);
  return t1_grad;
}


/** Compute derivative with respect to x of the Hermite type 2
    polynomial corresponding to interpolation point i. */
Real HermiteInterpPolynomial::type2_gradient(Real x, unsigned short i)
{ 
  int n2 = 2*interpPts.size(); Real x_copy = x, t2_val, t2_grad;
  webbur::hermite_interpolant_value(n2, &xValDiffTab[0], &yT2ValDiffTab[i][0],
    &xGradDiffTab[0], &yT2GradDiffTab[i][0], 1, &x_copy, &t2_val, &t2_grad);
  return t2_grad;
}


const RealArray& HermiteInterpPolynomial::
collocation_points(unsigned short order)
{
  if (order < 1) {
    PCerr << "Error: underflow in minimum order (1) in PiecewiseInterp"
	  << "Polynomial::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  bool rule_err = false;
  if (interpPts.size() != order) { // if not already computed
    interpPts.resize(order);
#ifdef HAVE_SPARSE_GRID
    switch (collocRule) {
    case GAUSS_PATTERSON:
      webbur::patterson_lookup_points(order, &interpPts[0]);           break;
    case CLENSHAW_CURTIS: 
      webbur::clenshaw_curtis_compute_points(order, &interpPts[0]);    break;
    case FEJER2: 
      webbur::fejer2_compute_points(order, &interpPts[0]);             break;
    case GAUSS_LEGENDRE:
      if (order <= 33) // retrieve full precision tabulated values
	webbur::legendre_lookup_points(order, &interpPts[0]);
      else { // sandia_rules calculates points/weights together
	RealArray legendre_wts(order);
	webbur::legendre_compute(order, &interpPts[0], &legendre_wts[0]);
      }
      break;
    default:
      rule_err = true; break;
    }
#else
    rule_err = true;
#endif
  }

  if (rule_err) {
    PCerr << "Error: unsupported collocation rule in HermiteInterpPolynomial"
	  << "::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  return interpPts;
}


const RealArray& HermiteInterpPolynomial::
type1_collocation_weights(unsigned short order)
{
  if (order < 1) {
    PCerr << "Error: underflow in minimum order (1) in HermiteInterpPolynomial"
	  << "::type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  bool rule_err = false;
  if (interpPts.size() != order)
    collocation_points(order);
  if (type1InterpWts.size() != order) { // if not already computed
    type1InterpWts.resize(order);
#ifdef HAVE_SPARSE_GRID
    // hermite_interpolant_rule returns type1 & type 2 weights
    RealArray wts(2*order);
    // Note: requires integration bounds
    // TO DO: consider adding support for Gaussian weighting
    //        (enabling STD_NORMAL_U transformations)
    webbur::hermite_interpolant_rule(order, -1., 1., &interpPts[0], &wts[0]);
    if (type2InterpWts.size() == order) // type2 already updated
      for (size_t i=0; i<order; ++i)
	type1InterpWts[i] = wtFactor * wts[2*i];
    else { // update both so subsequent type2 call returns immediately
      type2InterpWts.resize(order);
      for (size_t i=0; i<order; ++i) {
	type1InterpWts[i] = wtFactor * wts[2*i];
	type2InterpWts[i] = wtFactor * wts[2*i+1];
      }
    }
#else
    rule_err = true;
#endif
  }

  if (rule_err) {
    PCerr << "Error: unsupported type1 collocation weights in HermiteInterp"
	  << "Polynomial::type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  return type1InterpWts;
}


const RealArray& HermiteInterpPolynomial::
type2_collocation_weights(unsigned short order)
{
  if (order < 1) {
    PCerr << "Error: underflow in minimum order (1) in HermiteInterpPolynomial"
	  << "::type2_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  bool rule_err = false;
  if (interpPts.size() != order)
    collocation_points(order);
  if (type2InterpWts.size() != order) { // if not already computed
    type2InterpWts.resize(order);
#ifdef HAVE_SPARSE_GRID
    // hermite_interpolant_rule returns type1 & type2 weights
    RealArray wts(2*order);
    // Note: requires integration bounds
    // TO DO: consider adding support for Gaussian weighting
    //        (enabling STD_NORMAL_U transformations)
    webbur::hermite_interpolant_rule(order, -1., 1., &interpPts[0], &wts[0]);
    if (type1InterpWts.size() == order) // type1 already updated
      for (size_t i=0; i<order; ++i)
	type2InterpWts[i] = wtFactor * wts[2*i+1];
    else { // update both so subsequent type1 call returns immediately
      type1InterpWts.resize(order);
      for (size_t i=0; i<order; ++i) {
	type1InterpWts[i] = wtFactor * wts[2*i];
	type2InterpWts[i] = wtFactor * wts[2*i+1];
      }
    }
#else
    rule_err = true;
#endif
  }

  if (rule_err) {
    PCerr << "Error: unsupported type2 collocation weights in HermiteInterp"
	  << "Polynomial::type2_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  return type2InterpWts;
}

} // namespace Pecos
