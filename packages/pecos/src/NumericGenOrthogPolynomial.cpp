/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        NumericGenOrthogPolynomial
//- Description:  Implementation code for NumericGenOrthogPolynomial class
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "NumericGenOrthogPolynomial.hpp"
#include "ExponentialRandomVariable.hpp"
#include "Teuchos_LAPACK.hpp"
#ifdef HAVE_SPARSE_GRID
#include "sandia_rules.hpp"
#endif // HAVE_SPARSE_GRID

//#define DEBUG
//#define ADAPTIVE_DEBUG

namespace Pecos {


// TO DO:
// KDE: look at figTree --> if lightweight; else, code simple box kernel and
// Gaussian kernel
//
// sampling on stoch expansion within \xi --> kernel added for each sample -->
// approximate PDF for response --> accept a KDE PDF within NumGenOrthogPoly
//
// activate STOCHASTIC_EXPANSION allowing either moments (SC) or KDE (general)


/** Numbering conventions follow Gautschi. */
void NumericGenOrthogPolynomial::solve_eigenproblem(unsigned short m)
{
  // set up sets of flags to avoid unnecessary storage and computation
  bool pre_process_coeffs_norms = false, post_process_coeffs       = false,
    post_process_native_norm_sq = false, post_process_last_norm_sq = false;
  switch (distributionType) {
  // **************************
  // Analytic recursion coeffs:
  // **************************
  case LOGNORMAL:
    if (coeffsNormsFlag) // PCE
      pre_process_coeffs_norms = post_process_coeffs = true;
    break;
  // ********************************
  // Chebyshev moment-based approach:
  // ********************************
  //case LOGNORMAL:
  //case WEIBULL:
  case STOCHASTIC_EXPANSION:
    if (coeffsNormsFlag) // PCE
      pre_process_coeffs_norms = post_process_coeffs
        = post_process_native_norm_sq = post_process_last_norm_sq = true;
    break;
  // *******************************
  // Discretized Stieltjes approach:
  // *******************************
  default:
    pre_process_coeffs_norms = true;
    if (coeffsNormsFlag) // PCE
      post_process_last_norm_sq = true;
    break;
  }

  // size class-scope arrays
  if (pre_process_coeffs_norms) {
    polyCoeffs.resize(m+1);
    polyCoeffs[0].sizeUninitialized(1);
    polyCoeffs[0][0] = 1.; // constant term
    orthogPolyNormsSq.sizeUninitialized(m+1);
    orthogPolyNormsSq[0] = 1.; // integral of 1*1*PDF is 1
  }
  if (!m)
    return;

  int i, j;
  // TO DO: improve data persistence
  //if (alpha3TR.empty() || beta3TR.empty()) // || distParamsUpdate
    { alpha3TR.sizeUninitialized(m); beta3TR.sizeUninitialized(m); }
  //else if (m > alpha3TR.length() || m > beta3TR.length())
  //  { alpha3TR.resize(m); beta3TR.resize(m); } // preserve previous
  RealVector off_diag(m-1, false);
  switch (distributionType) {
  case LOGNORMAL: {
    // -----------------------------------------------------------------------
    // Recursion coefficients are available analytically:
    // (a) Simpson, I.C., "Numerical Integration over a Semi-Infinite Interval
    //     Using the Lognormal Distibution," Numer. Math. 31, pp.71-76, 1978.
    // (b) Wilck, M., "A general approximation method for solving integrals
    //     containing a lognormal weighting function," Aerosol Science 32,
    //     pp. 1111-1116, 2001.
    // -----------------------------------------------------------------------
    Real zeta_sq = distParams[1]*distParams[1], // denoted as "u" in Simpson
       e_zeta_sq = std::exp(zeta_sq), lambda = distParams[0],
       e_lambda  = std::exp(lambda);            // denoted as "m" in Simpson
    beta3TR[0] = 0.; // both Simpson and Wilck
    // compute analytic recursion coefficients alpha_i and beta_i
    for (i=0; i<m; ++i) {
      alpha3TR[i] = e_lambda * (std::pow(e_zeta_sq,i)*(e_zeta_sq + 1.) - 1.) *
	std::pow(e_zeta_sq,(2.*i-1.)/2.); // Simpson
      //alpha3TR[i] = ((e_zeta_sq + 1.) * std::pow(e_zeta_sq,i) - 1.) *
      //  std::pow(e_zeta_sq,(2.*i-1.)/2.); // Wilck
      if (i) {
	beta3TR[i] = e_lambda * e_lambda * (std::pow(e_zeta_sq,i) - 1.) *
	  std::pow(e_zeta_sq,3*i-2); // Simpson
	//beta3TR[i] = (std::pow(e_zeta_sq,i) - 1.) *
	//  std::pow(e_zeta_sq,3*i-2); // Wilck
	off_diag[i-1] = std::sqrt(beta3TR[i]);
      }
    }
    if (coeffsNormsFlag) {
      Real prod = 1.;
      for (i=1; i<=m; ++i) {
	prod *= std::pow(e_zeta_sq,i) - 1.;
	orthogPolyNormsSq[i] = std::pow(e_lambda,2*i) * prod * 
	  std::pow(e_zeta_sq,i*(3.*i-1.)/2.); // Simpson (not provided in Wilck)
      }
    }

    // Chebyshev option available below

    // Discretized option: use default approach below combined with discrete
    // inner_product using ~100k steps that sum pc1*pc2*pdf from 0->50.
    // --> this approach verified the Simpson results.
    // --> Wilck may work as well with different orthogPolyNormsSq scaling,
    //     (these norms are not provided), or there may be a typo (implied
    //     by dangling reference to "s" in Eqs. 10 and 11).
    break;
  }
  //case LOGNORMAL:
  //case WEIBULL:
  case STOCHASTIC_EXPANSION: {
    // -------------------------------------------------------------------------
    // Raw moments are available analytically: use moment-based approach
    // referred to as the (unmodified) Chebyshev algorithm by Gautschi -> avoids
    // need to approximate inner products for recursion coefficients, although
    // the highest order orthogPolyNormsSq[m] cannot be computed reliably with
    // native quadrature (must either be evaluated using another *_integral()
    // fn or another m+1 eigenproblem must be solved for m+1 Gauss pts/wts).
    // This method is unstable and loses accuracy for higher m -> therefore, it
    // is not the default; however, it is the _only_ method to not require a
    // full PDF of the random variable.
    // -------------------------------------------------------------------------
    beta3TR[0] = 1.;
    // compute raw moments
    RealVector raw_moments(2*m+1);
    raw_moments[0] = 1.; // integral of 1*PDF is 1

    // LOGNORMAL:
    //Real zeta_sq = distParams[1]*distParams[1],
    //     lambda  = distParams[0];
    //for (i=1; i<=2*m; ++i)
    //  raw_moments[i] = std::exp(i*lambda+i*i*zeta_sq/2.);

    // WEIBULL:
    //for (i=1; i<=2*m; ++i)
    //  raw_moments[i] = std::pow(distParams[1], i)
    //	* gamma_function(1. + i/distParams[0]);

    // STOCHASTIC_EXPANSION:
    // TO DO

    // form Hankel matrix
    RealSymMatrix H(m+1, false);
    for (i=0; i<=m; ++i)
      for (j=0; j<=i; ++j)
	H(i,j) = raw_moments[i+j];
    // perform Cholesky factorization
    RealSpdSolver chol_solver;
    chol_solver.setMatrix( Teuchos::rcp(&H, false) );
    chol_solver.factor(); // Cholesky factorization (LL^T) in place
    // compute recursion coefficients
    for (i=0; i<m; ++i) {
      alpha3TR[i] = H(i+1,i)/H(i,i); // lower triangle
      if (i) {
	Real Him1im1 = H(i-1,i-1);
	alpha3TR[i]  -= H(i,i-1) / Him1im1; // lower triangle
	off_diag[i-1] = H(i,i)   / Him1im1;
	beta3TR[i]    = std::pow(off_diag[i-1], 2);
      }
    }
    break;
  }
  default: {
    // ---------------------------------------------------------------------
    // The default approach is the discretized Stieltjes algorithm
    // (W. Gautschi, SIAM J. Sci. Stat. Comput., Vol. 3, No. 3, Sept. 1982)
    // ---------------------------------------------------------------------
    // set up the coefficients
    beta3TR[0] = orthogPolyNormsSq[0];
    RealVector xi_poly_coeffs_i;
    for (i=0; i<m; ++i) {
      const RealVector& poly_coeffs_i = polyCoeffs[i];
      if (i) { // compute norm^2, beta, and off_diag
	orthogPolyNormsSq[i] = inner_product(poly_coeffs_i, poly_coeffs_i);
	beta3TR[i] = orthogPolyNormsSq[i] / orthogPolyNormsSq[i-1];
	off_diag[i-1] = std::sqrt(beta3TR[i]);
      }
      // form xi_poly_coeffs_i to compute alpha_i
      int i_len = poly_coeffs_i.length();
      xi_poly_coeffs_i.sizeUninitialized(i_len+1);
      xi_poly_coeffs_i[0] = 0.;
      for (j=0; j<i_len; ++j)
	xi_poly_coeffs_i[j+1] = poly_coeffs_i[j];
      alpha3TR[i] = inner_product(xi_poly_coeffs_i, poly_coeffs_i)
	          / orthogPolyNormsSq[i];
      // update polyCoeffs[i+1]
      if (i == 0)
	polynomial_recursion(polyCoeffs[i+1], alpha3TR[i], poly_coeffs_i);
      else if (i < m-1 || coeffsNormsFlag)
	polynomial_recursion(polyCoeffs[i+1], alpha3TR[i], poly_coeffs_i,
			     beta3TR[i], polyCoeffs[i-1]);
    }
    break;
  }
  }

  if (post_process_coeffs) {
    polynomial_recursion(polyCoeffs[1], alpha3TR[0], polyCoeffs[0]);
    for (i=1; i<m; ++i)
      polynomial_recursion(polyCoeffs[i+1], alpha3TR[i], polyCoeffs[i],
			   beta3TR[i], polyCoeffs[i-1]);
  }

#ifdef DEBUG
  for (i=0; i<m; ++i)
    PCout << "alpha[" << i << "] = " << alpha3TR[i]
	  << " beta[" << i << "] = " << beta3TR[i] << '\n';
  if (coeffsNormsFlag)
    for (i=0; i<=m; ++i)
      PCout << "polyCoeffs[" << i << "] =\n" << polyCoeffs[i];
#endif

  // solve the symmetric tridiagonal eigenvalue problem using LAPACK
  int ldz = std::max((unsigned short)1,m);
  RealMatrix z_eigvec(ldz, m, false);
  Teuchos::LAPACK<int, Real> la;
  int info = 0;
  double* work = new double [std::max(1,2*m-2)]; // temporary work array
  // DSTEQR docs for 3rd field: (input/output)
  //   On entry, the diagonal elements of the tridiagonal matrix.
  //   On exit, if INFO = 0, the eigenvalues in ascending order.
  copy_data(alpha3TR, collocPoints); // eigenvalues are Gauss points
  la.STEQR('I', m, &collocPoints[0], off_diag.values(), z_eigvec.values(), ldz,
	   work, &info);
  if (info) {
    PCerr << "Error: nonzero return code (" << info << ") from LAPACK STEQR "
	  << "(symmetric tridiagonal eigensolution)\n       in "
	  << "NumericGenOrthogPolynomial::solve_eigenproblem()" << std::endl;
    abort_handler(-1);
  }
  delete [] work;

  // Gauss points are the eigenvalues which are updated in place by STEQR.
  // compute the Gauss weights from the eigenvectors:
  collocWeights.resize(m);
  //Real norm_sq_0 = orthogPolyNormsSq[0];
  for (i=0; i<m; ++i)
    collocWeights[i] = std::pow(z_eigvec(0, i), 2.);// * norm_sq_0;

  // orthogPolyNormsSq up to order m-1 are available using the just
  // computed Gauss points/weights
  if (post_process_native_norm_sq)
    for (i=0; i<m; ++i)
      orthogPolyNormsSq[i]
	= native_quadrature_integral(polyCoeffs[i], polyCoeffs[i]);
  // orthogPolyNormsSq of order m computed with (expensive) discrete integral
  if (post_process_last_norm_sq) // 2m integral exceeds native 2m-1 resolution
    orthogPolyNormsSq[m] = inner_product(polyCoeffs[m], polyCoeffs[m]);

#ifdef DEBUG
  if (coeffsNormsFlag)
    for (i=0; i<=m; ++i)
      PCout << "orthogPolyNormsSq[" << i << "] = " << orthogPolyNormsSq[i]
	    <<'\n';
  PCout << "collocPoints:\n"  << collocPoints
	<< "collocWeights:\n" << collocWeights;
#endif
}


Real NumericGenOrthogPolynomial::
inner_product(const RealVector& poly_coeffs1,
	      const RealVector& poly_coeffs2)
{
  // Use either classical Gaussian quadratures (Legendre/Jacobi,
  // Laguerre/Gen Laguerre, and Hermite) or discrete sums for bounded,
  // semi-bounded, and unbounded integration domains corresponding to
  // supported PDFs.
  Real dbl_inf = std::numeric_limits<Real>::infinity();
  switch (distributionType) {
  // *************************
  // * BOUNDED DISTRIBUTIONS *
  // *************************
  case BOUNDED_NORMAL: {
    // We must trap the case where only one finite bound has been specified:
    // Quick & dirty: trap +/-inf and replace with mu +/- x sigma
    //                with x sufficiently large [phi(15) ~ 5.53e-50].
    // More involved: replace bounded integral with semibounded integral;
    //                but we would need to add an option for [-inf, ub]:
    //                monotone mapping from [-inf, ub] to [-1,1] using
    //                x = ub + (1-z)/(1+z); dx/dz = -2/(1+z)^2
    Real lb = (real_compare(distParams[2], -dbl_inf)) ?
      distParams[0] - 15. * distParams[1] : distParams[2];
    Real ub = (real_compare(distParams[3],  dbl_inf)) ?
      distParams[0] + 15. * distParams[1] : distParams[3];

    // Alternate integrations:
    //return legendre_bounded_integral(poly_coeffs1, poly_coeffs2,
    //  BoundedNormalRandomVariable::pdf, lb, ub);
    return cc_bounded_integral(poly_coeffs1, poly_coeffs2, bounded_normal_pdf,
			       lb, ub, 500);
    //return riemann_bounded_integral(poly_coeffs1, poly_coeffs2,
    //  BoundedNormalRandomVariable::pdf, lb, ub);
    break;
  }
  case BOUNDED_LOGNORMAL: {
    // We must trap the case where a finite upper bound has not been specified:
    // Quick & dirty: trap +/-inf and replace with mu +/- x sigma
    //                with x sufficiently large [phi(15) ~ 5.53e-50].
    // More involved: replace bounded integral with semibounded integral;
    //                but we would need to add an option for [-inf, ub]:
    //                monotone mapping from [-inf, ub] to [-1,1] using
    //                x = ub + (1-z)/(1+z); dx/dz = -2/(1+z)^2
    Real ub = (real_compare(distParams[3], dbl_inf)) ?
      distParams[0] + 15. * distParams[1] : distParams[3];

    // Alternate integrations:
    //return legendre_bounded_integral(poly_coeffs1, poly_coeffs2,
    //  BoundedLognormalRandomVariable::pdf, distParams[2], ub);
    return cc_bounded_integral(poly_coeffs1, poly_coeffs2,
			       bounded_lognormal_pdf, distParams[2], ub, 500);
    //return riemann_bounded_integral(poly_coeffs1, poly_coeffs2,
    //  BoundedLognormalRandomVariable::pdf, distParams[2], ub);
    break;
  }
  case LOGUNIFORM:
    // Alternate integrations:
    //return legendre_bounded_integral(poly_coeffs1, poly_coeffs2,
    //  LoguniformRandomVariable::pdf, distParams[0], distParams[1]);
    return cc_bounded_integral(poly_coeffs1, poly_coeffs2, loguniform_pdf,
			       distParams[0], distParams[1], 500);
    //return riemann_bounded_integral(poly_coeffs1, poly_coeffs2,
    //  LoguniformRandomVariable::pdf, distParams[0], distParams[1]);
    break;
  case TRIANGULAR:
    // Alternate integrations:
    //return legendre_bounded_integral(poly_coeffs1, poly_coeffs2,
    //  TriangularRandomVariable::pdf, distParams[1], distParams[2]);
    return cc_bounded_integral(poly_coeffs1, poly_coeffs2, triangular_pdf,
			       distParams[0], distParams[2], 500);
    //return riemann_bounded_integral(poly_coeffs1, poly_coeffs2,
    //  TriangularRandomVariable::pdf, distParams[1], distParams[2]);
    break;
  case HISTOGRAM_BIN: {

    // BMA TODO: Use Gauss-Legendre integral on each uniform bin with
    // a number of points sufficient to exactly integrate p(x), q(x)
    // against a constant function (bin value).  Though see note about
    // Gautschi below...

    size_t dp_len = distParams.length(),
      u_bnd_index = (dp_len>=2) ? dp_len-2 : 0;
    // Alternate integrations:
    //return legendre_bounded_integral(poly_coeffs1, poly_coeffs2,
    //  HistogramBinRandomVariable::pdf,distParams[0],distParams[u_bnd_index]);
    return cc_bounded_integral(poly_coeffs1, poly_coeffs2,
      HistogramBinRandomVariable::pdf, distParams[0], distParams[u_bnd_index],
      50*dp_len); // 100 per bin

    // TO DO: consider FACTOR*ORDER1*ORDER2(*DP_LEN?); POLY ORDER IS MAIN ISSUE

    // Note: histograms provide a finite set of PDF data.
    // Gautschi: you cannot compute #Gauss pts > #discrete data points.
    // --> consider throwing an error (at a higher level) for order > #bins

    //return riemann_bounded_integral(poly_coeffs1, poly_coeffs2,
    //	HistogramBinRandomVariable::pdf,distParams[0],distParams[u_bnd_index]);
    break;
  }
  case HISTOGRAM_PT_INT: case HISTOGRAM_PT_STRING: case HISTOGRAM_PT_REAL: {
    // TOP: https://www.math.lsu.edu/~xlwan/papers/journal/megpc2.pdf
    // https://www.math.lsu.edu/~xlwan/papers/journal/beyondAskey.pdf
    // https://www.cs.purdue.edu/homes/wxg/Madrid.pdf
    // https://www.math.lsu.edu/~xlwan/publication.html

    // Discrete integral sum_i { w_i * p(x_i) * q(x_i) } for the
    // histogram points (x_i, w_i).  Assumes weights sum to 1.0.

    // distParams[even]: abscissas; distParams[odd]: weights (masses)
    Real sum = 0.0;
    for (size_t i=0; i<distParams.length(); i+=2) {
      Real v1 = type1_value(distParams[i], poly_coeffs1);
      sum += distParams[i+1] * v1 * type1_value(distParams[i], poly_coeffs2);
    }
    return sum;
    break;
  }
  // ******************************
  // * SEMI-BOUNDED DISTRIBUTIONS *
  // ******************************
  // Alternate integration:
  case LOGNORMAL:
    // Alternate integrations:
    //return laguerre_semibounded_integral(poly_coeffs1, poly_coeffs2,
    //					   LognormalRandomVariable::pdf);
    return fejer_semibounded_integral(poly_coeffs1, poly_coeffs2,
				      lognormal_pdf, 500);
    // medium left & right tail;
    // start is offset to avoid division by 0. in LognormalRandomVariable::pdf
    //return legendre_bounded_integral(poly_coeffs1, poly_coeffs2,
    //                                 LognormalRandomVariable::pdf,
    //				       1.e-5, mean+60.*stdev);
    //return cc_bounded_integral(poly_coeffs1, poly_coeffs2,
    //                           LognormalRandomVariable::pdf,
    //			         1.e-5, mean+60.*stdev, 500);
    //return riemann_bounded_integral(poly_coeffs1, poly_coeffs2,
    //                                LognormalRandomVariable::pdf,
    //                                1.e-5, mean+60.*stdev);
    break;
  case FRECHET: {
    // Alternate integration:
    // Better stability w/ Laguerre for Rosenbrock Frechet (infinite variance?)
    return laguerre_semibounded_integral(poly_coeffs1, poly_coeffs2,
					 frechet_pdf);
    //return fejer_semibounded_integral(poly_coeffs1, poly_coeffs2,
    //                                  FrechetRandomVariable::pdf, 500);
    //Real mean, stdev;
    //moments_from_frechet_params(distParams[0], distParams[1],
    //				         mean, stdev);
    // minimal left tail with heavy right tail;
    // start is offset to avoid division by 0. in FrechetRandomVariable::pdf
    //return legendre_bounded_integral(poly_coeffs1, poly_coeffs2,
    //                                 FrechetRandomVariable::pdf,
    //				       0.1, mean+300.*stdev);
    //return cc_bounded_integral(poly_coeffs1, poly_coeffs2,
    //                           FrechetRandomVariable::pdf, 0.1,
    //			         mean+300.*stdev, 500);
    //return riemann_bounded_integral(poly_coeffs1, poly_coeffs2,
    //                                FrechetRandomVariable::pdf,
    //                                0.1, mean+300.*stdev);
    break;
  }
  case WEIBULL: {
    // Alternate integration:
    //return laguerre_semibounded_integral(poly_coeffs1, poly_coeffs2,
    //                                     WeibullRandomVariable::pdf);
    return fejer_semibounded_integral(poly_coeffs1, poly_coeffs2,
				      weibull_pdf, 500);
    //Real mean, stdev;
    //moments_from_weibull_params(distParams[0], distParams[1], mean, stdev);
    // heavy left tail with minimal right tail;
    // start is offset to avoid division by negative power of 0. in
    // WeibullRandomVariable::pdf
    //return legendre_bounded_integral(poly_coeffs1, poly_coeffs2,
    //                                 WeibullRandomVariable::pdf,
    //				       1.e-50, mean+30.*stdev);
    //return cc_bounded_integral(poly_coeffs1, poly_coeffs2,
    //                           WeibullRandomVariable::pdf,
    //			         1.e-50, mean+30.*stdev, 500);
    //return riemann_bounded_integral(poly_coeffs1, poly_coeffs2,
    //                                WeibullRandomVariable::pdf,
    //                                1.e-50, mean+30.*stdev);
    break;
  }
  // ***************************
  // * UNBOUNDED DISTRIBUTIONS *
  // ***************************
  case GUMBEL: {
    // Alternate integration:
    //return hermite_unbounded_integral(poly_coeffs1, poly_coeffs2,
    //                                  GumbelRandomVariable::pdf);
    return fejer_unbounded_integral(poly_coeffs1, poly_coeffs2,
				    gumbel_pdf, 500);
    //Real mean, stdev;
    //moments_from_gumbel_params(distParams[0], distParams[1],
    //				        mean, stdev);
    // left and right tail treated equally
    //return legendre_bounded_integral(poly_coeffs1, poly_coeffs2,
    //                                 GumbelRandomVariable::pdf,
    //				       mean-25.*stdev, mean+25.*stdev);
    //return cc_bounded_integral(poly_coeffs1, poly_coeffs2,
    //                           GumbelRandomVariable::pdf,
    //			         mean-25.*stdev, mean+25.*stdev, 500);
    //return riemann_bounded_integral(poly_coeffs1, poly_coeffs2,
    //                                GumbelRandomVariable::pdf,
    //			              mean-25.*stdev, mean+25.*stdev);
    break;
  }
  default:
    PCerr << "Error: unsupported distribution type in NumericGenOrthog"
	  << "Polynomial::inner_product()." << std::endl;
    abort_handler(-1);
    break;
  }
  return 0.;
}


Real NumericGenOrthogPolynomial::
hermite_unbounded_integral(const RealVector& poly_coeffs1,
			   const RealVector& poly_coeffs2, NGFPType weight_fn)
{
  BasisPolynomial hermite_poly(HERMITE_ORTHOG);
  unsigned short quad_order = 170; // hardwired (could use adaptive loop)
  //quad_order = 175 has numerical problems -> nan's in recursion coeffs
  const RealArray& colloc_pts = hermite_poly.collocation_points(quad_order);
  const RealArray& colloc_wts
    = hermite_poly.type1_collocation_weights(quad_order);

  Real sum = 0., v1, gp_i;
  for (size_t i=0; i<quad_order; ++i) {
    gp_i = colloc_pts[i];
    v1 = type1_value(gp_i, poly_coeffs1); // cached due to update of const ref
    sum += colloc_wts[i] * v1 * type1_value(gp_i, poly_coeffs2)
        *  weight_fn(gp_i, distParams) / NormalRandomVariable::std_pdf(gp_i);
  }
  return sum;
}


Real NumericGenOrthogPolynomial::
fejer_unbounded_integral(const RealVector& poly_coeffs1,
			 const RealVector& poly_coeffs2, NGFPType weight_fn,
			 unsigned short quad_order)
{
  //unsigned short quad_order = 500; // hardwired (could use adaptive loop)
  RealVector fejer_pts(quad_order, false), fejer_wts(quad_order, false);
#ifdef HAVE_SPARSE_GRID
  webbur::fejer2_compute(quad_order, fejer_pts.values(), fejer_wts.values());
#else
  PCerr << "Error: quadrature package required for NumericGenOrthogPolynomial::"
	<< "fejer_unbounded_integral()." << std::endl;
  abort_handler(-1);
#endif // HAVE_SPARSE_GRID

  // monotone mapping from [-inf,inf] to [-1,1] using x = z/(1-z^2)
  // dx/dz = (1+z^2)/(1-z^2)^2
  Real sum = 0., v1, z_sq, one_m_z_sq, x_i, z_i;
  for (size_t i=0; i<quad_order; ++i) {
    z_i = fejer_pts[i];
    z_sq = z_i * z_i; one_m_z_sq = 1. - z_sq; x_i = z_i / one_m_z_sq;
    v1 = type1_value(x_i, poly_coeffs1); // cached due to update of const ref
    sum += fejer_wts[i] * v1 * type1_value(x_i, poly_coeffs2)
        *  weight_fn(x_i, distParams) * (1. + z_sq) / one_m_z_sq / one_m_z_sq;
  }
  // Notes on weighting: two different conceptual approaches
  // (1) fejer_wts assume a weighting fn = 1 (see http://people.sc.fsu.edu/
  //     ~burkardt/cpp_src/sandia_rules/sandia_rules.html) -> integrand is just
  //     poly1(x) * poly2(x) * weight_fn(x) * dx/dz where x ~ (-inf,inf) and
  //     z ~ (-1,1) on open intervals.
  // (2) could correct fejer_wts to PDF weighting by dividing by 2. and then
  //     divide sum by UniformRV::std_pdf() as in legendre_bounded_integral
  //     below --> cancels out.
  return sum;
}


Real NumericGenOrthogPolynomial::
laguerre_semibounded_integral(const RealVector& poly_coeffs1,
			      const RealVector& poly_coeffs2,
			      NGFPType weight_fn)
{
  BasisPolynomial laguerre_poly(LAGUERRE_ORTHOG);
  unsigned short quad_order = 95; // hardwired (could use adaptive loop)
  //quad_order = 100 has numerical problems: colloc_wts = inf
  const RealArray& colloc_pts = laguerre_poly.collocation_points(quad_order);
  const RealArray& colloc_wts
    = laguerre_poly.type1_collocation_weights(quad_order);

  Real sum = 0., v1, gp_i;
  for (size_t i=0; i<quad_order; ++i) {
    gp_i = colloc_pts[i];
    v1 = type1_value(gp_i, poly_coeffs1); // cached: update of const ref
    sum += colloc_wts[i] * v1 * type1_value(gp_i, poly_coeffs2)
      *  weight_fn(gp_i, distParams) / ExponentialRandomVariable::std_pdf(gp_i);
  }
  return sum;
}


Real NumericGenOrthogPolynomial::
fejer_semibounded_integral(const RealVector& poly_coeffs1,
			   const RealVector& poly_coeffs2,
			   NGFPType weight_fn, unsigned short quad_order)
{
  //unsigned short quad_order = 500; // hardwired (could use adaptive loop)
  RealVector fejer_pts(quad_order, false), fejer_wts(quad_order, false);
#ifdef HAVE_SPARSE_GRID
  webbur::fejer2_compute(quad_order, fejer_pts.values(), fejer_wts.values());
#else
  PCerr << "Error: quadrature package required for NumericGenOrthogPolynomial::"
	<< "fejer_semibounded_integral()." << std::endl;
  abort_handler(-1);
#endif // HAVE_SPARSE_GRID

  // monotone mapping from [0,inf] to [-1,1] using x = (1+z)/(1-z)
  // dx/dz = 2/(1-z)^2
  Real sum = 0., v1, one_m_z, x_i, z_i;
  for (size_t i=0; i<quad_order; ++i) {
    z_i = fejer_pts[i];
    one_m_z =  1. - z_i; x_i = (1. + z_i) / one_m_z;
    v1 = type1_value(x_i, poly_coeffs1); // cached due to update of const ref
    sum += fejer_wts[i] * v1 * type1_value(x_i, poly_coeffs2)
        *  weight_fn(x_i, distParams) * 2. / (one_m_z * one_m_z);
  }
  // Notes on weighting: two different conceptual approaches
  // (1) fejer_wts assume a weighting fn = 1 (see http://people.sc.fsu.edu/
  //     ~burkardt/cpp_src/sandia_rules/sandia_rules.html) -> integrand is just
  //     poly1(x) * poly2(x) * weight_fn(x) * dx/dz where x ~ (0,inf) and
  //     z ~ (-1,1) on open intervals.
  // (2) could correct fejer_wts to PDF weighting by dividing by 2. and then
  //     divide sum by UniformRV::std_pdf() as in legendre_bounded_integral
  //     below --> cancels out.
  return sum;
}


Real NumericGenOrthogPolynomial::
legendre_bounded_integral(const RealVector& poly_coeffs1,
			  const RealVector& poly_coeffs2, NGFPType weight_fn,
			  Real start, Real end)
{
  BasisPolynomial legendre_poly(LEGENDRE_ORTHOG);
  unsigned short quad_order = 50; // hardwired (could use adaptive loop)
  const RealArray& colloc_pts = legendre_poly.collocation_points(quad_order);
  const RealArray& colloc_wts
    = legendre_poly.type1_collocation_weights(quad_order);

  Real sum = 0., unscaled_gp_i, v1, half_range = (end - start) / 2.;
  for (size_t i=0; i<quad_order; ++i) {
    unscaled_gp_i = start + half_range * (colloc_pts[i]+1.);
    // change of variables: dx/dz = half_range
    v1 = type1_value(unscaled_gp_i, poly_coeffs1);// cached: update of const ref
    sum += colloc_wts[i] * v1 * type1_value(unscaled_gp_i, poly_coeffs2)
        *  weight_fn(unscaled_gp_i, distParams);
  }
  return sum / UniformRandomVariable::std_pdf() * half_range;
}


Real NumericGenOrthogPolynomial::
cc_bounded_integral(const RealVector& poly_coeffs1,
		    const RealVector& poly_coeffs2, NGFPType weight_fn,
		    Real start, Real end, unsigned short quad_order)
{
  //unsigned short quad_order = 500; // hardwired (could use adaptive loop)
  RealVector cc_pts(quad_order, false), cc_wts(quad_order, false);
#ifdef HAVE_SPARSE_GRID
  webbur::clenshaw_curtis_compute(quad_order, cc_pts.values(),
				  cc_wts.values());
#else
  PCerr << "Error: quadrature package required for NumericGenOrthogPolynomial::"
	<< "cc_bounded_integral()." << std::endl;
  abort_handler(-1);
#endif // HAVE_SPARSE_GRID

  Real sum = 0., unscaled_cc_i, v1, half_range = (end - start) / 2.;
  for (size_t i=0; i<quad_order; ++i) {
    unscaled_cc_i = start + half_range * (cc_pts[i]+1.);
    // change of variables: dx/dz = half_range
    v1 = type1_value(unscaled_cc_i, poly_coeffs1);// cached: update of const ref
    sum += cc_wts[i] * v1 * type1_value(unscaled_cc_i, poly_coeffs2)
        *  weight_fn(unscaled_cc_i, distParams);
  }
  // Notes on weighting: two different conceptual approaches
  // (1) cc_wts assume a weighting function = 1 (see http://people.sc.fsu.edu/
  //     ~burkardt/cpp_src/sandia_rules/sandia_rules.html) -> integrand is just
  //     poly1 * poly2 * weight_fn.
  // (2) could correct cc_wts to PDF weighting by dividing by 2. and then
  //     divide sum by UniformRV::std_pdf() as in legendre_bounded_integral
  //     above --> cancels out.
  return sum * half_range;
}


Real NumericGenOrthogPolynomial::
riemann_bounded_integral(const RealVector& poly_coeffs1,
			 const RealVector& poly_coeffs2, NGFPType weight_fn,
			 Real start, Real end)
{
  /*
  // fixed num_points:
  unsigned int num_points = num_partitions + 1;
  Real i_sum = 0., w, v1, w_sum = 0.,
    delta = (end-start)/num_partitions, x = start;
  for (unsigned int i=0; i<num_points; ++i) {
    w  = weight_fn(x, distParams);
    v1 = type1_value(x, poly_coeffs1); // cached due to update of const ref
    i_sum += v1 * type1_value(x, poly_coeffs2) * w;
    w_sum += w;
    x     += delta;
  }
  return i_sum / w_sum;
  */

  // adaptive num_points: 10 iters from 2000 results in a max of 10^6 points
  unsigned int i, num_partitions = 2000, num_points = num_partitions + 1,
    iter = 0, iter_max = 10;
  Real int_val, prev_int_val, i_sum = 0., w, v1, w_sum = 0., x = start,
    delta = (end-start)/num_partitions, rel_change = DBL_MAX, conv_tol = 1.e-6;
#ifdef ADAPTIVE_DEBUG
  PCout << "start = " << start << " end = " << end << std::endl;
#endif // ADAPTIVE_DEBUG
  while (rel_change > conv_tol && ++iter <= iter_max) {
    if (iter > 1) {
      if (iter == 2)
	num_points -= 1;
      else
	{ num_points *= 2; delta /= 2.; }
      x = start + delta/2.;
    }
    for (i=0; i<num_points; ++i) {
      w  = weight_fn(x, distParams);
#ifdef ADAPTIVE_DEBUG
      PCout << "x = " << x << " w = " << w << std::endl;
#endif // ADAPTIVE_DEBUG
      v1 = type1_value(x, poly_coeffs1); // cached due to update of const ref
      i_sum += v1 * type1_value(x, poly_coeffs2) * w;
      w_sum += w;
      x     += delta;
    }
    int_val = i_sum / w_sum;
    if (iter > 1)
      rel_change = std::fabs(int_val/prev_int_val - 1.);
    prev_int_val = int_val;
#ifdef ADAPTIVE_DEBUG
    PCout << "iter = " << iter << " delta = " << delta << " num_points = "
	  << num_points << " rel_change = " << rel_change << std::endl;
#endif // ADAPTIVE_DEBUG
  }
  return int_val;
}


Real NumericGenOrthogPolynomial::
native_quadrature_integral(const RealVector& poly_coeffs1,
			   const RealVector& poly_coeffs2)
{
  Real i_sum = 0., v1, gp_i;
  size_t i, num_colloc_pts = collocPoints.size();
  for (i=0; i<num_colloc_pts; ++i) {
    gp_i = collocPoints[i];
    v1 = type1_value(gp_i, poly_coeffs1); // cached due to update of const ref
    i_sum += v1 * type1_value(gp_i, poly_coeffs2) * collocWeights[i];
  }
  return i_sum;
}


Real NumericGenOrthogPolynomial::alpha_recursion(unsigned short order)
{
  if (alpha3TR.length() <= order)
    solve_eigenproblem(order+1);
  return alpha3TR[order];
}


Real NumericGenOrthogPolynomial::beta_recursion(unsigned short order)
{
  if (beta3TR.length() <= order)
    solve_eigenproblem(order+1);
  return beta3TR[order];
}


Real NumericGenOrthogPolynomial::type1_value(Real x, unsigned short order)
{
  if (polyCoeffs.size() <= order)
    solve_eigenproblem(order);
  return type1_value(x, polyCoeffs[order]);
}


Real NumericGenOrthogPolynomial::
type1_value(Real x, const RealVector& poly_coeffs)
{
  size_t num_terms = poly_coeffs.length();
  Real t1_val = poly_coeffs[0];
  for (int i=1; i<num_terms; ++i)
    t1_val += poly_coeffs[i]*std::pow(x,i);
  return t1_val;
}


Real NumericGenOrthogPolynomial::type1_gradient(Real x, unsigned short order)
{
  if (polyCoeffs.size() <= order)
    solve_eigenproblem(order);
  return type1_gradient(x, polyCoeffs[order]);
}


Real NumericGenOrthogPolynomial::
type1_gradient(Real x, const RealVector& poly_coeffs)
{
  // differentiate poly_coeffs with respect to x once
  size_t num_terms = poly_coeffs.length();
  Real t1_grad = (num_terms > 1) ? poly_coeffs[1] : 0.;
  for (int i=2; i<num_terms; ++i)
    t1_grad += i*poly_coeffs[i]*std::pow(x,i-1);
  return t1_grad;
}


Real NumericGenOrthogPolynomial::type1_hessian(Real x, unsigned short order)
{
  if (polyCoeffs.size() <= order)
    solve_eigenproblem(order);
  return type1_hessian(x, polyCoeffs[order]);
}


Real NumericGenOrthogPolynomial::
type1_hessian(Real x, const RealVector& poly_coeffs)
{
  // differentiate poly_coeffs with respect to x twice
  size_t num_terms = poly_coeffs.length();
  Real t1_hess = (num_terms > 2) ? poly_coeffs[2] : 0.;
  for (int i=3; i<num_terms; ++i)
    t1_hess += i*(i-1)*poly_coeffs[i]*std::pow(x,i-2);
  return t1_hess;
}


Real NumericGenOrthogPolynomial::norm_squared(unsigned short order)
{
  if (orthogPolyNormsSq.length() <= order)
    solve_eigenproblem(order);
  return orthogPolyNormsSq[order];
}


const RealArray& NumericGenOrthogPolynomial::
collocation_points(unsigned short order)
{
  // pull this out from default below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in "
	  << "NumericGenOrthogPolynomial::collocation_points()." << std::endl;
    abort_handler(-1);
  }

  if (collocPoints.size() != order) // if not already computed
    solve_eigenproblem(order);
  return collocPoints;
}


const RealArray& NumericGenOrthogPolynomial::
type1_collocation_weights(unsigned short order)
{
  // pull this out from default below since order=0 is initial colloc pts length
  if (order < 1) {
    PCerr << "Error: underflow in minimum quadrature order (1) in NumericGen"
	  << "OrthogPolynomial::type1_collocation_weights()." << std::endl;
    abort_handler(-1);
  }

  if (collocWeights.size() != order) // if not already computed
    solve_eigenproblem(order);
  return collocWeights;
}

} // namespace Pecos
