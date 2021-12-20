/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 CubatureDriver
//- Description: Implementation code for CubatureDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "CubatureDriver.hpp"
#include "sandia_cubature.hpp"
#include "SharedPolyApproxData.hpp"
#include "MarginalsCorrDistribution.hpp"
#include "NumericGenOrthogPolynomial.hpp"
//#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: CubatureDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void CubatureDriver::
initialize_grid(const MultivariateDistribution& mv_dist, unsigned short order,
		unsigned short rule)
{
  const ShortArray&  rv_types = mv_dist.random_variable_types();
  const BitArray& active_vars = mv_dist.active_variables();
  numVars = (active_vars.empty()) ? rv_types.size() : active_vars.count();

  integrand_order(order);
  collocation_rule(rule); // size collocRules and define first entry

  // check for isotropic rv_types
  if (verify_homogeneity(rv_types)) {
    PCerr << "Error: rv_types must be isotropic in CubatureDriver::"
	  << "initialize_grid(mv_dist)." << std::endl;
    abort_handler(-1);
  }

  ShortArray basis_types;
  // Cubature used for numerical integration of PCE
  // TO DO: require OPA/IPA switch? (see IntegrationDriver::initialize_grid())
  // TO DO: consider using a single BasisPolynomial for CubatureDriver (would
  // have to be expanded into array for PolynomialApproximation within NonDPCE).
  SharedPolyApproxData::initialize_polynomial_basis(basis_types, collocRules,
						    polynomialBasis);
}


void CubatureDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis)
{
  numVars         = poly_basis.size();
  polynomialBasis = poly_basis; // shallow copy

  // check for isotropic integration rules
  unsigned short rule0 = poly_basis[0].collocation_rule();
  for (size_t i=1; i<numVars; ++i)
    if (poly_basis[i].collocation_rule() != rule0) {
      PCerr << "Error: integration rule must be isotropic in CubatureDriver::"
	    << "initialize_grid(poly_basis)." << std::endl;
      abort_handler(-1);
    }
  collocation_rule(rule0);
}


void CubatureDriver::
initialize_grid_parameters(const MultivariateDistribution& mv_dist)
{
  // verify homogeneity in any polynomial parameterizations
  // (GAUSS_JACOBI, GEN_GAUSS_LAGUERRE, and GOLUB_WELSCH)
  bool err_flag = false;
  short rv_type0 = mv_dist.random_variable_type(0);
  std::shared_ptr<MarginalsCorrDistribution> mvd_rep =
    std::static_pointer_cast<MarginalsCorrDistribution>
    (mv_dist.multivar_dist_rep());
  switch (collocRules[0]) {
  case GAUSS_JACOBI: // STD_BETA: check only alpha/beta params
    err_flag =
      (verify_homogeneity(mvd_rep->pull_parameters<Real>(BETA, BE_ALPHA)) ||
       verify_homogeneity(mvd_rep->pull_parameters<Real>(BETA, BE_BETA)));
    break;
  case GEN_GAUSS_LAGUERRE: // STD_GAMMA: check only alpha params
    err_flag =
      verify_homogeneity(mvd_rep->pull_parameters<Real>(GAMMA, GA_ALPHA));
    break;
  case GOLUB_WELSCH: // numerically generated: check all params
    switch (rv_type0) { // rv_types verified in initialize_grid() above
    case BOUNDED_NORMAL:
      err_flag =
	(verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(BOUNDED_NORMAL, N_MEAN))    ||
	 verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(BOUNDED_NORMAL, N_STD_DEV)) ||
	 verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(BOUNDED_NORMAL, N_LWR_BND)) ||
	 verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(BOUNDED_NORMAL, N_UPR_BND)));
      break;
    case LOGNORMAL: // standardized on lambda, zeta
      err_flag =
	(verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(LOGNORMAL, LN_LAMBDA)) ||
	 verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(LOGNORMAL, LN_ZETA)));
      break;
    case BOUNDED_LOGNORMAL: // standardized on lambda, zeta
      err_flag =
	(verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(BOUNDED_LOGNORMAL, LN_LAMBDA))  ||
	 verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(BOUNDED_LOGNORMAL, LN_ZETA)) ||
	 verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(BOUNDED_LOGNORMAL, LN_LWR_BND)) ||
	 verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(BOUNDED_LOGNORMAL, LN_UPR_BND)));
       break;
    case LOGUNIFORM:
      err_flag =
	(verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(LOGUNIFORM, LU_LWR_BND)) ||
	 verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(LOGUNIFORM, LU_UPR_BND)));
      break;
    case TRIANGULAR:
      err_flag =
	(verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(TRIANGULAR, T_MODE)) ||
	 verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(TRIANGULAR, T_LWR_BND)) ||
	 verify_homogeneity(
	   mvd_rep->pull_parameters<Real>(TRIANGULAR, T_UPR_BND)));
      break;
    case GUMBEL:
      err_flag =
	(verify_homogeneity(mvd_rep->pull_parameters<Real>(GUMBEL, GU_ALPHA)) ||
	 verify_homogeneity(mvd_rep->pull_parameters<Real>(GUMBEL, GU_BETA)));
      break;
    case FRECHET:
      err_flag =
	(verify_homogeneity(mvd_rep->pull_parameters<Real>(FRECHET, F_ALPHA)) ||
	 verify_homogeneity(mvd_rep->pull_parameters<Real>(FRECHET, F_BETA)));
      break;
    case WEIBULL:
      err_flag =
	(verify_homogeneity(mvd_rep->pull_parameters<Real>(WEIBULL, W_ALPHA)) ||
	 verify_homogeneity(mvd_rep->pull_parameters<Real>(WEIBULL, W_BETA)));
      break;
    case HISTOGRAM_BIN:
      err_flag =
	verify_homogeneity(
	  mvd_rep->pull_parameters<RealRealMap>(HISTOGRAM_BIN, H_BIN_PAIRS));
      break;
    default: err_flag = true; break;
    }
    break;
  }

  if (err_flag) {
    PCerr << "Error: inhomogeneous distribution parameters in CubatureDriver::"
	  << "initialize_grid_parameters().\n       Consider using a variable "
	  << "transformation to standard form." << std::endl;
    abort_handler(-1);
  }

  // TO DO: consider using a single BasisPolynomial for CubatureDriver
  // (would have to be expanded into array for PolynomialApproximation
  // within NonDPCE).
  SharedPolyApproxData::
    update_basis_distribution_parameters(mv_dist, polynomialBasis);
}


int CubatureDriver::grid_size()
{
  if (numPts == 0) { // intial / special value indicates update required
    bool err_flag = false;
    switch(collocRules[0]) {
    case GAUSS_HERMITE:
      switch (integrandOrder) {
      case 1: numPts = webbur::en_her_01_1_size(numVars);    break; // 1
      case 2: numPts = webbur::en_her_02_xiu_size(numVars);  break; // n+1
    //case 3: numPts = webbur::en_her_03_1_size(numVars);    break; // 2n
      case 3: numPts = webbur::en_her_03_xiu_size(numVars);  break; // 2n
      case 5: numPts = (numVars >=2 && numVars <= 7) ?
	webbur::en_her_05_1_size(numVars) :                         // n^2+n+2
	webbur::en_her_05_2_size(numVars);                   break; // 2n^2+1
      default: err_flag = true;                              break;
      }
      break;
    case GAUSS_LEGENDRE:
      switch (integrandOrder) {
      case 1: numPts = webbur::cn_leg_01_1_size(numVars);    break; // 1
      case 2: numPts = webbur::cn_leg_02_xiu_size(numVars);  break; // n+1
    //case 3: numPts = webbur::cn_leg_03_1_size(numVars);    break; // 2n
      case 3: numPts = webbur::cn_leg_03_xiu_size(numVars);  break; // 2n
      case 5: numPts = (numVars >=4 && numVars <= 6) ?
	webbur::cn_leg_05_1_size(numVars) :                         // n^2+n+2
        webbur::cn_leg_05_2_size(numVars);                   break; // 2n^2+1
      default: err_flag = true;                              break;
      }
      break;
    case GAUSS_LAGUERRE:
      switch (integrandOrder) {
      case 1: numPts = webbur::epn_lag_01_1_size(numVars);   break;
      case 2: numPts = webbur::epn_lag_02_xiu_size(numVars); break;
      default: err_flag = true;                              break;
      }
      break;
    case GAUSS_JACOBI: {
      BasisPolynomial& poly0 = polynomialBasis[0];
      Real alpha_poly; poly0.pull_parameter(JACOBI_ALPHA, alpha_poly);
      Real  beta_poly; poly0.pull_parameter(JACOBI_BETA,   beta_poly);
      switch (integrandOrder) {
      case 1: numPts = webbur::cn_jac_01_1_size(numVars, alpha_poly, beta_poly);
	break;
      case 2: numPts = webbur::cn_jac_02_xiu_size(numVars,alpha_poly,beta_poly);
	break;
      default: err_flag = true; break;
      }
      break;
    }
    case GEN_GAUSS_LAGUERRE: {
      Real alpha_poly;
      polynomialBasis[0].pull_parameter(GENLAG_ALPHA, alpha_poly);
      switch (integrandOrder) {
      case 1: numPts = webbur::epn_glg_01_1_size(numVars,   alpha_poly); break;
      case 2: numPts = webbur::epn_glg_02_xiu_size(numVars, alpha_poly); break;
      default: err_flag = true;                                          break;
      }
      break;
    }
    case GOLUB_WELSCH:
      switch (integrandOrder) {
      case 2: numPts = webbur::gw_02_xiu_size(numVars); break;
      default: err_flag = true;                         break;
      }
      break;
    default:
      err_flag = true; break;
    }

    if (err_flag) {
      PCerr << "Error: unsupported rule in CubatureDriver::grid_size()."
	    << std::endl;
      abort_handler(-1);
    }
  }
  return numPts;
}


void CubatureDriver::compute_grid()
{
  // --------------------------------
  // Get number of collocation points
  // --------------------------------
  grid_size(); // ensure numPts is up to date
#ifdef DEBUG
  PCout << "Total number of cubature integration points: " << numPts << '\n';
#endif // DEBUG

  // ----------------------------------------------
  // Get collocation points and integration weights
  // ----------------------------------------------
  type1WeightSets.sizeUninitialized(numPts);
  variableSets.shapeUninitialized(numVars, numPts); // Teuchos: col major
  double *pts = variableSets.values(), *wts = type1WeightSets.values();
  bool err_flag = false, pt_scaling = false, wt_scaling = false;
  double pt_factor, wt_factor;
  BasisPolynomial& poly0 = polynomialBasis[0];
  switch(collocRules[0]) {
  case GAUSS_HERMITE: {
    switch (integrandOrder) {
    case 1: webbur::en_her_01_1(numVars,    numPts, pts, wts); break;
    case 2: webbur::en_her_02_xiu(numVars,  numPts, pts, wts); break;
  //case 3: webbur::en_her_03_1(numVars,    numPts, pts, wts); break;
    case 3: webbur::en_her_03_xiu(numVars,  numPts, pts, wts); break;
    case 5:
      if (numVars >=2 && numVars <= 7) {
	int option = 1; // two options for n=3,5,6
	webbur::en_her_05_1(numVars, option, numPts, pts, wts);
      }
      else
	webbur::en_her_05_2(numVars, numPts, pts, wts);        break;
    default: err_flag = true;                                  break;
    }
    pt_scaling = true; wt_scaling = true; break;
  }
  case GAUSS_LEGENDRE: {
    switch (integrandOrder) {
    case 1: webbur::cn_leg_01_1(numVars,    numPts, pts, wts); break;
    case 2: webbur::cn_leg_02_xiu(numVars,  numPts, pts, wts); break;
  //case 3: webbur::cn_leg_03_1(numVars,    numPts, pts, wts); break;
    case 3: webbur::cn_leg_03_xiu(numVars,  numPts, pts, wts); break;
    case 5:
      if (numVars >=4 && numVars <= 6) {
	int option = 1; // two options for n=5,6
	webbur::cn_leg_05_1(numVars, option, numPts, pts, wts);
      }
      else
	webbur::cn_leg_05_2(numVars, numPts, pts, wts);        break;
    default: err_flag = true;                                  break;
    }
    wt_scaling = true; break;
  }
  case GAUSS_LAGUERRE:
    switch (integrandOrder) {
    case 1: webbur::epn_lag_01_1(numVars,   numPts, pts, wts); break;
    case 2: webbur::epn_lag_02_xiu(numVars, numPts, pts, wts); break;
    default: err_flag = true;                                  break;
    } break;
  case GAUSS_JACOBI: {
    Real alpha_poly; poly0.pull_parameter(JACOBI_ALPHA, alpha_poly);
    Real  beta_poly; poly0.pull_parameter(JACOBI_BETA,   beta_poly);
    switch (integrandOrder) {
    case 1:
      webbur::cn_jac_01_1(numVars,   alpha_poly, beta_poly, numPts, pts, wts);
      break;
    case 2:
      webbur::cn_jac_02_xiu(numVars, alpha_poly, beta_poly, numPts, pts, wts);
      break;
    default: err_flag = true; break;
    }
    wt_scaling = true;        break;
  }
  case GEN_GAUSS_LAGUERRE: {
    Real alpha_poly; poly0.pull_parameter(GENLAG_ALPHA, alpha_poly);
    switch (integrandOrder) {
    case 1: webbur::epn_glg_01_1(numVars,   alpha_poly, numPts, pts, wts);
      break;
    case 2: webbur::epn_glg_02_xiu(numVars, alpha_poly, numPts, pts, wts);
      break;
    default: err_flag = true; break;
    }
    wt_scaling = true;        break;
  }
  case GOLUB_WELSCH:
    switch (integrandOrder) {
    case 2: {
      /*
      sandia_cubature (and "Numerical integration formulas of degree 2",
      Appl. Num. Math. (58), D. Xiu, 2008):
      ----------------
      x P(n,x) = An * P(n+1,x) + Bn * P(n,x) + Cn * P(n-1,x)  [from source]
      --> P(n+1,x) = x/An P(n,x) - Bn/An * P(n,x) - Cn/An * P(n-1,x)
      --> P(  1,x) = (x-B0)/A0 = GAMMA0x + DELTA0    [P(0,x) = 1, P(-1,x) = 0]
      --> P(  1,x) = GAMMA0x + DELTA0 --> GAMMA0 = 1./A0, DELTA0 = -B0/A0

      NumericGenOrthogPolynomial.hpp:
      -------------------------------
      poly_ip1 = x poly_i - alpha_i * poly_i - beta_i * poly_im1
      --> An = 1, alpha_i = Bn/An = Bn, beta_i = Cn/An = Cn
      --> consistent with monic polynomials (leading coefficient = 1)
      --> GAMMA0 = 1., DELTA0 = -alpha_0, C1 = beta_1
      */
      std::shared_ptr<NumericGenOrthogPolynomial> poly0_rep =
	std::static_pointer_cast<NumericGenOrthogPolynomial>
	(poly0.polynomial_rep());
      const Real& beta1  = poly0_rep->beta_recursion(1); // do order 1 first
      const Real& alpha0 = poly0_rep->alpha_recursion(0);
      webbur::gw_02_xiu(numVars, numPts, 1., -alpha0, beta1, 1., pts, wts);
      break;
    }
    default: err_flag = true; break;
    } break;
  default:
    err_flag = true; break;
  }

  if (err_flag) {
    PCerr << "Error: unsupported rule in CubatureDriver::compute_grid()."
	  << std::endl;
    abort_handler(-1);
  }

  // scale points and weights
  if (pt_scaling)
    variableSets.scale(poly0.point_factor());
  if (wt_scaling)
    type1WeightSets.scale(std::pow(poly0.weight_factor(), (int)numVars));
}

} // namespace Pecos
