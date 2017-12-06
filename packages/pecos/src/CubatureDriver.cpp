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
#include "DistributionParams.hpp"
#include "NumericGenOrthogPolynomial.hpp"
//#include "pecos_stat_util.hpp"

static const char rcsId[]="@(#) $Id: CubatureDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


void CubatureDriver::
initialize_grid(const ShortArray& u_types, unsigned short order,
		unsigned short rule)
{
  numVars = u_types.size();
  integrand_order(order);
  collocation_rule(rule); // size collocRules and define first entry

  // check for isotropic u_types
  short type0 = u_types[0];
  for (size_t i=1; i<numVars; ++i)
    if (u_types[i] != type0) {
      PCerr << "Error: u_types must be isotropic in CubatureDriver::"
	    << "initialize_grid(u_types)." << std::endl;
      abort_handler(-1);
    }

  ShortArray basis_types;
  // Cubature used for numerical integration of PCE
  // TO DO: require OPA/IPA switch? (see IntegrationDriver::initialize_grid())
  //BasisConfigOptions bc_options(false, false, false, false);
  //SharedPolyApproxData::initialize_basis_types(u_types, bc_options,
  //                                             basis_types);
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

  // check for isotropic u_types
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
initialize_grid_parameters(const ShortArray& u_types,
			   const AleatoryDistParams& adp)
{
  // verify homogeneity in any polynomial parameterizations
  // (GAUSS_JACOBI, GEN_GAUSS_LAGUERRE, and GOLUB_WELSCH)
  bool err_flag = false;
  switch (collocRules[0]) {
  case GAUSS_JACOBI: // STD_BETA: check only alpha/beta params
    err_flag = (verify_homogeneity(adp.beta_alphas()) ||
		verify_homogeneity(adp.beta_betas())); break;
  case GEN_GAUSS_LAGUERRE: // STD_GAMMA: check only alpha params
    err_flag = verify_homogeneity(adp.gamma_alphas()); break;
  case GOLUB_WELSCH: // numerically generated: check all params
    switch (u_types[0]) { // u_types verified in initialize_grid() above
    case BOUNDED_NORMAL:
      err_flag = (verify_homogeneity(adp.normal_means()) ||
		  verify_homogeneity(adp.normal_std_deviations()) ||
		  verify_homogeneity(adp.normal_lower_bounds()) ||
		  verify_homogeneity(adp.normal_upper_bounds())); break;
    case LOGNORMAL:
      err_flag = (verify_homogeneity(adp.lognormal_means()) ||
		  verify_homogeneity(adp.lognormal_std_deviations()) ||
		  verify_homogeneity(adp.lognormal_lambdas()) ||
		  verify_homogeneity(adp.lognormal_zetas()) ||
		  verify_homogeneity(adp.lognormal_error_factors())); break;
    case BOUNDED_LOGNORMAL:
      err_flag = (verify_homogeneity(adp.lognormal_means()) ||
		  verify_homogeneity(adp.lognormal_std_deviations()) ||
		  verify_homogeneity(adp.lognormal_lambdas()) ||
		  verify_homogeneity(adp.lognormal_zetas()) ||
		  verify_homogeneity(adp.lognormal_error_factors()) ||
		  verify_homogeneity(adp.lognormal_lower_bounds()) ||
		  verify_homogeneity(adp.lognormal_upper_bounds())); break;
    case LOGUNIFORM:
      err_flag = (verify_homogeneity(adp.loguniform_lower_bounds()) ||
		  verify_homogeneity(adp.loguniform_upper_bounds())); break;
    case TRIANGULAR:
      err_flag = (verify_homogeneity(adp.triangular_modes()) ||
		  verify_homogeneity(adp.triangular_lower_bounds()) ||
		  verify_homogeneity(adp.triangular_upper_bounds())); break;
    case GUMBEL:
      err_flag = (verify_homogeneity(adp.gumbel_alphas()) ||
		  verify_homogeneity(adp.gumbel_betas()));  break;
    case FRECHET:
      err_flag = (verify_homogeneity(adp.frechet_alphas()) ||
		  verify_homogeneity(adp.frechet_betas())); break;
    case WEIBULL:
      err_flag = (verify_homogeneity(adp.weibull_alphas()) ||
		  verify_homogeneity(adp.weibull_betas())); break;
    case HISTOGRAM_BIN:
      err_flag = verify_homogeneity(adp.histogram_bin_pairs()); break;
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
  SharedPolyApproxData::update_basis_distribution_parameters(u_types, adp,
							     polynomialBasis);
}


int CubatureDriver::grid_size()
{
  if (updateGridSize) {
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
      const Real& alpha_poly = poly0.alpha_polynomial();
      const Real& beta_poly  = poly0.beta_polynomial();
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
      const Real& alpha_poly = polynomialBasis[0].alpha_polynomial();
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
    else
      updateGridSize = false;
  }
  return numPts;
}


void CubatureDriver::compute_grid(RealMatrix& variable_sets)
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
  variable_sets.shapeUninitialized(numVars, numPts); // Teuchos: col major
  double *pts = variable_sets.values(), *wts = type1WeightSets.values();
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
    const Real& alpha_poly = poly0.alpha_polynomial();
    const Real& beta_poly  = poly0.beta_polynomial();
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
    const Real& alpha_poly = poly0.alpha_polynomial();
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
      NumericGenOrthogPolynomial* poly0_rep
	= (NumericGenOrthogPolynomial*)poly0.polynomial_rep();
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
    variable_sets.scale(poly0.point_factor());
  if (wt_scaling)
    type1WeightSets.scale(std::pow(poly0.weight_factor(), (int)numVars));
}

} // namespace Pecos
