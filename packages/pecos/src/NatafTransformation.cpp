/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "NatafTransformation.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
//#include "pecos_stat_util.hpp"
#include "NormalRandomVariable.hpp"
#include "UniformRandomVariable.hpp"

static const char rcsId[]="@(#) $Id: NatafTransformation.cpp 4768 2007-12-17 17:49:32Z mseldre $";

//#define DEBUG


namespace Pecos {


/** This procedure performs the transformation from u to x space.
    u_vars is the vector of random variables in uncorrelated standard
    normal space (u-space).  x_vars is the vector of random variables
    in the original user-defined x-space. */
void NatafTransformation::
trans_U_to_X(const RealVector& u_vars, RealVector& x_vars)
{
  if (xDist.correlation()) {
    RealVector z_vars;
    trans_U_to_Z(u_vars, z_vars);
    trans_Z_to_X(z_vars, x_vars);
  }
  else // z_vars = u_vars
    trans_Z_to_X(u_vars, x_vars);

#ifdef DEBUG
  PCout << "Transforming u_vars:\n" << u_vars << "to x_vars:\n" << x_vars;
#endif
}


/** This procedure computes the transformation from u to z space.
    u_vars is the vector of random variables in uncorrelated standard
    normal space (u-space).  z_vars is the vector of random variables
    in normal space with proper correlations (z-space). */
void NatafTransformation::
trans_U_to_Z(const RealVector& u_vars, RealVector& z_vars)
{
  // corrCholeskyFactorZ: the Cholesky factor of the modified correlation matrix

  int u_len = u_vars.length();
  if (z_vars.length() != u_len)
    z_vars.size(u_len);
  z_vars.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1., corrCholeskyFactorZ,
		  u_vars, 0.);
}


/** This procedure computes the transformation from z to x space.
    z_vars is the vector of random variables in normal space with
    proper correlations (z-space).  x_vars is the vector of random
    variables in the original user-defined x-space */
void NatafTransformation::
trans_Z_to_X(const RealVector& z_vars, RealVector& x_vars)
{
  int z_len = z_vars.length();
  if (x_vars.length() != z_len)
    x_vars.size(z_len);

  for (size_t i=0; i<z_len; ++i) {
    trans_Z_to_X(z_vars[i], x_vars[i], i);

#ifdef DEBUG
    PCout << "Z_to_X: z[" << i << "] = " << z_vars[i]
	  <<    " --> x[" << i << "] = " << x_vars[i] << std::endl;
#endif // DEBUG
  }
}


/** This procedure computes the transformation from z to x space.
    z_var is the random variable in standardized space with proper
    correlations (z-space).  x_var is the random variable in the
    original user-defined x-space */
void NatafTransformation::trans_Z_to_X(Real z, Real& x, size_t i)
{
  // This routine performs an inverse transformation based on CDF/CCDF
  // equivalence, e.g. F(X) = Phi(Z) in the case of a std normal z-space CDFs

  const RandomVariable& x_rv_i = xDist.active_random_variable(i);
  short x_type = x_rv_i.type(), u_type = uDist.active_random_variable_type(i);
  if (u_type == x_type)
    x = z;
  else if (u_type == STD_NORMAL) {
    switch (x_type) {
    case NORMAL:  x = x_rv_i.from_standard(z);  break;
    case LOGNORMAL: {
      Real lambda;  x_rv_i.pull_parameter(LN_LAMBDA, lambda);
      Real   zeta;  x_rv_i.pull_parameter(LN_ZETA,   zeta);
      x = std::exp(lambda + zeta * z);
      break;
    }
    /* log cdf offers no real benefit for normal target due to erf() cdf:
    case EXPONENTIAL: // Phi(z) = F(x) = 1 - e^(-x/beta)
    case WEIBULL:     // Phi(z) = F(x) = 1 - e^(-(x/beta)^alpha)
      x = x_rv_i.inverse_log_ccdf(NormalRandomVariable::log_std_ccdf(z)); break;
    case GUMBEL:  // Phi(z) = F(x) = e^(-e^(-alpha(x-beta)))
    case FRECHET: // Phi(z) = F(x) = e^(-(beta/x)^alpha)
      x = x_rv_i.inverse_log_cdf(NormalRandomVariable::log_std_cdf(z));   break;
    */
    default: // default mapping based on either CDF or CCDF equivalence
      x = (z > 0.) ? x_rv_i.inverse_ccdf(NormalRandomVariable::std_ccdf(z)) :
	x_rv_i.inverse_cdf(NormalRandomVariable::std_cdf(z));  break;
    }
  }
  else if (u_type == STD_UNIFORM)
    x = (z > 0.) ? x_rv_i.inverse_ccdf(UniformRandomVariable::std_ccdf(z)) :
      x_rv_i.inverse_cdf(UniformRandomVariable::std_cdf(z));
  else if ( (u_type == STD_EXPONENTIAL && x_type == EXPONENTIAL) ||
	    (u_type == STD_GAMMA       && x_type == GAMMA) ||
	    (u_type == STD_BETA        && x_type == BETA) )
    x = x_rv_i.from_standard(z);
  else {
    PCerr << "Error: unsupported variable mapping for variable " << i
	  << " in NatafTransformation::trans_Z_to_X()" << std::endl;
    abort_handler(-1);
  }
}


/** This procedure performs the transformation from x to u space
    u_vars is the vector of random variables in uncorrelated standard
    normal space (u-space).  x_vars is the vector of random variables
    in the original user-defined x-space. */
void NatafTransformation::
trans_X_to_U(const RealVector& x_vars, RealVector& u_vars) 
{ 
  if (xDist.correlation()) {
    RealVector z_vars;
    trans_X_to_Z(x_vars, z_vars);
    trans_Z_to_U(z_vars, u_vars);
  }
  else // z_vars = u_vars
    trans_X_to_Z(x_vars, u_vars);

#ifdef DEBUG
  PCout << "Transforming x_vars:\n" << x_vars << "to u_vars:\n" << u_vars;
#endif
}


/** This procedure performs the transformation from x to z space:
    z_vars is the vector of random variables in normal space with
    proper correlations (z-space).  x_vars is the vector of random
    variables in the original user-defined x-space. */
void NatafTransformation::
trans_X_to_Z(const RealVector& x_vars, RealVector& z_vars)
{
  int x_len = x_vars.length();
  if (z_vars.length() != x_len)
    z_vars.size(x_len);

  for (size_t i=0; i<x_len; ++i) {
    trans_X_to_Z(x_vars[i], z_vars[i], i);

#ifdef DEBUG
    PCout << "X_to_Z: x[" << i << "] = " << x_vars[i]
	  <<    " --> z[" << i << "] = " << z_vars[i] << std::endl;
#endif // DEBUG
  }
}


/** This procedure performs the transformation from x to z space:
    z_var is the random variable in standardized space with proper
    correlations (z-space).  x_var is the random variable in the
    original user-defined x-space. */
void NatafTransformation::trans_X_to_Z(Real x, Real& z, size_t i)
{
  // This routine performs an forward transformation based on CDF/CCDF
  // equivalence, e.g. F(X) = Phi(Z) in the case of a std normal z-space CDFs

  const RandomVariable& x_rv_i = xDist.active_random_variable(i);
  short x_type = x_rv_i.type(), u_type = uDist.active_random_variable_type(i);
  if (u_type == x_type)
    z = x;
  else if (u_type == STD_NORMAL) {
    switch (x_type) {
    case NORMAL:    z = x_rv_i.to_standard(x); break;
    case LOGNORMAL: {
      Real lambda;  x_rv_i.pull_parameter(LN_LAMBDA, lambda);
      Real   zeta;  x_rv_i.pull_parameter(LN_ZETA,   zeta);
      z = (std::log(x) - lambda) / zeta;
      break;
    }
    /* log cdf offers no real benefit for normal target due to erf() cdf:
    case EXPONENTIAL: // Phi(z) = F(x) = 1 - e^(-x/beta)
    case WEIBULL:     // Phi(z) = F(x) = 1 - e^(-(x/beta)^alpha)
      z = NormalRandomVariable::inverse_log_std_ccdf(x_rv_i.log_ccdf(x)); break;
    case GUMBEL:  // Phi(z) = F(x) = e^(-e^(-alpha(x-beta)))
    case FRECHET: // Phi(z) = F(x) = e^(-(beta/x)^alpha)
      z = NormalRandomVariable::inverse_log_std_cdf(x_rv_i.log_cdf(x));   break;
    */
    default: { // default mapping based on either CDF or CCDF equivalence
      Real xcdf = x_rv_i.cdf(x);
      z = (xcdf > .5) ? NormalRandomVariable::inverse_std_ccdf(x_rv_i.ccdf(x)) :
	NormalRandomVariable::inverse_std_cdf(xcdf);
      break;
    }
    }
  }
  else if (u_type == STD_UNIFORM) {
    Real xcdf = x_rv_i.cdf(x);
    z = (xcdf > .5) ? UniformRandomVariable::inverse_std_ccdf(x_rv_i.ccdf(x)) :
      UniformRandomVariable::inverse_std_cdf(xcdf);
  }
  else if ( (u_type == STD_EXPONENTIAL && x_type == EXPONENTIAL) ||
	    (u_type == STD_GAMMA       && x_type == GAMMA) ||
	    (u_type == STD_BETA        && x_type == BETA) )
    z = x_rv_i.to_standard(x);
  else {
    PCerr << "Error: unsupported variable mapping for variable " << i
	  << " in NatafTransformation::trans_X_to_Z()" << std::endl;
    abort_handler(-1);
  }
}


/** This procedure computes the transformation from z to u space.
    u_vars is the vector of random variables in uncorrelated standard
    normal space (u-space).  z_vars is the vector of random variables
    in normal space with proper correlations (z-space). */
void NatafTransformation::trans_Z_to_U(RealVector& z_vars, RealVector& u_vars)
{
  // corrCholeskyFactorZ: the Cholesky factor of the modified correlation matrix

  int z_len = z_vars.length();
  // The inbound u_vars object may be a Teuchos::View we want to
  // update with the solution.  However SerialDenseSolver::solve(),
  // line 647 will call (LHS_ = u_vars) = RHS_, disconnecting the
  // view.  To work around, we use a temporary.
  RealVector tmp_u_vars(z_len);
  // WARNING: this could also disconnect a view...
  if (u_vars.length() != z_len)
    u_vars.size(z_len);


  RealSolver corr_solver;
  corr_solver.setMatrix(  Teuchos::rcp(&corrCholeskyFactorZ, false) );
  corr_solver.setVectors( Teuchos::rcp(&tmp_u_vars, false),
			  Teuchos::rcp(&z_vars, false) );
  corr_solver.solveToRefinedSolution(true);
  corr_solver.solve();
  // Assign into u_vars, which may be a Teuchos::View
  u_vars.assign(tmp_u_vars);
}


/** This procedure modifies the correlation matrix input by the user for use in
    the Nataf distribution model (Der Kiureghian and Liu, ASCE JEM 112:1, 1986).
    It uses empirical expressionss derived from least-squares polynomial fits
    to numerical integration data.

    \li x_corr_matrix: the correlation coefficient matrix of the random
    variables provided by the user

    \li mod_corr_matrix: modified correlation matrix

    \li corrCholeskyFactorZ: Cholesky factor of the modified correlation matrix
    for use in Z_to_U and U_to_Z transformations.

    \li cf_var_{i,j}: coefficient of variation of Xi,Xj

    Note: The modification is exact for normal-normal, lognormal-lognormal, and
    normal-lognormal tranformations.  All other cases are approximations with
    some error as noted below. */
void NatafTransformation::transform_correlations()
{
  if (!xDist.correlation())
    return;

  const RealSymMatrix& x_corr_matrix = xDist.correlation_matrix();
  const BitArray&        active_vars = xDist.active_variables();
  const BitArray&        active_corr = xDist.active_correlations();
  bool no_v_mask = active_vars.empty(), no_c_mask = active_corr.empty();

  // Enumerate the active correlations (ignoring the active variable subset)
  const std::vector<RandomVariable>& x_rv = xDist.random_variables();
  const ShortArray&               u_types = uDist.random_variable_types();
  size_t rv_i, rv_j, c_i, c_j, v_i, v_j, num_rv = x_rv.size(),
    num_active_v = (no_v_mask) ? num_rv : active_vars.count(),
    num_active_c = (no_c_mask) ? num_rv : active_corr.count();

  // Loop over active variables using find_{first,next}:
  // RealSymMatrix mod_corr_matrix(x_corr_matrix); // copy
  // rv_i = rv_j = (no_c_mask) ? 0 : active_corr.find_first();
  // for (c_i=0; i<num_active_c; ++c_i) {
  //   //if (u_types[v_i] != STD_NORMAL) continue;
  //   if (c_i) rv_j = (no_c_mask) ? 0 : active_corr.find_first();
  //   for (c_j=0; c_j<c_i; ++c_j) {
  //     Real corr_ij = x_corr_matrix(c_i, c_j);
  //     if (/* u_types[v_j] == STD_NORMAL && */ std::abs(corr_ij) > 0.)
  // 	mod_corr_matrix(c_i, c_j) *=
  // 	  x_rv[rv_i].correlation_warping_factor(x_rv[rv_j], corr_ij);
  //     rv_j = (no_c_mask) ? c_j : active_corr.find_next(rv_j);
  //   }
  //   rv_i = (no_c_mask) ? c_i : active_corr.find_next(rv_i);
  // }

  // Loop over all random variables, checking active_corr bits:
  // RealSymMatrix mod_corr_matrix(x_corr_matrix); // copy
  // for (rv_i=0, c_i=0; rv_i<num_rv; ++rv_i)
  //   if (no_c_mask || active_corr[rv_i]) {
  //     //if (u_types[rv_i] == STD_NORMAL)
  // 	 for (rv_j=0, c_j=0; rv_j<rv_i; ++rv_j)
  // 	   if (no_c_mask || active_corr[rv_j]) {
  // 	     Real corr = x_corr_matrix(c_i, c_j);
  // 	     if (/* u_types[j] == STD_NORMAL && */ std::abs(corr) > 0.)
  // 	       mod_corr_matrix(c_i, c_j) *=
  // 		 x_rv[rv_i].correlation_warping_factor(x_rv[rv_j], corr);
  // 	     ++c_j;
  // 	   }
  //     ++c_i;
  //   }

  // Loop over all random variables, checking active_{vars,corr} bits:
  RealSymMatrix mod_corr_matrix(num_active_v); // init to 0
  bool active_c_i, active_v_i, active_c_j, active_v_j;
  for (rv_i=0, c_i=0, v_i=0; rv_i<num_rv; ++rv_i) {
    active_c_i = (no_c_mask || active_corr[rv_i]);
    active_v_i = (no_v_mask || active_vars[rv_i]);
    if (active_v_i) mod_corr_matrix(v_i, v_i) = 1.;
    if (active_c_i && active_v_i) {
      for (rv_j=0, c_j=0, v_j=0; rv_j<rv_i; ++rv_j) {
	active_c_j = (no_c_mask || active_corr[rv_j]);
        active_v_j = (no_v_mask || active_vars[rv_j]);
	if (active_c_j && active_v_j) {
	  Real corr = x_corr_matrix(c_i, c_j);
	  if (std::abs(corr) > 0.)
	    mod_corr_matrix(v_i, v_j) = corr *
	      x_rv[rv_i].correlation_warping_factor(x_rv[rv_j], corr);
	}
	if (active_c_j) ++c_j;
	if (active_v_j) ++v_j;
      }
    }
    if (active_c_i) ++c_i;
    if (active_v_i) ++v_i;
  }

  // Cholesky decomposition for modified correlation matrix
  RealSpdSolver corr_solver;
  corr_solver.setMatrix( Teuchos::rcp(&mod_corr_matrix, false) );
  corr_solver.factor(); // Cholesky factorization (LL^T) in place
  // Define corrCholeskyFactorZ to be L by assigning the lower triangle.
  // Inflate as needed if discrepancy between active_rv and active_corr
  if (corrCholeskyFactorZ.numRows() != num_active_v ||
      corrCholeskyFactorZ.numCols() != num_active_v)
    corrCholeskyFactorZ.shape(num_active_v, num_active_v); // init to 0

  // If mod_corr_matrix is kept in active_corr, need to map to active_vars
  //if (active_vars == active_corr)
  for (v_i=0; v_i<num_active_v; ++v_i)
    for (v_j=0; v_j<=v_i; ++v_j)
      corrCholeskyFactorZ(v_i, v_j) = mod_corr_matrix(v_i, v_j);
  //else { ...map to active_vars... }

#ifdef DEBUG
  PCout << "corrCholeskyFactorZ:\n" << corrCholeskyFactorZ;
#endif

  // could pre-compute L^-1 to avoid solving L u = z repeatedly for u
  //corrCholeskyFactorZInv.shape(num_active_vars, num_active_vars);
  //corrCholeskyFactorZInv = corrCholeskyFactorZ; // copy
  //RealSolver chol_solver;
  //chol_solver.setMatrix(corrCholeskyFactorZInv); 
  //chol_solver.invert();
  //PCout << "\ncorrCholeskyFactorZInv:" << corrCholeskyFactorZInv;
}


/** This procedure tranforms a gradient vector dg/dx from the original
    user-defined x-space (where evaluations are performed) to uncorrelated
    standard normal space (u-space) through application of the Jacobian dx/du.
    x_vars is the vector of random variables in x-space. */
void NatafTransformation::
trans_grad_X_to_U(const RealVector& fn_grad_x, RealVector& fn_grad_u,
		  const RealVector& x_vars, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  RealMatrix jacobian_xu;
  jacobian_dX_dU(x_vars, jacobian_xu);
  trans_grad_X_to_U(fn_grad_x, fn_grad_u, jacobian_xu, x_dvv, cv_ids);
}


/** This procedure tranforms a gradient vector dg/dx from the original
    user-defined x-space (where evaluations are performed) to uncorrelated
    standard normal space (u-space) through application of the Jacobian dx/du.
    This overloaded form allows for the separate calculation of jacobian_xu,
    as this matrix is independent of the response function index and can be
    pulled outside response function loops. */
void NatafTransformation::
trans_grad_X_to_U(const RealVector& fn_grad_x, RealVector& fn_grad_u,
		  const RealMatrix& jacobian_xu, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  // Jacobian dimensions = length of random variables = model.cv()
  int num_v = jacobian_xu.numRows();
  if (x_dvv == cv_ids) { // standard DVV
    if (fn_grad_x.length() != num_v) {
      PCerr << "Error: bad fn_grad_x dimension in NatafTransformation::"
	    << "trans_grad_X_to_U()." << std::endl;
      abort_handler(-1);
    }
    if (fn_grad_u.length() != num_v)
      fn_grad_u.size(num_v);
    fn_grad_u.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., jacobian_xu,
		       fn_grad_x, 0.);
  }
  else { // non-standard DVV
    RealVector fn_grad_x_trans(num_v), fn_grad_u_trans(num_v, false);
    size_t i, dvv_index, num_deriv_vars = x_dvv.size();
    SizetArray dvv_index_array(num_v);
    // extract relevant DVV components from fn_grad_x
    for (i=0; i<num_v; ++i) {
      dvv_index_array[i] = dvv_index = find_index(x_dvv, cv_ids[i]);
      if (dvv_index != _NPOS)
	fn_grad_x_trans[i] = fn_grad_x(dvv_index);
    }
    // perform transformation using full Jacobian
    fn_grad_u_trans.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., jacobian_xu,
			     fn_grad_x_trans, 0.);
    // copy relevant DVV components into fn_grad_u
    if (fn_grad_u.length() != num_deriv_vars)
      fn_grad_u.size(num_deriv_vars);
    for (i=0; i<num_v; ++i) {
      dvv_index = dvv_index_array[i];
      if (dvv_index != _NPOS)
	fn_grad_u(dvv_index) = fn_grad_u_trans[i];
    }
  }

#ifdef DEBUG
  PCout << "Transformed fn_grad_x:\n" << fn_grad_x
        << "to fn_grad_u:\n" << fn_grad_u;
#endif
}


/** This procedure tranforms a gradient vector dg/du from uncorrelated standard
    space (u-space) to the original user-defined x-space through application of
    the Jacobian du/dx.  x_vars is the vector of random variables in x-space. */
void NatafTransformation::
trans_grad_U_to_X(const RealVector& fn_grad_u, RealVector& fn_grad_x,
		  const RealVector& x_vars, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  RealMatrix jacobian_ux;
  jacobian_dU_dX(x_vars, jacobian_ux);
  trans_grad_U_to_X(fn_grad_u, fn_grad_x, jacobian_ux, x_dvv, cv_ids);
}


/** This procedure tranforms a gradient vector dg/du from uncorrelated standard
    space (u-space) to the original user-defined x-space through application of
    the Jacobian du/dx.  This overloaded form allows for the separate
    calculation of jacobian_ux, as this matrix is independent of the response
    function index and can be pulled outside response function loops. */
void NatafTransformation::
trans_grad_U_to_X(const RealVector& fn_grad_u, RealVector& fn_grad_x,
		  const RealMatrix& jacobian_ux, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  // Jacobian dimensions = length of random variables = model.cv()
  int u_len = jacobian_ux.numRows();
  if (x_dvv == cv_ids) { // standard DVV
    if (fn_grad_u.length() != u_len) {
      PCerr << "Error: bad fn_grad_u dimension in NatafTransformation::"
	    << "trans_grad_U_to_X()." << std::endl;
      abort_handler(-1);
    }
    if (fn_grad_x.length() != u_len)
      fn_grad_x.size(u_len);
    fn_grad_x.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., jacobian_ux,
		       fn_grad_u, 0.);
  }
  else { // non-standard DVV
    RealVector fn_grad_u_trans(u_len), fn_grad_x_trans(u_len, false);
    size_t dvv_index, num_deriv_vars = x_dvv.size();
    SizetArray dvv_index_array(u_len);
    // extract relevant DVV components from fn_grad_u
    for (int i=0; i<u_len; ++i) {
      dvv_index_array[i] = dvv_index = find_index(x_dvv, cv_ids[i]);
      if (dvv_index != _NPOS)
	fn_grad_u_trans[i] = fn_grad_u(dvv_index);
    }
    // perform transformation using full Jacobian
    fn_grad_x_trans.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., jacobian_ux,
			     fn_grad_u_trans, 0.);
    // copy relevant DVV components into fn_grad_x
    if (fn_grad_x.length() != num_deriv_vars)
      fn_grad_x.size(num_deriv_vars);
    for (int i=0; i<u_len; ++i) {
      dvv_index = dvv_index_array[i];
      if (dvv_index != _NPOS)
	fn_grad_x(dvv_index) = fn_grad_x_trans[i];
    }
  }

#ifdef DEBUG
  PCout << "Transformed fn_grad_u:\n" << fn_grad_u
        << "to fn_grad_x:\n" << fn_grad_x;
#endif
}


/** This procedure multiplies a gradient vector dg/dx from the
    original user-defined x-space (where evaluations are performed)
    with the design Jacobian dx/ds of the transformation x = x(u,s) to
    form the design gradient dg/ds.  x_vars is the vector of random
    variables in x-space. */
void NatafTransformation::
trans_grad_X_to_S(const RealVector& fn_grad_x, RealVector& fn_grad_s,
		  const RealVector& x_vars, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids,
		  SizetMultiArrayConstView acv_ids,
		  const SizetArray& acv_map1_indices,
		  const ShortArray& acv_map2_targets)
{
  RealMatrix jacobian_xs;
  jacobian_dX_dS(x_vars, jacobian_xs, cv_ids, acv_ids,
                 acv_map1_indices, acv_map2_targets);
  trans_grad_X_to_S(fn_grad_x, fn_grad_s, jacobian_xs, x_dvv, cv_ids, acv_ids,
		    acv_map1_indices, acv_map2_targets);
}


/** This procedure multiplies a gradient vector dg/dx from the
    original user-defined x-space (where evaluations are performed)
    with the design Jacobian dx/ds of the transformation x = x(u,s) to
    form the design gradient dg/ds.  This overloaded form allows
    for the separate calculation of jacobian_xs, as this matrix is
    independent of the response function index and can be pulled
    outside response function loops. */
void NatafTransformation::
trans_grad_X_to_S(const RealVector& fn_grad_x, RealVector& fn_grad_s,
		  const RealMatrix& jacobian_xs, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids,
		  SizetMultiArrayConstView acv_ids,
		  const SizetArray& acv_map1_indices,
		  const ShortArray& acv_map2_targets)
{
  // Jacobian dim is num_v by num_s, where
  // > num_v = length of random variables  = inner model.cv()
  // > num_s = acv_map1_indices.size() = outer model.cv()
  int num_v = jacobian_xs.numRows(), num_s = jacobian_xs.numCols();
  if (acv_map1_indices.empty() || acv_map2_targets.empty()) {
    PCerr << "Error: NatafTransformation::trans_grad_X_to_S() requires primary "
	  << "and secondary variable mappings to define S." << std::endl;
    abort_handler(-1);
  }

  bool std_dvv = (x_dvv == cv_ids);
  bool mixed_s = std::find(acv_map2_targets.begin(), acv_map2_targets.end(),
                 (short)NO_TARGET) != acv_map2_targets.end() ? true : false;
  RealVector fn_grad_x_std, fn_grad_s_std;

  // manage size of fn_grad_x input
  if (std_dvv) { // standard x-space DVV
    if (fn_grad_x.length() != num_v) {
      PCerr << "Error: bad fn_grad_x dimension in NatafTransformation::"
	    << "trans_grad_X_to_S()." << std::endl;
      abort_handler(-1);
    }
  }
  else { // non-standard x-space DVV
    fn_grad_x_std.size(num_v);
    // extract relevant DVV components from fn_grad_x
    size_t i, dvv_index;
    for (i=0; i<num_v; ++i) {
      dvv_index = find_index(x_dvv, cv_ids[i]);
      if (dvv_index != _NPOS)
	fn_grad_x_std[i] = fn_grad_x(dvv_index);
    }
  }

  // manage size of fn_grad_s output
  if (mixed_s || !std_dvv)
    fn_grad_s_std.size(num_s);
  size_t i, final_num_s;
  if (std_dvv)
    final_num_s = num_s;
  else {
    final_num_s = 0;
    for (i=0; i<num_s; ++i)
      if ( std::find(x_dvv.begin(), x_dvv.end(),
                     acv_ids[acv_map1_indices[i]]) != x_dvv.end() )
	++final_num_s;
  }
  if (fn_grad_s.length() != final_num_s)
    fn_grad_s.size(final_num_s);

  // perform transformation using full Jacobian
  const RealVector& fn_grad_x_trans
    = (std_dvv) ? fn_grad_x : fn_grad_x_std;
  RealVector& fn_grad_s_trans
    = (mixed_s || !std_dvv) ? fn_grad_s_std : fn_grad_s;
  fn_grad_s_trans.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., jacobian_xs,
			   fn_grad_x_trans, 0.);

  // reassemble final fn_grad_s
  if (mixed_s || !std_dvv) {
    size_t cntr = 0;
    for (i=0; i<num_s; ++i) {
      size_t acv_id = acv_ids[acv_map1_indices[i]],
	dvv_index = find_index(x_dvv, acv_id);
      if (dvv_index != _NPOS)
	fn_grad_s(cntr++) = (acv_map2_targets[i] == NO_TARGET) ?
	  fn_grad_x(dvv_index) : // no distribution parameter: if the missing
	  // fn_grad_s component is available in fn_grad_x, then use it; else it
	  // must be updated separately (as in NonDLocalReliability::dg_ds_eval)
	  fn_grad_s_trans[i];    // use the Jacobian-transformed component
      else if (std_dvv)
	fn_grad_s(cntr++) = 0.;
    }
  }

#ifdef DEBUG
  PCout << "\nfn_grad_x:\n"       << fn_grad_x
        << "\nfn_grad_x_trans:\n" << fn_grad_x_trans
        << "\njacobian_xs:\n"     << jacobian_xs
        << "\nfn_grad_s_trans:\n" << fn_grad_s_trans
        << "\nfn_grad_s:\n"       << fn_grad_s << std::endl;
#endif
}


/** This procedure tranforms a Hessian matrix from the original
    user-defined x-space (where evaluations are performed) to
    uncorrelated standard normal space (u-space).  x_vars is the
    vector of the random variables in x-space. */
void NatafTransformation::
trans_hess_X_to_U(const RealSymMatrix& fn_hess_x, RealSymMatrix& fn_hess_u,
		  const RealVector& x_vars, const RealVector& fn_grad_x,
		  const SizetArray& x_dvv, SizetMultiArrayConstView cv_ids)
{
  RealMatrix jacobian_xu;
  jacobian_dX_dU(x_vars, jacobian_xu);

  RealSymMatrixArray hessian_xu;
  bool nonlinear_vars_map = false;
  size_t i, num_v = x_vars.length(); short x_type, u_type;
  for (i=0; i<num_v; ++i) {
    x_type = xDist.active_random_variable_type(i);
    u_type = uDist.active_random_variable_type(i);
    if ( ( ( x_type == CONTINUOUS_RANGE || x_type == UNIFORM ||
	     x_type == CONTINUOUS_INTERVAL_UNCERTAIN ) &&
	   u_type   != STD_UNIFORM ) ||
	 ( x_type   == NORMAL      && u_type != STD_NORMAL ) ||
	 ( x_type   == EXPONENTIAL && u_type != STD_EXPONENTIAL ) ||
	 ( x_type   == BETA        && u_type != STD_BETA   ) ||
	 ( x_type   == GAMMA       && u_type != STD_GAMMA  ) ||
	 ( ( x_type == BOUNDED_NORMAL    || x_type == LOGNORMAL  ||
	     x_type == BOUNDED_LOGNORMAL || x_type == LOGUNIFORM ||
	     x_type == TRIANGULAR        || x_type == GUMBEL     ||
	     x_type == FRECHET           || x_type == WEIBULL ) &&
	   x_type   != u_type ) ||
	 ( x_type   == HISTOGRAM_BIN && u_type != STD_UNIFORM &&
	   x_type   != u_type ) )
      { nonlinear_vars_map = true; break; }
  }

  if (nonlinear_vars_map) // nonlinear transformation has Hessian
    hessian_d2X_dU2(x_vars, hessian_xu);

  trans_hess_X_to_U(fn_hess_x, fn_hess_u, jacobian_xu, hessian_xu,
		    fn_grad_x, x_dvv, cv_ids);
}


/** This procedure tranforms a Hessian matrix from the original
    user-defined x-space (where evaluations are performed) to
    uncorrelated standard normal space (u-space).  This overloaded
    form allows for the separate calculation of jacobian_xu and
    hessian_xu, since these are independent of the response function
    index and can be pulled outside response function loops. */
void NatafTransformation::
trans_hess_X_to_U(const RealSymMatrix& fn_hess_x, RealSymMatrix& fn_hess_u,
		  const RealMatrix& jacobian_xu,
		  const RealSymMatrixArray& hessian_xu,
		  const RealVector& fn_grad_x, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  // Jacobian dimensions = length of random variables = model.cv()
  int  num_v   = jacobian_xu.numRows();
  bool std_dvv = (x_dvv == cv_ids); // standard DVV
  bool nonlinear_vars_map = !hessian_xu.empty();

  RealSymMatrix fn_hess_x_std, fn_hess_u_std;
  RealVector fn_grad_x_std;
  SizetArray dvv_index_array;
  if (std_dvv) {
    if (fn_hess_x.numRows() != num_v) {
      PCerr << "Error: bad fn_hess_x dimension in NatafTransformation::"
	    << "trans_hess_X_to_U()." << std::endl;
      abort_handler(-1);
    }
    if ( nonlinear_vars_map &&
	 ( fn_grad_x.length() != num_v || hessian_xu.size() != num_v ) ) {
      PCerr << "Error: bad dimension in NatafTransformation::"
	    << "trans_hess_X_to_U()." << std::endl;
      abort_handler(-1);
    }
    if (fn_hess_u.numRows() != num_v)
      fn_hess_u.shape(num_v);
  }
  else { // extract relevant DVV components from fn_grad_x & fn_hess_x
    fn_hess_x_std.shape(num_v);
    fn_hess_u_std.shape(num_v);
    if (nonlinear_vars_map)
      fn_grad_x_std.size(num_v);
    size_t i, j, dvv_index_i, dvv_index_j, num_deriv_vars = x_dvv.size();
    dvv_index_array.resize(num_v);
    for (i=0; i<num_v; ++i)
      dvv_index_array[i] = dvv_index_i = find_index(x_dvv, cv_ids[i]);
    if (fn_hess_u.numRows() != num_deriv_vars)
      fn_hess_u.shape(num_deriv_vars);
    // extract relevant DVV components from fn_hess_x
    for (i=0; i<num_v; ++i) {
      dvv_index_i = dvv_index_array[i];
      if (dvv_index_i != _NPOS) {
	if (nonlinear_vars_map)
	  fn_grad_x_std[i] = fn_grad_x(dvv_index_i);
	for (j=0; j<num_v; ++j) {
	  dvv_index_j = dvv_index_array[j];
	  if (dvv_index_j != _NPOS)
	    fn_hess_x_std(i, j) = fn_hess_x(dvv_index_i, dvv_index_j);
	}
      }
    }
  }
  const RealVector&    fn_grad_x_trans = (std_dvv) ? fn_grad_x : fn_grad_x_std;
  const RealSymMatrix& fn_hess_x_trans = (std_dvv) ? fn_hess_x : fn_hess_x_std;
  RealSymMatrix&       fn_hess_u_trans = (std_dvv) ? fn_hess_u : fn_hess_u_std;

  // transform hess_x -> hess_u
  // d^2G/dU^2 = dG/dX^T d^2X/dU^2 + dX/dU^T d^2G/dX^2 dX/dU
  // Note: G(u) may have curvature even if g(x) is linear due to first term.
  Teuchos::symMatTripleProduct(Teuchos::TRANS, 1., fn_hess_x_trans,
                               jacobian_xu, fn_hess_u_trans);
#ifdef DEBUG
  PCout << "\nfnHessU 1st term J^T H_x J:" << fn_hess_u_trans;
#endif

  if (nonlinear_vars_map) { // nonlinear transformation has Hessian
    for (int i=0; i<num_v; ++i) {
      //hessian_xu[i].Scale(fn_grad_x[i]);
      //fn_hess_u += hessian_xu[i];
      const Real&          fn_grad_x_i  = fn_grad_x_trans[i];
      const RealSymMatrix& hessian_xu_i = hessian_xu[i];
      for (int j=0; j<num_v; ++j)
	for (int k=0; k<=j; k++)
	  fn_hess_u_trans(j,k) += fn_grad_x_i * hessian_xu_i(j,k);
#ifdef DEBUG
      PCout << "\nhessian_xu[" << i << "]:\n" << hessian_xu_i
	    << "fn_hess_u_trans increment:\n" << fn_hess_u_trans;
#endif
    }
  }

  if (!std_dvv) { // copy relevant DVV components back into fn_hess_u
    size_t i, j, dvv_index_i, dvv_index_j;
    for (i=0; i<num_v; ++i) {
      dvv_index_i = dvv_index_array[i];
      if (dvv_index_i != _NPOS) {
	for (j=0; j<num_v; ++j) {
	  dvv_index_j = dvv_index_array[j];
	  if (dvv_index_j != _NPOS)
	    fn_hess_u(dvv_index_i, dvv_index_j) = fn_hess_u_trans(i, j);
	}
      }
    }
  }

#ifdef DEBUG
  PCout << "Transformed fn_hess_x:\n" << fn_hess_x
        << "to fn_hess_u:\n" << fn_hess_u;
#endif
}


/** This procedure computes the Jacobian of the transformation x(u).
    x_vars is the vector of random variables in the original
    user-defined x-space. */
void NatafTransformation::
jacobian_dX_dU(const RealVector& x_vars, RealMatrix& jacobian_xu)
{
  if (xDist.correlation()) {
    // dX/dZ = diagonal
    RealMatrix jacobian_xz;
    jacobian_dX_dZ(x_vars, jacobian_xz);

    // dX/dU = dX/dZ dZ/dU = dX/dZ L = dense if variables are correlated
    int num_v = x_vars.length();
    if (jacobian_xu.numRows() != num_v || jacobian_xu.numCols() != num_v)
      jacobian_xu.shape(num_v, num_v);
    jacobian_xu.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1., jacobian_xz,
			 corrCholeskyFactorZ, 0.);
  }
  else // dX/dU = dX/dZ since dZ/dU = I
    jacobian_dX_dZ(x_vars, jacobian_xu);

#ifdef DEBUG
  PCout << "jacobian_dX_dU:\n" << jacobian_xu;
#endif
}


/** This procedure computes the Jacobian of the transformation x(z).
    x_vars is the vector of random variables in the original
    user-defined x-space. */
void NatafTransformation::
jacobian_dX_dZ(const RealVector& x_vars, RealMatrix& jacobian_xz)
{
  int num_v = x_vars.length();
  if (jacobian_xz.numRows() != num_v || jacobian_xz.numCols() != num_v)
    jacobian_xz.shape(num_v, num_v);

  // Rackwitz-Fiessler: Phi(z) = F(x)
  // d/dz -> phi(z) = f(x) dx/dz
  //         dx/dz = phi(z)/f(x)
  // dX/dZ is diagonal as defined by differentiation of trans_Z_to_X()

  Real z_var; short x_type, u_type;
  for (int i=0; i<num_v; ++i) {
    const RandomVariable& x_rv_i = xDist.active_random_variable(i);
    x_type = x_rv_i.type(); u_type = uDist.active_random_variable_type(i);
    if (u_type == x_type)
      jacobian_xz(i, i) = 1.;
    else if (u_type == STD_NORMAL)
      switch (x_type) {
      case NORMAL:      // z = (x - mean)/stdev
	x_rv_i.pull_parameter(N_STD_DEV, jacobian_xz(i, i));       break;
      case LOGNORMAL: { // z = (ln x - lambda)/zeta
	Real zeta;  x_rv_i.pull_parameter(LN_ZETA, zeta);
	jacobian_xz(i, i) = zeta * x_vars[i];
	break;
      }
      default:
	trans_X_to_Z(x_vars[i], z_var, i);
	jacobian_xz(i, i) = NormalRandomVariable::std_pdf(z_var)
	                  / x_rv_i.pdf(x_vars[i]);                        break;
      }
    else if (u_type == STD_UNIFORM) {
      //trans_X_to_Z(x_vars[i], z_var, i);
      //jacobian_xz(i, i) = UniformRandomVariable::std_pdf(z_var)
      //                  / x_rv_i.pdf(x_vars[i]);

      // don't bother to compute z for constant pdf:
      jacobian_xz(i, i)	= UniformRandomVariable::std_pdf(0.)
	                / x_rv_i.pdf(x_vars[i]);
    }
    else if (u_type == STD_EXPONENTIAL && x_type == EXPONENTIAL)
      x_rv_i.pull_parameter(E_BETA,  jacobian_xz(i, i));
    else if (u_type == STD_GAMMA       && x_type == GAMMA)
      x_rv_i.pull_parameter(GA_BETA, jacobian_xz(i, i));
    else if (u_type == STD_BETA        && x_type == BETA) {
      Real l_bnd;  x_rv_i.pull_parameter(BE_LWR_BND, l_bnd);
      Real u_bnd;  x_rv_i.pull_parameter(BE_UPR_BND, u_bnd);
      jacobian_xz(i, i) = (u_bnd - l_bnd) / 2.;
    }
    else {
      PCerr << "Error: unsupported variable mapping for variable " << i
	    << " in NatafTransformation::jacobian_dX_dZ()" << std::endl;
      abort_handler(-1);
    }
  }
}


/** This procedure computes the Jacobian of the transformation u(x).
    x_vars is the vector of random variables in the original
    user-defined x-space. */
void NatafTransformation::
jacobian_dU_dX(const RealVector& x_vars, RealMatrix& jacobian_ux)
{
  if (xDist.correlation()) {
    // dZ/dX = diagonal
    RealMatrix jacobian_zx;
    jacobian_dZ_dX(x_vars, jacobian_zx);

    // dU/dX = dU/dZ dZ/dX = L^-1 dZ/dX = dense if variables are correlated
    // Solve as L dU/dX = dZ/dX
    RealSolver corr_solver;
    corr_solver.setMatrix( Teuchos::rcp(&corrCholeskyFactorZ, false) );
    int num_v = x_vars.length();
    // The inbound jacobian_ux object may be a Teuchos::View we want
    // to update with the solution.  However
    // SerialDenseSolver::solve(), line 647 will call (LHS_ =
    // jacobian_ux) = RHS_, disconnecting the view.  To work around,
    // we use a temporary.  At this writing there are no such use
    // cases, but we don't want to have a latent bug.
    RealMatrix tmp_jac_ux(num_v, num_v);
    // WARNING: this could also disconnect a view...
    if (jacobian_ux.numRows() != num_v || jacobian_ux.numCols() != num_v)
      jacobian_ux.shape(num_v, num_v);
    corr_solver.setVectors( Teuchos::rcp(&tmp_jac_ux, false),
                            Teuchos::rcp(&jacobian_zx, false) );
    corr_solver.solveToRefinedSolution(true);
    corr_solver.solve();
    // Assign into jacobian_ux, which may be a Teuchos::View
    jacobian_ux.assign(tmp_jac_ux);
  }
  else // dU/dX = dZ/dX since dU/dZ = I
    jacobian_dZ_dX(x_vars, jacobian_ux);

#ifdef DEBUG
  PCout << "jacobian_dU_dX:\n" << jacobian_ux;
#endif
}


/** This procedure computes the Jacobian of the transformation z(x).
    x_vars is the vector of random variables in the original
    user-defined x-space. */
void NatafTransformation::
jacobian_dZ_dX(const RealVector& x_vars, RealMatrix& jacobian_zx) 
{
  int num_v = x_vars.length();
  if (jacobian_zx.numRows() != num_v || jacobian_zx.numCols() != num_v)
    jacobian_zx.shape(num_v, num_v);

  // Rackwitz-Fiessler: Phi(z) = F(x)
  // d/dx -> phi(z) dz/dx = f(x)
  //         dz/dx = f(x)/phi(z)
  // dZ/dX is diagonal as defined by differentiation of trans_X_to_Z()

  Real z_var; short x_type, u_type;
  for (int i=0; i<num_v; ++i) {
    const RandomVariable&   x_rv_i = xDist.active_random_variable(i);
    x_type = x_rv_i.type(); u_type = uDist.active_random_variable_type(i);
    if (u_type == x_type)
      jacobian_zx(i, i) = 1.;
    else if (u_type == STD_NORMAL)
      switch (x_type) {
      case NORMAL: {    // z = (x - mean)/stdev
	Real stdev;  x_rv_i.pull_parameter(N_STD_DEV, stdev);
	jacobian_zx(i, i) = 1. / stdev;  break;
      }
      case LOGNORMAL: { // z = (ln x - lambda)/zeta
	Real zeta;  x_rv_i.pull_parameter(LN_ZETA, zeta);
	jacobian_zx(i, i) = 1. / (zeta * x_vars[i]);  break;
      }
      default:
	trans_X_to_Z(x_vars[i], z_var, i);
	jacobian_zx(i, i)
	  = x_rv_i.pdf(x_vars[i]) / NormalRandomVariable::std_pdf(z_var); break;
      }
    else if (u_type == STD_UNIFORM) {
      //trans_X_to_Z(x_vars[i], z_var, i);
      //jacobian_zx(i, i) = x_rv_i.pdf(x_vars[i])
      //                  / UniformRandomVariable::std_pdf(z_var);

      // don't bother to compute z for constant pdf:
      jacobian_zx(i, i)	= x_rv_i.pdf(x_vars[i])
	                / UniformRandomVariable::std_pdf(0.);
    }
    else if (u_type == STD_EXPONENTIAL && x_type == EXPONENTIAL) {
      Real beta;  x_rv_i.pull_parameter(E_BETA, beta);
      jacobian_zx(i, i) = 1. / beta;
    }
    else if (u_type == STD_GAMMA       && x_type == GAMMA) {
      Real beta;  x_rv_i.pull_parameter(GA_BETA, beta);
      jacobian_zx(i, i) = 1. / beta;
    }
    else if (u_type == STD_BETA        && x_type == BETA) {
      Real l_bnd;  x_rv_i.pull_parameter(BE_LWR_BND, l_bnd);
      Real u_bnd;  x_rv_i.pull_parameter(BE_UPR_BND, u_bnd);
      jacobian_zx(i, i) = 2. / (u_bnd -	l_bnd);
    }
    else {
      PCerr << "Error: unsupported variable mapping for variable " << i
	    << " in NatafTransformation::jacobian_dZ_dX()" << std::endl;
      abort_handler(-1);
    }
  }
}


/** This procedure computes the derivative of the original variables x
    with respect to the random variable distribution parameters s.
    This provides the design Jacobian of the transformation for use in
    computing statistical design sensitivities for OUU.

    dX/dS is derived by differentiating trans_Z_to_X with respect to S.
    For the uncorrelated case, u and z are constants.  For the correlated
    case, u is a constant, but z(s) = L(s) u due to Nataf dependence on s
    and dz/ds = dL/ds u. */
void NatafTransformation::
jacobian_dX_dS(const RealVector& x_vars, RealMatrix& jacobian_xs,
	       SizetMultiArrayConstView cv_ids,
	       SizetMultiArrayConstView acv_ids,
	       const SizetArray& acv_map1_indices,
	       const ShortArray& acv_map2_targets)
{
  // Rectangular Jacobian = Gradient^T = num_X by num_S where num_S is the total
  // number of active continuous vars flowed down from a higher iteration level.
  // The number of distribution parameter insertions is <= num_S.
  size_t num_var_map_1c = acv_map1_indices.size();
  int num_v = x_vars.length();
  if (jacobian_xs.numRows() != num_v || jacobian_xs.numCols() != num_var_map_1c)
    jacobian_xs.shape(num_v, num_var_map_1c);

  RealVector z_vars;
  trans_X_to_Z(x_vars, z_vars);

  // For distributions without simple closed-form CDFs (beta, gamma), dx/ds is
  // computed numerically for nonlinear mappings.  If uncorrelated, then this
  // is only needed if the beta/gamma distribution parameters are design vars.
  // If correlated, then the beta/gamma distribution params do not have to be
  // design vars (dx/ds for beta/gamma x will include a dz/ds contribution).
  const RealSymMatrix& x_corr_matrix = xDist.correlation_matrix();
  const BitArray&      x_active_corr = xDist.active_correlations();
  short x_type, u_type;  size_t i, j, cntr_i, cntr_j;
  bool need_xs = false, non_std_beta_gamma_map = false,
       x_corr  = xDist.correlation(), no_mask = x_active_corr.empty();
  // non_std_beta_gamma_map detects unsupported X->U cases for dx_ds and
  // dz_ds_factor; this is augmented below with unsupported dist param targets.
  for (i=0, cntr_i=0; i<num_v; ++i) {
    x_type = xDist.active_random_variable_type(i);
    u_type = uDist.active_random_variable_type(i);
    if ( ( x_type == BETA  && u_type != STD_BETA  ) ||
	 ( x_type == GAMMA && u_type != STD_GAMMA ) ) {
      // non-std mapping invalidates support for BE_LWR_BND,BE_UPR_BND,GA_BETA
      // in {Beta,Gamma}RandomVariable::dx_ds()
      non_std_beta_gamma_map = true;
      // correlation with non-std mapped beta/gamma invalidates support in
      // {Beta,Gamma}RandomVariable::dz_ds_factor(), requiring numerical dx/ds
      if (x_corr)
	// since we don't check all rows, check *all* columns despite symmetry
	for (j=0, cntr_j=0; j<num_v; ++j)
	  if (no_mask || x_active_corr[j]) {
	    if (i != j && std::abs(x_corr_matrix(cntr_i,cntr_j)) > SMALL_NUMBER)
	      { need_xs = true; break; }
	    ++cntr_j;
	  }
    }
    if (no_mask || x_active_corr[i]) ++cntr_i;
  }
  // If numerical dx/ds not already reqd due to correlation-based contributions
  // from dz/ds, detect cases where dx_ds is not supported for particular s.
  if (!need_xs) {
    ShortArray::const_iterator cit;
    for (cit=acv_map2_targets.begin(); cit!=acv_map2_targets.end(); ++cit) {
      short tgt = *cit;
      if ( tgt == BE_ALPHA || tgt == BE_BETA || tgt == GA_ALPHA ||
	   ( non_std_beta_gamma_map &&
	     ( tgt == BE_LWR_BND || tgt == BE_UPR_BND || tgt == GA_BETA ) ) )
	{ need_xs = true; break; }
    }
  }

  // If need_xs, then we compute all of dx/ds and insert the necessary
  // components directly into jacobian_xs below.  If x-space correlations,
  // we compute all of dz/ds and utilize it within loop in combination with
  // dz_ds_factor().  Complete numerical jacobians are computed, even if
  // only certain portions (beta/gamma rows) are needed.
  RealMatrix num_dx_ds, num_dz_ds;
  if (need_xs || x_corr) // compute num_dx_ds and/or num_dz_ds
    numerical_design_jacobian(x_vars, need_xs, num_dx_ds, x_corr, num_dz_ds,
			      cv_ids, acv_ids, acv_map1_indices,
			      acv_map2_targets);
  if (need_xs)
    for (j=0; j<num_v; ++j)              // loop over X
      switch (xDist.active_random_variable_type(j)) {
      case BETA: case GAMMA:
	for (i=0; i<num_var_map_1c; ++i) // loop over S
	  jacobian_xs(j, i) = num_dx_ds(j, i);
	break;
      }

  Real x, z; bool numerical_xs;
  for (i=0; i<num_var_map_1c; ++i) { // loop over S
    size_t cv_index = find_index(cv_ids, acv_ids[acv_map1_indices[i]]);
    // If x_dvv were passed, it would be possible to distinguish different
    // fn_grad_x components, allowing passthrough for computing fn_grad_s for
    // augmented design variables.  For now, this has to be handled spearately
    // in NonDLocalReliability::dg_ds_eval() and NonD::trans_grad_X_to_S().
    //if (cv_index == _NPOS) // augmented variable: define identity mapping
    //  jacobian_xs(dvv_index, i) = 1.;
    //else {
    if (cv_index != _NPOS) {
      short target2 = acv_map2_targets[i];
      for (j=0; j<num_v; ++j) {      // loop over X
	// Jacobian row    = X value = j
	// Jacobian column = S value = i

	const RandomVariable& x_rv_j = xDist.active_random_variable(j);
	x_type = x_rv_j.type();
        numerical_xs = ( need_xs && (x_type == BETA || x_type == GAMMA) );
	if (!numerical_xs) { // else jacobian_xs already updated above
	  x = x_vars[j]; z = z_vars[j];
	  u_type = uDist.active_random_variable_type(j);
	  // corresponding variable has derivative w.r.t. its distribution param
	  if (j == cv_index)
	    jacobian_xs(j, i)  = x_rv_j.dx_ds(target2, u_type, x, z);
	  // dz/ds contributions (deriv of jth variable w.r.t. any variable's
	  // dist param) are included if correlated
	  // Note: BOUNDED_{NORMAL,LOGNORMAL},LOGUNIFORM,TRIANGULAR,BETA cases
	  // will currently be zero, but dL/ds can be nonzero in the future
	  // once correlation warping is more complete.
	  if (x_corr)
	    jacobian_xs(j, i) += x_rv_j.dz_ds_factor(u_type, x, z)
	                      *  num_dz_ds(j, i);
	}
      }
    }
  }
}


/** This procedure computes numerical derivatives of x and/or z with respect to
    distribution parameters s, and is used by jacobian_dX_dS() to provide data
    that is not available analytically.  Numerical dz/ds involves dL/ds
    (z(s) = L(s) u and dz/ds = dL/ds u) and is needed to evaluate dx/ds
    semi-analytically for correlated variables.  Numerical dx/ds is needed for
    distributions lacking simple closed-form CDF expressions (beta and gamma
    distributions). */
void NatafTransformation::
numerical_design_jacobian(const RealVector& x_vars,
			  bool xs, RealMatrix& num_jacobian_xs,
			  bool zs, RealMatrix& num_jacobian_zs,
			  SizetMultiArrayConstView cv_ids,
			  SizetMultiArrayConstView acv_ids,
			  const SizetArray& acv_map1_indices,
			  const ShortArray& acv_map2_targets)
{
  // For correlated vars, correlation matrix C = C(s) due to Nataf modifications
  //   z(s) = L(s) u  ->  dz/ds = dL/ds u  ->  need dL/ds
  //   C(s) = L(s) L(s)^T
  //   dC/ds (which could be derived analytically) = L dL/ds^T + dL/ds L^T
  // This takes the form dC/ds = A + A^T where A = L dL/ds^T
  // Unfortunately, solution of this equation for general A (which could
  // provide dL/ds) given symmetric dC/ds is not possible since it is nonunique.
  // Since we will not be differentiating the Cholesky solver, we will use
  // semi-analytic design sensitivities with numerical dL/ds.  Note that
  // numerical dz/ds is simpler and would likely be just as effective, but in
  // general, semi-analytic sensitivities should minimize the numerical portion.

  // Rectangular Jacobians = Gradient^T = num_Z x num_S where num_S is the total
  // number of active continuous vars flowed down from a higher iteration level.
  // The number of distribution parameter insertions is <= num_S.
  size_t i, j, k, num_var_map_1c = acv_map1_indices.size();
  int x_len = x_vars.length();
  if (xs && (num_jacobian_xs.numRows() != x_len ||
	     num_jacobian_xs.numCols() != num_var_map_1c) )
    num_jacobian_xs.shape(x_len, num_var_map_1c);
  if (zs && (num_jacobian_zs.numRows() != x_len ||
	     num_jacobian_zs.numCols() != num_var_map_1c) )
    num_jacobian_zs.shape(x_len, num_var_map_1c);

  RealMatrix L_s_plus_h, dL_dsi;
  RealVector dz_dsi;
  //RealVector z_vars_s_plus_h, z_vars_s_minus_h;
  RealVector x_vars_s_plus_h, x_vars_s_minus_h;
  if (zs) {
    L_s_plus_h.shape(x_len, x_len);
    dL_dsi.shape(x_len, x_len);
    dz_dsi.size(x_len);
  }

  RealVector u_vars;
  trans_X_to_U(x_vars, u_vars);

  Real fd_grad_ss = 1.e-4;
  RealMatrix chol_z0(corrCholeskyFactorZ);
  for (i=0; i<num_var_map_1c; i++) {

    size_t cv_index        = find_index(cv_ids, acv_ids[acv_map1_indices[i]]);
    short  acv_map2_target = acv_map2_targets[i];
    if (cv_index != _NPOS && acv_map2_target != NO_TARGET) {

      Pecos::RandomVariable& x_rv_i = xDist.active_random_variable(cv_index);

      Real s0;  x_rv_i.pull_parameter(acv_map2_target, s0);

      // Compute the offset for the ith gradient variable.
      // Enforce a minimum delta of fdgss*.01
      Real h_mag = fd_grad_ss * std::max(std::fabs(s0), .01);
      Real h = (s0 < 0.0) ? -h_mag : h_mag; // h has same sign as s0

      // -----------------------------------
      // Evaluate (L/z_vars/x_vars)_s_plus_h
      // -----------------------------------
      Real s1 = s0 + h;
      // update randomVars & corrCholeskyFactorZ:
      x_rv_i.push_parameter(acv_map2_target, s1);
      transform_correlations();
      if (zs) {
	L_s_plus_h = corrCholeskyFactorZ;        // L
	//trans_U_to_Z(u_vars, z_vars_s_plus_h); // z
      }
      if (xs)
	trans_U_to_X(u_vars, x_vars_s_plus_h);   // x

      // ------------------------------------
      // Evaluate (L/z_vars/x_vars)_s_minus_h
      // ------------------------------------
      s1 = s0 - h;
      // update randomVars & corrCholeskyFactorZ:
      x_rv_i.push_parameter(acv_map2_target, s1);
      transform_correlations();
      //if (zs) {
        // utilize corrCholeskyFactorZ below      // L
        //trans_U_to_Z(u_vars, z_vars_s_minus_h); // z
      //}
      if (xs)
	trans_U_to_X(u_vars, x_vars_s_minus_h);   // x

      // -------------------------------
      // Compute the central differences
      // -------------------------------
      if (zs) {
	for (j=0; j<x_len; j++)                            // dL/ds
	  for (k=0; k<=j; k++)
	    dL_dsi(j, k) = (L_s_plus_h(j, k) - corrCholeskyFactorZ(j, k))/2./h;
	dz_dsi.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1., dL_dsi,
			u_vars, 0.); // dz/ds
	for (j=0; j<x_len; j++)
	  num_jacobian_zs(j, i) = dz_dsi(j);
	//for (j=0; j<x_len; j++)                          // dz/ds (alt)
	//  num_jacobian_zs(j, i)=(z_vars_s_plus_h(j)-z_vars_s_minus_h(j))/2./h;
      }
      if (xs)
	for (j=0; j<x_len; j++)                            // dx/ds
	  num_jacobian_xs(j,i) = (x_vars_s_plus_h(j)-x_vars_s_minus_h(j))/2./h;

      // reset s0 (corrCholeskyFactorZ reset can be deferred):
      x_rv_i.push_parameter(acv_map2_target, s0);
    }
  }
  // reset corrCholeskyFactorZ:
  corrCholeskyFactorZ = chol_z0;
}


/** This procedure computes the Hessian of the transformation x(u).
    hessian_xu is a 3D tensor modeled as an array of matrices, where
    the i_th matrix is d^2X_i/dU^2.  x_vars is the vector of random
    variables in the original user-defined x-space. */
void NatafTransformation::
hessian_d2X_dU2(const RealVector& x_vars, RealSymMatrixArray& hessian_xu)
{
  if (xDist.correlation()) {
    // d^2X/dZ^2
    int num_v = x_vars.length();
    RealSymMatrixArray hessian_xz(num_v);
    hessian_d2X_dZ2(x_vars, hessian_xz);

    if (hessian_xu.size() != num_v)
      hessian_xu.resize(num_v);
    for (int i=0; i<num_v; ++i) {
      // d^2X/dU^2 = dX/dZ^T d^2Z/dU^2 + dZ/dU^T d^2X/dZ^2 dZ/dU
      //           = L^T d^2X/dZ^2 L
      if (hessian_xu[i].numRows() != num_v)
	hessian_xu[i].shape(num_v);
      Teuchos::symMatTripleProduct(Teuchos::TRANS, 1., hessian_xz[i],
                                   corrCholeskyFactorZ, hessian_xu[i]);
    }
  }
  else // d^2X/dU^2 = d^2X/dZ^2 since dZ/dU = I
    hessian_d2X_dZ2(x_vars, hessian_xu);
}


/** This procedure computes the Hessian of the transformation x(z).
    hessian_xz is a 3D tensor modeled as an array of matrices, where
    the i_th matrix is d^2X_i/dZ^2.  x_vars is the vector of random
    variables in the original user-defined x-space. */
void NatafTransformation::
hessian_d2X_dZ2(const RealVector& x_vars, RealSymMatrixArray& hessian_xz)
{
  // This routine calculates the Hessian of the transformation x(z):
  //
  // For z type is STD_NORMAL:
  // d^2x/dz^2 = d/dz (dx/dz) = d/dz ( phi(z)/f(x) )
  //           = (f(x) phi'(z) - phi(z) f'(x) dx/dz)/f(x)^2
  //           = -phi(z)/f(x)^2 (z f(x) + f'(x) dx/dz) [since phi'(z)=-z phi(z)]
  //           = -dx/dz (z + f'(x)/f(x) dx/dz)
  //
  // For z type is STD_UNIFORM:
  // d^2x/dz^2 = d/dz (dx/dz) = d/dz ( u(z)/f(x) )
  //           = (f(x) u'(z) - u(z) f'(x) dx/dz)/f(x)^2
  //           = -dx/dz u(z)/f(x) f'(x)/f(x)           [since u'(z)=0]
  //           = -dx/dz dx/dz f'(x)/f(x)
  //
  // This requires the additional calculation of f'(x), the derivative of
  // the PDF.  For cases with f'(x) = 0 (e.g., uniform), the expression can
  // be simplified to d^2x/dz^2 = -z dx/dz

  int num_v = x_vars.length();
  if (hessian_xz.size() != num_v)
    hessian_xz.resize(num_v);

  Real x, z, dx_dz, pdf; short x_type, u_type;
  for (int i=0; i<num_v; ++i) {
    const RandomVariable& x_rv_i = xDist.active_random_variable(i);
    x_type = x_rv_i.type(); u_type = uDist.active_random_variable_type(i);
    if (hessian_xz[i].numRows() != num_v)
      hessian_xz[i].shape(num_v);
    // each d^2X_i/dZ^2 has a single entry on the diagonal as defined by
    // differentiation of jacobian_dX_dZ()

    if (u_type == x_type)
      hessian_xz[i](i, i) = 0.;
    else if (u_type == STD_NORMAL)
      switch (x_type) {
      case NORMAL:
	hessian_xz[i](i, i) = 0.; break;
      case LOGNORMAL: {	// dx/dz = zeta x --> d^2x/dz^2 = zeta dx/dz = zeta^2 x
	Real zeta;  x_rv_i.pull_parameter(LN_ZETA, zeta);
	hessian_xz[i](i, i) = zeta * zeta * x_vars[i]; break;
      }
      case CONTINUOUS_RANGE: case CONTINUOUS_INTERVAL_UNCERTAIN:
      case UNIFORM:          case HISTOGRAM_BIN: // pdf_grad is zero
	x = x_vars[i]; trans_X_to_Z(x, z, i);
	dx_dz = NormalRandomVariable::std_pdf(z) / x_rv_i.pdf(x);
	hessian_xz[i](i, i) = -dx_dz * z;
	break;
      default:
	x = x_vars[i]; trans_X_to_Z(x, z, i); pdf = x_rv_i.pdf(x);
	dx_dz = NormalRandomVariable::std_pdf(z) / pdf;
	hessian_xz[i](i, i)
	  = -dx_dz * (z + x_rv_i.pdf_gradient(x) * dx_dz / pdf);
	break;
      }
    else if (u_type == STD_UNIFORM)
      switch (x_type) {
      case CONTINUOUS_RANGE: case CONTINUOUS_INTERVAL_UNCERTAIN:
      case UNIFORM:          case HISTOGRAM_BIN:
	hessian_xz[i](i, i) = 0.; break;
      case LOGUNIFORM: {
	Real l_bnd;  x_rv_i.pull_parameter(LU_LWR_BND, l_bnd);
	Real u_bnd;  x_rv_i.pull_parameter(LU_UPR_BND, u_bnd);
	Real log_range = std::log(u_bnd) - std::log(l_bnd);
	hessian_xz[i](i, i) = log_range * log_range * x_vars[i] / 4.; break;
      }
      default:
	x = x_vars[i]; pdf = x_rv_i.pdf(x);
	//trans_X_to_Z(x, z, i); 
	//dx_dz = UniformRandomVariable::std_pdf(z) / pdf;

	// don't bother to compute z constant pdf eval:
	dx_dz = UniformRandomVariable::std_pdf(0.) / pdf;
	hessian_xz[i](i, i) = -dx_dz * dx_dz * x_rv_i.pdf_gradient(x) / pdf;
	break;
      }
    else if ( (u_type == STD_EXPONENTIAL && x_type == EXPONENTIAL) ||
	      (u_type == STD_GAMMA       && x_type == GAMMA)       ||
	      (u_type == STD_BETA        && x_type == BETA) )
      hessian_xz[i](i, i) = 0.;
    else {
      PCerr << "Error: unsupported variable mapping for variable " << i
	    << " in NatafTransformation::hessian_d2X_dZ2()" << std::endl;
      abort_handler(-1);
    }
  }
}

} // namespace Pecos
