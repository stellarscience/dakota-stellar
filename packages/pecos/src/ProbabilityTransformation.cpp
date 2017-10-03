/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "pecos_stat_util.hpp"
#include "NatafTransformation.hpp"
#include "DistributionParams.hpp"

static const char rcsId[]="@(#) $Id: ProbabilityTransformation.cpp 4768 2007-12-17 17:49:32Z mseldre $";

namespace Pecos {


/** This constructor is the one which must build the base class data
    for all derived classes.  get_prob_trans() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_prob_trans() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~ProbabilityTransformation). */
ProbabilityTransformation::ProbabilityTransformation(BaseConstructor):
  correlationFlagX(false), probTransRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "ProbabilityTransformation::ProbabilityTransformation(Base"
        << "Constructor) called to build base class for letter." << std::endl;
#endif
}


/** The default constructor: probTransRep is NULL in this case.  This
    makes it necessary to check for NULL in the copy constructor,
    assignment operator, and destructor. */
ProbabilityTransformation::ProbabilityTransformation():
  probTransRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "ProbabilityTransformation::ProbabilityTransformation() called to "
        << "build empty envelope." << std::endl;
#endif
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_prob_trans, since ProbabilityTransformation(BaseConstructor)
    builds the actual base class data for the derived transformations. */
ProbabilityTransformation::
ProbabilityTransformation(const String& prob_trans_type):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "ProbabilityTransformation::ProbabilityTransformation(string&) "
        << "called to instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  probTransRep = get_prob_trans(prob_trans_type);
  if ( !probTransRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize probTransRep to the 
    appropriate derived type. */
ProbabilityTransformation* ProbabilityTransformation::
get_prob_trans(const String& prob_trans_type)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_prob_trans(string&)."
        << std::endl;
#endif

  if (prob_trans_type == "nataf")
    return new NatafTransformation();
  else {
    PCerr << "Error: ProbabilityTransformation type " << prob_trans_type
	  << " not available." << std::endl;
    return NULL;
  }
}


/** Copy constructor manages sharing of probTransRep and incrementing
    of referenceCount. */
ProbabilityTransformation::
ProbabilityTransformation(const ProbabilityTransformation& prob_trans)
{
  // Increment new (no old to decrement)
  probTransRep = prob_trans.probTransRep;
  if (probTransRep) // Check for an assignment of NULL
    probTransRep->referenceCount++;

#ifdef REFCOUNT_DEBUG
  PCout << "ProbabilityTransformation::ProbabilityTransformation("
        << "ProbabilityTransformation&)" << std::endl;
  if (probTransRep)
    PCout << "probTransRep referenceCount = " << probTransRep->referenceCount
	  << std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old probTransRep, assigns
    new probTransRep, and increments referenceCount for new probTransRep. */
ProbabilityTransformation ProbabilityTransformation::
operator=(const ProbabilityTransformation& prob_trans)
{
  if (probTransRep != prob_trans.probTransRep) { // normal case: old != new
    // Decrement old
    if (probTransRep) // Check for null pointer
      if (--probTransRep->referenceCount == 0) 
	delete probTransRep;
    // Assign and increment new
    probTransRep = prob_trans.probTransRep;
    if (probTransRep) // Check for an assignment of NULL
      probTransRep->referenceCount++;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  PCout << "ProbabilityTransformation::operator=(ProbabilityTransformation&)"
        << std::endl;
  if (probTransRep)
    PCout << "probTransRep referenceCount = " << probTransRep->referenceCount
	  << std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes probTransRep
    when referenceCount reaches zero. */
ProbabilityTransformation::~ProbabilityTransformation()
{ 
  // Check for NULL pointer 
  if (probTransRep) {
    --probTransRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    PCout << "probTransRep referenceCount decremented to " 
	  << probTransRep->referenceCount << std::endl;
#endif
    if (probTransRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      PCout << "deleting probTransRep" << std::endl;
#endif
      delete probTransRep;
    }
  }
}


/** This function provides a deep copy (the copy constructor supports
    shallow copies with shared reps) and is commonly used to publish
    tranformation data when the Model variables are in a transformed
    space (e.g., u-space) and x-space data may not be generated
    directly.  This allows for the use of inverse transformations to
    return the transformed space variables to their original states. */
void ProbabilityTransformation::
copy(const ProbabilityTransformation& prob_trans)
{
  if (probTransRep) // target is envelope
    probTransRep->copy(prob_trans);
  else { // target is letter
    if (prob_trans.probTransRep) { // source is envelope
      randomVarsX = prob_trans.probTransRep->randomVarsX;//[i].copy(); TO DO
      ranVarTypesU        = prob_trans.probTransRep->ranVarTypesU;
      correlationFlagX    = prob_trans.probTransRep->correlationFlagX;
      corrMatrixX         = prob_trans.probTransRep->corrMatrixX;
      corrCholeskyFactorZ = prob_trans.probTransRep->corrCholeskyFactorZ;
    }
    else { // source is letter
      randomVarsX         = prob_trans.randomVarsX;//[i].copy(); TO DO
      ranVarTypesU        = prob_trans.ranVarTypesU;
      correlationFlagX    = prob_trans.correlationFlagX;
      corrMatrixX         = prob_trans.corrMatrixX;
      corrCholeskyFactorZ = prob_trans.corrCholeskyFactorZ;
    }
  }
}


void ProbabilityTransformation::
initialize_random_variable_types(const ShortArray& x_types)
{
  if (probTransRep)
    probTransRep->initialize_random_variable_types(x_types);
  else {
    // (default) construction of x-space random variables occurs once
    // (updates to distribution parameters can occur repeatedly)
    size_t i, num_v = x_types.size();
    randomVarsX.resize(num_v);
    for (i=0; i<num_v; ++i)
      randomVarsX[i] = RandomVariable(x_types[i]);
  }
}


void ProbabilityTransformation::
initialize_random_variable_types(const ShortArray& x_types,
				 const ShortArray& u_types)
{
  if (probTransRep)
    probTransRep->initialize_random_variable_types(x_types, u_types);
  else {
    ranVarTypesU = u_types;
    initialize_random_variable_types(x_types);
  }
}


void ProbabilityTransformation::
initialize_random_variable_parameters(const RealVector& cd_l_bnds,
				      const RealVector& cd_u_bnds,
				      const AleatoryDistParams& adp,
				      const EpistemicDistParams& edp,
				      const RealVector& cs_l_bnds,
				      const RealVector& cs_u_bnds)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->
      initialize_random_variable_parameters(cd_l_bnds, cd_u_bnds, adp, edp,
					    cs_l_bnds, cs_u_bnds);
  else {
    size_t i, num_v = randomVarsX.size(), cd_cntr = 0, n_cntr = 0, ln_cntr = 0,
      u_cntr = 0, lu_cntr = 0, t_cntr = 0, e_cntr = 0, b_cntr = 0, ga_cntr = 0,
      gu_cntr = 0, f_cntr = 0, w_cntr = 0, h_cntr = 0, p_cntr = 0, bi_cntr = 0,
      nbi_cntr = 0, ge_cntr = 0, hge_cntr = 0, hpi_cntr = 0, hps_cntr = 0,
      hpr_cntr = 0, ci_cntr = 0, cs_cntr = 0;
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    RandomVariable* rv_rep_i;
    const RealVector& ln_means     = adp.lognormal_means();
    const RealVector& ln_std_devs  = adp.lognormal_std_deviations();
    const RealVector& ln_lambdas   = adp.lognormal_lambdas();
    const RealVector& ln_zetas     = adp.lognormal_zetas();
    const RealVector& ln_err_facts = adp.lognormal_error_factors();
    for (i=0; i<num_v; ++i) {
      rv_rep_i = randomVarsX[i].random_variable_rep();
      switch (rv_rep_i->type()) {
      case CONTINUOUS_DESIGN: {
	Real lwr = cd_l_bnds[cd_cntr], upr = cd_u_bnds[cd_cntr];
	if (lwr == -dbl_inf || upr == dbl_inf) {
	  PCerr << "Error: bounds specification required for design variable "
	       << "transformation to standard uniform." << std::endl;
	  abort_handler(-1);
	}
	((UniformRandomVariable*)rv_rep_i)->update(lwr, upr);
	++cd_cntr; break;
      }
      case NORMAL:
	((NormalRandomVariable*)rv_rep_i)->
	  update(adp.normal_mean(n_cntr), adp.normal_std_deviation(n_cntr));
	++n_cntr; break;
      case BOUNDED_NORMAL:
	((BoundedNormalRandomVariable*)rv_rep_i)->
	  update(adp.normal_mean(n_cntr), adp.normal_std_deviation(n_cntr),
		 adp.normal_lower_bound(n_cntr),adp.normal_upper_bound(n_cntr));
	++n_cntr; break;
      case LOGNORMAL: {
	Real lambda, zeta;
	params_from_lognormal_spec(ln_means, ln_std_devs, ln_lambdas, ln_zetas,
				   ln_err_facts, ln_cntr, lambda, zeta);
        ((LognormalRandomVariable*)rv_rep_i)->update(lambda, zeta);
	++ln_cntr; break;
      }
      case BOUNDED_LOGNORMAL: {
	Real lambda, zeta;
	params_from_lognormal_spec(ln_means, ln_std_devs, ln_lambdas, ln_zetas,
				   ln_err_facts, ln_cntr, lambda, zeta);
	((BoundedLognormalRandomVariable*)rv_rep_i)->
	  update(lambda, zeta, adp.lognormal_lower_bound(ln_cntr),
		 adp.lognormal_upper_bound(ln_cntr));
	++ln_cntr; break;
      }
      case UNIFORM:
	((UniformRandomVariable*)rv_rep_i)->
	  update(adp.uniform_lower_bound(u_cntr),
		 adp.uniform_upper_bound(u_cntr));
	++u_cntr; break;
      case LOGUNIFORM:
	((LoguniformRandomVariable*)rv_rep_i)->
	  update(adp.loguniform_lower_bound(lu_cntr),
		 adp.loguniform_upper_bound(lu_cntr));
	++lu_cntr; break;
      case TRIANGULAR:
	((TriangularRandomVariable*)rv_rep_i)->
	  update(adp.triangular_lower_bound(t_cntr),
		 adp.triangular_mode(t_cntr),
		 adp.triangular_upper_bound(t_cntr));
	++t_cntr; break;
      case EXPONENTIAL:
	((ExponentialRandomVariable*)rv_rep_i)->
	  update(adp.exponential_beta(e_cntr));
	++e_cntr; break;
      case BETA:
	((BetaRandomVariable*)rv_rep_i)->
	  update(adp.beta_alpha(b_cntr), adp.beta_beta(b_cntr),
		 adp.beta_lower_bound(b_cntr), adp.beta_upper_bound(b_cntr));
	++b_cntr; break;
      case GAMMA:
	((GammaRandomVariable*)rv_rep_i)->
	  update(adp.gamma_alpha(ga_cntr), adp.gamma_beta(ga_cntr));
	+ga_cntr; break;
      case GUMBEL:
	((GumbelRandomVariable*)rv_rep_i)->
	  update(adp.gumbel_alpha(gu_cntr), adp.gumbel_beta(gu_cntr));
	++gu_cntr; break;
      case FRECHET:
	((FrechetRandomVariable*)rv_rep_i)->
	  update(adp.frechet_alpha(f_cntr), adp.frechet_beta(f_cntr));
	++f_cntr; break;
      case WEIBULL:
	((WeibullRandomVariable*)rv_rep_i)->
	  update(adp.weibull_alpha(w_cntr), adp.weibull_beta(w_cntr));
	++w_cntr; break;
      case HISTOGRAM_BIN:
	((HistogramBinRandomVariable*)rv_rep_i)->
	  update(adp.histogram_bin_pairs(h_cntr));
	++h_cntr; break;

      // discrete int aleatory uncertain
      case POISSON:
	((PoissonRandomVariable*)rv_rep_i)->update(adp.poisson_lambda(p_cntr));
	++p_cntr; break;
      case BINOMIAL:
 	((BinomialRandomVariable*)rv_rep_i)->
	  update(adp.binomial_num_trials(bi_cntr),
		 adp.binomial_probability_per_trial(bi_cntr));
	++bi_cntr; break;
      case NEGATIVE_BINOMIAL:
 	((NegBinomialRandomVariable*)rv_rep_i)->
	  update(adp.negative_binomial_num_trials(nbi_cntr),
		 adp.negative_binomial_probability_per_trial(nbi_cntr));
	++nbi_cntr; break;
      case GEOMETRIC:
 	((GeometricRandomVariable*)rv_rep_i)->
	  update(adp.geometric_probability_per_trial(ge_cntr));
	++ge_cntr; break;
      case HYPERGEOMETRIC:
 	((HypergeometricRandomVariable*)rv_rep_i)->
	  update(adp.hypergeometric_total_population(hge_cntr),
		 adp.hypergeometric_selected_population(hge_cntr),
		 adp.hypergeometric_num_drawn(hge_cntr));
	++hge_cntr; break;
      case HISTOGRAM_PT_INT:
 	((HistogramPtRandomVariable*)rv_rep_i)->
	  update(adp.histogram_point_int_pairs(hpi_cntr));
	++hpi_cntr; break;

      // discrete string aleatory uncertain
      case HISTOGRAM_PT_STRING:
 	((HistogramPtRandomVariable*)rv_rep_i)->
	  update(adp.histogram_point_string_pairs(hps_cntr));
	++hps_cntr; break;

      // discrete real aleatory uncertain
      case HISTOGRAM_PT_REAL:
 	((HistogramPtRandomVariable*)rv_rep_i)->
	  update(adp.histogram_point_real_pairs(hpr_cntr));
	++hpr_cntr; break;

      case CONTINUOUS_INTERVAL: {
	const RealRealPairRealMap& ci_bp
	  = edp.continuous_interval_basic_probabilities(ci_cntr);
	// process interval sets for overall bounds; since intervals are sorted,
	// lb should be in 1st interval but go ahead and process lb same as ub
	Real lb = dbl_inf, ub = -dbl_inf;
	RealRealPairRealMap::const_iterator cit;
	for (cit=ci_bp.begin(); cit!=ci_bp.end(); ++cit) {
	  const RealRealPair& bnds = cit->first;
	  if (bnds.first  < lb) lb = bnds.first;
	  if (bnds.second > ub) ub = bnds.second;
	}
	((UniformRandomVariable*)rv_rep_i)->update(lb, ub);
	++ci_cntr; break;
      }

      // discrete int epistemic uncertain

      // discrete string epistemic uncertain

      // discrete real epistemic uncertain

      case CONTINUOUS_STATE: {
	Real lwr = cs_l_bnds[cs_cntr], upr = cs_u_bnds[cs_cntr];
	if (lwr == -dbl_inf || upr == dbl_inf) {
	  PCerr << "Error: bounds specification required for state variable "
	       << "transformation to standard uniform." << std::endl;
	  abort_handler(-1);
	}
	((UniformRandomVariable*)rv_rep_i)->update(lwr, upr);
	++cs_cntr; break;
      }
      case STD_NORMAL:      ++n_cntr; break;
      case STD_UNIFORM:     ++u_cntr; break;
      case STD_EXPONENTIAL: ++e_cntr; break;
      case STD_BETA:        ++b_cntr; break;
      case STD_GAMMA:      ++ga_cntr; break;
      //default:
      }
    }
  }
}


void ProbabilityTransformation::
initialize_random_variable_correlations(const RealSymMatrix& x_corr)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->initialize_random_variable_correlations(x_corr);
  else {
    corrMatrixX = x_corr;
    size_t num_ran_vars = x_corr.numRows();
    correlationFlagX = false;
    for (size_t i=1; i<num_ran_vars; i++)
      for (size_t j=0; j<i; j++)
	if (std::fabs(x_corr(i,j)) > SMALL_NUMBER)
	  correlationFlagX = true;
  }
}


void ProbabilityTransformation::
reshape_correlation_matrix(size_t num_leading_vars,
			   size_t num_probabilistic_vars,
			   size_t num_trailing_vars)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->reshape_correlation_matrix(num_leading_vars,
					     num_probabilistic_vars,
					     num_trailing_vars);
  else {
    if (!correlationFlagX)
      return;

    size_t i, j, offset, num_corr_vars = corrMatrixX.numRows(),
      num_active_vars = num_leading_vars + num_probabilistic_vars +
      num_trailing_vars;
    if (num_corr_vars != num_active_vars) {
      if (num_corr_vars != num_probabilistic_vars) {
	PCerr << "\nError: unknown symmetric matrix dim (" << num_corr_vars
	      << ") in ProbabilityTransformation::reshape_correlation_matrix()."
	      << std::endl;
	abort_handler(-1);
      }
      RealSymMatrix old_corr_matrix(corrMatrixX);
      corrMatrixX.shape(num_active_vars); // initializes to zero
      for (i=0; i<num_leading_vars; i++)
	corrMatrixX(i,i) = 1.;
      offset = num_leading_vars;
      for (i=0; i<num_probabilistic_vars; i++)
	for (j=0; j<num_probabilistic_vars; j++)
	  corrMatrixX(i+offset,j+offset) = old_corr_matrix(i,j);
      offset += num_probabilistic_vars;
      for (i=0; i<num_trailing_vars; i++)
	corrMatrixX(i+offset,i+offset) = 1.;
    }
  }
}


RealRealPairArray ProbabilityTransformation::u_moments() const
{
  if (probTransRep) return probTransRep->u_moments();
  else {
    size_t i, num_v = randomVarsX.size();
    RealRealPairArray u_mom(num_v);
    Real unif_stdev = 1./std::sqrt(3.);
    for (i=0; i<num_v; ++i)
      switch (ranVarTypesU[i]) {
      case STD_NORMAL:      u_mom[i] = RealRealPair(0., 1.);         break;
      case STD_UNIFORM:     u_mom[i] = RealRealPair(0., unif_stdev); break;
      case STD_EXPONENTIAL: u_mom[i] = RealRealPair(1., 1.);         break;
      case STD_BETA: {
	check_x_type(i, BETA);
	Real mean, stdev;
	BetaRandomVariable::
	  moments_from_params(randomVarsX[i].parameter(BE_ALPHA),
			      randomVarsX[i].parameter(BE_BETA), -1., 1.,
			      mean, stdev);
        u_mom[i] = RealRealPair(mean, stdev); break;
      }
      case STD_GAMMA: {
	check_x_type(i, GAMMA);
	Real mean, stdev;
	GammaRandomVariable::
	  moments_from_params(randomVarsX[i].parameter(GA_ALPHA), 1.,
			      mean, stdev);
        u_mom[i] = RealRealPair(mean, stdev); break;
      }
      default: // no transformation (e.g., PCE w/ numerically-generated basis)
	check_x_type(i, ranVarTypesU[i]);
	u_mom[i] = randomVarsX[i].moments(); break;
      }
    return u_mom;
  }
}


RealRealPairArray ProbabilityTransformation::u_bounds() const
{
  if (probTransRep) return probTransRep->u_bounds();
  else {
    size_t i, num_v = randomVarsX.size();
    RealRealPairArray u_bnds(num_v);
    Real dbl_inf = std::numeric_limits<Real>::infinity();
    for (i=0; i<num_v; ++i)
      switch (ranVarTypesU[i]) {
      case STD_NORMAL:
	u_bnds[i] = RealRealPair(-dbl_inf, dbl_inf); break;
      case STD_UNIFORM:     case STD_BETA:
	u_bnds[i] = RealRealPair(-1., 1.);           break;
      case STD_EXPONENTIAL: case STD_GAMMA:
	u_bnds[i] = RealRealPair(0., dbl_inf);       break;
      default: // no transformation (e.g., PCE w/ numerically-generated basis)
	check_x_type(i, ranVarTypesU[i]);
	u_bnds[i] = randomVarsX[i].bounds();         break;
      }
    return u_bnds;
  }
}


Real ProbabilityTransformation::u_pdf(Real u_val, size_t i) const
{
  if (probTransRep) return probTransRep->u_pdf(u_val, i);
  else {
    // can only use randomVarsX[i].standard_pdf() for cases where u_type is a
    // standardized form of the x_type.  For STD_NORMAL and STD_UNIFORM, many
    // x_types can be mapped to these u_types, so use global utility fns
    // whenever there are no auxilliary parameters to manage.
    switch (ranVarTypesU[i]) {
    case STD_NORMAL:  return  NormalRandomVariable::std_pdf(u_val); break;
    case STD_UNIFORM: return UniformRandomVariable::std_pdf();      break;
    case STD_EXPONENTIAL:
      return ExponentialRandomVariable::std_pdf(u_val);             break;
    case STD_BETA:
      check_x_type(i, BETA);  return randomVarsX[i].standard_pdf(u_val); break;
    case STD_GAMMA:
      check_x_type(i, GAMMA); return randomVarsX[i].standard_pdf(u_val); break;
    default: // no transformation (e.g., PCE with numerically-generated bases)
      check_x_type(i, ranVarTypesU[i]);
      return randomVarsX[i].pdf(u_val); break;
    }
  }
}


Real ProbabilityTransformation::u_log_pdf(Real u_val, size_t i) const
{
  if (probTransRep) return probTransRep->u_log_pdf(u_val, i);
  else {
    switch (ranVarTypesU[i]) {
    case STD_NORMAL:  return  NormalRandomVariable::log_std_pdf(u_val); break;
    case STD_UNIFORM: return UniformRandomVariable::log_std_pdf();      break;
    case STD_EXPONENTIAL:
      return ExponentialRandomVariable::log_std_pdf(u_val);             break;
    case STD_BETA:
      check_x_type(i, BETA);  // need alphaStat,betaStat for variable i
      return randomVarsX[i].log_standard_pdf(u_val); break;
    case STD_GAMMA:
      check_x_type(i, GAMMA); // need alphaStat for variable i
      return randomVarsX[i].log_standard_pdf(u_val); break;
    default: // no transformation (e.g., PCE with numerically-generated bases)
      check_x_type(i, ranVarTypesU[i]);
      return randomVarsX[i].log_pdf(u_val);          break;
    }
  }
}


Real ProbabilityTransformation::u_log_pdf_gradient(Real u_val, size_t i) const
{
  if (probTransRep) return probTransRep->u_log_pdf_gradient(u_val, i);
  else {
    switch (ranVarTypesU[i]) {
    case STD_NORMAL:
      return      NormalRandomVariable::log_std_pdf_gradient(u_val); break;
    case STD_UNIFORM:
      return     UniformRandomVariable::log_std_pdf_gradient(); break;
    case STD_EXPONENTIAL:
      return ExponentialRandomVariable::log_std_pdf_gradient(); break;
    case STD_BETA:
      check_x_type(i, BETA);  // need alphaStat,betaStat for variable i
      return randomVarsX[i].log_standard_pdf_gradient(u_val); break;
    case STD_GAMMA:
      check_x_type(i, GAMMA); // need alphaStat for variable i
      return randomVarsX[i].log_standard_pdf_gradient(u_val); break;
    default: // no transformation (e.g., PCE with numerically-generated bases)
      check_x_type(i, ranVarTypesU[i]);
      return randomVarsX[i].log_pdf_gradient(u_val);          break;
    }
  }
}


Real ProbabilityTransformation::u_log_pdf_hessian(Real u_val, size_t i) const
{
  if (probTransRep) return probTransRep->u_log_pdf_hessian(u_val, i);
  else {
    switch (ranVarTypesU[i]) {
    case STD_NORMAL:
      return      NormalRandomVariable::log_std_pdf_hessian(); break;
    case STD_UNIFORM:
      return     UniformRandomVariable::log_std_pdf_hessian(); break;
    case STD_EXPONENTIAL:
      return ExponentialRandomVariable::log_std_pdf_hessian(); break;
    case STD_BETA:
      check_x_type(i, BETA);  // need alphaStat,betaStat for variable i
      return randomVarsX[i].log_standard_pdf_hessian(u_val); break;
    case STD_GAMMA:
      check_x_type(i, GAMMA); // need alphaStat for variable i
      return randomVarsX[i].log_standard_pdf_hessian(u_val); break;
    default: // no transformation (e.g., PCE with numerically-generated bases)
      check_x_type(i, ranVarTypesU[i]);
      return randomVarsX[i].log_pdf_hessian(u_val);          break;
    }
  }
}


void ProbabilityTransformation::
trans_U_to_X(const RealVector& u_vars, RealVector& x_vars)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_U_to_X(u_vars, x_vars);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_U_to_X() virtual fn."
	  << "\nNo default defined at ProbabilityTransformation base class.\n"
	  << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_X_to_U(const RealVector& x_vars, RealVector& u_vars)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_X_to_U(x_vars, u_vars);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_X_to_U() virtual fn."
	  << "\nNo default defined at ProbabilityTransformation base class.\n"
	  << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::transform_correlations()
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->transform_correlations();
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine transform_correlations() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_grad_X_to_U(const RealVector& fn_grad_x, RealVector& fn_grad_u,
		  const RealVector& x_vars,    const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_grad_X_to_U(fn_grad_x, fn_grad_u, x_vars, x_dvv,
				    cv_ids);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_grad_X_to_U() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_grad_X_to_U(const RealVector& fn_grad_x,   RealVector& fn_grad_u,
		  const RealMatrix& jacobian_xu, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_grad_X_to_U(fn_grad_x, fn_grad_u, jacobian_xu, x_dvv,
				    cv_ids);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_grad_X_to_U() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_grad_X_to_S(const RealVector& fn_grad_x, RealVector& fn_grad_s,
		  const RealVector& x_vars, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids,
		  SizetMultiArrayConstView acv_ids,
		  const SizetArray& acv_map1_indices,
		  const ShortArray& acv_map2_targets)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_grad_X_to_S(fn_grad_x, fn_grad_s, x_vars, x_dvv, cv_ids,
				    acv_ids, acv_map1_indices,
				    acv_map2_targets);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_grad_X_to_S() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << "class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_grad_X_to_S(const RealVector& fn_grad_x, RealVector& fn_grad_s,
		  const RealMatrix& jacobian_xs, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids,
		  SizetMultiArrayConstView acv_ids,
		  const SizetArray& acv_map1_indices,
		  const ShortArray& acv_map2_targets)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_grad_X_to_S(fn_grad_x, fn_grad_s, jacobian_xs, x_dvv,
				    cv_ids, acv_ids, acv_map1_indices,
				    acv_map2_targets);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_grad_X_to_S() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_grad_U_to_X(const RealVector& fn_grad_u, RealVector& fn_grad_x,
		  const RealVector& x_vars,    const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_grad_U_to_X(fn_grad_u, fn_grad_x, x_vars, x_dvv,
				    cv_ids);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_grad_U_to_X() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_grad_U_to_X(const RealVector& fn_grad_u,   RealVector& fn_grad_x,
		  const RealMatrix& jacobian_ux, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_grad_U_to_X(fn_grad_u, fn_grad_x, jacobian_ux, x_dvv,
				    cv_ids);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_grad_U_to_X() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_hess_X_to_U(const RealSymMatrix& fn_hess_x, RealSymMatrix& fn_hess_u,
		  const RealVector& x_vars, const RealVector& fn_grad_x,
		  const SizetArray& x_dvv, SizetMultiArrayConstView cv_ids)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_hess_X_to_U(fn_hess_x, fn_hess_u, x_vars, fn_grad_x,
				    x_dvv, cv_ids);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_hess_X_to_U() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
trans_hess_X_to_U(const RealSymMatrix& fn_hess_x, RealSymMatrix& fn_hess_u,
		  const RealMatrix& jacobian_xu,
		  const RealSymMatrixArray& hessian_xu,
		  const RealVector& fn_grad_x, const SizetArray& x_dvv,
		  SizetMultiArrayConstView cv_ids)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->trans_hess_X_to_U(fn_hess_x, fn_hess_u, jacobian_xu,
				    hessian_xu, fn_grad_x, x_dvv, cv_ids);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine trans_hess_X_to_U() "
          << "virtual fn.\nNo default defined at ProbabilityTransformation base"
	  << " class.\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
jacobian_dX_dU(const RealVector& x_vars, RealMatrix& jacobian_xu)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->jacobian_dX_dU(x_vars, jacobian_xu);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine jacobian_dX_dU() virtual "
          << "fn.\nNo default defined at ProbabilityTransformation base class."
	  << "\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
jacobian_dU_dX(const RealVector& x_vars, RealMatrix& jacobian_ux)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->jacobian_dU_dX(x_vars, jacobian_ux);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine jacobian_dU_dX() virtual "
          << "fn.\nNo default defined at ProbabilityTransformation base class."
	  << "\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
jacobian_dX_dS(const RealVector& x_vars, RealMatrix& jacobian_xs,
	       SizetMultiArrayConstView cv_ids,
	       SizetMultiArrayConstView acv_ids,
	       const SizetArray& acv_map1_indices,
	       const ShortArray& acv_map2_targets)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->jacobian_dX_dS(x_vars, jacobian_xs, cv_ids, acv_ids,
				 acv_map1_indices, acv_map2_targets);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine jacobian_dX_dS() virtual "
          << "fn.\nNo default defined at ProbabilityTransformation base class."
	  << "\n" << std::endl;
    abort_handler(-1);
  }
}


void ProbabilityTransformation::
hessian_d2X_dU2(const RealVector& x_vars, RealSymMatrixArray& hessian_xu)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->hessian_d2X_dU2(x_vars, hessian_xu);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine hessian_d2X_dU2() virtual "
          << "fn.\nNo default defined at ProbabilityTransformation base class."
	  << "\n" << std::endl;
    abort_handler(-1);
  }
}


/** This procedure computes numerical derivatives of x and/or z with respect to
    distribution parameters s, and is used by jacobian_dX_dS() to provide data
    that is not available analytically.  Numerical dz/ds involves dL/ds
    (z(s) = L(s) u and dz/ds = dL/ds u) and is needed to evaluate dx/ds
    semi-analytically for correlated variables.  Numerical dx/ds is needed for
    distributions lacking simple closed-form CDF expressions (beta and gamma
    distributions). */
void ProbabilityTransformation::
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

      Real s0 = randomVarsX[cv_index].parameter(acv_map2_target);

      // Compute the offset for the ith gradient variable.
      // Enforce a minimum delta of fdgss*.01
      Real h_mag = fd_grad_ss * std::max(std::fabs(s0), .01);
      Real h = (s0 < 0.0) ? -h_mag : h_mag; // h has same sign as s0

      // -----------------------------------
      // Evaluate (L/z_vars/x_vars)_s_plus_h
      // -----------------------------------
      Real s1 = s0 + h;
      // update randomVars & corrCholeskyFactorZ:
      randomVarsX[cv_index].parameter(acv_map2_target, s1);
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
      randomVarsX[cv_index].parameter(acv_map2_target, s1);
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
      randomVarsX[cv_index].parameter(acv_map2_target, s0);
    }
  }
  // reset corrCholeskyFactorZ:
  corrCholeskyFactorZ = chol_z0;
}


#ifdef DERIV_DEBUG
void ProbabilityTransformation::
verify_trans_jacobian_hessian(const RealVector& v0)
{
  size_t i, j, k;
  bool fd_grad_flag = true, fd_hess_flag = true, fd_hess_by_fn_flag = false,
       fd_hess_by_grad_flag = true;

  Real fd_grad_ss = 1.e-8, fd_hess_by_fn_ss = 2.e-8, fd_hess_by_grad_ss = 1.e-8;

  RealVector trans_vars_v0;
  //trans_X_to_U(v0, trans_vars_v0); // v = x, trans_vars_v = u
  trans_U_to_X(v0, trans_vars_v0); // v = u, trans_vars_v = x
  int num_v = v0.Length(), num_tv = trans_vars_v0.Length();

  RealMatrix num_jac_dtv_dv(num_tv, num_v);
  RealSymMatrixArray num_hess_d2tv_dv2(num_tv);
  for (i=0; i<num_tv; i++)
    num_hess_d2tv_dv2[i].Shape(num_v);

  // ------------------------------
  // Estimate numerical derivatives
  // ------------------------------
  if (fd_grad_flag || fd_hess_flag) {
    RealVector v1 = v0; // for perturbed values
    // ---------------
    // Loop over num_v
    // ---------------
    for (j=0; j<num_v; j++) { // difference the 1st num_v vars
      if (fd_grad_flag) {

	// Compute the offset for the ith gradient variable.
	// Enforce a minimum delta of fdgss*.01
	Real h_mag = fd_grad_ss * std::max(std::fabs(v0(j)), .01);
	Real h = (v1(j) < 0.0) ? -h_mag : h_mag; // h has same sign as v1(j)

	// ----------------------------
	// Evaluate trans_vars_v_plus_h
	// ----------------------------
	RealVector trans_vars_v_plus_h;
	v1(j) = v0(j) + h;
	PCout << ">>>>> Pecos finite difference gradient evaluation for v["
	      << j+1 << "] + h:\n";
	//trans_X_to_U(v1, trans_vars_v_plus_h);
	trans_U_to_X(v1, trans_vars_v_plus_h);

	// -----------------------------
	// Evaluate trans_vars_v_minus_h
	// -----------------------------
	RealVector trans_vars_v_minus_h;
	v1(j) = v0(j) - h;
	PCout << ">>>>> Pecos finite difference gradient evaluation for v["
	      << j+1 << "] - h:\n";
	//trans_X_to_U(v1, trans_vars_v_minus_h);
	trans_U_to_X(v1, trans_vars_v_minus_h);

	// always use central diffs for verification purposes
	for (i=0; i<num_tv; i++)
	  num_jac_dtv_dv(i,j)
	    = (trans_vars_v_plus_h(i) - trans_vars_v_minus_h(i))/2./h;
      }

      if (fd_hess_flag) {

	if (fd_hess_by_fn_flag) {
	  RealVector trans_vars_v_plus_2h, trans_vars_v_minus_2h;

	  // Compute the 2nd-order Hessian offset for the ith variable.
	  // Enforce a minimum delta of fdhss*.01
	  Real h_mag = fd_hess_by_fn_ss * std::max(std::fabs(v0(j)), .01);
	  Real h = (v1(j) < 0.0) ? -h_mag : h_mag; // h has same sign as v1(j)

	  // evaluate diagonal term

	  // -----------------------------
	  // Evaluate trans_vars_v_plus_2h
	  // -----------------------------
	  v1(j) = v0(j) + 2.*h;
	  PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
	        << j+1 << "] + 2h:\n";
	  //trans_X_to_U(v1, trans_vars_v_plus_2h);
	  trans_U_to_X(v1, trans_vars_v_plus_2h);

	  // ------------------------------
	  // Evaluate trans_vars_v_minus_2h
	  // ------------------------------
	  v1(j) = v0(j) - 2.*h;
	  PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
	        << j+1 << "] - 2h:\n";
	  //trans_X_to_U(v1, trans_vars_v_minus_2h);
	  trans_U_to_X(v1, trans_vars_v_minus_2h);

	  for (i=0; i<num_tv; i++)
	    num_hess_d2tv_dv2[i](j,j)
	      = (trans_vars_v_plus_2h(i) - 2.*trans_vars_v0(i) +
		 trans_vars_v_minus_2h(i))/(4.*h*h);

	  // evaluate off-diagonal terms

	  for (k=j+1; k<num_v; k++) {
	    RealVector trans_vars_v_plus_h_plus_h,
	      trans_vars_v_plus_h_minus_h, trans_vars_v_minus_h_plus_h,
	      trans_vars_v_minus_h_minus_h;

	    // -----------------------------------
	    // Evaluate trans_vars_v_plus_h_plus_h
	    // -----------------------------------
	    v1(j) = v0(j) + h;
	    v1(k) = v0(k) + h;
	    PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
		  << j+1 << "] + h, v[" << k+1 << "] + h:\n";
	    //trans_X_to_U(v1, trans_vars_v_plus_h_plus_h);
	    trans_U_to_X(v1, trans_vars_v_plus_h_plus_h);
	    // ------------------------------------
	    // Evaluate trans_vars_v_plus_h_minus_h
	    // ------------------------------------
	    //v1(j) = v0(j) + h;
	    v1(k) = v0(k) - h;
	    PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
		  << j+1 << "] + h, v[" << k+1 << "] - h:\n";
	    //trans_X_to_U(v1, trans_vars_v_plus_h_minus_h);
	    trans_U_to_X(v1, trans_vars_v_plus_h_minus_h);
	    // ------------------------------------
	    // Evaluate trans_vars_v_minus_h_plus_h
	    // ------------------------------------
	    v1(j) = v0(j) - h;
	    v1(k) = v0(k) + h;
	    PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
		  << j+1 << "] - h, v[" << k+1 << "] + h:\n";
	    //trans_X_to_U(v1, trans_vars_v_minus_h_plus_h);
	    trans_U_to_X(v1, trans_vars_v_minus_h_plus_h);
	    // -------------------------------------
	    // Evaluate trans_vars_v_minus_h_minus_h
	    // -------------------------------------
	    //v1(j) = v0(j) - h;
	    v1(k) = v0(k) - h;
	    PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
		  << j+1 << "] - h, v[" << k+1 << "] - h:\n";
	    //trans_X_to_U(v1, trans_vars_v_minus_h_minus_h);
	    trans_U_to_X(v1, trans_vars_v_minus_h_minus_h);

	    for (i=0; i<num_tv; i++)
	      num_hess_d2tv_dv2[i](j,k) = num_hess_d2tv_dv2[i](k,j)
		= (trans_vars_v_plus_h_plus_h(i)
		- trans_vars_v_plus_h_minus_h(i)
		- trans_vars_v_minus_h_plus_h(i)
		+ trans_vars_v_minus_h_minus_h(i)) / (4.*h*h);

	    v1(k) = v0(k);
	  }
	}

	if (fd_hess_by_grad_flag) {

	  // Compute the 1st-order Hessian offset for the ith variable.
	  // Enforce a minimum delta of fdhss*.01
	  Real h_mag = fd_hess_by_grad_ss * std::max(std::fabs(v0(j)), .01);
	  Real h = (v1(j) < 0.0) ? -h_mag : h_mag; // h has same sign as v1(j)

	  // --------------------------
	  // Evaluate fn_grads_v_plus_h
	  // --------------------------
	  v1(j) = v0(j) + h;
	  PCout << ">>>>> Pecos finite difference Hessian evaluation for v["
	        << j+1 << "] + h:\n";
	  RealVector trans_vars_v_plus_h;
	  trans_U_to_X(v1, trans_vars_v_plus_h);
	  RealMatrix jac_v0, jac_v_plus_h;
	  // jacobian routines use x_vars:
	  jacobian_dX_dU(trans_vars_v0,       jac_v0);
	  jacobian_dX_dU(trans_vars_v_plus_h, jac_v_plus_h);
	  for (i=0; i<num_tv; i++)
	    for (k=0; k<num_v; k++)
	      num_hess_d2tv_dv2[i](j,k)	= (jac_v_plus_h(i,k) - jac_v0(i,k))/h;
	}
      }
      v1(j) = v0(j);
    }
  }

  // Enforce symmetry in the case of FD Hessians from 1st-order gradient
  // differences by averaging off-diagonal terms: H' = 1/2 (H + H^T)
  if (fd_hess_by_grad_flag)
    for (i=0; i<num_tv; i++)
      for (j=0; j<num_v; j++)
	for (k=j+1; k<num_v; k++)
	  num_hess_d2tv_dv2[i](j,k) = num_hess_d2tv_dv2[i](k,j)
	    = (num_hess_d2tv_dv2[i](j,k) + num_hess_d2tv_dv2[i](k,j))/2.;

  // Print out numerical and analytic:
  RealVector x0(num_tv);
  //x0 = v0;
  //RealMatrix jacobian_ux;
  //jacobian_dU_dX(x0, jacobian_ux);
  trans_U_to_X(v0, x0);
  RealMatrix jacobian_xu;
  jacobian_dX_dU(x0, jacobian_xu);
  PCout << "\nNumerical jacobian:" << num_jac_dtv_dv
        << "\nAnalytic jacobian:"  << jacobian_xu; //jacobian_ux;
  RealSymMatrixArray hessian_xu(num_tv);
  hessian_d2X_dU2(x0, hessian_xu);
  for (i=0; i<num_tv; i++)
    PCout << "\nNumerical Hessian:" << num_hess_d2tv_dv2[i]
	  << "\nAnalytic Hessian:"  << hessian_xu[i];
}


void ProbabilityTransformation::verify_design_jacobian(const RealVector& u0)
{
  RealVector x0;
  trans_U_to_X(u0, x0);

  RealMatrix num_jac_dx_ds, num_jac_dz_ds;
  numerical_design_jacobian(x0, true, num_jac_dx_ds, false, num_jac_dz_ds);

  RealMatrix jacobian_xs;
  jacobian_dX_dS(x0, jacobian_xs);

  // Print out numerical and analytic:
  PCout << "\nNumerical jacobian:" << num_jac_dx_ds
        << "\nAnalytic jacobian:"  << jacobian_xs;
}
#endif // DERIV_DEBUG

} // namespace Pecos
