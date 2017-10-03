/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "LHSDriver.hpp"
#include "BoostRNG_Monostate.hpp"
#include "LognormalRandomVariable.hpp"
#include <boost/lexical_cast.hpp>

static const char rcsId[]="@(#) $Id: LHSDriver.cpp 5248 2008-09-05 18:51:52Z wjbohnh $";

//#define DEBUG

#ifdef HAVE_LHS
#ifdef HAVE_CONFIG_H
// Use the classic, autotools Fortran name mangling macros in pecos_config.h
#define LHS_INIT_MEM_FC FC_FUNC_(lhs_init_mem,LHS_INIT_MEM)
#define LHS_PREP_FC     FC_FUNC_(lhs_prep,LHS_PREP)
#define LHS_RUN_FC      FC_FUNC_(lhs_run,LHS_RUN)
#define LHS_CLOSE_FC    FC_FUNC_(lhs_close,LHS_CLOSE)
#define LHS_OPTIONS2_FC FC_FUNC_(lhs_options2,LHS_OPTIONS2)
#define LHS_DIST2_FC    FC_FUNC_(lhs_dist2,LHS_DIST2)
#define LHS_UDIST2_FC   FC_FUNC_(lhs_udist2,LHS_UDIST2)
#define LHS_CONST2_FC   FC_FUNC_(lhs_const2,LHS_CONST2)
#define LHS_CORR2_FC    FC_FUNC_(lhs_corr2,LHS_CORR2)
#define LHS_FILES2_FC   FC_FUNC_(lhs_files2,LHS_FILES2)

#define defaultrnum1       FC_FUNC(defaultrnum1,DEFAULTRNUM1)
#define defaultrnum2       FC_FUNC(defaultrnum2,DEFAULTRNUM2)

#else
// Use the CMake-generated PREFIXED, fortran name mangling macros (no warnings)
#define LHS_INIT_MEM_FC LHS_GLOBAL_(lhs_init_mem,LHS_INIT_MEM)
#define LHS_PREP_FC     LHS_GLOBAL_(lhs_prep,LHS_PREP)
#define LHS_RUN_FC      LHS_GLOBAL_(lhs_run,LHS_RUN)
#define LHS_CLOSE_FC    LHS_GLOBAL_(lhs_close,LHS_CLOSE)
//#define LHS_OPTIONS2_FC LHS_GLOBAL_(lhs_options2,LHS_OPTIONS2)
//#define LHS_DIST2_FC    LHS_GLOBAL_(lhs_dist2,LHS_DIST2)
//#define LHS_UDIST2_FC   LHS_GLOBAL_(lhs_udist2,LHS_UDIST2)
//#define LHS_CONST2_FC   LHS_GLOBAL_(lhs_const2,LHS_CONST2)
//#define LHS_CORR2_FC    LHS_GLOBAL_(lhs_corr2,LHS_CORR2)
//#define LHS_FILES2_FC   LHS_GLOBAL_(lhs_files2,LHS_FILES2)
#define defaultrnum1    LHS_GLOBAL(defaultrnum1,DEFAULTRNUM1)
#define defaultrnum2    LHS_GLOBAL(defaultrnum2,DEFAULTRNUM2)

// BMA (20160315): Changed to use Fortran 2003 ISO C bindings.
// The Fortran symbol will be lowercase with same name as if in C
#define  LHS_OPTIONS2_FC lhs_options2
#define  LHS_DIST2_FC lhs_dist2
#define  LHS_UDIST2_FC lhs_udist2
#define  LHS_CONST2_FC lhs_const2
#define  LHS_CORR2_FC lhs_corr2
#define  LHS_FILES2_FC lhs_files2


#endif // HAVE_CONFIG_H

extern "C" {

// for these functions, a straight interface to F90 can be used
void LHS_INIT_MEM_FC( int& obs, int& seed, int& max_obs, int& max_samp_size,
		      int& max_var, int& max_interval, int& max_corr,
		      int& max_table, int& print_level, int& output_width,
		      int& ierror );

void LHS_PREP_FC( int& ierror, int& num_names, int& num_vars );

void LHS_RUN_FC( int& max_var, int& max_obs, int& max_names, int& ierror,
		 char* dist_names, int* name_order, Pecos::Real* ptvals,
		 int& num_names, Pecos::Real* sample_matrix, int& num_vars,
                 Pecos::Real* rank_matrix, int& rank_flag );

void LHS_CLOSE_FC( int& ierror );

// for these functions, bridge-function interfaces to F90 are required
void LHS_OPTIONS2_FC( int& num_replications, int& ptval_option,
		      const char* sampling_options, int& ierror );

void LHS_DIST2_FC( const char* label, int& ptval_flag, Pecos::Real& ptval,
		   const char* dist_type, Pecos::Real* dist_params,
		   int& num_params, int& ierror, int& dist_id, int& ptval_id );

void LHS_UDIST2_FC( const char* label, int& ptval_flag, Pecos::Real& ptval,
		    const char* dist_type, int& num_pts, Pecos::Real* x,
		    Pecos::Real* y, int& ierror, int& dist_id, int& ptval_id );

void LHS_CONST2_FC( const char* label, Pecos::Real& ptval, int& ierror,
		    int& ptval_id );

void LHS_CORR2_FC( char* label1, char* label2, Pecos::Real& corr, int& ierror );

void LHS_FILES2_FC( const char* lhsout, const char* lhsmsg, const char* lhstitl,
		    const char* lhsopts, int& ierror );

//void lhs_run2( int* max_var, int* max_obs, int* max_names, int* ierror,
//               char* dist_names, int* name_order, double* ptvals,
//               int* num_names, double* sample_matrix, int* num_vars );

Pecos::Real defaultrnum1( void );
Pecos::Real defaultrnum2( void );

//Pecos::Real rnumlhs10( void );
//Pecos::Real rnumlhs20( void );

}
#endif // HAVE_LHS


namespace Pecos {

/** Helper function to create labels that Fortran will see as
    character*16 values which are NOT null-terminated.  For convenience,
    the StringArray, lhs_names, is also populated. */
static const char* f77name16(const String& name, size_t index,
                             StringArray& lhs_names)
{
  String label = name + boost::lexical_cast<String>(index+1);
  label.resize(16, ' ');
  lhs_names[index] = label;
  return lhs_names[index].data();  // NOTE: no NULL terminator
}


void LHSDriver::seed(int seed)
{
  randomSeed = seed;
  // The Boost RNG is not set by LHS_INIT_MEM, so must be done here.
  if (BoostRNG_Monostate::randomNum == BoostRNG_Monostate::mt19937)
    BoostRNG_Monostate::seed(seed);
  // This would be redundant since the f77 ISeed is set in LHS_INIT_MEM:
  //else
  // lhs_setseed(&seed);
}


void LHSDriver::rng(String unif_gen)
{

  // check the environment (once only) for an RNG preference and cache it
  static bool first_entry = true;
  static const char *env_unifgen;
  if (first_entry) {
    env_unifgen = std::getenv("DAKOTA_LHS_UNIFGEN");
    first_entry = false;
  }

  // the environment overrides the passed rng specification
  if (env_unifgen) {
    unif_gen = env_unifgen;
    if (unif_gen != "rnum2" && unif_gen != "mt19937") {
      PCerr << "Error: LHSDriver::rng() expected $DAKOTA_LHS_UNIFGEN to be "
	    << "\"rnum2\" or \"mt19937\", not \"" << env_unifgen << "\".\n"
	    << std::endl;
      abort_handler(-1);
    }
  }

  // now point the monostate RNG to the desired generator function
  if (unif_gen == "mt19937" || unif_gen.empty()) {
    BoostRNG_Monostate::randomNum  = BoostRNG_Monostate::mt19937;
    BoostRNG_Monostate::randomNum2 = BoostRNG_Monostate::mt19937;
    allowSeedAdvance &= ~2; // drop 2 bit: disallow repeated seed update
  }
  else if (unif_gen == "rnum2") {
#ifdef HAVE_LHS
    BoostRNG_Monostate::randomNum  = (Rfunc)defaultrnum1;
    BoostRNG_Monostate::randomNum2 = (Rfunc)defaultrnum2;
    allowSeedAdvance |= 2;  // add 2 bit: allow repeated seed update
#else
    PCerr << "Error: LHSDriver::rng() Use of rnum2 for RNG selection is NOT "
	  << "supported in current (without-lhs) configuration" << std::endl;
    abort_handler(-1);
#endif
  }
  else {
    PCerr << "Error: LHSDriver::rng() expected string to be \"rnum2\" or "
	  << "\"mt19937\", not \"" << unif_gen << "\".\n" << std::endl;
    abort_handler(-1);
  }
}


void LHSDriver::abort_if_no_lhs()
{
#ifndef HAVE_LHS
  PCerr << "Error: LHSDriver not available as PECOS was configured without LHS."
        << std::endl;
  abort_handler(-1);
#endif
}


//String LHSDriver::rng()
//{
//  if (BoostRNG_Monostate::randomNum == BoostRNG_Monostate::random_num1)
//    return "mt19937";
//  if (BoostRNG_Monostate::randomNum == (Rfunc)rnumlhs10)
//    return "rnum2";
//  else
//    return "";
//}


/** While it would be desirable in some cases to carve this function
    into smaller parts and allow multiple invocations of LHS_RUN
    following a single initialization of types and arrays, the LHS
    code does not currently allow this: it will return an error if
    LHS_INIT/LHS_INIT_MEM, at least one distribution call (i.e.,
    LHS_DIST, LHS_UDIST or LHS_SDIST), and LHS_PREP are not called
    prior to each invocation of LHS_RUN.  Since LHS_INIT/LHS_INIT_MEM
    require input of a seed, the approach to computing multiple
    distinct sample sets must employ advance_seed_sequence() to
    re-seed multiple generate_samples() calls, rather than continuing
    an existing random number sequence. */
void LHSDriver::
generate_samples(const RealVector& cd_l_bnds,   const RealVector& cd_u_bnds,
		 const IntVector&  ddri_l_bnds, const IntVector&  ddri_u_bnds,
		 const IntSetArray& ddsi_values,
		 const StringSetArray& ddss_values,
		 const RealSetArray& ddsr_values,
		 const RealVector& cs_l_bnds,   const RealVector& cs_u_bnds,
		 const IntVector&  dsri_l_bnds, const IntVector&  dsri_u_bnds,
		 const IntSetArray& dssi_values,
		 const StringSetArray& dsss_values,
		 const RealSetArray& dssr_values, const AleatoryDistParams& adp,
		 const EpistemicDistParams& edp, int num_samples,
		 RealMatrix& samples, RealMatrix& sample_ranks)
{
#ifdef HAVE_LHS
  // generate samples within user-specified parameter distributions

  // error check on program parameters
  if (!num_samples) {
    PCerr << "\nError: number of samples in LHSDriver::generate_samples() must "
	  << "be nonzero." << std::endl;
    abort_handler(-1);
  }

  const RealVector& n_means = adp.normal_means();
  const RealVector& n_std_devs = adp.normal_std_deviations();
  const RealVector& n_l_bnds = adp.normal_lower_bounds();
  const RealVector& n_u_bnds = adp.normal_upper_bounds();
  const RealVector& ln_means = adp.lognormal_means();
  const RealVector& ln_std_devs = adp.lognormal_std_deviations();
  const RealVector& ln_lambdas = adp.lognormal_lambdas();
  const RealVector& ln_zetas = adp.lognormal_zetas();
  const RealVector& ln_err_facts = adp.lognormal_error_factors();
  const RealVector& ln_l_bnds = adp.lognormal_lower_bounds();
  const RealVector& ln_u_bnds = adp.lognormal_upper_bounds();
  const RealVector& u_l_bnds = adp.uniform_lower_bounds();
  const RealVector& u_u_bnds = adp.uniform_upper_bounds();
  const RealVector& lu_l_bnds = adp.loguniform_lower_bounds();
  const RealVector& lu_u_bnds = adp.loguniform_upper_bounds();
  const RealVector& t_modes = adp.triangular_modes();
  const RealVector& t_l_bnds = adp.triangular_lower_bounds();
  const RealVector& t_u_bnds = adp.triangular_upper_bounds();
  const RealVector& e_betas = adp.exponential_betas();
  const RealVector& b_alphas = adp.beta_alphas();
  const RealVector& b_betas = adp.beta_betas();
  const RealVector& b_l_bnds = adp.beta_lower_bounds();
  const RealVector& b_u_bnds = adp.beta_upper_bounds();
  const RealVector& ga_alphas = adp.gamma_alphas();
  const RealVector& ga_betas = adp.gamma_betas();
  const RealVector& gu_alphas = adp.gumbel_alphas();
  const RealVector& gu_betas = adp.gumbel_betas();
  const RealVector& f_alphas = adp.frechet_alphas();
  const RealVector& f_betas  = adp.frechet_betas();
  const RealVector& w_alphas = adp.weibull_alphas();
  const RealVector& w_betas = adp.weibull_betas();
  const RealRealMapArray& h_bin_prs = adp.histogram_bin_pairs();
  const RealVector& p_lambdas = adp.poisson_lambdas();
  const RealVector& bi_prob_per_tr = adp.binomial_probability_per_trial();
  const IntVector&  bi_num_tr = adp.binomial_num_trials();
  const RealVector& nb_prob_per_tr
    = adp.negative_binomial_probability_per_trial();
  const IntVector&  nb_num_tr = adp.negative_binomial_num_trials();
  const RealVector& ge_prob_per_tr = adp.geometric_probability_per_trial();
  const IntVector&  hg_total_pop = adp.hypergeometric_total_population();
  const IntVector&  hg_selected_pop = adp.hypergeometric_selected_population();
  const IntVector&  hg_num_drawn = adp.hypergeometric_num_drawn();
  const IntRealMapArray& h_pt_int_prs = adp.histogram_point_int_pairs();
  const StringRealMapArray& h_pt_string_prs
    = adp.histogram_point_string_pairs();
  const RealRealMapArray& h_pt_real_prs = adp.histogram_point_real_pairs();
  const RealSymMatrix& correlations = adp.uncertain_correlations();

  const RealRealPairRealMapArray& ci_bpa
    = edp.continuous_interval_basic_probabilities();
  const IntIntPairRealMapArray&   di_bpa
    = edp.discrete_interval_basic_probabilities();
  const IntRealMapArray& dusi_vals_probs
    = edp.discrete_set_int_values_probabilities();
  const StringRealMapArray& duss_vals_probs
    = edp.discrete_set_string_values_probabilities();
  const RealRealMapArray& dusr_vals_probs
    = edp.discrete_set_real_values_probabilities();

  bool correlation_flag = !correlations.empty();
  size_t i, j, num_cdv = cd_l_bnds.length(), num_ddriv = ddri_l_bnds.length(),
    num_ddsiv = ddsi_values.size(), num_ddssv = ddss_values.size(),
    num_ddsrv = ddsr_values.size(), num_nuv  = n_means.length(),
    num_lnuv  = std::max(ln_means.length(), ln_lambdas.length()),
    num_uuv   = u_l_bnds.length(),  num_luuv = lu_l_bnds.length(),
    num_tuv   = t_modes.length(),   num_exuv = e_betas.length(),
    num_beuv  = b_alphas.length(),  num_gauv = ga_alphas.length(),
    num_guuv  = gu_alphas.length(), num_fuv  = f_alphas.length(),
    num_wuv   = w_alphas.length(),  num_hbuv = h_bin_prs.size(),
    num_puv   = p_lambdas.length(),      num_biuv  = bi_prob_per_tr.length(),
    num_nbuv  = nb_prob_per_tr.length(), num_geuv  = ge_prob_per_tr.length(),
    num_hguv  = hg_num_drawn.length(),   num_hpuiv = h_pt_int_prs.size(),
    num_hpusv = h_pt_string_prs.size(),  num_hpurv = h_pt_real_prs.size(),
    num_ciuv  = ci_bpa.size(), num_diuv = di_bpa.size(),
    num_dusiv = dusi_vals_probs.size(),  num_dussv = duss_vals_probs.size(),
    num_dusrv = dusr_vals_probs.size(),  num_csv   = cs_l_bnds.length(),
    num_dsriv = dsri_l_bnds.length(),    num_dssiv = dssi_values.size(),
    num_dsssv = dsss_values.size(), num_dssrv = dssr_values.size(),
    num_dv   = num_cdv  + num_ddriv + num_ddsiv + num_ddssv + num_ddsrv,
    num_cauv = num_nuv  + num_lnuv + num_uuv  + num_luuv + num_tuv + num_exuv
             + num_beuv + num_gauv + num_guuv + num_fuv  + num_wuv + num_hbuv,
    num_dauiv = num_puv + num_biuv + num_nbuv + num_geuv + num_hguv + num_hpuiv,
    num_dausv = num_hpusv, num_daurv = num_hpurv,
    num_auv = num_cauv + num_dauiv + num_dausv + num_daurv,
    num_ceuv  = num_ciuv, num_deuiv = num_diuv + num_dusiv,
    num_deusv = num_dussv, num_deurv = num_dusrv, 
    num_euv = num_ceuv + num_deuiv + num_deusv + num_deurv,
    num_uv = num_auv + num_euv,
    num_sv = num_csv + num_dsriv + num_dssiv + num_dsssv + num_dssrv,
    num_av = num_dv  + num_uv    + num_sv;

  int err_code = 0, max_var = num_av, max_obs = num_samples,
      max_samp_size = num_av*num_samples, max_interval = -1,
      max_unc_corr = (num_uv*num_uv - num_uv)/2,
      max_table = -1, print_level = 0, output_width = 1;
  int max_corr = (num_uv > 1) ? max_unc_corr : -1;
  // randomSeed passed below propagates to ISeed in the f77 rnum2, but does
  // not propagate to Boost RNGs (LHSDriver::seed() must be used for that).
  LHS_INIT_MEM_FC(num_samples, randomSeed, max_obs, max_samp_size, max_var,
		  max_interval, max_corr, max_table, print_level, output_width,
		  err_code);
  check_error(err_code, "lhs_init_mem");

  // set sample type to either LHS (default) or random Monte Carlo (optional)
  bool call_lhs_option = false;
  String option_string("              ");
  if (sampleType == "random" || sampleType == "incremental_random") {
    option_string   = "RANDOM SAMPLE ";
    call_lhs_option = true;
  }
  // set mixing option to either restricted pairing (default) or random pairing
  // (optional).  For enforcing user-specified correlation, restricted pairing
  // is required.  And for uncorrelated variables, restricted pairing results
  // in lower correlation values than random pairing.  For these reasons, the
  // random pairing option is not currently active, although a specification
  // option for it could be added in the future if a use arises.
  bool random_pairing_flag = false; // this option hard-wired off for now
  if (!correlation_flag && random_pairing_flag) {
    option_string += "RANDOM PAIRING";
    call_lhs_option = true;
  }
  // else // use default of restricted pairing
  option_string.resize(32, ' ');
  if (call_lhs_option) {
    // Don't null-terminate the string since the '\0' is not used in Fortran
    int num_replic = 1, ptval_option = 1;
    LHS_OPTIONS2_FC(num_replic, ptval_option, option_string.data(), err_code);
    check_error(err_code, "lhs_options");
  }

  // Create files showing distributions and associated statistics.  Avoid
  // issues with null-terminated strings from C++ (which mess up the Fortran
  // output) by using std::string::data().
  String output_string("LHS_samples.out");
  output_string.resize(32, ' ');
  String message_string("LHS_distributions.out");
  message_string.resize(32, ' ');
  String title_string("Pecos::LHSDriver");
  title_string.resize(32, ' ');
  // From the LHS manual (p. 100): LHSRPTS is used to specify which reports LHS
  // will print in the message output file. If LHSRPTS is omitted, the message
  // file will contain only the title, run specifications, and descriptions of
  // the distributions sampled. If LHSRPTS is included, it must be followed by
  // one or more of the following three additional keywords:
  //   > CORR: Print both the achieved raw correlation matrix and the achieved
  //           rank correlation matrix.
  //   > HIST: Print a text-based histogram for each random variable.
  //   > DATA: Print the complete set of all data samples and their ranks.
  // Pecos::LHSDriver::reportFlag is set from Dakota::Iterator::subIteratorFlag,
  // which accomplishes two things: (1) it reduces some output when generating
  // multiple sample sets (the report files get overwritten anyway), and (2) it
  // avoids numerical problems with generating input variable histogram plots
  // as trust regions become small in SBO (mainly an issue before conversion of
  // f90 LHS to double precision).
  String options_string = (reportFlag) ? "LHSRPTS CORR HIST DATA" : " ";
  options_string.resize(32, ' ');
  LHS_FILES2_FC(output_string.data(), message_string.data(),
                title_string.data(), options_string.data(), err_code);
  check_error(err_code, "lhs_files");

  int num_params, cntr = 0, ptval_flag = 0;
  int dist_num, pv_num; // outputs (ignored)
  Real ptval = 0., dist_params[4];
  StringArray lhs_names(num_av);
  const char *name_string, *distname;
  Real dbl_inf = std::numeric_limits<Real>::infinity();

  // --------------------
  // CONTINUOUS VARIABLES
  // --------------------
  // continuous design (treated as uniform)
  for (i=0; i<num_cdv; ++i, ++cntr) {
    name_string = f77name16("ContDesign", cntr, lhs_names);
    String dist_string("uniform");
    dist_string.resize(32, ' ');
    num_params = 2;
    Real l_bnd = cd_l_bnds[i], u_bnd = cd_u_bnds[i];
    if (l_bnd > -dbl_inf && u_bnd < dbl_inf) {
      if (l_bnd >= u_bnd) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly "
	      << "less than upper bounds to\n       sample continuous design "
	      << "variables using uniform distributions." << std::endl;
	abort_handler(-1);
      }
      dist_params[0] = l_bnd;
      dist_params[1] = u_bnd;
    }
    else {
      PCerr << "\nError: Pecos::LHSDriver requires bounds to sample design "
	    << "variables using uniform\n       distributions." << std::endl;
      abort_handler(-1);
    }
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(continuous design)");
  }

  // normal uncertain
  bool n_bnd_spec = (!n_l_bnds.empty() && !n_u_bnds.empty());
  for (i=0; i<num_nuv; ++i, ++cntr) {
    name_string = f77name16("Normal", cntr, lhs_names);
    dist_params[0] = n_means[i];
    dist_params[1] = n_std_devs[i];
    // check for bounded normal
    if (n_bnd_spec && (n_l_bnds[i] > -dbl_inf || n_u_bnds[i] < dbl_inf) ) {
      if (n_l_bnds[i] >= n_u_bnds[i]) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly "
	      << "less than upper bounds to\n       sample using bounded "
	      << "normal distributions." << std::endl;
	abort_handler(-1);
      }
      num_params = 4;
      dist_params[2] = n_l_bnds[i];
      dist_params[3] = n_u_bnds[i];
      distname = "bounded normal";
    }
    else { // normal
      num_params = 2;
      distname = "normal";
    }
    String dist_string(distname);
    dist_string.resize(32, ' ');
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(normal)");
  }

  // lognormal uncertain
  bool ln_bnd_spec = (!ln_l_bnds.empty()  && !ln_u_bnds.empty());
  bool n_dist      = (!ln_lambdas.empty() || !ln_std_devs.empty());
  for (i=0; i<num_lnuv; ++i, ++cntr) {
    name_string = f77name16("Lognormal", cntr, lhs_names);
    if (n_dist) {
      // In the mean/std dev specification case, LHS expects the mean/std dev
      // of the underlying normal distribution (LHS manual, SAND#98-0210, p.39).
      // Therefore, map from the DAKOTA spec (mean/std_dev or lambda/zeta) to
      // the LHS input requirements (lambda/zeta), if required.
      if (!ln_lambdas.empty()) {
	dist_params[0] = ln_lambdas[i];
	dist_params[1] = ln_zetas[i];
      }
      else
	LognormalRandomVariable::params_from_moments(ln_means[i],
	  ln_std_devs[i], dist_params[0], dist_params[1]);
    }
    else {
      // In the mean/err factor specification case, DAKOTA and LHS are
      // consistent (LHS manual, SAND#98-0210, p.39) and no mapping is required.
      dist_params[0] = ln_means[i];
      dist_params[1] = ln_err_facts[i];
    }
    // check for bounded lognormal
    if (ln_bnd_spec && (ln_l_bnds[i] > 0. || ln_u_bnds[i] < dbl_inf) ) {
      if (ln_l_bnds[i] >= ln_u_bnds[i]) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly "
	      << "less than upper bounds to\n       sample using bounded "
	      << "lognormal distributions." << std::endl;
	abort_handler(-1);
      }
      num_params = 4;
      dist_params[2] = ln_l_bnds[i];
      dist_params[3] = ln_u_bnds[i];
      if (n_dist)
	distname = "bounded lognormal-n";
      else
	distname = "bounded lognormal";
    }
    else {
      num_params = 2;
      if (n_dist)
	distname = "lognormal-n";
      else
	distname = "lognormal";
    }
    String dist_string(distname);
    dist_string.resize(32, ' ');
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(lognormal)");
  }

  // uniform uncertain
  for (i=0; i<num_uuv; ++i, ++cntr) {
    name_string = f77name16("Uniform", cntr, lhs_names);
    String dist_string("uniform");
    dist_string.resize(32, ' ');
    num_params = 2;
    if (u_l_bnds[i] >= u_u_bnds[i]) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly "
	      << "less than upper bounds to\n       sample using uniform "
	      << "distributions." << std::endl;
	abort_handler(-1);
    }
    dist_params[0] = u_l_bnds[i];
    dist_params[1] = u_u_bnds[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(uniform)");
  }

  // loguniform uncertain
  for (i=0; i<num_luuv; ++i, ++cntr) {
    name_string = f77name16("Loguniform", cntr, lhs_names);
    String dist_string("loguniform");
    dist_string.resize(32, ' ');
    num_params = 2;
    if (lu_l_bnds[i] >= lu_u_bnds[i]) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly "
	      << "less than upper bounds to\n       sample using loguniform "
	      << "distributions." << std::endl;
	abort_handler(-1);
    }
    dist_params[0] = lu_l_bnds[i];
    dist_params[1] = lu_u_bnds[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(loguniform)");
  }

  // triangular uncertain
  for (i=0; i<num_tuv; ++i, ++cntr) {
    name_string = f77name16("Triangular", cntr, lhs_names);
    String dist_string("triangular");
    dist_string.resize(32, ' ');
    num_params = 3;
    if (t_l_bnds[i] >= t_u_bnds[i]) {
      PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly less "
	    << "than upper bounds to\n       sample using triangular "
	    << "distributions." << std::endl;
      abort_handler(-1);
    }
    dist_params[0] = t_l_bnds[i];
    dist_params[1] = t_modes[i];
    dist_params[2] = t_u_bnds[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(triangular)");
  }

  // exponential uncertain
  for (i=0; i<num_exuv; ++i, ++cntr) {
    name_string = f77name16("Exponential", cntr, lhs_names);
    String dist_string("exponential");
    dist_string.resize(32, ' ');
    num_params = 1;
    dist_params[0] = 1./e_betas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(exponential)");
  }

  // beta uncertain
  for (i=0; i<num_beuv; ++i, ++cntr) {
    name_string = f77name16("Beta", cntr, lhs_names);
    String dist_string("beta");
    dist_string.resize(32, ' ');
    num_params = 4;
    if (b_l_bnds[i] >= b_u_bnds[i]) {
      PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly less "
	    << "than upper bounds to\n       sample using beta "
	    << "distributions." << std::endl;
      abort_handler(-1);
    }
    dist_params[0] = b_l_bnds[i];
    dist_params[1] = b_u_bnds[i];
    dist_params[2] = b_alphas[i];
    dist_params[3] = b_betas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(beta)");
  }

  // gamma uncertain
  for (i=0; i<num_gauv; ++i, ++cntr) {
    name_string = f77name16("Gamma", cntr, lhs_names);
    String dist_string("gamma");
    dist_string.resize(32, ' ');
    num_params = 2;
    dist_params[0] = ga_alphas[i];
    dist_params[1] = 1./ga_betas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(gamma)");
  }

  // gumbel uncertain
  for (i=0; i<num_guuv; ++i, ++cntr) {
    name_string = f77name16("Gumbel", cntr, lhs_names);
    String dist_string("gumbel");
    dist_string.resize(32, ' ');
    num_params = 2;
    dist_params[0] = gu_alphas[i];
    dist_params[1] = gu_betas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(gumbel)");
  }

  // frechet uncertain
  for (i=0; i<num_fuv; ++i, ++cntr) {
    name_string = f77name16("Frechet", cntr, lhs_names);
    String dist_string("frechet");
    dist_string.resize(32, ' ');
    num_params = 2;
    dist_params[0] = f_alphas[i];
    dist_params[1] = f_betas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(frechet)");
  }

  // weibull uncertain
  for (i=0; i<num_wuv; ++i, ++cntr) {
    name_string = f77name16("Weibull", cntr, lhs_names);
    String dist_string("weibull");
    dist_string.resize(32, ' ');
    num_params = 2;
    dist_params[0] = w_alphas[i];
    dist_params[1] = w_betas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(weibull)");
  }

  // histogram bin uncertain: pairs are defined from an abscissa in the first
  // field and a count (not a density) in the second field.  The distinction
  // in the second field is only important for unequal bin widths.
  for (i=0; i<num_hbuv; ++i, ++cntr) {
    name_string = f77name16("HistBin", cntr, lhs_names);
    String dist_string("continuous linear");
    dist_string.resize(32, ' ');

    const RealRealMap& h_bin_prs_i = h_bin_prs[i];
    RRMCIter cit;

    num_params = h_bin_prs_i.size();
    Real* x_val = new Real [num_params];
    Real* y_val = new Real [num_params];
    // LHS requires accumulation of CDF with first y at 0 and last y at 1
    // Assume already normalized with sum = 1
    //Real sum = 0.;
    //RRMCIter end = --h_bin_prs_i.end(); // last y from DAKOTA must be zero
    //for (cit=h_bin_prs_i.begin(); cit!=end; ++cit)
    //  sum += cit->second;
    size_t end = num_params - 1;
    y_val[0] = 0.;
    for (j=0, cit=h_bin_prs_i.begin(); j<end; ++j, ++cit) {
      x_val[j]   = cit->first;
      y_val[j+1] = y_val[j] + cit->second/* /sum */;
    }
    x_val[end] = cit->first; // last prob value (cit->second) must be zero
    LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		  num_params, x_val, y_val, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_udist(histogram bin)");
    delete [] x_val;
    delete [] y_val;
  }

  // continuous interval uncertain: convert to histogram for sampling
  for (i=0; i<num_ciuv; ++i, ++cntr) {
    name_string = f77name16("ContInterval", cntr, lhs_names);
    String dist_string("continuous linear");
    dist_string.resize(32, ' ');

    const RealRealPairRealMap& ci_bpa_i = ci_bpa[i];
    RRPRMCIter cit;

    // x_sort_unique is a set with ALL of the interval bounds for this variable
    // in increasing order and unique.  For example, if there are 2 intervals
    // for a variable, and the bounds are (1,4) and (3,6), x_sorted will be
    // (1, 3, 4, 6).  If the intervals are contiguous, e.g. one interval is
    // (1,3) and the next is (3,5), x_sort_unique is (1,3,5).
    RealSet x_sort_unique;
    for (cit=ci_bpa_i.begin(); cit!=ci_bpa_i.end(); ++cit) {
      const RealRealPair& bounds = cit->first;
      x_sort_unique.insert(bounds.first);
      x_sort_unique.insert(bounds.second);
    }
    // convert sorted RealSet to x_val
    num_params = x_sort_unique.size();
    Real* x_val = new Real [num_params];
    RSIter it = x_sort_unique.begin();
    for (j=0; j<num_params; ++j, ++it)
      x_val[j] = *it;

    // Calculate the probability densities, and account for the cases where
    // there are intervals that are overlapping.  This section of code goes
    // through the original intervals and see where they fall relative to the
    // new, sorted intervals for the density calculation.
    RealVector prob_dens(num_params); // initialize to 0.
    for (cit=ci_bpa_i.begin(); cit!=ci_bpa_i.end(); ++cit) {
      const RealRealPair& bounds = cit->first;
      Real l_bnd = bounds.first, u_bnd = bounds.second;
      Real ci_density = cit->second / (u_bnd - l_bnd);
      int cum_int_index = 0;
      while (l_bnd > x_val[cum_int_index])
	++cum_int_index;
      ++cum_int_index;
      while (cum_int_index < num_params && x_val[cum_int_index] <= u_bnd)
	{ prob_dens[cum_int_index] += ci_density; ++cum_int_index; }
    }

    // put the densities in a cumulative format necessary for LHS histograms.
    // Note that x_val and y_val are defined as Real* for input to f77.
    Real* y_val = new Real [num_params];
    y_val[0] = 0.;
    for (j=1; j<num_params; ++j) {
      if (prob_dens[j] > 0.0)
	y_val[j] = y_val[j-1] + prob_dens[j] * (x_val[j] - x_val[j-1]);
      else // handle case where there is a gap
	y_val[j] = y_val[j-1] + 0.0001;
    }
    // normalize if necessary
    if (y_val[num_params-1] != 1.) {
      Real y_total = y_val[num_params-1];
      for (j=1; j<num_params; ++j)
	y_val[j] /= y_total;
    }

#ifdef DEBUG
    for (j=0; j<num_params; ++j)
      PCout << "ciuv[" << i << "]: x_val[" << j << "] is " << x_val[j]
	    << " y_val[" << j << "] is " << y_val[j] << '\n';
#endif // DEBUG
    LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		  num_params, x_val, y_val, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_udist(continuous interval)");
    delete [] x_val;
    delete [] y_val;
  }

  // continuous state (treated as uniform)
  for (i=0; i<num_csv; ++i, ++cntr) {
    name_string = f77name16("State", cntr, lhs_names);
    String dist_string("uniform");
    dist_string.resize(32, ' ');
    num_params = 2;
    if (cs_l_bnds[i] > -dbl_inf && cs_u_bnds[i] < dbl_inf) {
      if (cs_l_bnds[i] >= cs_u_bnds[i]) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds strictly "
	      << "less than upper bounds to\n       sample state variables "
	      << "using uniform distributions." << std::endl;
	abort_handler(-1);
      }
      dist_params[0] = cs_l_bnds[i];
      dist_params[1] = cs_u_bnds[i];
    }
    else {
      PCerr << "\nError: Pecos::LHSDriver requires bounds to sample state "
	    << "variables using uniform\n       distributions." << std::endl;
      abort_handler(-1);
    }
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(state)");
  }

  // --------------------------
  // DISCRETE INTEGER VARIABLES
  // --------------------------
  // discrete design range (treated as discrete histogram)
  for (i=0; i<num_ddriv; ++i, ++cntr) {
    name_string = f77name16("DiscDesRange", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');
    int lb_i = ddri_l_bnds[i], ub_i = ddri_u_bnds[i];
    if (lb_i > INT_MIN && ub_i < INT_MAX) {
      if (lb_i > ub_i) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds <= upper "
	      << "bounds to\n       sample discrete design variables using "
	      << "discrete histogram distributions." << std::endl;
	abort_handler(-1);
      }
      num_params = ub_i - lb_i + 1;
      if (num_params == 1) {
	Real pt_val = (Real)lb_i;
	LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
	check_error(err_code, "lhs_const(discrete design range)");
      }
      else {
	Real* x_val = new Real [num_params];
	Real* y_val = new Real [num_params];
	for (j=0; j<num_params; ++j) {
	  x_val[j] = (Real)(lb_i+j);
	  y_val[j] = 1.;
	}
	LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		      num_params, x_val, y_val, err_code, dist_num, pv_num);
	check_error(err_code, "lhs_udist(discrete design range)");
	delete [] x_val;
	delete [] y_val;
      }
    }
    else {
      PCerr << "\nError: Pecos::LHSDriver requires bounds to sample discrete "
	    << "design variables\n       using discrete histogram "
	    << "distributions." << std::endl;
      abort_handler(-1);
    }
  }

  // discrete design set integer (treated as discrete histogram)
  for (i=0; i<num_ddsiv; ++i, ++cntr) {
    name_string = f77name16("DiscDesSetI", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');
    const IntSet& ddsi_vals_i = ddsi_values[i];
    num_params  = ddsi_vals_i.size();
    ISCIter cit = ddsi_vals_i.begin();
    if (num_params == 1) {
      Real pt_val = (Real)(*cit);
      LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
      check_error(err_code, "lhs_const(discrete design set int)");
    }
    else {
      Real* x_val = new Real [num_params];
      Real* y_val = new Real [num_params];
      for (j=0; j<num_params; ++j, ++cit) {
	x_val[j] = (Real)(*cit);
	y_val[j] = 1.;
      }
      LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		    num_params, x_val, y_val, err_code, dist_num, pv_num);
      check_error(err_code, "lhs_udist(discrete design set int)");
      delete [] x_val;
      delete [] y_val;
    }
  }

  // poisson uncertain
  for (i=0; i<num_puv; ++i, ++cntr) {
    name_string = f77name16("Poisson", cntr, lhs_names);
    String dist_string("poisson");
    dist_string.resize(32, ' ');
    num_params = 1;
    dist_params[0] = p_lambdas[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(poisson)");
  }

  // binomial uncertain
  for (i=0; i<num_biuv; ++i, ++cntr) {
    name_string = f77name16("Binomial", cntr, lhs_names);
    String dist_string("binomial");
    dist_string.resize(32, ' ');
    num_params = 2;
    dist_params[0] = bi_prob_per_tr[i];
    dist_params[1] = bi_num_tr[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(binomial)");
  }

  // negative binomial uncertain
  for (i=0; i<num_nbuv; ++i, ++cntr) {
    name_string = f77name16("NegBinomial", cntr, lhs_names);
    String dist_string("negative binomial");
    dist_string.resize(32, ' ');
    num_params = 2;
    dist_params[0] = nb_prob_per_tr[i];
    dist_params[1] = nb_num_tr[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(negative binomial)");
  }

  // geometric uncertain
  for (i=0; i<num_geuv; ++i, ++cntr) {
    name_string = f77name16("Geometric", cntr, lhs_names);
    String dist_string("geometric");
    dist_string.resize(32, ' ');
    num_params = 1;
    dist_params[0] = ge_prob_per_tr[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(geometric)");
  }

  // hypergeometric uncertain
  for (i=0; i<num_hguv; ++i, ++cntr) {
    name_string = f77name16("Hypergeom", cntr, lhs_names);
    String dist_string("hypergeometric");
    dist_string.resize(32, ' ');
    num_params = 3;
    dist_params[0] = hg_total_pop[i];
    dist_params[1] = hg_num_drawn[i];
    dist_params[2] = hg_selected_pop[i];
    LHS_DIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		 dist_params, num_params, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_dist(hypergeometric)");
  }

  // histogram point int: map from an integer abscissa to a probability
  for (i=0; i<num_hpuiv; ++i, ++cntr) {
    name_string = f77name16("HistPtInt", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');

    const IntRealMap& h_pt_int_prs_i = h_pt_int_prs[i];
    num_params = h_pt_int_prs_i.size();
    IRMCIter cit = h_pt_int_prs_i.begin();

    if (num_params == 1) {
      Real pt_val = (Real)cit->first; // frequency is 1
      LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
      check_error(err_code, "lhs_const(histogram pt int)");
    }
    else {
      Real* x_val = new Real [num_params];
      Real* y_val = new Real [num_params];
      // LHS can use discrete frequency information directly
      for (j=0; j<num_params; ++j, ++cit) {
	x_val[j] = (Real)cit->first;
	y_val[j] = cit->second;
      }
      LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		    num_params, x_val, y_val, err_code, dist_num, pv_num);
      check_error(err_code, "lhs_udist(histogram pt int)");
      delete [] x_val;
      delete [] y_val;
    }
  }

  // discrete interval uncertain
  for (i=0; i<num_diuv; ++i, ++cntr) {
    name_string = f77name16("DiscInterval", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');

    const IntIntPairRealMap& di_bpa_i = di_bpa[i];
    IIPRMCIter cit;

    // x_sort_unique contains ALL of the unique integer values for this
    // discrete interval variable in increasing order.  For example, if
    // there are 3 intervals for a variable and the bounds are (1,4),
    // (3,6), and [9,10], x_sorted will be (1,2,3,4,5,6,9,10).
    IntSet x_sort_unique;
    for (cit=di_bpa_i.begin(); cit!=di_bpa_i.end(); ++cit) {
      const RealRealPair& bounds = cit->first;
      int val, u_bnd = bounds.second;
      for (val=bounds.first; val<=u_bnd; ++val)
	x_sort_unique.insert(val);
    }
    // copy sorted IntSet to x_val
    num_params = x_sort_unique.size();
    Real* x_val = new Real [num_params];
    ISIter it = x_sort_unique.begin();
    for (j=0; j<num_params; ++j, ++it)
      x_val[j] = *it;

    // Calculate probability densities and account for overlapping intervals.
    // Loop over the original intervals and see where they fall relative to
    // the new, sorted intervals for the density calculation.
    Real* y_val = new Real [num_params];
    for (j=0; j<num_params; ++j) y_val[j] = 0.;
    int l_bnd, u_bnd; size_t index;
    for (cit=di_bpa_i.begin(); cit!=di_bpa_i.end(); ++cit) {
      const RealRealPair& bounds = cit->first;
      int val, l_bnd = bounds.first, u_bnd = bounds.second;
      Real di_density = cit->second / (u_bnd - l_bnd + 1); // prob/#integers
      it = x_sort_unique.find(l_bnd);
      if (it == x_sort_unique.end()) {
	PCerr << "Error: lower bound not found in sorted set within LHSDriver "
	      << "mapping of discrete interval uncertain variable."<< std::endl;
	abort_handler(-1);
      }
      index = std::distance(x_sort_unique.begin(), it);
      for (val=l_bnd; val<=u_bnd; ++val, ++index)
	y_val[index] += di_density;
    }

#ifdef DEBUG
    for (j=0; j<num_params; ++j)
      PCout << "diuv[" << i << "]: x_val[" << j << "] is " << x_val[j]
	    << " y_val[" << j << "] is " << y_val[j] << '\n';
#endif // DEBUG
    LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		  num_params, x_val, y_val, err_code, dist_num, pv_num);
    check_error(err_code, "lhs_udist(discrete interval)");
    delete [] x_val;
    delete [] y_val;
  }

  // discrete uncertain set integer (treated as discrete histogram)
  for (i=0; i<num_dusiv; ++i, ++cntr) {
    name_string = f77name16("DiscUncSetI", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');

    const IntRealMap& dusi_v_p_i = dusi_vals_probs[i];
    num_params = dusi_v_p_i.size();
    IRMCIter cit = dusi_v_p_i.begin();

    if (num_params == 1) {
      Real pt_val = (Real)(cit->first);
      LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
      check_error(err_code, "lhs_const(discrete uncertain set int)");
    }
    else {
      Real* x_val = new Real [num_params];
      Real* y_val = new Real [num_params];
      for (j=0; j<num_params; ++j, ++cit) {
	x_val[j] = (Real)(cit->first); // discrete uncertain set value
	y_val[j] = cit->second;        // basic probability
      }
      LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		    num_params, x_val, y_val, err_code, dist_num, pv_num);
      check_error(err_code, "lhs_udist(discrete uncertain set int)");
      delete [] x_val;
      delete [] y_val;
    }
  }

  // discrete state range (treated as discrete histogram)
  for (i=0; i<num_dsriv; ++i, ++cntr) {
    name_string = f77name16("DiscStateRange", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');
    int l_bnd = dsri_l_bnds[i], u_bnd = dsri_u_bnds[i];
    if (l_bnd > INT_MIN && u_bnd < INT_MAX) {
      if (l_bnd > u_bnd) {
	PCerr << "\nError: Pecos::LHSDriver requires lower bounds <= upper "
	      << "bounds to\n       sample discrete state variables using "
	      << "discrete histogram distributions." << std::endl;
	abort_handler(-1);
      }
      num_params = u_bnd - l_bnd + 1;
      if (num_params == 1) {
	Real pt_val = (Real)l_bnd;
	LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
	check_error(err_code, "lhs_const(discrete state range)");
      }
      else {
	Real* x_val = new Real [num_params];
	Real* y_val = new Real [num_params];
	for (j=0; j<num_params; ++j) {
	  x_val[j] = (Real)(l_bnd+j);
	  y_val[j] = 1.;
	}
	LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		      num_params, x_val, y_val, err_code, dist_num, pv_num);
	check_error(err_code, "lhs_udist(discrete state range)");
	delete [] x_val;
	delete [] y_val;
      }
    }
    else {
      PCerr << "\nError: Pecos::LHSDriver requires bounds to sample discrete "
	    << "state variables\n       using discrete histogram distributions."
	    << std::endl;
      abort_handler(-1);
    }
  }

  // discrete state set integer (treated as discrete histogram)
  for (i=0; i<num_dssiv; ++i, ++cntr) {
    name_string = f77name16("DiscStateSetI", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');

    const IntSet& dssi_vals_i = dssi_values[i];
    num_params = dssi_vals_i.size();
    ISCIter cit = dssi_vals_i.begin();

    if (num_params == 1) {
      Real pt_val = (Real)(*cit);
      LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
      check_error(err_code, "lhs_const(discrete state set int)");
    }
    else {
      Real* x_val = new Real [num_params];
      Real* y_val = new Real [num_params];
      for (j=0; j<num_params; ++j, ++cit) {
	x_val[j] = (Real)(*cit);
	y_val[j] = 1.;
      }
      LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		    num_params, x_val, y_val, err_code, dist_num, pv_num);
      check_error(err_code, "lhs_udist(discrete state set int)");
      delete [] x_val;
      delete [] y_val;
    }
  }

  // -------------------------
  // DISCRETE STRING VARIABLES
  // -------------------------
  // discrete design set string (treated as discrete histogram)
  for (i=0; i<num_ddssv; ++i, ++cntr) {
    name_string = f77name16("DiscDesSetS", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');

    const StringSet& ddss_vals_i = ddss_values[i];
    num_params = ddss_vals_i.size();
    //SSCIter cit = ddss_vals_i.begin();

    if (num_params == 1) {
      Real pt_val = 0.;//*cit; // index value used to define string
      LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
      check_error(err_code, "lhs_const(discrete design set string)");
    }
    else {
      Real* x_val = new Real [num_params];
      Real* y_val = new Real [num_params];
      for (j=0; j<num_params; ++j) {//, ++cit) {
	x_val[j] = (Real)j;//*cit; // index value used to define string
	y_val[j] = 1.;             // equal probability
      }
      LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		    num_params, x_val, y_val, err_code, dist_num, pv_num);
      check_error(err_code, "lhs_udist(discrete design set string)");
      delete [] x_val;
      delete [] y_val;
    }
  }

  // histogram point string: map from a string abscissa to a probability
  for (i=0; i<num_hpusv; ++i, ++cntr) {
    name_string = f77name16("HistPtString", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');

    const StringRealMap& h_pt_string_prs_i = h_pt_string_prs[i];
    num_params = h_pt_string_prs_i.size();
    SRMCIter cit = h_pt_string_prs_i.begin();

    if (num_params == 1) {
      Real pt_val = 0.;//cit->first; // index value used to define string
      // probability information in cit->second is discarded
      LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
      check_error(err_code, "lhs_const(histogram pt string)");
    }
    else {
      Real* x_val = new Real [num_params];
      Real* y_val = new Real [num_params];
      // LHS can use discrete frequency information directly
      for (j=0; j<num_params; ++j, ++cit) {
	x_val[j] = (Real)j;//cit->first; // index value used to define string
	y_val[j] = cit->second;
      }
      LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		    num_params, x_val, y_val, err_code, dist_num, pv_num);
      check_error(err_code, "lhs_udist(histogram pt string)");
      delete [] x_val;
      delete [] y_val;
    }
  }

  // discrete uncertain set string (treated as discrete histogram)
  for (i=0; i<num_dussv; ++i, ++cntr) {
    name_string = f77name16("DiscUncSetR", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');

    const StringRealMap& duss_v_p_i = duss_vals_probs[i];
    num_params = duss_v_p_i.size();
    SRMCIter cit = duss_v_p_i.begin();

    if (num_params == 1) {
      Real pt_val = 0.;//cit->first; // index value used to define string
      // probability information in cit->second is discarded
      LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
      check_error(err_code, "lhs_const(discrete uncertain set string)");
    }
    else {
      Real* x_val = new Real [num_params];
      Real* y_val = new Real [num_params];
      for (j=0; j<num_params; ++j, ++cit) {
	x_val[j] = (Real)j;//cit->first; // index value used to define string
	y_val[j] = cit->second; // basic probability
      }
      LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		    num_params, x_val, y_val, err_code, dist_num, pv_num);
      check_error(err_code, "lhs_udist(discrete uncertain set string)");
      delete [] x_val;
      delete [] y_val;
    }
  }

  // discrete state set string (treated as discrete histogram)
  for (i=0; i<num_dsssv; ++i, ++cntr) {
    name_string = f77name16("DiscStateSetR", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');

    const StringSet& dsss_vals_i = dsss_values[i];
    num_params = dsss_vals_i.size();
    //SSCIter cit = dsss_vals_i.begin();

    if (num_params == 1) {
      Real pt_val = 0.;//*cit; // index value used to define string
      LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
      check_error(err_code, "lhs_const(discrete state set string)");
    }
    else {
      Real* x_val = new Real [num_params];
      Real* y_val = new Real [num_params];
      for (j=0; j<num_params; ++j) {//, ++cit) {
	x_val[j] = (Real)j;//*cit; // index value used to define string
	y_val[j] = 1.;             // equal probability
      }
      LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		    num_params, x_val, y_val, err_code, dist_num, pv_num);
      check_error(err_code, "lhs_udist(discrete state set string)");
      delete [] x_val;
      delete [] y_val;
    }
  }

  // -----------------------
  // DISCRETE REAL VARIABLES
  // -----------------------
  // discrete design set real (treated as discrete histogram)
  for (i=0; i<num_ddsrv; ++i, ++cntr) {
    name_string = f77name16("DiscDesSetR", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');

    const RealSet& ddsr_vals_i = ddsr_values[i];
    num_params = ddsr_vals_i.size();
    RSCIter cit = ddsr_vals_i.begin();

    if (num_params == 1) {
      Real pt_val = *cit;
      LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
      check_error(err_code, "lhs_const(discrete design set real)");
    }
    else {
      Real* x_val = new Real [num_params];
      Real* y_val = new Real [num_params];
      for (j=0; j<num_params; ++j, ++cit) {
	x_val[j] = *cit;
	y_val[j] = 1.;
      }
      LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		    num_params, x_val, y_val, err_code, dist_num, pv_num);
      check_error(err_code, "lhs_udist(discrete design set real)");
      delete [] x_val;
      delete [] y_val;
    }
  }

  // histogram point real: map from a real abscissa to a probability
  for (i=0; i<num_hpurv; ++i, ++cntr) {
    name_string = f77name16("HistPtReal", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');

    const RealRealMap& h_pt_real_prs_i = h_pt_real_prs[i];
    num_params = h_pt_real_prs_i.size();
    RRMCIter cit = h_pt_real_prs_i.begin();

    if (num_params == 1) {
      Real pt_val = cit->first; // frequency is 1
      LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
      check_error(err_code, "lhs_const(histogram pt real)");
    }
    else {
      Real* x_val = new Real [num_params];
      Real* y_val = new Real [num_params];
      // LHS can use discrete frequency information directly
      for (j=0; j<num_params; ++j, ++cit) {
	x_val[j] = cit->first;
	y_val[j] = cit->second;
      }
      LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		    num_params, x_val, y_val, err_code, dist_num, pv_num);
      check_error(err_code, "lhs_udist(histogram pt real)");
      delete [] x_val;
      delete [] y_val;
    }
  }

  // discrete uncertain set real (treated as discrete histogram)
  for (i=0; i<num_dusrv; ++i, ++cntr) {
    name_string = f77name16("DiscUncSetR", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');

    const RealRealMap& dusr_v_p_i = dusr_vals_probs[i];
    num_params = dusr_v_p_i.size();
    RRMCIter cit = dusr_v_p_i.begin();

    if (num_params == 1) {
      Real pt_val = cit->first; // basic probability is 1
      LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
      check_error(err_code, "lhs_const(discrete uncertain set real)");
    }
    else {
      Real* x_val = new Real [num_params];
      Real* y_val = new Real [num_params];
      for (j=0; j<num_params; ++j, ++cit) {
	x_val[j] = cit->first;  // discrete uncertain set value
	y_val[j] = cit->second; // basic probability
      }
      LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		    num_params, x_val, y_val, err_code, dist_num, pv_num);
      check_error(err_code, "lhs_udist(discrete uncertain set real)");
      delete [] x_val;
      delete [] y_val;
    }
  }

  // discrete state set real (treated as discrete histogram)
  for (i=0; i<num_dssrv; ++i, ++cntr) {
    name_string = f77name16("DiscStateSetR", cntr, lhs_names);
    String dist_string("discrete histogram");
    dist_string.resize(32, ' ');

    const RealSet& dssr_vals_i = dssr_values[i];
    num_params = dssr_vals_i.size();
    RSCIter cit = dssr_vals_i.begin();

    if (num_params == 1) {
      Real pt_val = *cit;
      LHS_CONST2_FC(name_string, pt_val, err_code, pv_num);
      check_error(err_code, "lhs_const(discrete state set real)");
    }
    else {
      Real* x_val = new Real [num_params];
      Real* y_val = new Real [num_params];
      for (j=0; j<num_params; ++j, ++cit) {
	x_val[j] = *cit;
	y_val[j] = 1.;
      }
      LHS_UDIST2_FC(name_string, ptval_flag, ptval, dist_string.data(),
		    num_params, x_val, y_val, err_code, dist_num, pv_num);
      check_error(err_code, "lhs_udist(discrete state set real)");
      delete [] x_val;
      delete [] y_val;
    }
  }

  // specify the rank correlations among the uncertain vars (no correlation
  // currently supported for design and state vars in allVars mode).  Only
  // non-zero values in the lower triangular portion of the rank correlation
  // matrix are specified.
  if (correlation_flag) {
    for (i=1; i<num_auv; ++i) {
      for (j=0; j<i; ++j) {
	Real corr_val = correlations(i,j);
	if (fabs(corr_val) > SMALL_NUMBER) {
	  // jump over cdv, ceuv, csv, ddv int, dsv int, and ddv real as needed:
	  size_t offset_i = num_cdv, offset_j = num_cdv;
	  if (i>=num_cauv)
	    offset_i += num_ceuv + num_csv + num_ddriv + num_ddsiv;
	  if (i>=num_cauv+num_dauiv)
	    offset_i += num_dsriv + num_dssiv + num_ddsrv;
	  if (j>=num_cauv)
	    offset_j += num_ceuv + num_csv + num_ddriv + num_ddsiv;
	  if (j>=num_cauv+num_dauiv)
	    offset_j += num_dsriv + num_dssiv + num_ddsrv;
	  LHS_CORR2_FC(const_cast<char*>(lhs_names[i+offset_i].data()),
		       const_cast<char*>(lhs_names[j+offset_j].data()),
		       corr_val, err_code);
	  check_error(err_code, "lhs_corr");
	}
      }
    }
  }

  // perform internal checks on input to LHS
  int num_nam = num_av, num_var = num_av;
  LHS_PREP_FC(err_code, num_nam, num_var);
  check_error(err_code, "lhs_prep");

  // allocate the memory to hold samples, pt values, variable names, etc.
  int   max_nam        = num_av;
  int*  index_list     = new int    [max_nam];       // output
  Real* ptval_list     = new Real   [max_nam];       // output
  char* dist_name_list = new char   [16*max_nam];    // output
  // dist_name_list is a bit tricky since the f90 array is declared as
  // CHARACTER(LEN=16) :: LSTNAME(MAXNAM), which would seem to be a
  // noncontiguous memory model.  However, a char** does not map correctly to
  // f90.  Rather, f90 takes the contiguous memory block from the C++ char*
  // allocation and indexes into it as if it were a vector of 16 char arrays
  // arranged head to tail.

  // The matrix of parameter samples from Fortran 90 is arranged in column-major
  // order with all variables for sample 1, followed by all variables for
  // sample 2, etc.  Teuchos::SerialDenseMatrix using column-major memory layout
  // as well, so use samples(var#,sample#) or samples[sample#][var#] for access.
  if (samples.numRows() != num_av || samples.numCols() != num_samples)
    samples.shapeUninitialized(num_av, num_samples);
  if (sampleRanksMode && sample_ranks.empty()) {
    if (sampleRanksMode == SET_RANKS || sampleRanksMode == SET_GET_RANKS) {
      PCerr << "Error: empty sample ranks array cannot be set in Pecos::"
	    << "LHSDriver::get_parameter_sets()" << std::endl;
      abort_handler(-1);
    }
    else if (sampleRanksMode == GET_RANKS)
      sample_ranks.shapeUninitialized(num_av, num_samples);
  }

  // generate the samples
  int rflag = sampleRanksMode; // short -> int
  LHS_RUN_FC(max_var, max_obs, max_nam, err_code, dist_name_list,
	     index_list, ptval_list, num_nam, samples.values(), num_var,
	     sample_ranks.values(), rflag);
  check_error(err_code, "lhs_run");

  // deallocate LHS memory
  LHS_CLOSE_FC(err_code);
  check_error(err_code, "lhs_close");

  // clean up memory
  delete [] index_list;
  delete [] ptval_list;
  delete [] dist_name_list;
#endif // HAVE_LHS
}

void LHSDriver::
generate_unique_samples( const RealVector& cd_l_bnds,
			 const RealVector& cd_u_bnds,
			 const IntVector&  ddri_l_bnds, 
			 const IntVector&  ddri_u_bnds,
			 const IntSetArray& ddsi_values,
			 const StringSetArray& ddss_values,
			 const RealSetArray& ddsr_values,
			 const RealVector& cs_l_bnds,
			 const RealVector& cs_u_bnds,
			 const IntVector&  dsri_l_bnds,
			 const IntVector&  dsri_u_bnds,
			 const IntSetArray& dssi_values,
			 const StringSetArray& dsss_values,
			 const RealSetArray& dssr_values, 
			 const AleatoryDistParams& adp,
			 const EpistemicDistParams& edp, int num_samples,
			 RealMatrix& samples, RealMatrix& sample_ranks )
{
  // NonDSampling ordering of variables
  // Design    - continuous, discrete range, discrete set integer, 
  //             discrete set string, discrete set real
  // Aleatory  - continuous, discrete range, discrete set integer, 
  //             discrete set string, discrete set real
  // Epistemic - continuous, discrete range, discrete set integer, 
  //             discrete set string, discrete set real
  // State     - continuous, discrete range, discrete set integer, 
  //             discrete set string, discrete set real
 
  // LHSDriver ordering of variables returned in allSamples matrix:
  // Continuous - design, aleatory, epistemic, state
  // Discrete range - design,    Discrete set int - design
  // Discrete range - aleatory,  Discrete set int - aleatory
  // Discrete range - epistemic, Discrete set int - epistemic
  // Discrete range - state,     Discrete set int - state
  // Discrete set string - design, aleatory, epistemic, state
  // Discrete set real   - design, aleatory, epistemic, state

  // LHSDriver ordering of discrete uncertain variables returned in 
  // allSamples matrix
  // poisson uncertain, binomial uncertain, negative binomial uncertain, 
  // geometric uncertain, hypergeometric uncertain, histogram point int, 
  // discrete interval uncertain, discrete uncertain set integer

  // allSamples is a num_dims x num_samples RealMatrix
  // the values are stored in the following manner:
  // continuous: value
  // discrete set string: integer (mapped to Real) 
  //   I think an alphabetical ordering is applied internally by dakota so
  //   string values input with order 'b' 'a' are assigned int values of 1 0
  // discrete set int: value (mapped to Real)
  // discrete integer range: value (mapped to Real)
  // discrete set real: value

  // get total number of variables and number of discrete variables
  size_t num_cd_vars = cd_l_bnds.length(), num_cs_vars = cs_l_bnds.length(),
    num_cau_vars = adp.cauv(), num_ceu_vars = edp.ceuv(), 
    num_dau_vars = adp.dauv(), num_deu_vars = edp.deuv(),
    //num_daui_vars=adp.dauiv(),num_dausv_vars=adp.dausv(),num_daur_vars=daurv(),
    num_ddri_vars = ddri_l_bnds.length(), num_ddsi_vars = ddsi_values.size(),
    num_ddss_vars = ddss_values.size(), num_ddsr_vars = ddsr_values.size(),
    num_dsri_vars = dsri_l_bnds.length(), num_dssi_vars = dssi_values.size(),
    num_dsss_vars = dsss_values.size(), num_dssr_vars = dssr_values.size();
  size_t num_continuous_vars = num_cd_vars + num_cs_vars + num_cau_vars + 
    num_ceu_vars;
  size_t num_discrete_vars = num_dau_vars + num_deu_vars + num_ddri_vars + 
    num_ddsi_vars + num_ddss_vars + num_ddsr_vars + num_dsri_vars + 
    num_dssi_vars + num_dsss_vars + num_dssr_vars;
  size_t num_vars = num_continuous_vars + num_discrete_vars;

  // compute the number of total possible combinations of discrete variables
  // TODO Must look over all data structures passed into function
  // for range variables use ub-lb+1
  // for set variables use sum_i (ddsi_values[i].size())
  // if num_values**d < num_samples then call generate_samples
  int k=0;
  IntVector num_discrete_strata_1d( num_discrete_vars, false );
  // Discrete design variables
  for (int i=0; i<num_ddri_vars; i++)
    {  num_discrete_strata_1d[k] = ddri_u_bnds[i]-ddri_l_bnds[i]+1; k++; }
  for (int i=0; i<num_ddsi_vars; i++)
    {  num_discrete_strata_1d[k] = ddsi_values[i].size(); k++; }
  for (int i=0; i<num_ddss_vars; i++)
    {  num_discrete_strata_1d[k] = ddss_values[i].size(); k++; }
  for (int i=0; i<num_ddsr_vars; i++)
    {  num_discrete_strata_1d[k] = ddsr_values[i].size(); k++; }
  for (int i=0; i<num_dsri_vars; i++)

  // Discrete state variables
    {  num_discrete_strata_1d[k] = dsri_u_bnds[i]-dsri_l_bnds[i]+1; k++; }
  for (int i=0; i<num_dssi_vars; i++)
    {  num_discrete_strata_1d[k] = dssi_values.size(); k++; }
  for (int i=0; i<num_dsss_vars; i++)
    {  num_discrete_strata_1d[k] = dsss_values.size(); k++; }
  for (int i=0; i<num_dssr_vars; i++)
    {  num_discrete_strata_1d[k] = dssr_values.size(); k++; }

  // Epistemic discrete variables
  IntIntPairRealMapArray di_bpa = edp.discrete_interval_basic_probabilities();
  int num_diuv = di_bpa.size();
  for (int i=0; i<num_diuv; ++i){
    const IntIntPairRealMap& di_bpa_i = di_bpa[i];
    IIPRMCIter cit;
    // x_sort_unique contains ALL of the unique integer values for this
    // discrete interval variable in increasing order.  For example, if
    // there are 3 intervals for a variable and the bounds are (1,4),
    // (3,6), and [9,10], x_sorted will be (1,2,3,4,5,6,9,10).
    IntSet x_sort_unique;
    for (cit=di_bpa_i.begin(); cit!=di_bpa_i.end(); ++cit) {
      const RealRealPair& bounds = cit->first;
      int val, u_bnd = bounds.second;
      for (val=bounds.first; val<=u_bnd; ++val)
	x_sort_unique.insert(val);
    }
    int num_vals = x_sort_unique.size();
    num_discrete_strata_1d[k] = num_vals;
    k++;
  }
  IntRealMapArray dusi_vals_probs = edp.discrete_set_int_values_probabilities();
  int num_dusiv = dusi_vals_probs.size();
  for (int i=0; i<num_dusiv; ++i){
    const IntRealMap& dusi_v_p_i = dusi_vals_probs[i];
    int num_vals = dusi_v_p_i.size();
    num_discrete_strata_1d[k] = num_vals;
    k++;
  }
  StringRealMapArray duss_vals_probs = 
    edp.discrete_set_string_values_probabilities();
  int num_dussv = duss_vals_probs.size();
  for (int i=0; i<num_dussv; ++i){
    const StringRealMap& duss_v_p_i = duss_vals_probs[i];
    int num_vals = duss_v_p_i.size();
    num_discrete_strata_1d[k] = num_vals;
    k++;
  }
  RealRealMapArray dusr_vals_probs = 
    edp.discrete_set_real_values_probabilities();
  int num_dusrv = dusr_vals_probs.size();
  for (int i=0; i<num_dusrv; ++i){
    const RealRealMap& dusr_v_p_i = dusr_vals_probs[i];
    int num_vals = dusr_v_p_i.size();
    num_discrete_strata_1d[k] = num_vals;
    k++;
  }

  // Aleatory discrete variables with finite support
  IntVector bnt = adp.binomial_num_trials();
  for ( int i=0; i < bnt.length(); i++ )
    { num_discrete_strata_1d[k] = bnt[i]+1; k++;}
  IntVector htp = adp.hypergeometric_total_population();
  IntVector hsp = adp.hypergeometric_selected_population();
  IntVector hnd = adp.hypergeometric_num_drawn();
  for ( int i=0; i < htp.length(); i++ )
    { int num_fail=hnd[i], num_total_pop=htp[i], num_sel_pop=hsp[i]; 
      // Todo confirm this
      num_discrete_strata_1d[k] = std::min(num_fail,num_sel_pop)-
	std::max(0,num_sel_pop+num_fail-num_total_pop)+1;
      k++;
    }

  // Aleatory discrete variables with infinite support.  If any of these
  // variables are present then backfill can always be used.
  int max_num_unique_discrete_samples = 1;
  RealVector pl = adp.poisson_lambdas();
  RealVector nbppt = adp.negative_binomial_probability_per_trial();
  RealVector gppt = adp.geometric_probability_per_trial();
  if ( pl.length() > 0 ) max_num_unique_discrete_samples = INT_MAX;
  else if ( nbppt.length() > 0 ) max_num_unique_discrete_samples = INT_MAX;
  else if ( gppt.length() > 0 ) max_num_unique_discrete_samples = INT_MAX;
  else
    {
      num_discrete_strata_1d.resize( num_discrete_vars-pl.length()-
				     nbppt.length()- gppt.length());
      for (int k=0; k < num_discrete_strata_1d.length(); k++)
	  max_num_unique_discrete_samples *= num_discrete_strata_1d[k];
    }
  
  if ( max_num_unique_discrete_samples >= num_samples )
    // If the number of samples requested is greater than the maximum possible
    // number of discrete samples then we must allow replicates of the 
    // disscrete variables to obtain the desired number of variables. 
    // If not then we can proceed with generating a unique set of discrete 
    // samples.
    {
      if (samples.numRows() != num_vars || samples.numCols() != num_samples)
	samples.shapeUninitialized(num_vars, num_samples);
      // Currently sample_ranks will always be returned empty. It should only be
      // filled when NonDSampling.sampleRanksMode>0. But I cannot see anywhere
      // in the code where this is true.
      //sample_ranks.shapeUninitialized( num_vars, num_samples );

      RealMatrix sample_ranks_rm, samples_rm;

      // unique index of all discrete variables if any
      std::set<RealArray>::iterator it;
      std::set<RealArray> sorted_discrete_samples; 
      RealArray discrete_sample( num_discrete_vars );

      // Determine the columns in samples_rm that contain discrete variables
      IntVector discrete_samples_map( num_discrete_vars, false );
      for (int i=0; i<num_discrete_vars; i++)
	discrete_samples_map[i]=num_continuous_vars+i;
  
      // Eliminate redundant samples by resampling if necessary.  Could pad
      // num_samples in anticipation of duplicates, but this would alter LHS
      // stratification that could be intended, so use num_samples for now.
      bool complete = false, initial = true;
      int num_unique_samples = 0;
      while (!complete) {
	generate_samples(cd_l_bnds, cd_u_bnds, ddri_l_bnds, ddri_u_bnds, 
			 ddsi_values, ddss_values, ddsr_values, 
			 cs_l_bnds, cs_u_bnds, dsri_l_bnds, dsri_u_bnds,
			 dssi_values, dsss_values, dssr_values, adp, edp,
			 num_samples, samples_rm, sample_ranks_rm);

	if (initial) { // pack initial sample set
	  for (int i=0; i<num_samples; ++i) { // or matrix->set<vector> ?
	    //std::cout << "[";
	    for (int j=0; j<num_discrete_vars; j++) {
	      int index = discrete_samples_map[j];
	      discrete_sample[j] = samples_rm(index,i);
	      //std::cout << discrete_sample[j] << ",";
	    }
	    //std::cout << "]\n";
	    sorted_discrete_samples.insert(discrete_sample);
	    if ( sorted_discrete_samples.size() > num_unique_samples ){
	      // copy sample into samples matrix
	      for (int j=0; j<num_vars; j++)
		samples(j,num_unique_samples) = samples_rm(j,i);
	      num_unique_samples++;
	    }
	  }
	  if (num_unique_samples == num_samples) complete = true;
	  else initial = false;
	}
	else { // backfill duplicates with new samples
	  //std::cout << num_unique_samples << "," << sorted_discrete_samples.size() << "," << num_discrete_vars << "," << num_vars << "," << num_continuous_vars << std::endl;

	  for (int i=0; i<num_samples; ++i) {
	    if (num_unique_samples < num_samples) {
	      //std::cout << "[";
	      for (int j=0; j<num_discrete_vars; j++) {
		int index = discrete_samples_map[j];
		discrete_sample[j] = samples_rm(index,i);
		//std::cout << discrete_sample[j] << ",";
	      }
	      //std::cout << "]\n";
	      sorted_discrete_samples.insert(discrete_sample);
	      if ( sorted_discrete_samples.size() > num_unique_samples ){
		// copy sample into samples matrix
		for (int j=0; j<num_vars; j++)
		  samples(j,num_unique_samples) = samples_rm(j,i);
		num_unique_samples++;
	      }
	    }
	    else
	      { complete = true; break; }
	  }
	}
      }
    }
  else
    {
      PCout << "LHS backfill was requested, but the discrete variables "
	    << "specified do not have enough unique values ("
	    << max_num_unique_discrete_samples
	    << ") to obtain the number of samples requested so replicated "
	    << "discrete samples have been allowed.\n";
      generate_samples( cd_l_bnds, cd_u_bnds, ddri_l_bnds, ddri_u_bnds,
			ddsi_values, ddss_values, ddsr_values, cs_l_bnds,
			cs_u_bnds, dsri_l_bnds, dsri_u_bnds, dssi_values,
			dsss_values,dssr_values, adp, edp, num_samples,
			samples, sample_ranks );
    }
}


void LHSDriver::
generate_unique_index_samples(const IntVector& index_l_bnds,
			      const IntVector& index_u_bnds, int num_samples,
			      IntMatrix& index_samples )
{
  // For    uniform probability, model as discrete design range (this fn).
  // For nonuniform probability, model as discrete uncertain set integer.

  RealVector  empty_rv;  RealMatrix empty_rm, samples_rm;
  IntVector   empty_iv;  IntArray sample;
  IntSetArray empty_isa; StringSetArray empty_ssa; RealSetArray empty_rsa;
  AleatoryDistParams adp; EpistemicDistParams edp;
  generate_unique_samples( empty_rv, empty_rv, index_l_bnds, index_u_bnds,
			   empty_isa, empty_ssa, empty_rsa, empty_rv, empty_rv,
			   empty_iv, empty_iv, empty_isa, empty_ssa, empty_rsa,
			   adp, edp, num_samples, samples_rm, empty_rm);
  
  index_samples.shapeUninitialized( samples_rm.numRows(), samples_rm.numCols() );
  for (int j=0; j<samples_rm.numCols(); j++){
    for (int i=0; i<samples_rm.numRows(); i++)
      index_samples(i,j) = (int)samples_rm(i,j);
  }
}

} // namespace Pecos
