/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "LHSDriver.hpp"
#include "BoostRNG_Monostate.hpp"
#include "RangeVariable.hpp"
#include "SetVariable.hpp"
#include "DiscreteSetRandomVariable.hpp"
#include "IntervalRandomVariable.hpp"
#include "HistogramBinRandomVariable.hpp"

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

#define defaultrnum1    FC_FUNC(defaultrnum1,DEFAULTRNUM1)
#define defaultrnum2    FC_FUNC(defaultrnum2,DEFAULTRNUM2)

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
#define LHS_OPTIONS2_FC lhs_options2
#define LHS_DIST2_FC    lhs_dist2
#define LHS_UDIST2_FC   lhs_udist2
#define LHS_CONST2_FC   lhs_const2
#define LHS_CORR2_FC    lhs_corr2
#define LHS_FILES2_FC   lhs_files2

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


void LHSDriver::abort_if_no_lhs()
{
#ifndef HAVE_LHS
  PCerr << "Error: LHSDriver not available as PECOS was configured without LHS."
        << std::endl;
  abort_handler(-1);
#endif
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


//String LHSDriver::rng()
//{
//  if (BoostRNG_Monostate::randomNum == BoostRNG_Monostate::random_num1)
//    return "mt19937";
//  if (BoostRNG_Monostate::randomNum == (Rfunc)rnumlhs10)
//    return "rnum2";
//  else
//    return "";
//}


void LHSDriver::
lhs_dist_register(const char* var_name, const char* dist_name, size_t rv,
		  const RealArray& dist_params)
{
  String dist_string;                 f77name32(dist_name,   dist_string);
  String& var_string = lhsNames[rv];  f77name16(var_name, rv, var_string);
  int num_params = dist_params.size(), err_code = 0, ptval_flag = 0, // inputs
      dist_num, pv_num; // outputs (not used)
  Real ptval = 0.;

  LHS_DIST2_FC(var_string.data(), ptval_flag, ptval, dist_string.data(),
	       const_cast<Real*>(&dist_params[0]), num_params, err_code,
	       dist_num, pv_num);
  check_error(err_code, "lhs_dist()", var_string.data());
}


void LHSDriver::
lhs_udist_register(const char* var_name, const char* dist_name, size_t rv,
		   const RealArray& x_val, const RealArray& y_val)
{
  String dist_string;                 f77name32(dist_name,   dist_string);
  String& var_string = lhsNames[rv];  f77name16(var_name, rv, var_string);
  int num_params = std::min(x_val.size(), y_val.size()), err_code = 0,
    ptval_flag = 0, dist_num, pv_num;
  Real ptval = 0.;

  LHS_UDIST2_FC(var_string.data(), ptval_flag, ptval, dist_string.data(),
		num_params, const_cast<Real*>(&x_val[0]),
		const_cast<Real*>(&y_val[0]), err_code, dist_num, pv_num);
  check_error(err_code, "lhs_udist()", var_string.data());
}


void LHSDriver::lhs_const_register(const char* var_name, size_t rv, Real pt_val)
{
  String& var_string = lhsNames[rv];  f77name16(var_name, rv, var_string);
  int err_code = 0, pv_num;

  LHS_CONST2_FC(var_string.data(), pt_val, err_code, pv_num);
  check_error(err_code, "lhs_const()", var_string.data());
}


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
generate_samples(const std::vector<RandomVariable>& random_vars,
		 const RealSymMatrix& corr, int num_samples,
		 RealMatrix& samples, RealMatrix& sample_ranks,
		 const BitArray& active_vars, const BitArray& active_corr)
{
#ifdef HAVE_LHS
  // generate samples within user-specified parameter distributions

  // error check on program parameters
  if (!num_samples) {
    PCerr << "\nError: number of samples in LHSDriver::generate_samples() must "
	  << "be nonzero." << std::endl;
    abort_handler(-1);
  }

  // active_vars identifies the active subset of random_vars for which we will
  // generate samples; it will vary based on sampling context.
  size_t i, num_rv = random_vars.size(), num_active_rv, av_cntr;
  bool subset_rv = false;
  if (active_vars.empty())
    num_active_rv = num_rv;
  else {
    num_active_rv = active_vars.count();
    if (num_active_rv < num_rv) subset_rv = true;
  }
  lhsNames.resize(num_active_rv);

  // active_corr identifies the subset of random_vars for which we specify
  // correlations; it need not vary based on sampling context, potentially
  // indicating a static subset that admits correlations in general (e.g.,
  // cont/disc aleatory RV in the current Dakota XML spec).
  // > Both active subsets are relative to the full RV vector (i.e., not
  //   nested subsets) such that they overlay when specifying correlations
  //   to LHS, e.g., uncertain vars are (currently) active for sampling and
  //   correlations are (always) bound to aleatory uncertain vars, so both
  //   bits must be on to call LHS_CORR2.
  // > avoid resizing the correlation matrix (inflating with 1's on diagonal)
  //   for active sampling context: the size of the corr matrix is defined by
  //   the (static) RV subset that admits correlations (corr.numRows()
  //   == active_corr.count()).  Since this is independent of the active
  //   sampling context, a subset of corr matrix could be passed to LHS_CORR2.
  bool correlation_flag = !corr.empty(), subset_corr = false;
  int  num_corr = 0, num_active_corr = 0;
  if (correlation_flag) {
    if (active_corr.empty())
      num_corr = corr.numRows();
    else {
      for (i=0; i<num_rv; ++i)
	if (active_corr[i])
	  { ++num_active_corr; if (active_vars[i]) ++num_corr; }
      if (num_active_corr < num_rv) subset_corr = true;
    }
  }

  int max_corr_size = (num_corr > 1) ? num_corr * (num_corr - 1) / 2 : -1,
      max_var = num_active_rv, max_samp_size = max_var * num_samples,
      max_interval = -1, max_table = -1, err_code = 0, print_level = 0,
      output_width = 1;
  // randomSeed passed below propagates to ISeed in the f77 rnum2, but does
  // not propagate to Boost RNGs (LHSDriver::seed() must be used for that).
  LHS_INIT_MEM_FC(num_samples, randomSeed, num_samples, max_samp_size, max_var,
		  max_interval, max_corr_size, max_table, print_level,
		  output_width, err_code);
  check_error(err_code, "lhs_init_mem");

  // set sample type to either LHS (default) or random Monte Carlo (optional)
  bool call_lhs_option = false;
  String option_string("              ");
  if (sampleType == "random" || sampleType == "incremental_random")
    { option_string = "RANDOM SAMPLE "; call_lhs_option = true; }
  // set mixing option to either restricted pairing (default) or random pairing
  // (optional).  For enforcing user-specified correlation, restricted pairing
  // is required.  And for uncorrelated variables, restricted pairing results
  // in lower correlation values than random pairing.  For these reasons, the
  // random pairing option is not currently active, although a specification
  // option for it could be added in the future if a use arises.
  bool random_pairing_flag = false; // this option hard-wired off for now
  if (!correlation_flag && random_pairing_flag)
    { option_string += "RANDOM PAIRING"; call_lhs_option = true; }
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

  //////////////////////////////////////////////////////////
  // Register RandomVariables with lhs_{dist,udist,const} //
  //////////////////////////////////////////////////////////
  RealArray dist_params;
  for (i=0, av_cntr=0; i<num_rv; ++i) {
    if (subset_rv && !active_vars[i]) continue; // skip this RV if not active

    const RandomVariable& rv_i = random_vars[i];
    switch (rv_i.type()) {
    case CONTINUOUS_RANGE: {
      std::shared_ptr<RangeVariable<Real>> rv_rep =
	std::static_pointer_cast<RangeVariable<Real>>(rv_i.random_variable_rep());
      Real l_bnd;  rv_rep->pull_parameter(CR_LWR_BND, l_bnd);
      Real u_bnd;  rv_rep->pull_parameter(CR_UPR_BND, u_bnd);
      check_finite(l_bnd, u_bnd);
      if (u_bnd > l_bnd) {
	dist_params.resize(2);
	dist_params[0] = l_bnd;  dist_params[1] = u_bnd;
	lhs_dist_register("ContRange", "uniform", av_cntr, dist_params);
      }
      else {
	check_range(l_bnd, u_bnd, true); // allow equal
	lhs_const_register("ContRange", av_cntr, (Real)l_bnd);
      }
      break;
    }
    case DISCRETE_RANGE: {
      std::shared_ptr<RangeVariable<int>> rv_rep =
	std::static_pointer_cast<RangeVariable<int>>(rv_i.random_variable_rep());
      int l_bnd;  rv_rep->pull_parameter(DR_LWR_BND, l_bnd);
      int u_bnd;  rv_rep->pull_parameter(DR_UPR_BND, u_bnd);
      check_finite(l_bnd, u_bnd);
      if (u_bnd > l_bnd) {
	RealArray x_val, y_val;
	int_range_to_xy_pdf(l_bnd, u_bnd, x_val, y_val);
	lhs_udist_register("DiscRange", "discrete histogram", av_cntr,
			   x_val, y_val);
      }
      else {
	check_range(l_bnd, u_bnd, true); // allow equal
	lhs_const_register("DiscRange", av_cntr, (Real)l_bnd);
      }
      break;
    }
    case DISCRETE_SET_INT: {
      std::shared_ptr<SetVariable<int>> rv_rep =
	std::static_pointer_cast<SetVariable<int>>(rv_i.random_variable_rep());
      IntSet int_set;  rv_rep->pull_parameter(DSI_VALUES, int_set);
      size_t set_size = int_set.size();
      if (set_size > 1) {
	RealArray x_val, y_val;
	set_to_xy_pdf(int_set, x_val, y_val); // set values
	lhs_udist_register("DiscSetI", "discrete histogram", av_cntr,
			   x_val, y_val);
      }
      else if (set_size)
	lhs_const_register("DiscSetI", av_cntr, (Real)(*int_set.begin()));
      break;
    }
    case DISCRETE_SET_STRING: {
      std::shared_ptr<SetVariable<String>> rv_rep =
	std::static_pointer_cast<SetVariable<String>>(rv_i.random_variable_rep());
      StringSet str_set;  rv_rep->pull_parameter(DSS_VALUES, str_set);
      int set_size = str_set.size();
      if (set_size > 1) {
	RealArray x_val, y_val;
	int_range_to_xy_pdf(0, set_size - 1, x_val, y_val); // indices
	lhs_udist_register("DiscSetS", "discrete histogram", av_cntr,
			   x_val, y_val);
      }
      else if (set_size)
	lhs_const_register("DiscSetS", av_cntr, 0.);
      break;
    }
    case DISCRETE_SET_REAL: {
      std::shared_ptr<SetVariable<Real>> rv_rep =
	std::static_pointer_cast<SetVariable<Real>>(rv_i.random_variable_rep());
      RealSet real_set;  rv_rep->pull_parameter(DSR_VALUES, real_set);
      size_t set_size = real_set.size();
      if (set_size > 1) {
	RealArray x_val, y_val;
	set_to_xy_pdf(real_set, x_val, y_val); // set values
	lhs_udist_register("DiscSetR", "discrete histogram", av_cntr,
			   x_val, y_val);
      }
      else if (set_size)
	lhs_const_register("DiscSetR", av_cntr, *real_set.begin());
      break;
    }
    case NORMAL: case STD_NORMAL:
      dist_params.resize(2);
      rv_i.pull_parameter(N_MEAN,    dist_params[0]);
      rv_i.pull_parameter(N_STD_DEV, dist_params[1]);
      lhs_dist_register("Normal", "normal", av_cntr, dist_params);
      break;
    case BOUNDED_NORMAL:
      dist_params.resize(4);
      rv_i.pull_parameter(N_MEAN,    dist_params[0]);
      rv_i.pull_parameter(N_STD_DEV, dist_params[1]);
      rv_i.pull_parameter(N_LWR_BND, dist_params[2]);
      rv_i.pull_parameter(N_UPR_BND, dist_params[3]);
      check_range(dist_params[2], dist_params[3], false);
      lhs_dist_register("Normal", "bounded normal", av_cntr, dist_params);
      break;
    case LOGNORMAL:     // LognormalRandomVariable standardizes on Lambda/Zeta
      dist_params.resize(2);
      rv_i.pull_parameter(LN_LAMBDA,  dist_params[0]);
      rv_i.pull_parameter(LN_ZETA,    dist_params[1]);
      lhs_dist_register("Lognormal", "lognormal-n", av_cntr, dist_params);
      break;
    case BOUNDED_LOGNORMAL: // BoundedLognormalRandomVariable uses Lambda/Zeta
      dist_params.resize(4);
      rv_i.pull_parameter(LN_LAMBDA,  dist_params[0]);
      rv_i.pull_parameter(LN_ZETA,    dist_params[1]);
      rv_i.pull_parameter(LN_LWR_BND, dist_params[2]);
      rv_i.pull_parameter(LN_UPR_BND, dist_params[3]);
      check_range(dist_params[2], dist_params[3], false);
      lhs_dist_register("Lognormal", "bounded lognormal-n",av_cntr,dist_params);
      break;
    case UNIFORM: case STD_UNIFORM:
      dist_params.resize(2);
      rv_i.pull_parameter(U_LWR_BND, dist_params[0]);
      rv_i.pull_parameter(U_UPR_BND, dist_params[1]);
      check_range(dist_params[0],  dist_params[1], false);
      check_finite(dist_params[0], dist_params[1]);
      lhs_dist_register("Uniform", "uniform", av_cntr, dist_params);
      break;
    case LOGUNIFORM:
      dist_params.resize(2);
      rv_i.pull_parameter(LU_LWR_BND, dist_params[0]);
      rv_i.pull_parameter(LU_UPR_BND, dist_params[1]);
      check_range(dist_params[0],  dist_params[1], false);
      check_finite(dist_params[0], dist_params[1]);
      lhs_dist_register("Loguniform", "loguniform", av_cntr, dist_params);
      break;
    case TRIANGULAR:
      dist_params.resize(3);
      rv_i.pull_parameter(T_LWR_BND, dist_params[0]);
      rv_i.pull_parameter(T_MODE,    dist_params[1]);
      rv_i.pull_parameter(T_UPR_BND, dist_params[2]);
      check_range(dist_params[0],  dist_params[2], false);
      check_finite(dist_params[0], dist_params[2]);
      lhs_dist_register("Triangular", "triangular", av_cntr, dist_params);
      break;
    case EXPONENTIAL: case STD_EXPONENTIAL: {
      dist_params.resize(1);
      Real e_beta; rv_i.pull_parameter(E_BETA, e_beta);
      dist_params[0] = 1./e_beta; // convert to LHS convention
      lhs_dist_register("Exponential", "exponential", av_cntr, dist_params);
      break;
    }
    case BETA: case STD_BETA:
      dist_params.resize(4);
      rv_i.pull_parameter(BE_LWR_BND, dist_params[0]);
      rv_i.pull_parameter(BE_UPR_BND, dist_params[1]);
      rv_i.pull_parameter(BE_ALPHA,   dist_params[2]);
      rv_i.pull_parameter(BE_BETA,    dist_params[3]);
      check_range(dist_params[0],  dist_params[1], false);
      check_finite(dist_params[0], dist_params[1]);
      lhs_dist_register("Beta", "beta", av_cntr, dist_params);
      break;
    case GAMMA: case STD_GAMMA: {
      dist_params.resize(2);
      rv_i.pull_parameter(GA_ALPHA, dist_params[0]);
      Real ga_beta; rv_i.pull_parameter(GA_BETA, ga_beta);
      dist_params[1] = 1./ga_beta; // convert to LHS convention
      lhs_dist_register("Gamma", "gamma", av_cntr, dist_params);
      break;
    }
    case GUMBEL:
      dist_params.resize(2);
      rv_i.pull_parameter(GU_ALPHA, dist_params[0]);
      rv_i.pull_parameter(GU_BETA,  dist_params[1]);
      lhs_dist_register("Gumbel", "gumbel", av_cntr, dist_params);
      break;
    case FRECHET:
      dist_params.resize(2);
      rv_i.pull_parameter(F_ALPHA, dist_params[0]);
      rv_i.pull_parameter(F_BETA,  dist_params[1]);
      lhs_dist_register("Frechet", "frechet", av_cntr, dist_params);
      break;
    case WEIBULL:
      dist_params.resize(2);
      rv_i.pull_parameter(W_ALPHA, dist_params[0]);
      rv_i.pull_parameter(W_BETA,  dist_params[1]);
      lhs_dist_register("Weibull", "weibull", av_cntr, dist_params);
      break;
    case HISTOGRAM_BIN: {
      std::shared_ptr<HistogramBinRandomVariable> rv_rep =
	std::static_pointer_cast<HistogramBinRandomVariable>
	(rv_i.random_variable_rep());
      RealRealMap h_bin_pairs; rv_rep->pull_parameter(H_BIN_PAIRS, h_bin_pairs);
      RealArray x_val, y_val;  bins_to_xy_cdf(h_bin_pairs, x_val, y_val);
      // Note: continuous linear accumulates CDF with first y=0 and last y=1
      lhs_udist_register("HistBin", "continuous linear", av_cntr, x_val, y_val);
      break;
    }
    case POISSON:
      dist_params.resize(1);
      rv_i.pull_parameter(P_LAMBDA, dist_params[0]);
      lhs_dist_register("Poisson", "poisson", av_cntr, dist_params);
      break;
    case BINOMIAL: {
      dist_params.resize(2);  unsigned int num_tr;
      rv_i.pull_parameter(BI_P_PER_TRIAL, dist_params[0]);
      rv_i.pull_parameter(BI_TRIALS, num_tr); dist_params[1] = (Real)num_tr;
      lhs_dist_register("Binomial", "binomial", av_cntr, dist_params);
      break;
    }
    case NEGATIVE_BINOMIAL: {
      dist_params.resize(2);  unsigned int num_tr;
      rv_i.pull_parameter(NBI_P_PER_TRIAL, dist_params[0]);
      rv_i.pull_parameter(NBI_TRIALS, num_tr); dist_params[1] = (Real)num_tr;
      lhs_dist_register("NegBinomial", "negative binomial",av_cntr,dist_params);
      break;
    }
    case GEOMETRIC:
      dist_params.resize(1);
      rv_i.pull_parameter(GE_P_PER_TRIAL, dist_params[0]);
      lhs_dist_register("Geometric", "geometric", av_cntr, dist_params);
      break;
    case HYPERGEOMETRIC: {
      dist_params.resize(3);  unsigned int tot_pop, num_drw, sel_pop;
      rv_i.pull_parameter(HGE_TOT_POP, tot_pop); dist_params[0] = (Real)tot_pop;
      rv_i.pull_parameter(HGE_DRAWN,   num_drw); dist_params[1] = (Real)num_drw;
      rv_i.pull_parameter(HGE_SEL_POP, sel_pop); dist_params[2] = (Real)sel_pop;
      lhs_dist_register("Hypergeom", "hypergeometric", av_cntr, dist_params);
      break;
    }
    case HISTOGRAM_PT_INT: {
      std::shared_ptr<DiscreteSetRandomVariable<int>> rv_rep =
	std::static_pointer_cast<DiscreteSetRandomVariable<int>>
	(rv_i.random_variable_rep());
      IntRealMap ir_map;  rv_rep->pull_parameter(H_PT_INT_PAIRS, ir_map);
      size_t map_size = ir_map.size();
      if (map_size > 1) {
	RealArray x_val, y_val;  map_to_xy_pdf(ir_map, x_val, y_val);
	lhs_udist_register("HistPtInt", "discrete histogram", av_cntr,
			   x_val, y_val);
      }
      else if (map_size)
	lhs_const_register("HistPtInt", av_cntr, (Real)ir_map.begin()->first);
      break;
    }
    case HISTOGRAM_PT_STRING: {
      std::shared_ptr<DiscreteSetRandomVariable<String>> rv_rep =
	std::static_pointer_cast<DiscreteSetRandomVariable<String>>
	(rv_i.random_variable_rep());
      StringRealMap sr_map;  rv_rep->pull_parameter(H_PT_STR_PAIRS, sr_map);
      int map_size = sr_map.size();
      if (map_size > 1) {
	RealArray x_val, y_val;  map_indices_to_xy_pdf(sr_map, x_val, y_val);
	lhs_udist_register("HistPtString","discrete histogram",av_cntr,
			   x_val, y_val);
      }
      else if (map_size)
	lhs_const_register("HistPtString", av_cntr, 0.);
      break;
    }
    case HISTOGRAM_PT_REAL: {
      std::shared_ptr<DiscreteSetRandomVariable<Real>> rv_rep =
	std::static_pointer_cast<DiscreteSetRandomVariable<Real>>
	(rv_i.random_variable_rep());
      RealRealMap rr_map;  rv_rep->pull_parameter(H_PT_REAL_PAIRS, rr_map);
      size_t map_size = rr_map.size();
      if (map_size > 1) {
	RealArray x_val, y_val;  map_to_xy_pdf(rr_map, x_val, y_val);
	lhs_udist_register("HistPtReal", "discrete histogram", av_cntr,
			   x_val, y_val);
      }
      else if (map_size)
	lhs_const_register("HistPtReal", av_cntr, rr_map.begin()->first);
      break;
    }
    case CONTINUOUS_INTERVAL_UNCERTAIN: {
      std::shared_ptr<IntervalRandomVariable<Real>> rv_rep =
	std::static_pointer_cast<IntervalRandomVariable<Real>>
	(rv_i.random_variable_rep());
      RealRealPairRealMap ci_bpa;  rv_rep->pull_parameter(CIU_BPA, ci_bpa);
      // Note: continuous linear accumulates CDF with first y=0 and last y=1
      RealArray x_val, y_val;  intervals_to_xy_cdf(ci_bpa, x_val, y_val);
      lhs_udist_register("ContInterval", "continuous linear", av_cntr,
			 x_val, y_val);
      break;
    }
    case DISCRETE_INTERVAL_UNCERTAIN: {
      std::shared_ptr<IntervalRandomVariable<int>> rv_rep =
	std::static_pointer_cast<IntervalRandomVariable<int>>
	(rv_i.random_variable_rep());
      IntIntPairRealMap di_bpa;  rv_rep->pull_parameter(DIU_BPA, di_bpa);
      IntArray i_val;  RealArray x_val, y_val;
      intervals_to_xy_pdf(di_bpa, i_val, y_val);  cast_data(i_val, x_val);
      lhs_udist_register("DiscInterval", "discrete histogram", av_cntr,
			 x_val, y_val);
      break;
    }
    case DISCRETE_UNCERTAIN_SET_INT: {
      std::shared_ptr<DiscreteSetRandomVariable<int>> rv_rep =
	std::static_pointer_cast<DiscreteSetRandomVariable<int>>
	(rv_i.random_variable_rep());
      IntRealMap ir_map;  rv_rep->pull_parameter(DUSI_VALUES_PROBS, ir_map);
      size_t map_size = ir_map.size();
      if (map_size > 1) {
	RealArray x_val, y_val;  map_to_xy_pdf(ir_map, x_val, y_val);
	lhs_udist_register("DiscUncSetI","discrete histogram", av_cntr,
			   x_val, y_val);
      }
      else if (map_size)
	lhs_const_register("DiscUncSetI", av_cntr, (Real)ir_map.begin()->first);
      break;
    }
    case DISCRETE_UNCERTAIN_SET_STRING: {
      std::shared_ptr<DiscreteSetRandomVariable<String>> rv_rep =
	std::static_pointer_cast<DiscreteSetRandomVariable<String>>
	(rv_i.random_variable_rep());
      StringRealMap sr_map;  rv_rep->pull_parameter(DUSS_VALUES_PROBS, sr_map);
      int map_size = sr_map.size();
      if (map_size > 1) {
	RealArray x_val, y_val;  map_indices_to_xy_pdf(sr_map, x_val, y_val);
	lhs_udist_register("DiscUncSetS","discrete histogram", av_cntr,
			   x_val, y_val);
      }
      else if (map_size)
	lhs_const_register("DiscUncSetS", av_cntr, 0.);
      break;
    }
    case DISCRETE_UNCERTAIN_SET_REAL: {
      std::shared_ptr<DiscreteSetRandomVariable<Real>> rv_rep =
	std::static_pointer_cast<DiscreteSetRandomVariable<Real>>
	(rv_i.random_variable_rep());
      RealRealMap rr_map;  rv_rep->pull_parameter(DUSR_VALUES_PROBS, rr_map);
      size_t map_size = rr_map.size();
      if (map_size > 1) {
	RealArray x_val, y_val;  map_to_xy_pdf(rr_map, x_val, y_val);
	lhs_udist_register("DiscUncSetR","discrete histogram", av_cntr,
			   x_val, y_val);
      }
      else if (map_size)
	lhs_const_register("DiscUncSetR", av_cntr, rr_map.begin()->first);
      break;
    }
    }
    ++av_cntr;
  }

  /////////////////////////////////////////
  // Register correlations with lhs_corr //
  /////////////////////////////////////////
  // specify the rank correlations among the RV.  Only non-zero values in the
  // lower triangular portion of the rank correlation matrix are specified.
  if (correlation_flag) {
    // Spec order: {cdv, ddv}, {cauv, dauv, corr}, {ceuv, deuv}, {csv, dsv}
    // > pass in bit array for active RV's to sample + another for active corr's
    // > Default empty arrays --> all RVs active; corr matrix applies to all RVs
    size_t j, ac_cntr_i, ac_cntr_j, av_cntr_i, av_cntr_j;
    bool av_i, av_j, ac_i, ac_j;
    for (i=0, ac_cntr_i=0, av_cntr_i=0; i<num_rv; ++i) {
      av_i = (!subset_rv   || active_vars[i]);
      ac_i = (!subset_corr || active_corr[i]);
      if (av_i && ac_i) {
	for (j=0, ac_cntr_j=0, av_cntr_j=0; j<i; ++j) {
	  av_j = (!subset_rv   || active_vars[j]);
	  ac_j = (!subset_corr || active_corr[j]);
	  if (av_j && ac_j) {
	    Real corr_val = corr(ac_cntr_i, ac_cntr_j);
	    if (std::abs(corr_val) > SMALL_NUMBER) {
	      LHS_CORR2_FC(const_cast<char*>(lhsNames[av_cntr_i].data()),
			   const_cast<char*>(lhsNames[av_cntr_j].data()),
			   corr_val, err_code);
	      check_error(err_code, "lhs_corr");
	    }
	  }
	  if (av_j) ++av_cntr_j;
	  if (ac_j) ++ac_cntr_j;
	}
      }
      if (av_i) ++av_cntr_i;
      if (ac_i) ++ac_cntr_i;
    }
  }

  /////////////////////
  // RUN THE SAMPLER //
  /////////////////////
  // perform internal checks on input to LHS
  int num_nam = num_rv, num_var = num_rv;
  LHS_PREP_FC(err_code, num_nam, num_var);
  check_error(err_code, "lhs_prep");

  // allocate the memory to hold samples, pt values, variable names, etc.
  int*      index_list = new  int    [num_nam];  // output
  Real*     ptval_list = new Real    [num_nam];  // output
  char* dist_name_list = new char [16*num_nam];  // output
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
  if (samples.numRows() != num_var || samples.numCols() != num_samples)
    samples.shapeUninitialized(num_var, num_samples);
  if (sampleRanksMode && sample_ranks.empty()) {
    if (sampleRanksMode == SET_RANKS || sampleRanksMode == SET_GET_RANKS) {
      PCerr << "Error: empty sample ranks array cannot be set in Pecos::"
	    << "LHSDriver::get_parameter_sets()" << std::endl;
      abort_handler(-1);
    }
    else if (sampleRanksMode == GET_RANKS)
      sample_ranks.shapeUninitialized(num_var, num_samples);
  }

  // generate the samples
  int rflag = sampleRanksMode; // short -> int
  LHS_RUN_FC(max_var, num_samples, num_nam, err_code, dist_name_list,
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
generate_unique_samples(const std::vector<RandomVariable>& random_vars,
			const RealSymMatrix& corr, int num_samples,
			RealMatrix& samples, RealMatrix& sample_ranks,
			const BitArray& active_vars,
			const BitArray& active_corr)
{
  // ***************************************************************************
  // IMPORTANT NOTE: this heuristic approach emphasizes having unique discrete
  // samples irregardless of their relative probability. As such, ANY STATISTICS
  // GENERATED FROM THESE SAMPLES WILL BE INVALID --> routine should only be
  // used for generating space-filling designs, e.g. for building surrogates
  // or approximating epistemic bounds.
  // ***************************************************************************
  // Approach 1 (current): accept a new sample if unique in cont+disc space.
  // > if mixed, first aggregate sample will be accepted with strata intact
  //   but discrete stratification will not be enhanced.
  // > if discrete only, iteration will occur irregardless of strata/probs.
  // Approach 2 (original): accept a new sample if unique in discrete subset.
  // > if mixed or discrete, iteration will occur irregardless of strata/probs
  //   and accepted samples are only from newly generated sets.
  // Approach 3 (TO DO?): preserve cont strata from a single cont+disc sample
  // --> backfill only the disc subset with generated sets of discrete samples.
  // > This approach preserves cont strata; but ignores disc strata/probs
  // ***************************************************************************

  // determine if the RV set has a finite set of possible sample points
  size_t i, num_rv = random_vars.size();//, num_finite_dv = 0;
  bool finite_combinations = true, no_mask = active_vars.empty();
  size_t num_av = (no_mask) ? num_rv : active_vars.count();
  // track discrete variables that have finite support
  for (i=0; i<num_rv; ++i) {
    if (no_mask || active_vars[i]) {
      switch (random_vars[i].type()) {
      case DISCRETE_RANGE:
      case DISCRETE_SET_INT: case DISCRETE_SET_STRING: case DISCRETE_SET_REAL:
      case BINOMIAL:         case HYPERGEOMETRIC:
      case HISTOGRAM_PT_INT: case HISTOGRAM_PT_STRING: case HISTOGRAM_PT_REAL:
      case DISCRETE_INTERVAL_UNCERTAIN:   case DISCRETE_UNCERTAIN_SET_INT:
      case DISCRETE_UNCERTAIN_SET_STRING: case DISCRETE_UNCERTAIN_SET_REAL:
	//++num_finite_dv; // if count needed, don't break out of for loop
	break;
      default: // any RV with countably/uncountably infinite support
	finite_combinations = false; break;
      }
    }
    if (!finite_combinations) break;
  }

  // if finite, compute number of possible discrete combinations
  // > for range variables, use ub-lb+1
  // > for BPA variables, overlay unique ranges
  // > for {set,map} variables, use {set,map}.size()
  // > for finite support in other discrete, compute #support pts
  size_t max_unique = _NPOS;
  //IntArray discrete_strata_1d; discrete_strata_1d.reserve(num_finite_dv);
  if (finite_combinations) {
    max_unique = 1;
    for (i=0; i<num_rv; ++i)
      if (no_mask || active_vars[i]) {
	const RandomVariable& rv_i = random_vars[i];
	switch (rv_i.type()) {
	// discrete design, state
	case DISCRETE_RANGE: {
	  int l_bnd;  rv_i.pull_parameter(DR_LWR_BND, l_bnd);
	  int u_bnd;  rv_i.pull_parameter(DR_UPR_BND, u_bnd);
	  //discrete_strata_1d.push_back(u_bnd - l_bnd + 1);
	  max_unique *= u_bnd - l_bnd + 1; break;
	}
	case DISCRETE_SET_INT: {
	  std::shared_ptr<SetVariable<int>> rv_rep =
	    std::static_pointer_cast<SetVariable<int>>
	    (rv_i.random_variable_rep());
	  IntSet i_set;  rv_rep->pull_parameter(DSI_VALUES, i_set);
	  //discrete_strata_1d.push_back(i_set.size());
	  max_unique *= i_set.size(); break;
	}
	case DISCRETE_SET_STRING: {
	  std::shared_ptr<SetVariable<String>> rv_rep =
	    std::static_pointer_cast<SetVariable<String>>
	    (rv_i.random_variable_rep());
	  StringSet s_set;  rv_rep->pull_parameter(DSS_VALUES, s_set);
	  //discrete_strata_1d.push_back(s_set.size());
	  max_unique *= s_set.size();  break;
	}
	case DISCRETE_SET_REAL: {
	  std::shared_ptr<SetVariable<Real>> rv_rep =
	    std::static_pointer_cast<SetVariable<Real>>
	    (rv_i.random_variable_rep());
	  RealSet r_set;  rv_rep->pull_parameter(DSR_VALUES, r_set);
	  //discrete_strata_1d.push_back(r_set.size());
	  max_unique *= r_set.size();  break;
	}
	// discrete aleatory uncertain
	case BINOMIAL: { // finite support
	  unsigned int num_tr;  rv_i.pull_parameter(BI_TRIALS, num_tr);
	  //discrete_strata_1d.push_back(1 + num_tr);
	  max_unique *= 1 + num_tr;  break;
	}
	case HYPERGEOMETRIC: { // finite support
	  unsigned int tot_p;  rv_i.pull_parameter(HGE_TOT_POP, tot_p);
	  unsigned int sel_p;  rv_i.pull_parameter(HGE_SEL_POP, sel_p);
	  unsigned int num_d;  rv_i.pull_parameter(HGE_DRAWN,   num_d);
	  //discrete_strata_1d.push_back(1 + std::min(num_d, sel_p) -
	  //			        std::max(0, sel_p + num_d - tot_p));
	  max_unique *= 1 + std::min(num_d, sel_p) + tot_p
	    - std::max(tot_p, sel_p + num_d); // care with unsigned
	  break;
	}
	case HISTOGRAM_PT_INT: {
	  std::shared_ptr<DiscreteSetRandomVariable<int>> rv_rep =
	    std::static_pointer_cast<DiscreteSetRandomVariable<int>>
	    (rv_i.random_variable_rep());
	  IntRealMap ir_map;  rv_rep->pull_parameter(H_PT_INT_PAIRS, ir_map);
	  //discrete_strata_1d.push_back(ir_map.size());
	  max_unique *= ir_map.size();   break;
	}
	case HISTOGRAM_PT_STRING: {
	  std::shared_ptr<DiscreteSetRandomVariable<String>> rv_rep =
	    std::static_pointer_cast<DiscreteSetRandomVariable<String>>
	    (rv_i.random_variable_rep());
	  StringRealMap sr_map;  rv_rep->pull_parameter(H_PT_STR_PAIRS, sr_map);
	  //discrete_strata_1d.push_back(sr_map.size());
	  max_unique *= sr_map.size();  break;
	}
	case HISTOGRAM_PT_REAL: {
	  std::shared_ptr<DiscreteSetRandomVariable<Real>> rv_rep =
	    std::static_pointer_cast<DiscreteSetRandomVariable<Real>>
	    (rv_i.random_variable_rep());
	  RealRealMap rr_map;  rv_rep->pull_parameter(H_PT_REAL_PAIRS, rr_map);
	  //discrete_strata_1d.push_back(rr_map.size());
	  max_unique *= rr_map.size();  break;
	}
	// discrete epistemic uncertain
	case DISCRETE_INTERVAL_UNCERTAIN: {
	  std::shared_ptr<IntervalRandomVariable<int>> rv_rep =
	    std::static_pointer_cast<IntervalRandomVariable<int>>
	    (rv_i.random_variable_rep());
	  IntIntPairRealMap di_bpa;  rv_rep->pull_parameter(DIU_BPA, di_bpa);
	  // x_sort_unique contains ALL of the unique integer values for this
	  // discrete interval variable in increasing order.  For example, if
	  // there are 3 intervals for a variable and the bounds are (1,4),
	  // (3,6), and [9,10], x_sorted will be (1,2,3,4,5,6,9,10).
	  IntSet x_sort_unique;
	  for (IIPRMCIter cit=di_bpa.begin(); cit!=di_bpa.end(); ++cit) {
	    const RealRealPair& bounds = cit->first;
	    int val, u_bnd = bounds.second;
	    for (val=bounds.first; val<=u_bnd; ++val)
	      x_sort_unique.insert(val);
	  }
	  //discrete_strata_1d.push_back(x_sort_unique.size());
	  max_unique *= x_sort_unique.size();  break;
	}
	case DISCRETE_UNCERTAIN_SET_INT: {
	  std::shared_ptr<DiscreteSetRandomVariable<int>> rv_rep =
	    std::static_pointer_cast<DiscreteSetRandomVariable<int>>
	    (rv_i.random_variable_rep());
	  IntRealMap ir_map;  rv_rep->pull_parameter(DUSI_VALUES_PROBS, ir_map);
	  //discrete_strata_1d.push_back(ir_map.size());
	  max_unique *= ir_map.size();  break;
	}
	case DISCRETE_UNCERTAIN_SET_STRING: {
	  std::shared_ptr<DiscreteSetRandomVariable<String>> rv_rep =
	    std::static_pointer_cast<DiscreteSetRandomVariable<String>>
	    (rv_i.random_variable_rep());
	  StringRealMap sr_map;
	  rv_rep->pull_parameter(DUSS_VALUES_PROBS, sr_map);
	  //discrete_strata_1d.push_back(sr_map.size());
	  max_unique *= sr_map.size();    break;
	}
	case DISCRETE_UNCERTAIN_SET_REAL: {
	  std::shared_ptr<DiscreteSetRandomVariable<Real>> rv_rep =
	    std::static_pointer_cast<DiscreteSetRandomVariable<Real>>
	    (rv_i.random_variable_rep());
	  RealRealMap rr_map; rv_rep->pull_parameter(DUSR_VALUES_PROBS, rr_map);
	  //discrete_strata_1d.push_back(rr_map.size());
	  max_unique *= rr_map.size();    break;
	}
	}
      }
  }

  bool complete = false;
  RealMatrix new_samples;  RealArray new_samp_i(num_av);
  std::set<RealArray> unique_samples;
  size_t unique_cntr = 0, iter = 0, max_iter = 1000;

  // If the number of samples requested is greater than the maximum possible
  // number of discrete samples then we must allow replicates of the discrete
  // variables to obtain the desired number of variables.  If not then we can
  // proceed with generating a unique set of discrete samples.
  if (max_unique < num_samples) {

    // Allow iteration to continue if new > old
    PCout << "Warning: LHS backfill was requested, but the discrete variables "
	  << "provided\n         do not have enough unique values ("
	  << max_unique << ") to obtain the number of\n         samples "
	  << "requested.  Backfill iterations will be attempted to increase\n"
	  << "         the number of unique samples until no further progress "
	  << "is detected." << std::endl;

    size_t unique_cntr_prev = 0;
    while (!complete && iter < max_iter) {
      generate_samples(random_vars, corr, num_samples, new_samples,
		       sample_ranks, active_vars, active_corr);

      // Sort real-valued samples and replace until # unique >= # requested.
      // > For now, don't try to replace only the discrete subset.
      // > Preserve original sample ordering --> don't truncate an ordered set
      //   (omitting tail of ordered set could bias coverage).
      for (i=0; i<num_samples; ++i) {
	const Real* samp_i = new_samples[i];
	if (test_unique(random_vars, active_vars, samp_i, unique_samples)) {
	  copy_data(samp_i, num_av, samples[unique_cntr++]); // append
	  if (unique_cntr >= num_samples)
	    { complete = true; break; }
	}
      }
      ++iter;
      if (unique_cntr == unique_cntr_prev)  complete = true;
      else               unique_cntr_prev = unique_cntr;
    }
  }
  else {
    if (samples.numRows() != num_av || samples.numCols() != num_samples)
      samples.shapeUninitialized(num_av, num_samples);

    // Currently sample_ranks will always be returned empty. It should only be
    // filled when NonDSampling::sampleRanksMode > 0. But I cannot see anywhere
    // in the code where this is true.
    //sample_ranks.shapeUninitialized( num_vars, num_samples );

    /*
    // unique index of all discrete variables if any
    std::set<RealArray>::iterator it;
    std::set<RealArray> sorted_discrete_samples; 
    RealArray discrete_sample( num_discrete_vars );

    // Determine the columns in samples_rm that contain discrete variables
    IntVector discrete_samples_map( num_discrete_vars, false );
    for (int i=0; i<num_discrete_vars; i++)
      discrete_samples_map[i]=num_continuous_vars+i;

    int num_unique_samples = 0;
    */

    // Eliminate redundant samples by resampling if necessary.  Could pad
    // num_samples in anticipation of duplicates, but this would alter LHS
    // stratification that could be intended, so use num_samples for now.
    while (!complete && iter < max_iter) {
      generate_samples(random_vars, corr, num_samples, new_samples,
		       sample_ranks, active_vars, active_corr);

      // Sort real-valued samples and replace until # unique >= # requested.
      // > For now, don't try to replace only the discrete subset.
      // > Preserve original sample ordering --> don't truncate an ordered set
      //   (omitting tail of ordered set could bias coverage).
      for (i=0; i<num_samples; ++i) {
	const Real* samp_i = new_samples[i];
	if (test_unique(random_vars, active_vars, samp_i, unique_samples)) {
	  copy_data(samp_i, num_av, samples[unique_cntr++]); // append
	  if (unique_cntr >= num_samples)
	    { complete = true; break; }
	}
      }
      ++iter;

      ////////////////////////////////////////
      /*
      if (initial) { // pack initial sample set
	for (i=0; i<num_samples; ++i) { // or matrix->set<vector> ?
	  //PCout << "[";
	  for (int j=0; j<num_discrete_vars; j++) {
	    int index = discrete_samples_map[j];
	    discrete_sample[j] = samples_rm(index,i);
	    //PCout << discrete_sample[j] << ",";
	  }
	  //PCout << "]\n";
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
	//PCout << num_unique_samples << "," << sorted_discrete_samples.size()
	//      << "," << num_discrete_vars << "," << num_vars << ","
	//      << num_continuous_vars << std::endl;

	for (i=0; i<num_samples; ++i) {
	  if (num_unique_samples < num_samples) {
	    //PCout << "[";
	    for (int j=0; j<num_discrete_vars; j++) {
	      int index = discrete_samples_map[j];
	      discrete_sample[j] = samples_rm(index,i);
	      //PCout << discrete_sample[j] << ",";
	    }
	    //PCout << "]\n";
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
      */
      ////////////////////////////////////////
    }
  }
  if (unique_cntr < num_samples) {
    PCerr << "Warning: iterations completed with number of unique samples ("
	  << unique_cntr << ")\n         less than target (" << num_samples
	  << ")." << std::endl;
    samples.reshape(num_av, unique_cntr);
  }
}


bool LHSDriver::
test_unique(const std::vector<RandomVariable>& random_vars,
	    const BitArray& active_vars, const Real* new_samp,
	    std::set<RealArray>& unique_samples)
{
  bool full_match = false; // hardwire this for right now

  size_t num_rv = random_vars.size();
  bool  no_mask = active_vars.empty();
  RealArray new_samp_v;
  if (full_match) {
    size_t num_av = (no_mask) ? num_rv : active_vars.count();
    new_samp_v.resize(num_av);
    copy_data(new_samp, num_av, new_samp_v);  // vector for sorting
  }
  else { // check for discrete match
    size_t i, cntr = 0;
    for (i=0; i<num_rv; ++i)
      if (no_mask || active_vars[i]) {
	switch (random_vars[i].type()) {
	case DISCRETE_RANGE:
	case DISCRETE_SET_INT: case DISCRETE_SET_STRING: case DISCRETE_SET_REAL:
	case POISSON: case BINOMIAL: case NEGATIVE_BINOMIAL:
	case GEOMETRIC: case HYPERGEOMETRIC: case HISTOGRAM_PT_INT:
	case HISTOGRAM_PT_STRING:            case HISTOGRAM_PT_REAL:
	case DISCRETE_INTERVAL_UNCERTAIN:    case DISCRETE_UNCERTAIN_SET_INT:
	case DISCRETE_UNCERTAIN_SET_STRING:  case DISCRETE_UNCERTAIN_SET_REAL:
	  new_samp_v.push_back(new_samp[cntr]);  break;
	}
	++cntr;
      }
  }
  // returns pair<iterator,bool>
  return unique_samples.insert(new_samp_v).second;
}

} // namespace Pecos
