/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_GLOBAL_DEFS_H
#define PECOS_GLOBAL_DEFS_H

#include <iostream>
#include <cfloat>  // for DBL_MIN, DBL_MAX
#include <cmath>
#include <cstdlib>

#include <boost/math/constants/constants.hpp>

namespace Pecos {

// --------------
// Special values
// --------------
/// the value for PI used in various numerical routines
const double PI = boost::math::constants::pi<double>();

/// special value returned by index() when entry not found
//const size_t _NPOS  = ~(size_t)0; // one's complement
const size_t SZ_MAX = std::numeric_limits<size_t>::max();
const size_t _NPOS  = SZ_MAX; // alias old definition

/// used in ostream data output functions
const int WRITE_PRECISION = 10;

/// small value used for protecting division by zero, etc.; an alternative
/// to DBL_MIN that is less likely to cause underflow/overflow when numbers
/// larger than it are used in calculations
const double SMALL_NUMBER = 1.e-25;
/// large value used as a surrogate for infinity in error traps; an alternative
/// to DBL_MAX or inf that is less likely to cause underflow/overflow when used
/// in subsequent calculations
const double LARGE_NUMBER = 1.e+50;

// define special values for vector/matrix data copying modes
enum { DEFAULT_COPY=0, SHALLOW_COPY, DEEP_COPY };

// define special values for ExpansionConfigOptions::outputLevel
enum { SILENT_OUTPUT, QUIET_OUTPUT, NORMAL_OUTPUT, VERBOSE_OUTPUT,
       DEBUG_OUTPUT };

// define special values for multivariate distribution types
enum { NO_DIST=0, MARGINALS_CORRELATIONS, MULTIVARIATE_NORMAL, JOINT_KDE };
     //GAUSSIAN_COPULA, ...

// define special values for type of moment results
enum { NO_MOMENTS=0, STANDARD_MOMENTS, CENTRAL_MOMENTS };

// define special values for univariate RandomVariable types
enum { NO_TYPE=0,
       // continuous + discrete design / state / other
       CONTINUOUS_RANGE, DISCRETE_RANGE,
       DISCRETE_SET_INT, DISCRETE_SET_STRING, DISCRETE_SET_REAL,
       STOCHASTIC_EXPANSION,
       // continuous + discrete aleatory
       STD_NORMAL, NORMAL, BOUNDED_NORMAL, LOGNORMAL, BOUNDED_LOGNORMAL,
       STD_UNIFORM, UNIFORM, LOGUNIFORM, TRIANGULAR, STD_EXPONENTIAL,
       EXPONENTIAL, STD_BETA, BETA, STD_GAMMA, GAMMA, INV_GAMMA,
       GUMBEL, FRECHET, WEIBULL, HISTOGRAM_BIN,
       POISSON, BINOMIAL, NEGATIVE_BINOMIAL, GEOMETRIC, HYPERGEOMETRIC,
       HISTOGRAM_PT_INT, HISTOGRAM_PT_STRING, HISTOGRAM_PT_REAL,
       // continuous + discrete epistemic
       // (Note: INTERVAL has multiple BPAs --> distinct from RANGE)
       CONTINUOUS_INTERVAL_UNCERTAIN, DISCRETE_INTERVAL_UNCERTAIN,
       DISCRETE_UNCERTAIN_SET_INT, DISCRETE_UNCERTAIN_SET_STRING,
       DISCRETE_UNCERTAIN_SET_REAL };

// define special values for secondaryACVarMapTargets/secondaryADVarMapTargets
enum { NO_TARGET=0,
       CR_LWR_BND, CR_UPR_BND, DR_LWR_BND, DR_UPR_BND,
       DSI_VALUES, DSS_VALUES, DSR_VALUES,
       N_MEAN, N_STD_DEV, N_LWR_BND, N_UPR_BND, N_LOCATION, N_SCALE, N_VARIANCE,
       LN_MEAN, LN_STD_DEV, LN_LAMBDA, LN_ZETA, LN_ERR_FACT,
       LN_LWR_BND, LN_UPR_BND,
       U_LWR_BND, U_UPR_BND, U_LOCATION, U_SCALE, LU_LWR_BND, LU_UPR_BND,
       T_MODE, T_LWR_BND, T_UPR_BND, T_LOCATION, T_SCALE, E_BETA, E_SCALE,
       BE_ALPHA, BE_BETA, BE_LWR_BND, BE_UPR_BND, JACOBI_ALPHA, JACOBI_BETA,
       GA_ALPHA, GA_BETA, GA_SHAPE, GA_SCALE, GENLAG_ALPHA, IGA_ALPHA, IGA_BETA,
       GU_ALPHA, GU_BETA, F_ALPHA, F_BETA, W_ALPHA, W_BETA, H_BIN_PAIRS,
       P_LAMBDA, BI_P_PER_TRIAL, BI_TRIALS, NBI_P_PER_TRIAL, NBI_TRIALS,
       GE_P_PER_TRIAL, HGE_TOT_POP, HGE_SEL_POP, HGE_DRAWN,
       H_PT_INT_PAIRS, H_PT_STR_PAIRS, H_PT_REAL_PAIRS, CIU_BPA, DIU_BPA,
       DUSI_VALUES_PROBS, DUSS_VALUES_PROBS, DUSR_VALUES_PROBS };

/// derived basis approximation types
enum { NO_BASIS=0, //FOURIER_BASIS, EIGEN_BASIS,
       GLOBAL_NODAL_INTERPOLATION_POLYNOMIAL,
       PIECEWISE_NODAL_INTERPOLATION_POLYNOMIAL,
       GLOBAL_HIERARCHICAL_INTERPOLATION_POLYNOMIAL,
       PIECEWISE_HIERARCHICAL_INTERPOLATION_POLYNOMIAL,
       GLOBAL_REGRESSION_ORTHOGONAL_POLYNOMIAL,
       GLOBAL_PROJECTION_ORTHOGONAL_POLYNOMIAL,
       GLOBAL_ORTHOGONAL_POLYNOMIAL };
       //PIECEWISE_REGRESSION_ORTHOGONAL_POLYNOMIAL,
       //PIECEWISE_PROJECTION_ORTHOGONAL_POLYNOMIAL,
       //PIECEWISE_ORTHOGONAL_POLYNOMIAL };

/// derived basis polynomial types (orthogonal polynomial order follows
/// uncertain variable spec order of normal, uniform, exponential, beta, gamma)
enum { NO_POLY=0, HERMITE_ORTHOG, LEGENDRE_ORTHOG, LAGUERRE_ORTHOG,
       JACOBI_ORTHOG, GEN_LAGUERRE_ORTHOG, CHEBYSHEV_ORTHOG, NUM_GEN_ORTHOG,
       LAGRANGE_INTERP, HERMITE_INTERP, PIECEWISE_LINEAR_INTERP,
       PIECEWISE_QUADRATIC_INTERP, PIECEWISE_CUBIC_INTERP,
       KRAWTCHOUK_DISCRETE, MEIXNER_DISCRETE, CHARLIER_DISCRETE,
       HAHN_DISCRETE };

/// integration rules within VPISparseGrid (1-12: CC through User-closed)
/// and beyond (GOLUB_WELSCH, NEWTON_COTES)
enum { NO_RULE=0, CLENSHAW_CURTIS, FEJER2, GAUSS_PATTERSON, GAUSS_LEGENDRE,
       GAUSS_HERMITE, GEN_GAUSS_HERMITE, GAUSS_LAGUERRE, GEN_GAUSS_LAGUERRE,
       GAUSS_JACOBI, GENZ_KEISTER, /*USER_OPEN, USER_CLOSED,*/ GOLUB_WELSCH,
       NEWTON_COTES, GAUSS_KRAWTCHOUK, GAUSS_MEIXNER, GAUSS_CHARLIER,
       GAUSS_HAHN };

// growth rules within VPISparseGrid
//enum { DEFAULT_GROWTH=0, SLOW_LINEAR, SLOW_LINEAR_ODD, MODERATE_LINEAR,
//       SLOW_EXPONENTIAL, MODERATE_EXPONENTIAL, FULL_EXPONENTIAL };

/// options for synchronizing linear and exponential growth rule settings
/// (consistent with slow/moderate/full growth for new level_to_growth_*
/// functions in sandia_rules.cpp)
enum { SLOW_RESTRICTED_GROWTH, MODERATE_RESTRICTED_GROWTH,
       UNRESTRICTED_GROWTH };

/// solution approaches for calculating the polynomial basis coefficients
/// (options for ExpansionConfigOptions::expCoeffsSolnApproach)
enum { QUADRATURE, CUBATURE, LIGHTWEIGHT_SPARSE_GRID, COMBINED_SPARSE_GRID,
       INCREMENTAL_SPARSE_GRID, HIERARCHICAL_SPARSE_GRID, SAMPLING,
       DEFAULT_REGRESSION, DEFAULT_LEAST_SQ_REGRESSION, SVD_LEAST_SQ_REGRESSION,
       EQ_CON_LEAST_SQ_REGRESSION, BASIS_PURSUIT, BASIS_PURSUIT_DENOISING,
       ORTHOG_MATCH_PURSUIT, LASSO_REGRESSION, LEAST_ANGLE_REGRESSION,
       ORTHOG_LEAST_INTERPOLATION };
/// options for BasisConfigOptions::nestingOverride (inactive)
enum { NO_NESTING_OVERRIDE=0, NESTED, NON_NESTED };
/// options for overriding the default growth restriction policy
enum { NO_GROWTH_OVERRIDE=0, RESTRICTED, UNRESTRICTED };
/// options for ExpansionConfigOptions::refineType (inactive)
enum { NO_REFINEMENT=0, P_REFINEMENT, H_REFINEMENT };
/// options for ExpansionConfigOptions::refineControl
// Note: C3 augments UNIFORM_CONTROL with a c3AdvancementType;
//       could consider a similar breakout for DIMENSION_ADAPTIVE_CONTROL
enum { NO_CONTROL=0, UNIFORM_CONTROL, LOCAL_ADAPTIVE_CONTROL,
       DIMENSION_ADAPTIVE_CONTROL_SOBOL, DIMENSION_ADAPTIVE_CONTROL_DECAY,
       DIMENSION_ADAPTIVE_CONTROL_GENERALIZED };
/// options for ExpansionConfigOptions::refineMetric
enum { NO_METRIC=0, COVARIANCE_METRIC, LEVEL_STATS_METRIC, MIXED_STATS_METRIC };
/// options for ExpansionConfigOptions::refineStatsType
enum { NO_EXPANSION_STATS=0, DEFAULT_EXPANSION_STATS, ACTIVE_EXPANSION_STATS,
       COMBINED_EXPANSION_STATS };

/// options for expansion basis type
enum { DEFAULT_BASIS=0, TENSOR_PRODUCT_BASIS, TOTAL_ORDER_BASIS,
       ADAPTED_BASIS_GENERALIZED, ADAPTED_BASIS_EXPANDING_FRONT,
       NODAL_INTERPOLANT, HIERARCHICAL_INTERPOLANT };

/// mode of integration driver: integration versus interpolation
enum { DEFAULT_MODE=0, INTEGRATION_MODE, INTERPOLATION_MODE };

/// options for local basis functions within PiecewiseInterpPolynomial
enum { LINEAR_EQUIDISTANT, LINEAR, QUADRATIC_EQUIDISTANT, QUADRATIC,
       CUBIC_EQUIDISTANT, CUBIC };

/// special values for nodal interpolation of variance and variance gradient
enum { INTERPOLATION_OF_PRODUCTS, REINTERPOLATION_OF_PRODUCTS,
       PRODUCT_OF_INTERPOLANTS_FAST, PRODUCT_OF_INTERPOLANTS_FULL };

/// special values for polynomial expansion combination
enum { NO_COMBINE=0, ADD_COMBINE, MULT_COMBINE, ADD_MULT_COMBINE };

/// special values for type of sequence across a set of model forms/resolutions
enum { NO_SEQUENCE=0, MODEL_FORM_SEQUENCE, RESOLUTION_LEVEL_SEQUENCE };

/// special values for type of discrepancy emulation
enum { NO_DISCREPANCY=0, DISTINCT_DISCREPANCY, RECURSIVE_DISCREPANCY };

/// special values for key type indicating augmentation of single data sets
/// with reduction data (e.g., for augmenting raw data with discrepancy data);
/// these are arranged to allow bit detection (raw = 1, reduction = 2, both = 3)
enum { NO_DATA=0, RAW_DATA, SINGLE_REDUCTION, RAW_WITH_REDUCTION };
       //, MULTIPLE_REDUCTION = 4 bit, ...

/// special values for filtering keyed SurrogateData maps
enum { NO_FILTER = 0, SINGLETON_FILTER, AGGREGATED_FILTER, RAW_DATA_FILTER,
       REDUCTION_DATA_FILTER, RAW_WITH_REDUCTION_DATA_FILTER };
       //, SINGLE_REDUCTION_DATA_FILTER, MULTIPLE_REDUCTION_DATA_FILTER, ...

// ----------------
// Standard streams
// ----------------
#define PCout std::cout
#define PCerr std::cerr


// --------------
// Global objects
// --------------
/// Dummy struct for overloading letter-envelope constructors.
/** BaseConstructor is used to overload the constructor for the base class
    portion of letter objects.  It avoids infinite recursion (Coplien p.139)
    in the letter-envelope idiom by preventing the letter from instantiating
    another envelope.  Putting this struct here avoids circular dependencies. */
struct BaseConstructor {
  BaseConstructor(int = 0) {} ///< C++ structs can have constructors
};


// ----------------
// Global functions
// ----------------

/// global function which handles serial or parallel aborts
void abort_handler(int code);


inline void abort_handler(int code)
{ std::exit(code); } // for now, prior to use of MPI


/** Templatized abort_handler_t method that allows for convenient return from 
    methods that otherwise have no sensible return from error clauses.  Usage:
    MyType& method() { return abort_handler<MyType&>(-1); } */
template <typename T>
T abort_handler_t(int code)
{
  abort_handler(code);
  throw code;
}

} // namespace Pecos

#endif // PECOS_GLOBAL_DEFS_H
