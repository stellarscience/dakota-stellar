/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 RandomVariable
//- Description: Implementation code for RandomVariable class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "RandomVariable.hpp"
#include "RangeVariable.hpp"
#include "SetVariable.hpp"
#include "BoundedNormalRandomVariable.hpp"
#include "BoundedLognormalRandomVariable.hpp"
#include "LoguniformRandomVariable.hpp"
#include "TriangularRandomVariable.hpp"
#include "BetaRandomVariable.hpp"
#include "GammaRandomVariable.hpp"
#include "InvGammaRandomVariable.hpp"
#include "GumbelRandomVariable.hpp"
#include "FrechetRandomVariable.hpp"
#include "WeibullRandomVariable.hpp"
#include "HistogramBinRandomVariable.hpp"
#include "PoissonRandomVariable.hpp"
#include "BinomialRandomVariable.hpp"
#include "NegBinomialRandomVariable.hpp"
#include "GeometricRandomVariable.hpp"
#include "HypergeometricRandomVariable.hpp"
#include "DiscreteSetRandomVariable.hpp"
#include "IntervalRandomVariable.hpp"

static const char rcsId[]="@(#) $Id: RandomVariable.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {


/** This constructor is the one which must build the base class data
    for all derived classes.  get_random_variable() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_random_variable() again).  Since the
    letter IS the representation, its rep pointer is set to NULL. */
RandomVariable::RandomVariable(BaseConstructor)
{ /* empty ctor */ }


/** The default constructor: ranVarRep is NULL in this case. */
RandomVariable::RandomVariable()
{ /* empty ctor */ }


/** Envelope constructor only needs to extract enough data to properly
    execute get_random_variable, since RandomVariable(BaseConstructor)
    builds the actual base class data for the derived basis functions. */
RandomVariable::RandomVariable(short ran_var_type):
  // Set the rep pointer to the appropriate derived type
  ranVarRep(get_random_variable(ran_var_type))
{
  if ( !ranVarRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize ranVarRep to the 
    appropriate derived type. */
std::shared_ptr<RandomVariable>
RandomVariable::get_random_variable(short ran_var_type)
{
  std::shared_ptr<RandomVariable> ran_var_rep;
  switch (ran_var_type) {
  // continuous / discrete range / set variables
  case CONTINUOUS_RANGE:
    ran_var_rep = std::make_shared<RangeVariable<Real>>();                break;
  case DISCRETE_RANGE:
    ran_var_rep = std::make_shared<RangeVariable<int>>();                 break;
  case DISCRETE_SET_INT:
    ran_var_rep = std::make_shared<SetVariable<int>>();                   break;
  case DISCRETE_SET_STRING:
    ran_var_rep = std::make_shared<SetVariable<String>>();                break;
  case DISCRETE_SET_REAL:
    ran_var_rep = std::make_shared<SetVariable<Real>>();                  break;
  // continuous aleatory random variables:
  case STD_NORMAL: case NORMAL:
    ran_var_rep = std::make_shared<NormalRandomVariable>();               break;
  case BOUNDED_NORMAL:
    ran_var_rep = std::make_shared<BoundedNormalRandomVariable>();        break;
  case LOGNORMAL:
    ran_var_rep = std::make_shared<LognormalRandomVariable>();            break;
  case BOUNDED_LOGNORMAL:
    ran_var_rep = std::make_shared<BoundedLognormalRandomVariable>();     break;
  case STD_UNIFORM: case UNIFORM:
    ran_var_rep = std::make_shared<UniformRandomVariable>();              break;
  case LOGUNIFORM:
    ran_var_rep = std::make_shared<LoguniformRandomVariable>();           break;
  case TRIANGULAR:
    ran_var_rep = std::make_shared<TriangularRandomVariable>();           break;
  case STD_EXPONENTIAL: case EXPONENTIAL:
    ran_var_rep = std::make_shared<ExponentialRandomVariable>();          break;
  case STD_BETA: case BETA:
    ran_var_rep = std::make_shared<BetaRandomVariable>();                 break;
  case STD_GAMMA: case GAMMA:
    ran_var_rep = std::make_shared<GammaRandomVariable>();                break;
  case GUMBEL:
    ran_var_rep = std::make_shared<GumbelRandomVariable>();               break;
  case FRECHET: 
    ran_var_rep = std::make_shared<FrechetRandomVariable>();              break;
  case WEIBULL:
    ran_var_rep = std::make_shared<WeibullRandomVariable>();              break;
  case HISTOGRAM_BIN:
    ran_var_rep = std::make_shared<HistogramBinRandomVariable>();         break;
  // hyper-parameter distributions:
  case INV_GAMMA:
    ran_var_rep = std::make_shared<InvGammaRandomVariable>();             break;
  // discrete aleatory random variables:
  case POISSON:
    ran_var_rep = std::make_shared<PoissonRandomVariable>();              break;
  case BINOMIAL:
    ran_var_rep = std::make_shared<BinomialRandomVariable>();             break;
  case NEGATIVE_BINOMIAL:
    ran_var_rep = std::make_shared<NegBinomialRandomVariable>();          break;
  case GEOMETRIC:
    ran_var_rep = std::make_shared<GeometricRandomVariable>();            break;
  case HYPERGEOMETRIC:
    ran_var_rep = std::make_shared<HypergeometricRandomVariable>();       break;
  // continuous / discrete epistemic intervals: distinct from HistogramBin
  // in ability to support overlapping/disjoint intervals
  case CONTINUOUS_INTERVAL_UNCERTAIN:
    ran_var_rep = std::make_shared<IntervalRandomVariable<Real>>();       break;
  case DISCRETE_INTERVAL_UNCERTAIN:
    ran_var_rep = std::make_shared<IntervalRandomVariable<int>>();        break;
  // aleatory / epistemic sets: distinct only in interpretation of set probs
  // (statistical expectations should not be used in epistemic case)
  case HISTOGRAM_PT_INT:    case DISCRETE_UNCERTAIN_SET_INT:
    ran_var_rep = std::make_shared<DiscreteSetRandomVariable<int>>();     break;
  case HISTOGRAM_PT_STRING: case DISCRETE_UNCERTAIN_SET_STRING:
    ran_var_rep = std::make_shared<DiscreteSetRandomVariable<String>>();  break;
  case HISTOGRAM_PT_REAL:   case DISCRETE_UNCERTAIN_SET_REAL:
    ran_var_rep = std::make_shared<DiscreteSetRandomVariable<Real>>();    break;
  default:
    PCerr << "Error: RandomVariable type " << ran_var_type << " not available."
	  << std::endl;
  }

  // some derived classes (especially template classes) cover multiple
  // ranVarTypes, so override ctor assignments for those cases:
  if (ran_var_rep)
    ran_var_rep->ranVarType = ran_var_type;

  return ran_var_rep;
}


/** Copy constructor manages sharing of ranVarRep. */
RandomVariable::RandomVariable(const RandomVariable& ran_var):
  ranVarRep(ran_var.ranVarRep)
{ /* empty ctor */ }


RandomVariable RandomVariable::operator=(const RandomVariable& ran_var)
{
  ranVarRep = ran_var.ranVarRep;
  return *this; // calls copy constructor since returned by value
}


RandomVariable::~RandomVariable()
{ /* empty dtor */ }


Real RandomVariable::cdf(Real x) const
{
  if (!ranVarRep) {
    PCerr << "Error: cdf() not supported for this random variable type ("
	  << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->cdf(x); // forward to letter
}


Real RandomVariable::ccdf(Real x) const
{
  if (ranVarRep)
    return ranVarRep->ccdf(x); // forward to letter
  else
    return 1. - cdf(x); // default (overriden in most cases)
}


Real RandomVariable::inverse_cdf(Real p_cdf) const
{
  if (!ranVarRep) {
    PCerr << "Error: inverse_cdf() not supported for this random variable "
	  << "type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->inverse_cdf(p_cdf); // forward to letter
}


Real RandomVariable::inverse_ccdf(Real p_ccdf) const
{
  if (ranVarRep)
    return ranVarRep->inverse_ccdf(p_ccdf); // forward to letter
  else
    return inverse_cdf(1. - p_ccdf); // default (overriden in most cases)
}


Real RandomVariable::pdf(Real x) const
{
  if (!ranVarRep) {
    PCerr << "Error: pdf() not supported for this random variable type ("
	  << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->pdf(x); // forward to letter
}


Real RandomVariable::pdf_gradient(Real x) const
{
  if (!ranVarRep) {
    PCerr << "Error: pdf_gradient() not supported for this random variable "
	  << "type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->pdf_gradient(x); // forward to letter
}


Real RandomVariable::pdf_hessian(Real x) const
{
  if (!ranVarRep) {
    PCerr << "Error: pdf_hessian() not supported for this random variable "
	  << "type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->pdf_hessian(x); // forward to letter
}


Real RandomVariable::log_pdf(Real x) const
{
  if (ranVarRep)
    return ranVarRep->log_pdf(x); // forward to letter
  else
    return std::log(pdf(x)); // default (overriden for exponential pdf cases)
}


Real RandomVariable::log_pdf_gradient(Real x) const
{
  if (ranVarRep)
    return ranVarRep->log_pdf_gradient(x); // forward to letter
  else // nabla [log pi] = nabla pi / pi
    return pdf_gradient(x) / pdf(x); // default
}


Real RandomVariable::log_pdf_hessian(Real x) const
{
  if (ranVarRep)
    return ranVarRep->log_pdf_hessian(x); // forward to letter
  else {
    // nabla^2 [log pi] = nabla [ nabla pi / pi ]
    //   = nabla^2 pi / pi - ( nabla pi / pi )^2
    Real val = pdf(x), grad_val_ratio = pdf_gradient(x) / val;
    return pdf_hessian(x) / val - grad_val_ratio * grad_val_ratio; // default
  }
}


Real RandomVariable::inverse_standard_cdf(Real p_cdf) const
{
  if (!ranVarRep) {
    PCerr << "Error: inverse_standard_cdf() not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->inverse_standard_cdf(p_cdf); // forward to letter
}


Real RandomVariable::standard_pdf(Real z) const
{
  if (!ranVarRep) {
    PCerr << "Error: standard_pdf() not supported for this random variable "
	  << "type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->standard_pdf(z); // forward to letter
}


Real RandomVariable::log_standard_pdf(Real z) const
{
  if (ranVarRep) // forward to letter
    return ranVarRep->log_standard_pdf(z);
  else // default (overriden for exponential standard_pdf)
    return std::log(standard_pdf(z));
}


Real RandomVariable::log_standard_pdf_gradient(Real z) const
{
  if (!ranVarRep) {
    PCerr << "Error: log_standard_pdf_gradient() not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  // default usage of standard_pdf() and standard_pdf_gradient() does
  // not need to be supported since all standard PDF cases implement
  // log_standard_pdf_gradient()
  return ranVarRep->log_standard_pdf_gradient(z); // forward to letter
}


Real RandomVariable::log_standard_pdf_hessian(Real z) const
{
  if (!ranVarRep) {
    PCerr << "Error: log_standard_pdf_hessian() not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  // default usage of standard_pdf{,_gradient,_hessian}() does not
  // need to be supported since all standard PDF cases implement
  // log_standard_pdf_hessian()
  return ranVarRep->log_standard_pdf_hessian(z); // forward to letter
}


Real RandomVariable::to_standard(Real x) const
{
  if (!ranVarRep) {
    PCerr << "Error: to_standard() not supported for this random variable "
	  << "type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->to_standard(x); // forward to letter
}


Real RandomVariable::from_standard(Real x) const
{
  if (!ranVarRep) {
    PCerr << "Error: from_standard() not supported for this random variable "
	  << "type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->from_standard(x); // forward to letter
}


void RandomVariable::pull_parameter(short dist_param, Real& val) const
{
  if (ranVarRep)
    ranVarRep->pull_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: pull_parameter(Real) not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::pull_parameter(short dist_param, int& val) const
{
  if (ranVarRep)
    ranVarRep->pull_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: pull_parameter(int) not supported for this "
	  << "random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::pull_parameter(short dist_param, unsigned int& val) const
{
  if (ranVarRep)
    ranVarRep->pull_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: pull_parameter(unsigned int) not supported for this "
	  << "random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::pull_parameter(short dist_param, IntSet& val) const
{
  if (ranVarRep)
    ranVarRep->pull_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: pull_parameter(IntSet) not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::pull_parameter(short dist_param, StringSet& val) const
{
  if (ranVarRep)
    ranVarRep->pull_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: pull_parameter(StringSet) not supported for this "
	  << "random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::pull_parameter(short dist_param, RealSet& val) const
{
  if (ranVarRep)
    ranVarRep->pull_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: pull_parameter(RealSet) not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::pull_parameter(short dist_param, IntRealMap& val) const
{
  if (ranVarRep)
    ranVarRep->pull_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: pull_parameter(IntRealMap) not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::pull_parameter(short dist_param, StringRealMap& val) const
{
  if (ranVarRep)
    ranVarRep->pull_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: pull_parameter(StringRealMap) not supported for this "
	  << "random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::pull_parameter(short dist_param, RealRealMap& val) const
{
  if (ranVarRep)
    ranVarRep->pull_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: pull_parameter(RealRealMap) not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::
pull_parameter(short dist_param, IntIntPairRealMap& val) const
{
  if (ranVarRep)
    ranVarRep->pull_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: pull_parameter(IntIntPairRealMap) not supported for this "
	  << "random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::
pull_parameter(short dist_param, RealRealPairRealMap& val) const
{
  if (ranVarRep)
    ranVarRep->pull_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: pull_parameter(RealRealPairRealMap) not supported for "
	  << "this random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::push_parameter(short dist_param, Real val)
{
  if (ranVarRep)
    ranVarRep->push_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: push_parameter(Real) not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::push_parameter(short dist_param, int val)
{
  if (ranVarRep)
    ranVarRep->push_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: push_parameter(int) not supported for this "
	  << "random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::push_parameter(short dist_param, unsigned int val)
{
  if (ranVarRep)
    ranVarRep->push_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: push_parameter(unsigned int) not supported for this "
	  << "random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::push_parameter(short dist_param, const IntSet& val)
{
  if (ranVarRep)
    ranVarRep->push_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: push_parameter(IntSet) not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::push_parameter(short dist_param, const StringSet& val)
{
  if (ranVarRep)
    ranVarRep->push_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: push_parameter(StringSet) not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::push_parameter(short dist_param, const RealSet& val)
{
  if (ranVarRep)
    ranVarRep->push_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: push_parameter(RealSet) not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::push_parameter(short dist_param, const IntRealMap& val)
{
  if (ranVarRep)
    ranVarRep->push_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: push_parameter(IntRealMap) not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::push_parameter(short dist_param, const StringRealMap& val)
{
  if (ranVarRep)
    ranVarRep->push_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: push_parameter(StringRealMap) not supported for this "
	  << "random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::push_parameter(short dist_param, const RealRealMap& val)
{
  if (ranVarRep)
    ranVarRep->push_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: push_parameter(RealRealMap) not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::
push_parameter(short dist_param, const IntIntPairRealMap& val)
{
  if (ranVarRep)
    ranVarRep->push_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: push_parameter(IntIntPairRealMap) not supported for this "
	  << "random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::
push_parameter(short dist_param, const RealRealPairRealMap& val)
{
  if (ranVarRep)
    ranVarRep->push_parameter(dist_param, val); // forward to letter
  else {
    PCerr << "Error: push_parameter(RealRealPairRealMap) not supported for "
	  << "this random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


void RandomVariable::copy_parameters(const RandomVariable& rv)
{
  if (ranVarRep)
    ranVarRep->copy_parameters(rv);
  else {
    PCerr << "Error: copy_parameters(RandomVariable) not supported for this "
	  << "random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
}


Real RandomVariable::mean() const
{
  if (!ranVarRep) {
    PCerr << "Error: mean() not supported for this random variable type ("
	  << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->mean(); // forward to letter
}


Real RandomVariable::median() const
{
  if (ranVarRep)
    return ranVarRep->median(); // forward to letter
  else
    return inverse_cdf(.5);
}


Real RandomVariable::mode() const
{
  if (!ranVarRep) {
    PCerr << "Error: mode() not supported for this random variable type ("
	  << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->mode(); // forward to letter
}


Real RandomVariable::standard_deviation() const
{
  if (ranVarRep)
    return ranVarRep->standard_deviation(); // forward to letter
  else
    return std::sqrt(variance());
}


Real RandomVariable::variance() const
{
  if (!ranVarRep) {
    PCerr << "Error: variance() not supported for this random variable type ("
	  << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->variance(); // forward to letter
}


RealRealPair RandomVariable::moments() const
{
  if (ranVarRep)
    return ranVarRep->moments(); // forward to letter
  else
    return RealRealPair(mean(), standard_deviation()); // default used by most
}


RealRealPair RandomVariable::distribution_bounds() const
{
  if (!ranVarRep) {
    PCerr << "Error: distribution_bounds() not supported for this random "
	  << "variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->distribution_bounds(); // forward to letter
}


void RandomVariable::lower_bound(Real l_bnd)
{
  if (ranVarRep)
    ranVarRep->lower_bound(l_bnd); // forward to letter
  // else {
  //   PCerr << "Error: lower_bound(Real) not supported for this random "
  // 	  << "variable type (" << ranVarType << ")." << std::endl;
  //   abort_handler(-1);
  // }
}


void RandomVariable::lower_bound(int l_bnd)
{
  if (ranVarRep)
    ranVarRep->lower_bound(l_bnd); // forward to letter
  // else {
  //   PCerr << "Error: lower_bound(int) not supported for this random "
  // 	  << "variable type (" << ranVarType << ")." << std::endl;
  //   abort_handler(-1);
  // }
}


void RandomVariable::upper_bound(Real u_bnd)
{
  if (ranVarRep)
    ranVarRep->upper_bound(u_bnd); // forward to letter
  // else {
  //   PCerr << "Error: upper_bound(Real) not supported for this random "
  // 	  << "variable type (" << ranVarType << ")." << std::endl;
  //   abort_handler(-1);
  // }
}


void RandomVariable::upper_bound(int u_bnd)
{
  if (ranVarRep)
    ranVarRep->upper_bound(u_bnd); // forward to letter
  // else {
  //   PCerr << "Error: upper_bound(int) not supported for this random "
  // 	  << "variable type (" << ranVarType << ")." << std::endl;
  //   abort_handler(-1);
  // }
}


Real RandomVariable::coefficient_of_variation() const
{
  if (ranVarRep)
    return ranVarRep->coefficient_of_variation(); // forward to letter
  else {
    RealRealPair moms = moments();
    return moms.second / moms.first; // default used by most
  }
}


Real RandomVariable::
correlation_warping_factor(const RandomVariable& rv, Real corr) const
{
  if (!ranVarRep) {
    PCerr << "Error: correlation_warping_factor() not supported for this "
	  << "random variable type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->correlation_warping_factor(rv, corr); // forward to letter
}


Real RandomVariable::dx_ds(short dist_param, short u_type, Real x, Real z) const
{
  if (!ranVarRep) {
    PCerr << "Error: dx_ds() not supported for this random variable type ("
	  << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->dx_ds(dist_param, u_type, x, z); // forward to letter
}


Real RandomVariable::dz_ds_factor(short u_type, Real x, Real z) const
{
  if (!ranVarRep) {
    PCerr << "Error: dz_ds_factor() not supported for this random variable "
	  << "type (" << ranVarType << ")." << std::endl;
    abort_handler(-1);
  }
  return ranVarRep->dz_ds_factor(u_type, x, z); // forward to letter
}

} // namespace Pecos
