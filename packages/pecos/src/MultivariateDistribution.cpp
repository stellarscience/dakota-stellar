/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "pecos_stat_util.hpp"
#include "MarginalsCorrDistribution.hpp"
#include "MultivariateNormalDistribution.hpp"

static const char rcsId[]="@(#) $Id: MultivariateDistribution.cpp 4768 2007-12-17 17:49:32Z mseldre $";

namespace Pecos {


/** This constructor is the one which must build the base class data
    for all derived classes.  get_mv_dist() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_mv_dist() again).  Since the
    letter IS the representation, its rep pointer is set to NULL. */
MultivariateDistribution::MultivariateDistribution(BaseConstructor):
  correlationFlag(false)
{ /* empty ctor */ }


/** The default constructor: mvDistRep is NULL in this case. */
MultivariateDistribution::MultivariateDistribution()
{ /* empty ctor */ }


/** Envelope constructor only needs to extract enough data to properly
    execute get_mv_dist, since MultivariateDistribution(BaseConstructor)
    builds the actual base class data for the derived transformations. */
MultivariateDistribution::
MultivariateDistribution(short mv_dist_type):
  // Set the rep pointer to the appropriate derived type
  mvDistRep(get_distribution(mv_dist_type))
{
  if ( !mvDistRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize mvDistRep to the 
    appropriate derived type. */
std::shared_ptr<MultivariateDistribution>
MultivariateDistribution::get_distribution(short mv_dist_type) const
{
  std::shared_ptr<MultivariateDistribution> mvd_rep;
  switch (mv_dist_type) {
  case MARGINALS_CORRELATIONS:
    mvd_rep = std::make_shared<MarginalsCorrDistribution>();      break;
  case MULTIVARIATE_NORMAL:
    mvd_rep = std::make_shared<MultivariateNormalDistribution>(); break;
  //case JOINT_KDE:
  //  mvd_rep = new JointKDEDistribution();         break;
  //case GAUSSIAN_COPULA:
  //  mvd_rep = new CopulaDistribution<Gaussian>(); break; // if templated...
  //etc.
  default:
    PCerr << "Error: MultivariateDistribution type " << mv_dist_type
	  << " not available." << std::endl;
  }

  // some derived classes (especially template classes) cover multiple
  // ranVarTypes, so override ctor assignments for those cases:
  if (mvd_rep)
    mvd_rep->mvDistType = mv_dist_type;

  return mvd_rep;
}


/** Copy constructor manages sharing of mvDistRep. */
MultivariateDistribution::
MultivariateDistribution(const MultivariateDistribution& mv_dist):
  mvDistRep(mv_dist.mvDistRep)
{ /* empty rep */ }


MultivariateDistribution MultivariateDistribution::
operator=(const MultivariateDistribution& mv_dist)
{
  mvDistRep = mv_dist.mvDistRep;
  return *this; // calls copy constructor since returned by value
}


MultivariateDistribution::~MultivariateDistribution()
{ /* empty dtor */ }


const RandomVariable& MultivariateDistribution::random_variable(size_t i) const
{
  if (mvDistRep)
    return mvDistRep->random_variable(i);
  else // default definition
    return random_variables()[i];
}


RandomVariable& MultivariateDistribution::random_variable(size_t i)
{
  if (mvDistRep)
    return mvDistRep->random_variable(i);
  else // default definition
    return random_variables()[i];
}


const std::vector<RandomVariable>& MultivariateDistribution::
random_variables() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: random_variables() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->random_variables();
}


std::vector<RandomVariable>& MultivariateDistribution::
random_variables()
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: random_variables() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->random_variables();
}


const ShortArray& MultivariateDistribution::random_variable_types() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: random_variable_types() not supported for this "
	  << "multivariate distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->random_variable_types();
}


void MultivariateDistribution::random_variable_types(const ShortArray& rv_types)
{
  if (mvDistRep)
    mvDistRep->random_variable_types(rv_types);
  else { // forward to letter
    PCerr << "Error: random_variable_types(ShortArray) not supported for this "
	  << "multivariate distribution type." << std::endl;
    abort_handler(-1);
  }
}


short MultivariateDistribution::random_variable_type(size_t i) const
{
  if (mvDistRep)
    return mvDistRep->random_variable_type(i);
  else // default definition
    return random_variable_types()[i];
}


void MultivariateDistribution::random_variable_type(short rv_type, size_t i)
{
  if (mvDistRep)
    return mvDistRep->random_variable_type(rv_type, i);
  else { // forward to letter
    PCerr << "Error: random_variable_type(short, size_t) not supported for "
	  << "this multivariate distribution type." << std::endl;
    abort_handler(-1);
  }
}


size_t MultivariateDistribution::active_variable_index(size_t i) const
{
  if (mvDistRep)
    return mvDistRep->active_variable_index(i);
  else { // forward to letter
    const BitArray& active_vars = active_variables();
    if (active_vars.empty()) // no subset, all variables are active
      return i;

    /* This approach may need to store an index mapping for fast lookup
    size_t v, index = _NPOS, num_v = active_vars.size(), count = 0, id = i+1;
    for (v=0; v<num_v; ++v) {
      if (active_vars[v]) ++count;
      if (count == id) { index = v; break; }
    }
    */
    // Should be more efficient (presuming based on bit shifting)
    size_t index = active_vars.find_first(), cntr = 0;
    while (cntr < i)
      { index = active_vars.find_next(index); ++cntr; }
    return index;
  }
}


const RandomVariable& MultivariateDistribution::
active_random_variable(size_t i) const
{
  if (mvDistRep)
    return mvDistRep->active_random_variable(i);
  else // not virtual
    return random_variable(active_variable_index(i));
}


RandomVariable& MultivariateDistribution::active_random_variable(size_t i)
{
  if (mvDistRep)
    return mvDistRep->active_random_variable(i);
  else // not virtual
    return random_variable(active_variable_index(i));
}


short MultivariateDistribution::active_random_variable_type(size_t i) const
{
  if (mvDistRep)
    return mvDistRep->active_random_variable_type(i);
  else // not virtual
    return random_variable_type(active_variable_index(i));
}


const BitArray& MultivariateDistribution::active_variables() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: active_variables() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->active_variables();
}


const BitArray& MultivariateDistribution::active_correlations() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: active_correlations() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->active_correlations();
}


const RealSymMatrix& MultivariateDistribution::correlation_matrix() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: correlation_matrix() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->correlation_matrix();
}


void MultivariateDistribution::correlation_matrix(const RealSymMatrix& corr)
{
  if (mvDistRep)
    mvDistRep->correlation_matrix(corr);
  else { // forward to letter
    PCerr << "Error: correlation_matrix(RealSymMatrix) not supported for this "
	  << "multivariate distribution type." << std::endl;
    abort_handler(-1);
  }
}


void MultivariateDistribution::
pull_distribution_parameters(const MultivariateDistribution& mv_dist)
{
  if (mvDistRep)
    mvDistRep->pull_distribution_parameters(mv_dist);
  else { // forward to letter
    PCerr << "Error: pull_distribution_parameters(MultivariateDistribution) "
	  << "not supported for this multivariate distribution type."
	  << std::endl;
    abort_handler(-1);
  }
}


void MultivariateDistribution::
pull_distribution_parameters(const MultivariateDistribution& mv_dist,
			     const StringArray& pull_labels,
			     const StringArray& push_labels)
{
  if (mvDistRep)
    mvDistRep->pull_distribution_parameters(mv_dist, pull_labels, push_labels);
  else { // forward to letter
    PCerr << "Error: pull_distribution_parameters(MultivariateDistribution, "
	  << "StringArray, StringArray)\n       not supported for this "
	  <<  "multivariate distribution type." << std::endl;
    abort_handler(-1);
  }
}


RealRealPairArray MultivariateDistribution::moments() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: moments() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->moments();
}


RealVector MultivariateDistribution::means() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: means() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->means();
}


RealVector MultivariateDistribution::variances() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: variances() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->variances();
}


RealVector MultivariateDistribution::std_deviations() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: std_deviations() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->std_deviations();
}


RealRealPairArray MultivariateDistribution::distribution_bounds() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: distribution_bounds() not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->distribution_bounds();
}


RealVector MultivariateDistribution::distribution_lower_bounds() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: distribution_lower_bounds() not supported for this "
	  << "multivariate distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->distribution_lower_bounds();
}


RealVector MultivariateDistribution::distribution_upper_bounds() const
{
  if (!mvDistRep) { // forward to letter
    PCerr << "Error: distribution_upper_bounds() not supported for this "
	  << "multivariate distribution type." << std::endl;
    abort_handler(-1);
  }

  return mvDistRep->distribution_upper_bounds();
}


bool MultivariateDistribution::global_bounds() const
{ return (mvDistRep) ? mvDistRep->global_bounds() : false; }


void MultivariateDistribution::
lower_bounds(const RealVector& l_bnds, const BitArray& mask)
{
  if (mvDistRep) // forward to letter
    mvDistRep->lower_bounds(l_bnds, mask);
  else {
    PCerr << "Error: lower_bounds(RealVector, BitArray)\n       not "
	  << "supported for this multivariate distribution type." << std::endl;
    abort_handler(-1);
  }
}


void MultivariateDistribution::
lower_bounds(const IntVector& l_bnds, const BitArray& mask)
{
  if (mvDistRep) // forward to letter
    mvDistRep->lower_bounds(l_bnds, mask);
  else {
    PCerr << "Error: lower_bounds(IntVector, BitArray)\n       not "
	  << "supported for this multivariate distribution type." << std::endl;
    abort_handler(-1);
  }
}


void MultivariateDistribution::lower_bound(Real l_bnd, size_t rv_index)
{
  if (mvDistRep) // forward to letter
    mvDistRep->lower_bound(l_bnd, rv_index);
  else {
    PCerr << "Error: lower_bound(Real, size_t)\n       not supported for "
	  << "this multivariate distribution type." << std::endl;
    abort_handler(-1);
  }
}


void MultivariateDistribution::lower_bound(int l_bnd, size_t rv_index)
{
  if (mvDistRep) // forward to letter
    mvDistRep->lower_bound(l_bnd, rv_index);
  else {
    PCerr << "Error: lower_bound(int, size_t)\n       not supported for "
	  << "this multivariate distribution type." << std::endl;
    abort_handler(-1);
  }
}


void MultivariateDistribution::
upper_bounds(const RealVector& u_bnds, const BitArray& mask)
{
  if (mvDistRep) // forward to letter
    mvDistRep->upper_bounds(u_bnds, mask);
  else {
    PCerr << "Error: upper_bounds(RealVector, BitArray)\n       not "
	  << "supported for this multivariate distribution type." << std::endl;
    abort_handler(-1);
  }
}


void MultivariateDistribution::
upper_bounds(const IntVector& u_bnds, const BitArray& mask)
{
  if (mvDistRep) // forward to letter
    mvDistRep->upper_bounds(u_bnds, mask);
  else {
    PCerr << "Error: upper_bounds(IntVector, BitArray)\n       not "
	  << "supported for this multivariate distribution type." << std::endl;
    abort_handler(-1);
  }
}


void MultivariateDistribution::upper_bound(Real u_bnd, size_t rv_index)
{
  if (mvDistRep) // forward to letter
    mvDistRep->upper_bound(u_bnd, rv_index);
  else {
    PCerr << "Error: upper_bound(Real, size_t)\n       not supported for "
	  << "this multivariate distribution type." << std::endl;
    abort_handler(-1);
  }
}


void MultivariateDistribution::upper_bound(int u_bnd, size_t rv_index)
{
  if (mvDistRep) // forward to letter
    mvDistRep->upper_bound(u_bnd, rv_index);
  else {
    PCerr << "Error: upper_bound(int, size_t)\n       not supported for "
	  << "this multivariate distribution type." << std::endl;
    abort_handler(-1);
  }
}

// Vector input / Joint:

Real MultivariateDistribution::pdf(const RealVector& pt) const
{
  if (mvDistRep)
    return mvDistRep->pdf(pt);
  else { // forward to letter
    PCerr << "Error: pdf(RealVector) not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
    return 0.;
  }
}


Real MultivariateDistribution::log_pdf(const RealVector& pt) const
{
  if (mvDistRep)
    return mvDistRep->log_pdf(pt);
  else // default implementation (exponential-based distribs will override)
    return std::log(pdf(pt));
}

// Scalar input / Marginals:

Real MultivariateDistribution::pdf(Real val, size_t i) const
{
  if (mvDistRep)
    return mvDistRep->pdf(val, i);
  else { // default implementation (exponential-based distribs will override)
    PCerr << "Error: pdf(Real, size_t) not supported for this multivariate "
	  << "distribution type." << std::endl;
    abort_handler(-1);
    return 0.;
  }
}


Real MultivariateDistribution::pdf_gradient(Real val, size_t i) const
{
  if (mvDistRep)
    return mvDistRep->pdf_gradient(val, i);
  else { // default implementation (exponential-based distribs will override)
    PCerr << "Error: pdf_gradient(Real, size_t) not supported for this "
	  << "multivariate distribution type." << std::endl;
    abort_handler(-1);
    return 0.;
  }
}


Real MultivariateDistribution::pdf_hessian(Real val, size_t i) const
{
  if (mvDistRep)
    return mvDistRep->pdf_hessian(val, i);
  else { // default implementation (exponential-based distribs will override)
    PCerr << "Error: pdf_hessian(Real, size_t) not supported for this "
	  << "multivariate distribution type." << std::endl;
    abort_handler(-1);
    return 0.;
  }
}


Real MultivariateDistribution::log_pdf(Real val, size_t i) const
{
  if (mvDistRep)
    return mvDistRep->log_pdf(val, i);
  else // default implementation (exponential-based distribs will override)
    return std::log(pdf(val, i));
}


Real MultivariateDistribution::log_pdf_gradient(Real val, size_t i) const
{
  if (mvDistRep)
    return mvDistRep->log_pdf_gradient(val, i);
  else // default implementation (exponential-based distribs will override)
    return std::log(pdf_gradient(val, i));
}


Real MultivariateDistribution::log_pdf_hessian(Real val, size_t i) const
{
  if (mvDistRep)
    return mvDistRep->log_pdf_hessian(val, i);
  else // default implementation (exponential-based distribs will override)
    return std::log(pdf_hessian(val, i));
}
   

/** This function provides a deep copy, creating a envelope and copying
    data from the current letter (if defined) into the new envelope. */
MultivariateDistribution MultivariateDistribution::copy() const
{
  MultivariateDistribution mvd; // target

  if (mvDistRep) { // source has letter
    // allocate a responseRep letter, copy data attributes, share the srd
    mvd.mvDistRep = get_distribution(mvDistRep->mvDistType);
    // allow derived classes to specialize copy_rep if they augment
    // the base class data
    mvd.mvDistRep->copy_rep(mvDistRep);
  }

  return mvd;
}


/** Default overridden by derived classes */
void MultivariateDistribution::
copy_rep(std::shared_ptr<MultivariateDistribution> source_rep)
{ correlationFlag = source_rep->correlationFlag; }

} // namespace Pecos
