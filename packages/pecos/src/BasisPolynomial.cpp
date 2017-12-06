/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        BasisPolynomial
//- Description:  Class implementation of base class for basis polynomials
//-               
//- Owner:        Mike Eldred, Sandia National Laboratories

#include "BasisPolynomial.hpp"
#include "HermiteOrthogPolynomial.hpp"
#include "LegendreOrthogPolynomial.hpp"
#include "JacobiOrthogPolynomial.hpp"
#include "LaguerreOrthogPolynomial.hpp"
#include "GenLaguerreOrthogPolynomial.hpp"
#include "ChebyshevOrthogPolynomial.hpp"
#include "NumericGenOrthogPolynomial.hpp"
#include "LagrangeInterpPolynomial.hpp"
#include "HermiteInterpPolynomial.hpp"
#include "PiecewiseInterpPolynomial.hpp"
#include "KrawtchoukOrthogPolynomial.hpp"
#include "MeixnerOrthogPolynomial.hpp"
#include "CharlierOrthogPolynomial.hpp"
#include "HahnOrthogPolynomial.hpp"


namespace Pecos {

/** This constructor is the one which must build the base class data
    for all derived classes.  get_polynomial() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_polynomial() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~BasisPolynomial). */
BasisPolynomial::BasisPolynomial(BaseConstructor):// basisPolyType(-1),
  parametricUpdate(false), wtFactor(1.), ptFactor(1.), polyRep(NULL),
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "BasisPolynomial::BasisPolynomial(BaseConstructor) called "
	<< "to build base class for letter." << std::endl;
#endif
}


/** The default constructor is used in Array<BasisPolynomial>
    instantiations and by the alternate envelope constructor.  polyRep
    is NULL in this case (problem_db is needed to build a meaningful
    instance).  This makes it necessary to check for NULL in the copy
    constructor, assignment operator, and destructor. */
BasisPolynomial::BasisPolynomial(): polyRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "BasisPolynomial::BasisPolynomial() called to build empty "
	<< "basis polynomial object." << std::endl;
#endif
}


/** Envelope constructor which does not require access to problem_db.
    This constructor executes get_polynomial(type), which invokes the
    default constructor of the derived letter class, which in turn
    invokes the BaseConstructor of the base class. */
BasisPolynomial::BasisPolynomial(short poly_type, short rule):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "BasisPolynomial::BasisPolynomial(short) called to "
	<< "instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  polyRep = get_polynomial(poly_type, rule);
  if (poly_type && !polyRep) // bad type, insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize polyRep to the 
    appropriate derived type. */
BasisPolynomial* BasisPolynomial::get_polynomial(short poly_type, short rule)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_polynomial(short, short)."
	<< std::endl;
#endif

  BasisPolynomial* polynomial;
  // In orthogonal polynomial and global interpolation polynomial cases,
  // basisPolyType is not available at construct time, but is thereafter.
  switch (poly_type) {
  case NO_POLY:
    polynomial = NULL;                                                    break;
  case HERMITE_ORTHOG:  // var_type == "normal"
    polynomial = (rule) ? new HermiteOrthogPolynomial(rule)
                        : new HermiteOrthogPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case LEGENDRE_ORTHOG: // var_type == "uniform"
    polynomial = (rule) ? new LegendreOrthogPolynomial(rule)
                        : new LegendreOrthogPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case LAGUERRE_ORTHOG: // var_type == "exponential"
    polynomial = new LaguerreOrthogPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case JACOBI_ORTHOG:   // var_type == "beta"
    polynomial = new JacobiOrthogPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case GEN_LAGUERRE_ORTHOG: // var_type == "gamma"
    polynomial = new GenLaguerreOrthogPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case CHEBYSHEV_ORTHOG: // for Clenshaw-Curtis and Fejer
    polynomial = (rule) ? new ChebyshevOrthogPolynomial(rule)
                        : new ChebyshevOrthogPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case NUM_GEN_ORTHOG:
    polynomial = new NumericGenOrthogPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case LAGRANGE_INTERP:
    polynomial = new LagrangeInterpPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case HERMITE_INTERP:
    polynomial = (rule) ? new HermiteInterpPolynomial(rule)
                        : new HermiteInterpPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;              break;
  // PIECEWISE options include poly order, point type, and point data order:
  // LINEAR/QUADRATIC/CUBIC covers poly order, rule covers EQUIDISTANT/GENERAL
  // point type, and data order is inferred from poly order (grads for CUBIC).
  case PIECEWISE_LINEAR_INTERP: case PIECEWISE_QUADRATIC_INTERP:
  case PIECEWISE_CUBIC_INTERP:
    polynomial = (rule) ? new PiecewiseInterpPolynomial(poly_type, rule)
                        : new PiecewiseInterpPolynomial(poly_type);
    break;
  // Some Discrete orthogonal polynomials
  case KRAWTCHOUK_DISCRETE:   // var_type == "binomial"
    polynomial = new KrawtchoukOrthogPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case MEIXNER_DISCRETE:   // var_type == "negative binomial"
    polynomial = new MeixnerOrthogPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case CHARLIER_DISCRETE:   // var_type == "poisson"
    polynomial = new CharlierOrthogPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case HAHN_DISCRETE:   // var_type == "hygergeometric"
    polynomial = new HahnOrthogPolynomial();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  default:
    PCerr << "Error: BasisPolynomial type " << poly_type << " not available."
	 << std::endl;
    polynomial = NULL;                                                    break;
  }
  return polynomial;
}


/** Copy constructor manages sharing of polyRep and incrementing of
    referenceCount. */
BasisPolynomial::BasisPolynomial(const BasisPolynomial& polynomial)
{
  // Increment new (no old to decrement)
  polyRep = polynomial.polyRep;
  if (polyRep) // Check for an assignment of NULL
    polyRep->referenceCount++;

#ifdef REFCOUNT_DEBUG
  PCout << "BasisPolynomial::BasisPolynomial(BasisPolynomial&)" << std::endl;
  if (polyRep)
    PCout << "polyRep referenceCount = " << polyRep->referenceCount <<std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old polyRep,
    assigns new polyRep, and increments referenceCount for new polyRep. */
BasisPolynomial BasisPolynomial::operator=(const BasisPolynomial& polynomial)
{
  if (polyRep != polynomial.polyRep) { // normal case: old != new
    // Decrement old
    if (polyRep) // Check for NULL
      if ( --polyRep->referenceCount == 0 ) 
	delete polyRep;
    // Assign and increment new
    polyRep = polynomial.polyRep;
    if (polyRep) // Check for NULL
      polyRep->referenceCount++;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  PCout << "BasisPolynomial::operator=(BasisPolynomial&)" << std::endl;
  if (polyRep)
    PCout << "polyRep referenceCount = " << polyRep->referenceCount
	  << std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes polyRep when
    referenceCount reaches zero. */
BasisPolynomial::~BasisPolynomial()
{ 
  // Check for NULL pointer 
  if (polyRep) {
    --polyRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    PCout << "polyRep referenceCount decremented to " << polyRep->referenceCount
	  << std::endl;
#endif
    if (polyRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      PCout << "deleting polyRep" << std::endl;
#endif
      delete polyRep;
    }
  }
}


Real BasisPolynomial::type1_value(unsigned short n)
{
  if (!polyRep) {
    PCerr << "Error: type1_value(unsigned short) not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->type1_value(n);
}


Real BasisPolynomial::type1_value(Real x, unsigned short n)
{
  if (!polyRep) {
    PCerr << "Error: type1_value(Real, unsigned short) not available for this "
	  << "basis polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->type1_value(x, n);
}


Real BasisPolynomial::type2_value(Real x, unsigned short n)
{
  if (!polyRep) {
    PCerr << "Error: type2_value(Real, unsigned short) not available for this "
	  << "basis polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->type2_value(x, n);
}


Real BasisPolynomial::type1_gradient(unsigned short n)
{
  if (!polyRep) {
    PCerr << "Error: type1_gradient(unsigned short) not available for this "
	  << "basis polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->type1_gradient(n);
}


Real BasisPolynomial::type1_gradient(Real x, unsigned short n)
{
  if (!polyRep) {
    PCerr << "Error: type1_gradient(Real, unsigned short) not available for "
	  << "this basis polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->type1_gradient(x, n);
}


Real BasisPolynomial::type2_gradient(Real x, unsigned short n)
{
  if (!polyRep) {
    PCerr << "Error: type2_gradient(Real, unsigned short) not available for "
	  << "this basis polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->type2_gradient(x, n);
}


Real BasisPolynomial::type1_hessian(Real x, unsigned short n)
{
  if (!polyRep) {
    PCerr << "Error: type1_hessian(Real, unsigned short) not available for "
	  << "this basis polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->type1_hessian(x, n);
}


Real BasisPolynomial::norm_squared(unsigned short n)
{
  if (!polyRep) {
    PCerr << "Error: norm_squared(unsigned short) not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->norm_squared(n);
}


const RealArray& BasisPolynomial::collocation_points(unsigned short n)
{
  if (!polyRep) {
    PCerr << "Error: collocation_points() not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->collocation_points(n);
}


const RealArray& BasisPolynomial::type1_collocation_weights(unsigned short n)
{
  if (!polyRep) {
    PCerr << "Error: type1_collocation_weights() not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->type1_collocation_weights(n);
}


const RealArray& BasisPolynomial::type2_collocation_weights(unsigned short n)
{
  if (!polyRep) {
    PCerr << "Error: type2_collocation_weights() not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->type2_collocation_weights(n);
}


void BasisPolynomial::set_new_point(Real x, short order)
{
  if (polyRep)
    polyRep->set_new_point(x, order);
  else {
    PCerr << "Error: set_new_point(Real, short) not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
}


void BasisPolynomial::
set_new_point(Real x, short order, const UShortArray& delta_key)
{
  if (polyRep)
    polyRep->set_new_point(x, order, delta_key);
  else {
    PCerr << "Error: set_new_point(Real, short, UShortArray) not available for "
	  << "this basis polynomial type." << std::endl;
    abort_handler(-1);
  }
}


size_t BasisPolynomial::exact_index() const
{
  if (!polyRep) {
    PCerr << "Error: exact_index() not available for this basis polynomial "
	  << "type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->exact_index();
}


size_t BasisPolynomial::exact_delta_index() const
{
  if (!polyRep) {
    PCerr << "Error: exact_delta_index() not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->exact_delta_index();
}


const RealVector& BasisPolynomial::barycentric_value_factors() const
{
  if (!polyRep) {
    PCerr << "Error: barycentric_value_factors() not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->barycentric_value_factors();
}


Real BasisPolynomial::barycentric_value_factor(unsigned short i) const
{
  if (!polyRep) {
    PCerr << "Error: barycentric_value_factor() not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->barycentric_value_factor(i);
}


const RealVector& BasisPolynomial::barycentric_gradient_factors() const
{
  if (!polyRep) {
    PCerr << "Error: barycentric_gradient_factors() not available for this "
	  << "basis polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->barycentric_gradient_factors();
}


Real BasisPolynomial::barycentric_gradient_factor(unsigned short i) const
{
  if (!polyRep) {
    PCerr << "Error: barycentric_gradient_factor() not available for this "
	  << "basis polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->barycentric_gradient_factor(i);
}


Real BasisPolynomial::barycentric_value_factor_sum() const
{
  if (!polyRep) {
    PCerr << "Error: barycentric_value_factor_sum() not available for this "
	  << "basis polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->barycentric_value_factor_sum();
}


Real BasisPolynomial::barycentric_difference_product() const
{
  if (!polyRep) {
    PCerr << "Error: barycentric_difference_product() not available for this "
	  << "basis polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->barycentric_difference_product();
}


void BasisPolynomial::reset_gauss()
{
  if (polyRep)
    polyRep->reset_gauss();
  else {
    PCerr << "Error: reset_gauss() not available for this basis polynomial "
	  << "type." << std::endl;
    abort_handler(-1);
  }
}


Real BasisPolynomial::point_factor()
{
  if (polyRep)
    return polyRep->point_factor();
  else // default is used whenever ptFactor does not need to be updated
    return ptFactor;
}


Real BasisPolynomial::weight_factor()
{
  if (polyRep)
    return polyRep->weight_factor();
  else // default is used whenever wtFactor does not need to be updated
    return wtFactor;
}


Real BasisPolynomial::alpha_polynomial() const
{
  if (!polyRep) {
    PCerr << "Error: alpha_polynomial() not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->alpha_polynomial();
}


Real BasisPolynomial::beta_polynomial() const
{
  if (!polyRep) {
    PCerr << "Error: beta_polynomial() not available for this basis polynomial "
	  << "type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->beta_polynomial();
}


void BasisPolynomial::alpha_stat(Real alpha)
{
  if (polyRep)
    polyRep->alpha_stat(alpha);
  else {
    PCerr << "Error: alpha_stat() not available for this basis polynomial type."
	  << std::endl;
    abort_handler(-1);
  }
}


void BasisPolynomial::beta_stat(Real beta)
{
  if (polyRep)
    polyRep->beta_stat(beta);
  else {
    PCerr << "Error: beta_stat() not available for this basis polynomial type."
	  << std::endl;
    abort_handler(-1);
  }
}


void BasisPolynomial::collocation_rule(short rule)
{
  if (polyRep)
    polyRep->collocation_rule(rule);
  else {
    PCerr << "Error: collocation_rule(short) not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
}


short BasisPolynomial::collocation_rule() const
{
  if (!polyRep) {
    PCerr << "Error: collocation_rule() not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->collocation_rule();
}


size_t BasisPolynomial::interpolation_size() const
{
  if (!polyRep) {
    PCerr << "Error: interpolation_size() not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->interpolation_size();
}


void BasisPolynomial::interpolation_points(const RealArray& interpolation_pts)
{
  if (polyRep)
    polyRep->interpolation_points(interpolation_pts);
  else {
    PCerr << "Error: interpolation_points() not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
}


const RealArray& BasisPolynomial::interpolation_points() const
{
  if (!polyRep) {
    PCerr << "Error: interpolation_points() not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->interpolation_points();
}


bool BasisPolynomial::parameterized() const
{
  if (polyRep)
    return polyRep->parameterized();
  else
    return false; // default if not overridden
}


Real BasisPolynomial::length_scale() const
{
  if (!polyRep) {
    PCerr << "Error: length_scale() not available for this basis polynomial "
	  << "type." << std::endl;
    abort_handler(-1);
  }
  return polyRep->length_scale();
}


void BasisPolynomial::precompute_rules(unsigned short order)
{
  if (polyRep)
    polyRep->precompute_rules(order); // else no-op
  //else {
  //  PCerr << "Error: precompute_rules() not available for this basis "
  //	    << "polynomial type." << std::endl;
  //  abort_handler(-1);
  //}
}

} // namespace Pecos
