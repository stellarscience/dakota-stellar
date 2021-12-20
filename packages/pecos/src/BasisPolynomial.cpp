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
    letter IS the representation, its rep pointer is set to NULL. */
BasisPolynomial::BasisPolynomial(BaseConstructor):// basisPolyType(-1),
  wtFactor(1.), ptFactor(1.)
{ /* empty ctor */ }


/** The default constructor is used in Array<BasisPolynomial>
    instantiations and by the alternate envelope constructor.  polyRep
    is NULL in this case (problem_db is needed to build a meaningful
    instance). */
BasisPolynomial::BasisPolynomial()
{ /* empty ctor */ }


/** Envelope constructor which does not require access to problem_db.
    This constructor executes get_polynomial(type), which invokes the
    default constructor of the derived letter class, which in turn
    invokes the BaseConstructor of the base class. */
BasisPolynomial::BasisPolynomial(short poly_type, short rule):
  // Set the rep pointer to the appropriate derived type
  polyRep(get_polynomial(poly_type, rule))
{
  if (poly_type && !polyRep) // bad type, insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize polyRep to the 
    appropriate derived type. */
std::shared_ptr<BasisPolynomial>
BasisPolynomial::get_polynomial(short poly_type, short rule)
{
  std::shared_ptr<BasisPolynomial> polynomial;
  // In orthogonal polynomial and global interpolation polynomial cases,
  // basisPolyType is not available at construct time, but is thereafter.
  switch (poly_type) {
  case NO_POLY:
    break;
  case HERMITE_ORTHOG:  // var_type == "normal"
    polynomial = (rule) ? std::make_shared<HermiteOrthogPolynomial>(rule) :
      std::make_shared<HermiteOrthogPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case LEGENDRE_ORTHOG: // var_type == "uniform"
    polynomial = (rule) ? std::make_shared<LegendreOrthogPolynomial>(rule) :
      std::make_shared<LegendreOrthogPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case LAGUERRE_ORTHOG: // var_type == "exponential"
    polynomial = std::make_shared<LaguerreOrthogPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case JACOBI_ORTHOG:   // var_type == "beta"
    polynomial = std::make_shared<JacobiOrthogPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case GEN_LAGUERRE_ORTHOG: // var_type == "gamma"
    polynomial = std::make_shared<GenLaguerreOrthogPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case CHEBYSHEV_ORTHOG: // for Clenshaw-Curtis and Fejer
    polynomial = (rule) ? std::make_shared<ChebyshevOrthogPolynomial>(rule) :
      std::make_shared<ChebyshevOrthogPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case NUM_GEN_ORTHOG:
    polynomial = std::make_shared<NumericGenOrthogPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case LAGRANGE_INTERP:
    polynomial = std::make_shared<LagrangeInterpPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case HERMITE_INTERP:
    polynomial = (rule) ? std::make_shared<HermiteInterpPolynomial>(rule) :
      std::make_shared<HermiteInterpPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  // PIECEWISE options include poly order, point type, and point data order:
  // LINEAR/QUADRATIC/CUBIC covers poly order, rule covers EQUIDISTANT/GENERAL
  // point type, and data order is inferred from poly order (grads for CUBIC).
  case PIECEWISE_LINEAR_INTERP: case PIECEWISE_QUADRATIC_INTERP:
  case PIECEWISE_CUBIC_INTERP:
    polynomial = (rule) ?
      std::make_shared<PiecewiseInterpPolynomial>(poly_type, rule) :
      std::make_shared<PiecewiseInterpPolynomial>(poly_type);
    break;
  // Some Discrete orthogonal polynomials
  case KRAWTCHOUK_DISCRETE:   // var_type == "binomial"
    polynomial = std::make_shared<KrawtchoukOrthogPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case MEIXNER_DISCRETE:   // var_type == "negative binomial"
    polynomial = std::make_shared<MeixnerOrthogPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case CHARLIER_DISCRETE:   // var_type == "poisson"
    polynomial = std::make_shared<CharlierOrthogPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  case HAHN_DISCRETE:   // var_type == "hygergeometric"
    polynomial = std::make_shared<HahnOrthogPolynomial>();
    if (polynomial) polynomial->basisPolyType = poly_type;                break;
  default:
    PCerr << "Error: BasisPolynomial type " << poly_type << " not available."
	 << std::endl;
    break;
  }
  return polynomial;
}


/** Copy constructor manages sharing of polyRep. */
BasisPolynomial::BasisPolynomial(const BasisPolynomial& polynomial):
  polyRep(polynomial.polyRep)
{ /* empty ctor */ }


BasisPolynomial BasisPolynomial::operator=(const BasisPolynomial& polynomial)
{
  polyRep = polynomial.polyRep;
  return *this; // calls copy constructor since returned by value
}


BasisPolynomial::~BasisPolynomial()
{ /* empty dtor */ }


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


bool BasisPolynomial::parameter_update() const
{
  if (polyRep)
    return polyRep->parameter_update();
  else // default for non-parametric polynomials (e.g., PiecewiseInterp)
    return false;
}


bool BasisPolynomial::points_defined(unsigned short order) const
{
  if (polyRep)
    return polyRep->points_defined(order);
  else // default
    return false;
}


bool BasisPolynomial::type1_weights_defined(unsigned short order) const
{
  if (polyRep)
    return polyRep->type1_weights_defined(order);
  else // default
    return false;
}


bool BasisPolynomial::type2_weights_defined(unsigned short order) const
{
  if (polyRep)
    return polyRep->type2_weights_defined(order);
  else // default
    return false;
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


void BasisPolynomial::pull_parameter(short dist_param, Real& param) const
{
  if (polyRep)
    polyRep->pull_parameter(dist_param, param);
  else {
    PCerr << "Error: pull_parameter(Real) not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
}


void BasisPolynomial::
pull_parameter(short dist_param, unsigned int& param) const
{
  if (polyRep)
    polyRep->pull_parameter(dist_param, param);
  else {
    PCerr << "Error: pull_parameter(unsigned int) not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
}


void BasisPolynomial::push_parameter(short dist_param, Real param)
{
  if (polyRep)
    polyRep->push_parameter(dist_param, param);
  else {
    PCerr << "Error: push_parameter(Real) not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
}


void BasisPolynomial::push_parameter(short dist_param, unsigned int param)
{
  if (polyRep)
    polyRep->push_parameter(dist_param, param);
  else {
    PCerr << "Error: push_parameter(unsigned int) not available for this basis "
	  << "polynomial type." << std::endl;
    abort_handler(-1);
  }
}


bool BasisPolynomial::parameterized() const
{
  if (polyRep)
    return polyRep->parameterized();
  else
    return false; // default if not overridden
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
