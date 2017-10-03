/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 IntegrationDriver
//- Description: Implementation code for IntegrationDriver class
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#include "IntegrationDriver.hpp"
#include "CubatureDriver.hpp"
#include "TensorProductDriver.hpp"
#include "LightweightSparseGridDriver.hpp"
#include "CombinedSparseGridDriver.hpp"
#include "HierarchSparseGridDriver.hpp"
#include "SharedInterpPolyApproxData.hpp"
#include "SharedOrthogPolyApproxData.hpp"

static const char rcsId[]="@(#) $Id: IntegrationDriver.C,v 1.57 2004/06/21 19:57:32 mseldre Exp $";

//#define DEBUG

namespace Pecos {

UShortArray IntegrationDriver::orderGenzKeister;
UShortArray IntegrationDriver::precGenzKeister;


/** This constructor is the one which must build the base class data
    for all derived classes.  get_driver() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_driver() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~IntegrationDriver). */
IntegrationDriver::IntegrationDriver(BaseConstructor):
  driverMode(DEFAULT_MODE), computeType2Weights(false),
  activeReinterpIndex(0), driverRep(NULL), referenceCount(1)
{
  /* Standard 5 step sequence is fully nested (1+2+6+10+16 = 1, 3, 9, 19, 35)
  if (orderGenzKeister.empty()) {
    orderGenzKeister.resize(5); //orderGenzKeister = { 1, 3, 9, 19, 35 };
    orderGenzKeister[0] =  1; orderGenzKeister[1] =  3; orderGenzKeister[2] = 9;
    orderGenzKeister[3] = 19; orderGenzKeister[4] = 35;
  }
  if (precGenzKeister.empty()) {
    precGenzKeister.resize(5); //precGenzKeister = { 1, 5, 15, 29, 51 }; 
    precGenzKeister[0] =  1; precGenzKeister[1] =  5; precGenzKeister[2] = 15;
    precGenzKeister[3] = 29; precGenzKeister[4] = 51;
  }
  */

  // To maximize the available precision, we augment this fully nested 5 step
  // sequence with a 6th step that reuses the 1+2+6+10 sequence portion, but
  // replaces the +16 with +24.
  if (orderGenzKeister.empty()) {
    orderGenzKeister.resize(6); //orderGenzKeister = { 1, 3, 9, 19, 35, 43 };
    orderGenzKeister[0] =  1; orderGenzKeister[1] =  3;
    orderGenzKeister[2] =  9; orderGenzKeister[3] = 19;
    orderGenzKeister[4] = 35; orderGenzKeister[5] = 43;
  }
  if (precGenzKeister.empty()) {
    precGenzKeister.resize(6); //precGenzKeister = { 1, 5, 15, 29, 51, 67 }; 
    precGenzKeister[0] =  1; precGenzKeister[1] =  5; precGenzKeister[2] = 15;
    precGenzKeister[3] = 29; precGenzKeister[4] = 51; precGenzKeister[5] = 67;
  }

#ifdef REFCOUNT_DEBUG
  PCout << "IntegrationDriver::IntegrationDriver(BaseConstructor) called to "
        << "build base class for letter." << std::endl;
#endif
}


/** The default constructor: driverRep is NULL in this case.  This
    makes it necessary to check for NULL in the copy constructor,
    assignment operator, and destructor. */
IntegrationDriver::IntegrationDriver(): driverRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "IntegrationDriver::IntegrationDriver() called to build empty "
        << "envelope." << std::endl;
#endif
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_driver, since IntegrationDriver(BaseConstructor)
    builds the actual base class data for the derived basis functions. */
IntegrationDriver::IntegrationDriver(short driver_type):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "IntegrationDriver::IntegrationDriver(short) called to "
        << "instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  driverRep = get_driver(driver_type);
  if ( !driverRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize driverRep to the 
    appropriate derived type. */
IntegrationDriver* IntegrationDriver::get_driver(short driver_type)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_driver(short)." << std::endl;
#endif

  switch (driver_type) {
  case QUADRATURE:              return new TensorProductDriver();         break;
  case CUBATURE:                return new CubatureDriver();              break;
  case LIGHTWEIGHT_SPARSE_GRID: return new LightweightSparseGridDriver(); break;
  case COMBINED_SPARSE_GRID:    return new CombinedSparseGridDriver();    break;
  case HIERARCHICAL_SPARSE_GRID: return new HierarchSparseGridDriver();   break;
  default:
    PCerr << "Error: IntegrationDriver type " << driver_type
	  << " not available." << std::endl;
    return NULL;                                           break;
  }
}


/** Copy constructor manages sharing of driverRep and incrementing
    of referenceCount. */
IntegrationDriver::IntegrationDriver(const IntegrationDriver& driver)
{
  // Increment new (no old to decrement)
  driverRep = driver.driverRep;
  if (driverRep) // Check for an assignment of NULL
    driverRep->referenceCount++;

#ifdef REFCOUNT_DEBUG
  PCout << "IntegrationDriver::IntegrationDriver(IntegrationDriver&)"
	<< std::endl;
  if (driverRep)
    PCout << "driverRep referenceCount = " << driverRep->referenceCount
	  << std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old driverRep,
    assigns new driverRep, and increments referenceCount for new
    driverRep. */
IntegrationDriver IntegrationDriver::operator=(const IntegrationDriver& driver)
{
  if (driverRep != driver.driverRep) { // std case: old != new
    // Decrement old
    if (driverRep) // Check for null pointer
      if (--driverRep->referenceCount == 0) 
	delete driverRep;
    // Assign and increment new
    driverRep = driver.driverRep;
    if (driverRep) // Check for an assignment of NULL
      driverRep->referenceCount++;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  PCout << "IntegrationDriver::operator=(IntegrationDriver&)" << std::endl;
  if (driverRep)
    PCout << "driverRep referenceCount = " << driverRep->referenceCount
	  << std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes driverRep
    when referenceCount reaches zero. */
IntegrationDriver::~IntegrationDriver()
{ 
  // Check for NULL pointer 
  if (driverRep) {
    --driverRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    PCout << "driverRep referenceCount decremented to "
	  << driverRep->referenceCount << std::endl;
#endif
    if (driverRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      PCout << "deleting driverRep" << std::endl;
#endif
      delete driverRep;
    }
  }
}


void IntegrationDriver::
assign_rep(IntegrationDriver* driver_rep, bool ref_count_incr)
{ 
  if (driverRep == driver_rep) {
    // if ref_count_incr = true (rep from another envelope), do nothing as
    // referenceCount should already be correct (see also operator= logic).
    // if ref_count_incr = false (rep from on the fly), then this is an error.
    if (!ref_count_incr) {
      PCerr << "Error: duplicated driver_rep pointer assignment without "
	    << "reference count increment in IntegrationDriver::assign_rep()."
	    << std::endl;
      abort_handler(-1);
    }
  }
  else { // normal case: old != new
    // Decrement old
    if (driverRep) // Check for NULL
      if ( --driverRep->referenceCount == 0 ) 
	delete driverRep;
    // Assign new
    driverRep = driver_rep;
    // Increment new
    if (driverRep && ref_count_incr) // Check for NULL & honor ref_count_incr
      driverRep->referenceCount++;
  }

#ifdef REFCOUNT_DEBUG
  PCout << "IntegrationDriver::assign_rep(IntegrationDriver*)" << std::endl;
  if (driverRep)
    PCout << "driverRep referenceCount = " << driverRep->referenceCount
	  << std::endl;
#endif
}


void IntegrationDriver::
initialize_grid_parameters(const ShortArray& u_types, 
			   const AleatoryDistParams& adp)
{
  if (driverRep)
    driverRep->initialize_grid_parameters(u_types, adp); // forward to letter
  else // default implementation
    SharedPolyApproxData::update_basis_distribution_parameters(u_types, adp,
							       polynomialBasis);
}


void IntegrationDriver::compute_grid(RealMatrix& var_sets)
{
  if (driverRep)
    driverRep->compute_grid(var_sets); // forward to letter
  else {
    PCerr << "Error: compute_grid(RealMatrix&) not available for this "
	  << "driver type." << std::endl;
    abort_handler(-1);
  }
}


int IntegrationDriver::grid_size()
{
  if (!driverRep) {
    PCerr << "Error: grid_size() not available for this driver type."
	  << std::endl;
    abort_handler(-1);
  }
  return driverRep->grid_size(); // forward to letter
}


void IntegrationDriver::
reinterpolated_tensor_grid(const UShortArray& lev_index,
			   const SizetList& reinterp_indices)
{
  if (driverRep) // forward to letter
    driverRep->reinterpolated_tensor_grid(lev_index, reinterp_indices);
  else {
    PCerr << "Error: reinterpolated_tensor_grid() not available for this "
	  << "driver type." << std::endl;
    abort_handler(-1);
  }
}


const RealVector& IntegrationDriver::type1_weight_sets() const
{
  if (!driverRep) {
    PCerr << "Error: type1_weight_sets() not available for this driver type."
	  << std::endl;
    abort_handler(-1);
  }
  return driverRep->type1_weight_sets();
}


const RealMatrix& IntegrationDriver::type2_weight_sets() const
{
  if (!driverRep) {
    PCerr << "Error: type2_weight_sets() not available for this driver type."
	  << std::endl;
    abort_handler(-1);
  }
  return driverRep->type2_weight_sets();
}


/** protected function called only from derived class letters. */
void IntegrationDriver::
initialize_grid(const ShortArray& u_types,
		const ExpansionConfigOptions& ec_options,
		const BasisConfigOptions& bc_options)
{
  numVars = u_types.size();
  ShortArray basis_types; bool dist_params;
  if (ec_options.expBasisType == NODAL_INTERPOLANT ||
      ec_options.expBasisType == HIERARCHICAL_INTERPOLANT) {
    driverMode = INTERPOLATION_MODE;
    dist_params = SharedInterpPolyApproxData::
      initialize_driver_types_rules(u_types, bc_options, basis_types,
				    collocRules);
  }
  else {
    driverMode = INTEGRATION_MODE;
    dist_params = SharedPolyApproxData::
      initialize_orthogonal_basis_types_rules(u_types, bc_options, basis_types,
					      collocRules);
  }

  SharedPolyApproxData::
    initialize_polynomial_basis(basis_types, collocRules, polynomialBasis);
  // TO DO: need AleatoryDistParams instance
  //if (dist_params)
  //  SharedPolyApproxData::
  //    update_basis_distribution_parameters(u_types, adp, polynomialBasis);

  for (size_t i=0; i<numVars; i++)
    if (basis_types[i] == HERMITE_INTERP ||
	basis_types[i] == PIECEWISE_CUBIC_INTERP)
      { computeType2Weights = true; break; }
}


/** protected function called only from derived class letters. */
void IntegrationDriver::
initialize_grid(const std::vector<BasisPolynomial>& poly_basis)
{
  if (driverRep)
    driverRep->initialize_grid(poly_basis);
  else {
    numVars         = poly_basis.size();
    polynomialBasis = poly_basis; // shallow copy

    // For setting driverMode, basis_type is insufficient since incoming basis
    // is the driver basis (not the interp poly basis).  collocRules is also
    // insufficient since we want to sync on #pts (not integrand precision) for
    // nested Gauss as well.  Currently, it is set separately via driver.mode().
    //driverMode = (ec_options.expBasisType == NODAL_INTERPOLANT ||
    //		    ec_options.expBasisType == HIERARCHICAL_INTERPOLANT)
    //           ? INTERPOLATION_MODE : INTEGRATION_MODE;

    collocRules.resize(numVars);
    for (size_t i=0; i<numVars; i++) {
      // update collocRules
      collocRules[i] = poly_basis[i].collocation_rule();
      // define computeType2Weights
      short basis_type = poly_basis[i].basis_type();
      if (basis_type == HERMITE_INTERP || basis_type == PIECEWISE_CUBIC_INTERP)
	computeType2Weights = true;
    }
  }
}


void IntegrationDriver::store_grid(size_t index)
{ } // default is no-op


void IntegrationDriver::restore_grid(size_t index)
{ } // default is no-op


void IntegrationDriver::remove_stored_grid(size_t index)
{ } // default is no-op


void IntegrationDriver::clear_stored()
{ } // default is no-op


size_t IntegrationDriver::maximal_grid() const
{
  if (!driverRep) {
    PCerr << "Error: maximal_grid() not available for this driver type."
	  << std::endl;
    abort_handler(-1);
  }
  return driverRep->maximal_grid();
}


void IntegrationDriver::swap_grid(size_t index)
{
  if (driverRep)
    driverRep->swap_grid(index);
  else {
    PCerr << "Error: swap_grid() not available for this driver type."
	  << std::endl;
    abort_handler(-1);
  }
}


void IntegrationDriver::
compute_tensor_grid(const UShortArray& quad_order, const UShortArray& lev_index,
		    RealMatrix& variable_sets,  RealVector& t1_weight_sets,
		    RealMatrix& t2_weight_sets, UShort2DArray& colloc_key)
{
  size_t i, j, k, num_colloc_pts = 1;
  for (i=0; i<numVars; ++i)
    num_colloc_pts *= quad_order[i];

  // update collocPts1D, type1CollocWts1D, and type2CollocWts1D
  update_1d_collocation_points_weights(quad_order, lev_index);

  // Tensor-product quadrature: Integral of f approximated by
  // Sum_i1 Sum_i2 ... Sum_in (w_i1 w_i2 ... w_in) f(x_i1, x_i2, ..., x_in)
  // > project 1-D colloc point arrays (of potentially different type and order)
  //   into an n-dimensional stencil
  // > compute and store products of 1-D colloc weights at each point in stencil
  t1_weight_sets.sizeUninitialized(num_colloc_pts);
  if (computeType2Weights)
    t2_weight_sets.shapeUninitialized(numVars, num_colloc_pts);
  variable_sets.shapeUninitialized(numVars, num_colloc_pts);//Teuchos: col major
  colloc_key.resize(num_colloc_pts);
  UShortArray colloc_indices(numVars, 0);
  for (i=0; i<num_colloc_pts; ++i) {
    Real& t1_wt_i = t1_weight_sets[i]; t1_wt_i = 1.;
    Real*    pt_i = variable_sets[i]; // column vector i
    for (j=0; j<numVars; ++j) {
      pt_i[j]  =      collocPts1D[lev_index[j]][j][colloc_indices[j]];
      t1_wt_i *= type1CollocWts1D[lev_index[j]][j][colloc_indices[j]];
    }
    if (computeType2Weights) {
      Real* t2_wt_i = t2_weight_sets[i]; // column vector i
      for (j=0; j<numVars; ++j) {
	Real& t2_wt_ij = t2_wt_i[j]; t2_wt_ij = 1.;
	for (k=0; k<numVars; ++k)
	  t2_wt_ij *= (k == j) ?
	    type2CollocWts1D[lev_index[k]][k][colloc_indices[k]] :
	    type1CollocWts1D[lev_index[k]][k][colloc_indices[k]];
      }
    }
    colloc_key[i] = colloc_indices;
    // increment the n-dimensional collocation point index set
    if (i != num_colloc_pts-1)
      SharedPolyApproxData::increment_indices(colloc_indices, quad_order,false);
  }

#ifdef DEBUG
  PCout << "\nvariable_sets:\n";
  write_data(PCout, variable_sets, false, true, true);
  PCout << "\nt1_weight_sets:\n"; write_data(PCout, t1_weight_sets);
  if (computeType2Weights) {
    PCout << "\nt2_weight_sets:\n";
    write_data(PCout, t2_weight_sets, false, true, true);
  }
#endif
}


void IntegrationDriver::
compute_tensor_grid(const UShortArray& quad_order, const UShortArray& lev_index,
		    const SizetList& subset_indices, RealMatrix& variable_sets,
		    UShort2DArray& colloc_key)
{
  size_t i, j, k, num_colloc_pts = 1;
  for (i=0; i<numVars; ++i)
    num_colloc_pts *= quad_order[i];

  // update collocPts1D only for the subset variables
  update_1d_collocation_points_weights(quad_order, lev_index, subset_indices);

  // Tensor-product quadrature: Integral of f approximated by
  // Sum_i1 Sum_i2 ... Sum_in (w_i1 w_i2 ... w_in) f(x_i1, x_i2, ..., x_in)
  // > project 1-D colloc point arrays (of potentially different type and order)
  //   into an n-dimensional stencil
  variable_sets.shapeUninitialized(numVars, num_colloc_pts);//Teuchos: col major
  colloc_key.resize(num_colloc_pts);
  UShortArray colloc_indices(numVars, 0);
  for (i=0; i<num_colloc_pts; ++i) {
    Real*    pt_i = variable_sets[i]; // column vector i
    // assign pts for all of the variables (previous collocPts1D is sufficient
    // for non-subset variables)
    for (j=0; j<numVars; ++j)
      pt_i[j] = collocPts1D[lev_index[j]][j][colloc_indices[j]];
    colloc_key[i] = colloc_indices;
    // increment the n-dimensional collocation point index set
    if (i != num_colloc_pts-1)
      SharedPolyApproxData::increment_indices(colloc_indices, quad_order,false);
  }

#ifdef DEBUG
  PCout << "\nvariable_sets:\n";
  write_data(PCout, variable_sets, false, true, true);
#endif
}


void IntegrationDriver::
update_1d_collocation_points_weights(const UShortArray& quad_order,
				     const UShortArray& lev_index)
{
  // resize arrays
  size_t i, size_1d = collocPts1D.size(), max_index = lev_index[0];
  for (i=1; i<numVars; ++i)
    if (lev_index[i] > max_index)
      max_index = lev_index[i];
  if (max_index >= size_1d) {
    collocPts1D.resize(max_index+1); type1CollocWts1D.resize(max_index+1);
    for (i=size_1d; i<=max_index; ++i)
      { collocPts1D[i].resize(numVars); type1CollocWts1D[i].resize(numVars); }
    if (computeType2Weights) {
      type2CollocWts1D.resize(max_index+1);
      for (i=size_1d; i<=max_index; ++i)
	type2CollocWts1D[i].resize(numVars);
    }
  }
  // assign values
  for (i=0; i<numVars; ++i)
    assign_1d_collocation_points_weights(i, quad_order[i], lev_index[i]);
}


void IntegrationDriver::
update_1d_collocation_points_weights(const UShortArray& quad_order,
				     const UShortArray& lev_index,
				     const SizetList& subset_indices)
{
  // resize arrays (all variables for simplicity)
  size_t i, size_1d = collocPts1D.size(), max_index = lev_index[0];
  for (i=1; i<numVars; ++i)
    if (lev_index[i] > max_index)
      max_index = lev_index[i];
  if (max_index >= size_1d) {
    collocPts1D.resize(max_index+1); type1CollocWts1D.resize(max_index+1);
    for (i=size_1d; i<=max_index; ++i)
      { collocPts1D[i].resize(numVars); type1CollocWts1D[i].resize(numVars); }
    if (computeType2Weights) {
      type2CollocWts1D.resize(max_index+1);
      for (i=size_1d; i<=max_index; ++i)
	type2CollocWts1D[i].resize(numVars);
    }
  }
  // assign values for subset variables (for memory efficiency)
  SizetList::const_iterator cit;
  for (cit=subset_indices.begin(); cit!=subset_indices.end(); ++cit) {
    i = *cit;
    assign_1d_collocation_points_weights(i, quad_order[i], lev_index[i]);
  }
}


void IntegrationDriver::
assign_1d_collocation_points_weights(size_t i, unsigned short quad_order,
				     unsigned short lev_index)
{
  BasisPolynomial& poly_i =             polynomialBasis[i];
  RealArray&       pts_1d =      collocPts1D[lev_index][i];
  RealArray&    t1_wts_1d = type1CollocWts1D[lev_index][i];
  if (poly_i.parametric_update() || pts_1d.empty() || t1_wts_1d.empty()) {
    pts_1d    = poly_i.collocation_points(quad_order);
    t1_wts_1d = poly_i.type1_collocation_weights(quad_order);
  }
#ifdef DEBUG
  PCout << "collocPts1D[" << lev_index << "][" << i << "]:\n" << pts_1d
	<< "type1CollocWts1D[" << lev_index << "][" << i << "]:\n" << t1_wts_1d;
#endif // DEBUG
  if (computeType2Weights) {
    RealArray& t2_wts_1d = type2CollocWts1D[lev_index][i];
    if (t2_wts_1d.empty())
      t2_wts_1d = poly_i.type2_collocation_weights(quad_order);
#ifdef DEBUG
    PCout << "type2CollocWts1D[" << lev_index << "][" << i << "]:\n"
	  << t2_wts_1d;
#endif // DEBUG
  }
}

} // namespace Pecos
