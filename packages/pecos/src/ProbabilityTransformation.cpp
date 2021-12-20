/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "pecos_stat_util.hpp"
#include "NatafTransformation.hpp"

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
  probTransRep(NULL), referenceCount(1)
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
    return the transformed space variables to their original states.
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
*/


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
numerical_design_jacobian(const RealVector& x_vars, bool xs,
			  RealMatrix& num_jacobian_xs, bool zs,
			  RealMatrix& num_jacobian_zs,
			  SizetMultiArrayConstView cv_ids,
			  SizetMultiArrayConstView acv_ids,
			  const SizetArray& acv_map1_indices,
			  const ShortArray& acv_map2_targets)
{
  if (probTransRep) // envelope fwd to letter
    probTransRep->numerical_design_jacobian(x_vars, xs, num_jacobian_xs, zs,
					    num_jacobian_zs, cv_ids, acv_ids,
					    acv_map1_indices, acv_map2_targets);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine numerical_design_jacobian"
          << "() virtual fn.\nNo default defined at ProbabilityTransformation "
	  << "base class.\n" << std::endl;
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
