#ifndef FAULT_TOLERANCE_HPP
#define FAULT_TOLERANCE_HPP

#include "LinearAlgebra.hpp"
#include "MathTools.hpp"
#include "pecos_data_types.hpp"

namespace Pecos {

struct FaultInfo
{
  size_t constr_eqns, anchor_fn, anchor_grad, num_data_pts_fn,
    num_data_pts_grad, total_eqns, num_surr_data_pts, num_vars, num_grad_rhs;
  bool under_determined, reuse_solver_data, use_derivatives;
  
  void set_info( size_t constr_eqns_in, size_t anchor_fn_in, 
		 size_t anchor_grad_in, 
		 bool under_determined_in, size_t num_data_pts_fn_in, 
		 size_t num_data_pts_grad_in, bool reuse_solver_data_in,
		 size_t total_eqns_in, size_t num_surr_data_pts_in,
		 size_t num_vars_in, bool use_derivatives_in, 
		 size_t num_grad_rhs_in )
  {
    constr_eqns = constr_eqns_in; anchor_fn = anchor_fn_in; 
    anchor_grad = anchor_grad_in; under_determined = under_determined_in; 
    num_data_pts_fn = num_data_pts_fn_in; 
    num_data_pts_grad = num_data_pts_grad_in; 
    reuse_solver_data = reuse_solver_data_in;
    total_eqns = total_eqns_in; num_surr_data_pts = num_surr_data_pts_in;
    num_vars = num_vars_in; use_derivatives = use_derivatives_in;
    num_grad_rhs = num_grad_rhs_in;
  }
};

// \todo modified from OrthogPolyApproximation. Perhaps consider removing from
// OrthogPolyApproximation
void fail_booleans(SizetShortMap::const_iterator& fit, size_t j,
		   bool& add_val, bool& add_grad,
		   const SizetShortMap& failed_response_data );

void remove_faulty_data( RealMatrix &A, RealMatrix &B, 
			 RealMatrix &points,
			 IntVector &index_mapping,
			 FaultInfo fault_info,
			 const SizetShortMap& failed_resp_data );

} //namespace Pecos

#endif // FAULT_TOLERANCE_HPP
