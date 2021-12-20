/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:       DiscrepancyCalculator
//- Description: Implementation code for the DiscrepancyCalculator class
//- Owner:       Mike Eldred
//- Checked by:

#include "DiscrepancyCalculator.hpp"
#include "SurrogateData.hpp"

static const char rcsId[]="@(#) $Id: DiscrepancyCalculator.cpp 7024 2010-10-16 01:24:42Z mseldre $";


namespace Pecos {

bool DiscrepancyCalculator::
check_multiplicative(Real truth_fn, Real approx_fn, short corr_order)
{
  // Multiplicative will fail if response functions are near zero.
  //   0th order:     a truth_val == 0 causes a zero scaling which will cause
  //                  optimization failure; an approx_val == 0 will cause a
  //                  division by zero FPE.
  //   1st/2nd order: a truth_val == 0 is OK (so long as the total scaling
  //                  function != 0); an approx_val == 0 will cause a division
  //                  by zero FPE.
  return
    ( std::fabs(approx_fn) < Pecos::SMALL_NUMBER ||
      ( corr_order == 0 && std::fabs(truth_fn) < Pecos::SMALL_NUMBER ) );
}


bool DiscrepancyCalculator::
check_multiplicative(const RealVector& truth_fns, const RealVector& approx_fns,
		     short corr_order)
{
  bool bad_scaling = false;
  size_t i, num_fns = std::min(truth_fns.length(), approx_fns.length());
  for (i=0; i<num_fns; ++i)
    if ( std::fabs(approx_fns[i]) < Pecos::SMALL_NUMBER ||
	 ( corr_order == 0 && std::fabs(truth_fns[i]) < Pecos::SMALL_NUMBER ) )
      { bad_scaling = true; break; }

  return bad_scaling;
}


void DiscrepancyCalculator::
compute_additive(Real truth_fn, Real approx_fn, Real& discrep_fn)
{
  // -----------------------------
  // Additive 0th order correction
  // -----------------------------
  discrep_fn = truth_fn - approx_fn;
}


void DiscrepancyCalculator::
compute_additive(const RealVector& truth_grad, const RealVector& approx_grad,
		 RealVector& discrep_grad)
{
  // -----------------------------
  // Additive 1st order correction
  // -----------------------------
  size_t v, num_v = std::min(truth_grad.length(), approx_grad.length());
  if (discrep_grad.length() != num_v) discrep_grad.sizeUninitialized(num_v);
  for (v=0; v<num_v; ++v)
    discrep_grad[v] = truth_grad[v] - approx_grad[v];
}


void DiscrepancyCalculator::
compute_additive(const RealSymMatrix& truth_hess,
		 const RealSymMatrix& approx_hess, RealSymMatrix& discrep_hess)
{
  // -----------------------------
  // Additive 2nd order correction
  // -----------------------------
  size_t r, c, num_v = std::min(truth_hess.numRows(), approx_hess.numRows());
  if (discrep_hess.numRows() != num_v)
    discrep_hess.shapeUninitialized(num_v);
  for (r=0; r<num_v; ++r)
    for (c=0; c<=r; ++c) // lower half
      discrep_hess(r,c) = truth_hess(r,c) - approx_hess(r,c);
}


void DiscrepancyCalculator::
compute_multiplicative(Real truth_fn, Real approx_fn, Real& discrep_fn)
{
  // -----------------------------------
  // Multiplicative 0th order correction
  // -----------------------------------
  discrep_fn = truth_fn / approx_fn;
}


void DiscrepancyCalculator::
compute_multiplicative(Real truth_fn, const RealVector& truth_grad,
		       Real approx_fn, const RealVector& approx_grad,
		       RealVector& discrep_grad)
{
  // -----------------------------------
  // Multiplicative 1st order correction
  // -----------------------------------
  // The beta-correction method is based on the work of Chang and Haftka,
  // and Alexandrov.  It is a multiplicative correction like the "scaled"
  // correction method, but it uses gradient information to achieve
  // 1st-order consistency (matches the high-fidelity function values and
  // the high-fidelity gradients at the center of the approximation region).
  size_t v, num_v = std::min(truth_grad.length(), approx_grad.length());
  if (discrep_grad.length() != num_v) discrep_grad.sizeUninitialized(num_v);
  Real ratio = truth_fn / approx_fn;
  for (v=0; v<num_v; ++v)
    discrep_grad[v] = (truth_grad[v] - approx_grad[v] * ratio) / approx_fn;
}


void DiscrepancyCalculator::
compute_multiplicative(Real truth_fn, const RealVector& truth_grad,
		       const RealSymMatrix& truth_hess,
		       Real approx_fn, const RealVector& approx_grad,
		       const RealSymMatrix& approx_hess,
		       RealSymMatrix& discrep_hess)
{
  // -----------------------------------
  // Multiplicative 2nd order correction
  // -----------------------------------
  size_t r, c, num_v = std::min(truth_hess.numRows(), approx_hess.numRows());
  if (discrep_hess.numRows() != num_v)
    discrep_hess.shapeUninitialized(num_v);
  // consider use of Teuchos assign and operator-=
  Real ratio2 = 2. * truth_fn / approx_fn, approx_sq = approx_fn * approx_fn;
  for (r=0; r<num_v; ++r)
    for (c=0; c<=r; ++c) // lower half
      discrep_hess(r,c) = ( truth_hess(r,c) * approx_fn - approx_hess(r,c) *
	truth_fn + ratio2 * approx_grad[r] * approx_grad[c] - truth_grad[r] *
	approx_grad[c] - approx_grad[r] * truth_grad[c] ) / approx_sq;
}


void DiscrepancyCalculator::
compute(const SDRArray& hf_sdr_array, const SDRArray& lf_sdr_array,
	SDRArray& delta_sdr_array, short combine_type)
{
  size_t i, num_pts = std::min(hf_sdr_array.size(), lf_sdr_array.size());

  // Note: SurrogateData::size_active_sdr() is used to allocate delta_sdr_array

  switch (combine_type) {
  case MULT_COMBINE:
    for (i=0; i<num_pts; ++i) {
      const SurrogateDataResp& lf_sdr  =    lf_sdr_array[i];
      const SurrogateDataResp& hf_sdr  =    hf_sdr_array[i];
      SurrogateDataResp&    delta_sdr  = delta_sdr_array[i];
      short                 delta_bits = delta_sdr.active_bits();
      short                 corr_order = (delta_bits & 2) ? 1 : 0;
      if (check_multiplicative(hf_sdr.response_function(),
			       lf_sdr.response_function(), corr_order)) {
	PCerr << "Error: numerical FPE in computing multiplicative discrepancy."
	      << "\n       Please change to additive discrepancy." << std::endl;
	abort_handler(-1);
      }
      if (delta_bits & 1)
	compute_multiplicative(hf_sdr.response_function(),
			       lf_sdr.response_function(),
			       delta_sdr.response_function_view());
      if (delta_bits & 2) {
	RealVector delta_grad(delta_sdr.response_gradient_view());
	compute_multiplicative(hf_sdr.response_function(),
			       hf_sdr.response_gradient(),
			       lf_sdr.response_function(),
			       lf_sdr.response_gradient(), delta_grad);
      }
      if (delta_bits & 4) {
	RealSymMatrix delta_hess(delta_sdr.response_hessian_view());
	compute_multiplicative(hf_sdr.response_function(),
			       hf_sdr.response_gradient(),
			       hf_sdr.response_hessian(),
			       lf_sdr.response_function(),
			       lf_sdr.response_gradient(),
			       lf_sdr.response_hessian(), delta_hess);
      }
    }
    break;
  default: //case ADD_COMBINE: (correction specification not required)
    for (i=0; i<num_pts; ++i) {
      const SurrogateDataResp& lf_sdr  =    lf_sdr_array[i];
      const SurrogateDataResp& hf_sdr  =    hf_sdr_array[i];
      SurrogateDataResp&    delta_sdr  = delta_sdr_array[i];
      short                 delta_bits = delta_sdr.active_bits();
      if (delta_bits & 1)
	compute_additive(hf_sdr.response_function(), lf_sdr.response_function(),
			 delta_sdr.response_function_view());
      if (delta_bits & 2) {
	RealVector delta_grad(delta_sdr.response_gradient_view());
	compute_additive(hf_sdr.response_gradient(), lf_sdr.response_gradient(),
			 delta_grad);
      }
      if (delta_bits & 4) {
	RealSymMatrix delta_hess(delta_sdr.response_hessian_view());
	compute_additive(hf_sdr.response_hessian(), lf_sdr.response_hessian(),
			 delta_hess);
      }
    }
    break;
  }
}


void DiscrepancyCalculator::
compute(SurrogateData& surr_data, const ActiveKey& delta_key,
	short combine_type)
{
  ActiveKey hf_key, lf_key;  delta_key.extract_keys(hf_key, lf_key);

  const std::map<ActiveKey, SDRArray>& resp_map
    = surr_data.response_data_map();
  std::map<ActiveKey, SDRArray>::const_iterator
    lf_cit = resp_map.find(lf_key), hf_cit = resp_map.find(hf_key);
  if (lf_cit == resp_map.end() || hf_cit == resp_map.end()) {
    PCerr << "Error: key lookup failure for individual fidelity in Discrepancy"
	  << "Calculator::compute()" << std::endl;
    abort_handler(-1);
  }

  surr_data.active_key(delta_key);
  surr_data.variables_data(surr_data.variables_data(hf_key)); // shallow copies
  surr_data.anchor_index(surr_data.anchor_index(hf_key));
  surr_data.pop_count_stack(surr_data.pop_count_stack(hf_key));

  // TO DO: do this more incrementally (based on curr state of delta sdr_array)
  const SDRArray& hf_sdr_array = hf_cit->second;
  surr_data.size_active_sdr(hf_sdr_array);
  SDRArray& delta_sdr_array = surr_data.response_data();
  compute(hf_sdr_array, lf_cit->second, delta_sdr_array, combine_type);

  // compute discrepancy faults from scratch for delta_sdr_array
  // (aggregates LF,HF failures)
  surr_data.data_checks();
}

} // namespace Pecos
