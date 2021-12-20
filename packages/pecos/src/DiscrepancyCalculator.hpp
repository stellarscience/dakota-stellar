/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:       DiscrepancyCalculator
//- Description: Utility for computing discrepancies between a truth model
//-              and an approximation.
//- Owner:       Mike Eldred
//- Checked by:
//- Version: $Id: DiscrepancyCalculator.hpp 7024 2010-10-16 01:24:42Z mseldre $

#ifndef DISCREPANCY_CALCULATOR_H
#define DISCREPANCY_CALCULATOR_H

#include "pecos_data_types.hpp"


namespace Pecos {

class SurrogateData; // forward declaration
class ActiveKey;


/// Class for discrepancy calculations

/** The DiscrepancyCalculator class provides common functions for
    computing additive and multiplicative discrepancies. */

class DiscrepancyCalculator
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  DiscrepancyCalculator();
  /// destructor
  ~DiscrepancyCalculator();

  //
  //- Heading: Member functions
  //

  /// check for numerical issues with multiplicative discrepancy calculations
  static bool check_multiplicative(Real truth_fn, Real approx_fn,
				   short corr_order);
  /// check for numerical issues with multiplicative discrepancy calculations
  static bool check_multiplicative(const RealVector& truth_fns,
				   const RealVector& approx_fns,
				   short corr_order);

  /// compute additive 0th order correction between truth and approximate values
  static void compute_additive(Real truth_fn, Real approx_fn, Real& discrep_fn);
  /// compute additive 1st order correction between truth and
  /// approximate gradients
  static void compute_additive(const RealVector& truth_grad,
			       const RealVector& approx_grad,
			       RealVector& discrep_grad);
  /// compute additive 2nd order correction between truth and
  /// approximate Hessians
  static void compute_additive(const RealSymMatrix& truth_hess,
			       const RealSymMatrix& approx_hess,
			       RealSymMatrix& discrep_hess);
  /// compute additive corrections between truth and approximate responses
  /// for all requested orders
  static void compute_additive(Real truth_fn, const RealVector& truth_grad,
			       const RealSymMatrix& truth_hess,
			       Real approx_fn, const RealVector& approx_grad,
			       const RealSymMatrix& approx_hess,
			       Real& discrep_fn, RealVector& discrep_grad,
			       RealSymMatrix& discrep_hess,
			       short data_bits = 7);

  /// compute multiplicative 0th order correction between truth and
  /// approximate values
  static void compute_multiplicative(Real truth_fn, Real approx_fn,
				     Real& discrep_fn);
  /// compute multiplicative 1st order correction between truth and
  /// approximate gradients
  static void compute_multiplicative(Real truth_fn,
				     const RealVector& truth_grad,
				     Real approx_fn,
				     const RealVector& approx_grad,
				     RealVector& discrep_grad);
  /// compute multiplicative 2nd order correction between truth and
  /// approximate Hessians
  static void compute_multiplicative(Real truth_fn,
				     const RealVector& truth_grad,
				     const RealSymMatrix& truth_hess,
				     Real approx_fn,
				     const RealVector& approx_grad,
				     const RealSymMatrix& approx_hess,
				     RealSymMatrix& discrep_hess);
  /// compute multiplicative corrections between truth and approximate
  /// responses for all required orders
  static void compute_multiplicative(Real truth_fn,
				     const RealVector& truth_grad,
				     const RealSymMatrix& truth_hess,
				     Real approx_fn,
				     const RealVector& approx_grad,
				     const RealSymMatrix& approx_hess,
				     Real& discrep_fn, RealVector& discrep_grad,
				     RealSymMatrix& discrep_hess,
				     short data_bits = 7);

  /// compute discrepancies for arrays of SurrogateDataResp
  static void compute(const SDRArray& hf_sdr_array,
		      const SDRArray& lf_sdr_array,
		      SDRArray& delta_sdr_array, short combine_type);
  /// compute discrepancies between two model keys in surr_data and store
  /// results in the discrepancy key
  static void compute(SurrogateData& surr_data, const ActiveKey& delta_key,
		      short combine_type);

  /*
  /// function for applying additive correction to an approximate response
  void apply_additive(const Variables& vars, Response& approx_response);
  /// function for applying multiplicative correction to an approximate response
  void apply_multiplicative(const Variables& vars,
			    Response& approx_response);

  /// update correctionType
  void correction_type(short corr_type);
  /// return correctionType
  short correction_type() const;
  /// update correctionOrder
  void correction_order(short order);
  /// return correctionOrder
  short correction_order() const;
  /// update dataOrder
  void data_order(short order);
  /// return dataOrder
  short data_order() const;
  */

private:

  //
  //- Heading: Convenience functions
  //

  //
  //- Heading: Data
  //

  /*
  /// approximation correction approach to be used: NO_CORRECTION,
  /// ADDITIVE_CORRECTION, MULTIPLICATIVE_CORRECTION, or COMBINED_CORRECTION.
  short correctionType;
  /// approximation correction order to be used: 0, 1, or 2
  short correctionOrder;
  /// order of correction data in 3-bit format: overlay of 1 (value),
  /// 2 (gradient), and 4 (Hessian)
  short dataOrder;
  */
};


inline DiscrepancyCalculator::DiscrepancyCalculator()
  //correctionType(NO_CORRECTION), correctionOrder(0), dataOrder(1)
{ }


inline DiscrepancyCalculator::~DiscrepancyCalculator()
{ }


inline void DiscrepancyCalculator::
compute_additive(Real truth_fn, const RealVector& truth_grad,
		 const RealSymMatrix& truth_hess,
		 Real approx_fn, const RealVector& approx_grad,
		 const RealSymMatrix& approx_hess,
		 Real& discrep_fn, RealVector& discrep_grad,
		 RealSymMatrix& discrep_hess, short data_bits)
{
  if (data_bits & 1)
    compute_additive(truth_fn, approx_fn, discrep_fn);
  if (data_bits & 2)
    compute_additive(truth_grad, approx_grad, discrep_grad);
  if (data_bits & 4)
    compute_additive(truth_hess, approx_hess, discrep_hess);
}


inline void DiscrepancyCalculator::
compute_multiplicative(Real truth_fn, const RealVector& truth_grad,
		       const RealSymMatrix& truth_hess,
		       Real approx_fn, const RealVector& approx_grad,
		       const RealSymMatrix& approx_hess,
		       Real& discrep_fn, RealVector& discrep_grad,
		       RealSymMatrix& discrep_hess, short data_bits)
{
  if (data_bits & 1)
    compute_multiplicative(truth_fn, approx_fn, discrep_fn);
  if (data_bits & 2)
    compute_multiplicative(truth_fn, truth_grad, approx_fn, approx_grad,
			   discrep_grad);
  if (data_bits & 4)
    compute_multiplicative(truth_fn, truth_grad, truth_hess, approx_fn,
			   approx_grad, approx_hess, discrep_hess);
}

} // namespace Pecos

#endif
