/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       NonDMultilevelStochCollocation
//- Description: Iterator for the stochastic collocation approach to UQ
//- Owner:       Mike Eldred, Sandia National Laboratories

#ifndef NOND_MULTILEVEL_STOCH_COLLOCATION_H
#define NOND_MULTILEVEL_STOCH_COLLOCATION_H

#include "NonDStochCollocation.hpp"


namespace Dakota {


/// Nonintrusive stochastic collocation approaches to uncertainty
/// quantification

/** The NonDMultilevelStochCollocation class uses a stochastic
    collocation (SC) approach to approximate the effect of parameter
    uncertainties on response functions of interest.  It utilizes the
    InterpPolyApproximation class to manage multidimensional Lagrange
    polynomial interpolants. */

class NonDMultilevelStochCollocation: public NonDStochCollocation
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// standard constructor
  NonDMultilevelStochCollocation(ProblemDescDB& problem_db, Model& model);
  /// alternate constructor
  NonDMultilevelStochCollocation(Model& model, short exp_coeffs_approach,
    const UShortArray& num_int_seq, const RealVector& dim_pref,
    short u_space_type, short refine_type, short refine_control,
    short covar_control, short ml_alloc_cntl, short ml_discrep,
    short rule_nest, short rule_growth, bool piecewise_basis, bool use_derivs);
  /// destructor
  ~NonDMultilevelStochCollocation();

  //
  //- Heading: Virtual function redefinitions
  //

  bool resize();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void initialize_u_space_model();
  void core_run();
  void assign_specification_sequence();
  void increment_specification_sequence();
  //void combined_to_active();
  void print_results(std::ostream& s, short results_state = FINAL_RESULTS);

private:

  //
  //- Heading: Utility functions
  //

  /// verify supported and define default discrepancy emulation mode
  void assign_discrepancy_mode();
  /// define the surrogate response mode for a hierarchical model in 
  /// multilevel/multifidelity expansions
  void assign_hierarchical_response_mode();

  //
  //- Heading: Data
  //

  /// type of sample allocation scheme for discretization levels / model forms
  /// within multilevel / multifidelity methods
  short mlmfAllocControl;

  /// user request of quadrature order
  UShortArray quadOrderSeqSpec;
  /// user request of sparse grid level
  UShortArray ssgLevelSeqSpec;
  /// sequence index for quadOrderSeqSpec and ssgLevelSeqSpec
  size_t sequenceIndex;
};

} // namespace Dakota

#endif
