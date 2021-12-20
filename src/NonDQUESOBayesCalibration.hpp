/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	 NonDBayesCalibration
//- Description: Base class for generic Bayesian inference
//- Owner:       Laura Swiler, Brian Adams
//- Checked by:
//- Version:

#ifndef NOND_QUESO_BAYES_CALIBRATION_H
#define NOND_QUESO_BAYES_CALIBRATION_H

#include "NonDBayesCalibration.hpp"

// forward declare here to isolate QUESO includes to Dakota .cpp files
namespace QUESO {
  class GslVector;
  class GslMatrix;
  class EnvOptionsValues;
  class FullEnvironment;
  template<class V, class M> class VectorSpace;
  template<class V, class M> class BoxSubset;
  template<class V, class M> class GenericScalarFunction;
  template<class V, class M> class BaseVectorRV;
  template<class V, class M> class GenericVectorRV;
  template<class V, class M> class StatisticalInverseProblem;
  class SipOptionsValues;
  class MhOptionsValues;
}

namespace Dakota {

// forward declarations of Dakota specializations
template <class V, class M> class DerivInformedPropCovTK;
template <class V, class M> class DerivInformedPropCovLogitTK;


/// Bayesian inference using the QUESO library from UT Austin

/** This class wraps the Quantification of Uncertainty for Estimation,
    Simulation, and Optimization (QUESO) library, developed as part of
    the Predictive Science Academic Alliance Program (PSAAP)-funded
    Predictive Engineering and Computational Sciences (PECOS) Center
    at UT Austin. */
class NonDQUESOBayesCalibration: public NonDBayesCalibration
{
  /// Random walk transition kernel needs callback access to QUESO details
  friend class DerivInformedPropCovTK<QUESO::GslVector, QUESO::GslMatrix>;
  /// Logit random walk transition kernel needs callback access to QUESO details
  friend class DerivInformedPropCovLogitTK<QUESO::GslVector, QUESO::GslMatrix>;

public:

  //
  //- Heading: Constructors and destructor
  //

  /// standard constructor
  NonDQUESOBayesCalibration(ProblemDescDB& problem_db, Model& model);
  /// destructor
  ~NonDQUESOBayesCalibration();

  //
  //- Heading: public member functions
  //

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void calibrate();
  void print_results(std::ostream& s, short results_state = FINAL_RESULTS);

  /// convenience function to print calibration parameters, e.g., for
  /// MAP / best parameters
  void print_variables(std::ostream& s, const RealVector& c_vars);

  /// initialize the QUESO FullEnvironment on the Dakota MPIComm
  void init_queso_environment();

  /// initialize the ASV value for preconditioned cases
  void init_precond_request_value();
  
  /// define solver options, likelihood callback, posterior RV, and
  /// inverse problem
  virtual void init_queso_solver();

  /// use derivative information from the emulator to define the proposal
  /// covariance (inverse of misfit Hessian)
  void precondition_proposal(unsigned int chain_index);

  /// perform the MCMC process
  void run_queso_solver();

  void map_pre_solve();

  /// short term option to restart the MCMC chain with updated proposal
  /// density computed from the emulator at a new starting point
  void run_chain();

  /// cache the chain to acceptanceChain and acceptedFnVals
  void cache_chain();

  /// log at most batchSize best chain points into bestSamples
  void log_best();

  /// extract batchSize points from the MCMC chain and store final
  /// aggregated set within allSamples; unique points with best
  /// conditioning are selected, as determined by pivoted LU
  void filter_chain_by_conditioning();

  /// copy bestSamples to allSamples to use in surrogate update
  void best_to_all();

  /// evaluates allSamples on iteratedModel and update the mcmcModel emulator
  /// with all{Samples,Responses}
  void update_model();

  /// compute the L2 norm of the change in emulator coefficients
  Real assess_emulator_convergence();

  /// intialize the QUESO parameter space, min, max, initial, and domain
  void init_parameter_domain();

  void init_proposal_covariance();

  /// use covariance of prior distribution for setting proposal covariance
  void prior_proposal_covariance();

  /// set proposal covariance from user-provided diagonal or matrix
  void user_proposal_covariance(const String& input_fmt, 
				const RealVector& cov_data, 
				const String& cov_filename);

  // perform sanity checks on proposalCovMatrix
  void validate_proposal();

  /// set inverse problem options calIpOptionsValues common to all solvers
  void set_ip_options(); 
  /// set MH-specific inverse problem options calIpMhOptionsValues
  void set_mh_options();
  /// update MH-specific inverse problem options calIpMhOptionsValues
  void update_chain_size(unsigned int size);

  //The likelihood routine is in the format that QUESO requires, 
  //with a particular argument list that QUESO expects. 
  //We are not using all of these arguments but may in the future.
  /// Log Likelihood function for call-back from QUESO to DAKOTA for evaluation
  static double dakotaLogLikelihood(const QUESO::GslVector& paramValues,
				    const QUESO::GslVector* paramDirection,
				    const void*             functionDataPtr,
				    QUESO::GslVector*       gradVector,
				    QUESO::GslMatrix*       hessianMatrix,
				    QUESO::GslVector*       hessianEffect);

  /// local copy_data utility from GslVector to RealVector
  void copy_gsl(const QUESO::GslVector& qv, RealVector& rv);
  /// local copy_data utility from RealVector to GslVector
  void copy_gsl(const RealVector& rv, QUESO::GslVector& qv);

  /// local copy_data utility from portion of GslVector to RealVector
  void copy_gsl_partial(const QUESO::GslVector& qv, size_t start, 
			RealVector& rv);
  /// local copy_data utility from RealVector to portion of GslVector
  void copy_gsl_partial(const RealVector& rv, QUESO::GslVector& qv,
			size_t start);

  /// local copy_data utility from GslVector to column in RealMatrix
  void copy_gsl(const QUESO::GslVector& qv, RealMatrix& rm, int i);

  /// equality tester for two GslVectors
  bool equal_gsl(const QUESO::GslVector& qv1, const QUESO::GslVector& qv2);

  //
  //- Heading: Data
  //

  /// MCMC type ("dram" or "delayed_rejection" or "adaptive_metropolis" 
  /// or "metropolis_hastings" or "multilevel",  within QUESO) 
  String mcmcType;
  /// period (number of accepted chain samples) for proposal covariance update
  int propCovUpdatePeriod;
  /// number of points to add to surrogate at each iteration
  unsigned int batchSize;
  /// the active set request value to use in proposal preconditioning
  short precondRequestValue;
  /// flag indicating user activation of logit transform option
  /** this option is useful for preventing rejection or resampling for
      out-of-bounds samples by transforming bounded domains to [-inf,inf]. */
  bool logitTransform;

  // the following QUESO objects listed in order of construction;
  // scoped_ptr more appropriate, but don't want to include QUESO
  // headers here (would be needed for checked delete on scoped_ptr)

  // TODO: see if this can be a local withing init function
  /// options for setting up the QUESO Environment
  boost::shared_ptr<QUESO::EnvOptionsValues> envOptionsValues;

  /// top-level QUESO Environment
  boost::shared_ptr<QUESO::FullEnvironment> quesoEnv;
  
  /// QUESO parameter space based on number of calibrated parameters
  boost::shared_ptr<QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> > 
  paramSpace;

  /// QUESO parameter domain: hypercube based on min/max values
  boost::shared_ptr<QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix> >
  paramDomain;

  /// initial parameter values at which to start chain
  boost::shared_ptr<QUESO::GslVector> paramInitials;

  /// random variable for the prior
  boost::shared_ptr<QUESO::BaseVectorRV<QUESO::GslVector,QUESO::GslMatrix> >
    priorRv;

  /// proposal covariance for DRAM
  boost::shared_ptr<QUESO::GslMatrix> proposalCovMatrix;

  /// optional multiplier to scale prior-based proposal covariance 
  double priorPropCovMult;

  /// general inverse problem options
  boost::shared_ptr<QUESO::SipOptionsValues> calIpOptionsValues;

  /// MH-specific inverse problem options
  boost::shared_ptr<QUESO::MhOptionsValues> calIpMhOptionsValues;

  boost::shared_ptr<QUESO::GenericScalarFunction<QUESO::GslVector,
    QUESO::GslMatrix> > likelihoodFunctionObj;

  /// random variable for the posterior
  boost::shared_ptr<QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix> >
    postRv;

  /// QUESO inverse problem solver
  boost::shared_ptr<QUESO::StatisticalInverseProblem<QUESO::GslVector,
    QUESO::GslMatrix> > inverseProb;

  /// advanced options file name (GPMSA only); settings from this file
  /// override any C++ / Dakota input file settings
  String advancedOptionsFile;

  /// Pointer to current class instance for use in static callback functions
  static NonDQUESOBayesCalibration* nonDQUESOInstance;

private:

  //
  // - Heading: Data
  // 
  
  /// cache previous expansion coefficients for assessing convergence of
  /// emulator refinement process
  RealVectorArray prevCoeffs;
  /// store MAP estimate for warm-starting pre-solves
  RealVector mapSoln;
  
};

} // namespace Dakota

#endif
