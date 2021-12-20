#ifndef AMPROPOSAL_H_
#define AMPROPOSAL_H_

#include "MUQ/SamplingAlgorithms/MHProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @ingroup MCMCProposals
        @class AMProposal
        @brief An implemental of the Adaptive Metropolis algorithm
        @details <B>Configuration Parameters:</B>

        Parameter Key | Type | Default Value | Description |
        ------------- | ------------- | ------------- | ------------- |
        "ProposalVariance" | Double | 1.0 | Initial proposal variance before adaptation begins. |
        "AdaptSteps"  | Integer | - | How often the proposal covariance should be updated. |
        "AdaptStart"  | Integer | - | How many steps should be taken before beginning to adapt the covariance. |
        "AdaptScale"  | Double  | - | Scaling of the sample covariance used to define the proposal covariance: \f$\Sigma_{prop} = \alpha \hat{\Sigma}\f$ |
    */
    class AMProposal : public MCMCProposal {
    public:

      AMProposal(boost::property_tree::ptree                     pt,
                 std::shared_ptr<AbstractSamplingProblem> const& prob);

      /**
      Construct the adaptive Metropolis proposal with a specific estimate of the
      target covariance matrix.  Note that initailCov will be scaled by the adaptScale
      parameter set in the ptree and should therefore be an estiamte of the
      target covariance, not a specification of the proposal covariance.
      */
      AMProposal(boost::property_tree::ptree                     pt,
                 std::shared_ptr<AbstractSamplingProblem> const& prob,
                 Eigen::MatrixXd                          const& initialCov);


      virtual ~AMProposal() = default;

      /// Adapt the proposal after each step
      /**
	 Adapt the proposal covariance.
	 @param[in] t The current step
	 @param[in] state The current state
       */
      virtual void Adapt(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& states) override;

      /** Returns the current proposal covariance. */
      virtual Eigen::MatrixXd ProposalCovariance() const;

    private:

      static Eigen::MatrixXd ConstructCovariance(boost::property_tree::ptree                                const& pt ,
                                                 std::shared_ptr<AbstractSamplingProblem> const& prob);


      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> const& currentState) override;

      virtual double LogDensity(std::shared_ptr<SamplingState> const& currState,
                                std::shared_ptr<SamplingState> const& propState) override;

      /// Samples seen between visits
      Eigen::MatrixXd newSamps;
      unsigned int numNewSamps=0;

      /// Total number of samples that have been used so far to compute mean and covariance
      unsigned int numAdaptSamps=0;

      /// The current mean
      Eigen::VectorXd mean;

      /// The current proposal covariance
      Eigen::MatrixXd propCov;

      /// The Cholesky factorization of the current proposal covariance
      Eigen::LLT<Eigen::MatrixXd, Eigen::Lower> propChol;

      /// How frequently should we update the adaption?
      const unsigned int adaptSteps;

      /// When should we start adapting?
      const unsigned int adaptStart;

      /// After how many iterations should we stop adapting?
      const unsigned int adaptEnd;

      // Multiplicative scaling of the sample covariance
      const double adaptScale;

    };

  } // namespace SamplingAlgorithms
} // namespace muq

#endif
