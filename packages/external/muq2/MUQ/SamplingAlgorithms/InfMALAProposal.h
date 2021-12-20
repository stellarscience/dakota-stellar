#ifndef INFMALAPROPOSAL_H_
#define INFMALAPROPOSAL_H_

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/ModPiece.h"

#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
    @ingroup MCMCProposals
    @class InfMALAProposal
    @brief An implement of the dimension-independent MALA (or Inf-MALA) proposal.

    @details
    This class implements the a modification of the \f$\infty-\f$MALA proposal of
    Beskos et al., 2017.   The proposal takes the form
    \f[
    u^\prime = \bar{u} + \sqrt{1-\beta^2} (u_c-\bar{u}) + \beta\left( z + \tau \Sigma\, \nabla \log \pi(u_c) \right)
    \f]
    where \f$\bar{u}\f$ is an estimate of the posterior mean (e.g., the prior mean or posterior MAP), \f$u_c\f$ is the current state of the chain, \f$z\sim N(0,C)\f$ is a normal
    random variable with a strategically chosen covariance \f$C\f$ (often the prior covariance), and \f$u^\prime\f$
    is the propsed point.  The parameters \f$\beta\f$ and \f$\tau\f$ control the size and offset or the proposal.
    When \f$\tau=0\f$, this proposal is identical to the standard CrankNicolsonProposal.  Also note that unlike the usual MALA
    proposal, the variance of the \f$\infty-\f$MALA proposal does not depend on the step size \f$\tau\f$.

    <B>Configuration Parameters:</B>
    Parameter Key | Type | Default Value | Description |
    ------------- | ------------- | ------------- | ------------- |
    "Beta"  | double | 0.5 | The proposal scaling \f$\beta\f$ defined above. |
    "StepSize" | double | 1.0 | The parameter \f$\tau\f$ governing the MALA shift in the direction of the gradient. |
    "PriorNode" | string | - | (Optional)  If specified, this class assumes the target density was constructed from a WorkGraph and will set the value of the covariance \f$C\f$ to the covariance of the Gaussian density at the specified node. |
    */
    class InfMALAProposal : public CrankNicolsonProposal {
    public:

      InfMALAProposal(boost::property_tree::ptree           const& pt,
                      std::shared_ptr<AbstractSamplingProblem>     prob,
                      std::shared_ptr<muq::Modeling::GaussianBase> prior);

      InfMALAProposal(boost::property_tree::ptree       const& pt,
                      std::shared_ptr<AbstractSamplingProblem> prob);

      virtual ~InfMALAProposal() = default;

    protected:

      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> const& currentState) override;

      virtual double LogDensity(std::shared_ptr<SamplingState> const& currState,
                                std::shared_ptr<SamplingState> const& propState) override;


      /** Returns the product of the proposal covariance times gradient of the
          log target density at the current state.  Checks the metadata of the
          currentstate first to see if this has already been computed.
      */
      Eigen::VectorXd GetSigmaGrad(std::shared_ptr<SamplingState> const& currentState) const;

      double stepSize;

    }; // class InfMALAProposal

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
