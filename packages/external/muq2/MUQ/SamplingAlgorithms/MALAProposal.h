#ifndef MALAPROPOSAL_H_
#define MALAPROPOSAL_H_

#include "MUQ/Modeling/Distributions/GaussianBase.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @ingroup MCMCProposals
        @class MALAProposal
        @brief Implementation preconditioned Langevin proposal used in the MALA algorithm.
        @details
        Defines a proposal of the form
        \f[
        x^{\prime} = x^{(k)} + \tau \Sigma \nabla \log \pi(x^{(k)}) + \sqrt(2\tau) \Sigma^{1/2} \tilde{z},
        \f]
        where \f$x^{(k)}\f$ is the current position of the chain, \f$\Sigma\f$
         is the proposal covariance, \f$\tau\f$ is a stepsize parameter,
        \f$z\f$ is a standard normal random variable, and \f$x^{\prime}\f$ is the
        proposed move.  Note that this proposal can also be defined in terms of
        a Gaussian random variable \f$z=\Sigma^{1/2} \tilde{z}\f$, which has
        mean zero and covariance \f$\Sigma\f$.  The Gaussian distribution for
        \f$z\f$ can be specified in the constructor of this proposal.  If not
        specified, a default identity covariance is employed, \f$\Sigma=I\f$.

        <B>Configuration Parameters:</B>

        Parameter Key | Type | Default Value | Description |
        ------------- | ------------- | ------------- | ------------- |
        "StepSize"  | Double | - | The step size parameter \f$\tau\f$. |
    */
    class MALAProposal : public MCMCProposal {
    public:

      MALAProposal(boost::property_tree::ptree                   pt,
                 std::shared_ptr<AbstractSamplingProblem> const& probIn);

      MALAProposal(boost::property_tree::ptree                   pt,
                 std::shared_ptr<AbstractSamplingProblem> const& probIn,
                 std::shared_ptr<muq::Modeling::GaussianBase>    zDistIn);

      virtual ~MALAProposal() = default;

    protected:

      double stepSize;

      /// The proposal distribution
      std::shared_ptr<muq::Modeling::GaussianBase> zDist;

      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> const& currentState) override;

      virtual double
      LogDensity(std::shared_ptr<SamplingState> const& currState,
                 std::shared_ptr<SamplingState> const& propState) override;


    };

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
