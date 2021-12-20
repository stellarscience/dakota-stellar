#ifndef INVERSEGAMMAPROPOSAL_H_
#define INVERSEGAMMAPROPOSAL_H_

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/InverseGamma.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
      @ingroup MCMCProposals
      @class InverseGammaProposal
      @brief Defines a proposal using the analytic conditional Inverse Gamma distribution for the variance of a Gaussian distribution
      @details Consider a Metropolis-Within-Gibbs sampler for a problem where
      the inverse Gamma distribution is used to model the variance of a Gaussian
      distribution.  In that setting, the distribution of the variance given the
      state of the Gaussian random variable is known analytically and can be sampled
      from directly.  This proposal leverages that fact. It assumes a prior distribution
      over the variance is given by \f$\sigma^2 \sim IG(\alpha,\beta)\f$ and that
      the Gaussian random variable has zero mean and covariance \f$\simga^2 I\f$.  Then,
      given an observation of the Gaussian random variable \f$x=[x_1,x_2, \ldots, x_N]\f$,
      we know that
      \f[
      \sigma^2 | x \sim IG(\alpha + \frac{N}{2}, \beta + \frac{1}{2}\sum_{i=1}^Nx_i^2.
      \f]
      This is the proposal density defined by this class.

      The class assumes that the provided AbstractSamplingProblem is an instance
      of the SamplingProblem class, which contains a WorkGraph defining the
      relationship between the inverse gamma and Gaussian distributions.

      <B>Configuration Parameters:</B>

      Parameter Key      | Type   | Default Value | Description |
      ------------------ | ------ | ------------- | ------------- |
      "InverseGammaNode" | string | -             |  The name of the node in the WorkGraph that contains the InverseGamma density. |
      "GaussianNode"     | string | -             |  The name of the node in the WorkGraph that contains the Gaussian distribution with a variance input. |

    */
    class InverseGammaProposal : public MCMCProposal {
    public:

      InverseGammaProposal(boost::property_tree::ptree       const& pt,
                           std::shared_ptr<AbstractSamplingProblem> prob);

      virtual ~InverseGammaProposal() = default;

    protected:

      /// The prior value of alpha
      const Eigen::VectorXd alpha;

      /// The prior value of beta
      const Eigen::VectorXd beta;

      /// The index of the Gaussian block
      const unsigned int gaussBlock;

      /// The mean of the Gaussian distribution
      const Eigen::VectorXd gaussMean;

      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> const& currentState) override;

      virtual double LogDensity(std::shared_ptr<SamplingState> const& currState,
                                std::shared_ptr<SamplingState> const& propState) override;


      static Eigen::VectorXd ExtractMean(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& gaussNode);

      static std::shared_ptr<muq::Modeling::InverseGamma> ExtractInverseGamma(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& igNode);
      static Eigen::VectorXd ExtractAlpha(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& igNode);
      static Eigen::VectorXd ExtractBeta(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& igNode);
      static unsigned int ExtractGaussInd(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& gaussNode);
    };

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
