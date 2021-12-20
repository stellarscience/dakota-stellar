#ifndef MIXTUREPROPOSAL_H_
#define MIXTUREPROPOSAL_H_

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @ingroup MCMCProposals
        @class MixtureProposal
        @brief This class implements a weighted mixture of other proposals.
        @details This class implements proposals of the form
        \f[
        p(x^\prime | x) = \sum_{i=1}^N w_i p_i(x^\prime | x),
        \f]
        where \f$N\f$ is the total number of proposal components and \f$w_i\f$
        is the weight of the \f$i^{th}\f$ component.  Note that the weights must
        sum to one and are thus normalized after being passed to this class.

        The components \f$p_i(x^\prime | x)\f$ can be any other MCMC proposal.

        <B>Configuration Parameters:</B>

        Parameter Key | Type | Default Value | Description |
        ------------- | ------------- | ------------- | ------------- |
        "Components"  | string | - | A comma separated list with the names of other blocks that define each component in the proposal mixture. |
        "Weights"   | string | 1.0,1.0,... | A comma separated list of floats defining the weights (potentially unnormalized) for each mixture component.  If not specified, even weights are prescribed. |
    */
    class MixtureProposal : public MCMCProposal {
    public:

      MixtureProposal(boost::property_tree::ptree                pt,
                 std::shared_ptr<AbstractSamplingProblem> const& prob);

      MixtureProposal(boost::property_tree::ptree                       pt,
                      std::shared_ptr<AbstractSamplingProblem>   const& prob,
                      std::vector<std::shared_ptr<MCMCProposal>> const& proposals,
                      std::vector<double>                        const& weights);

      virtual ~MixtureProposal() = default;

    protected:

      /** Extracts the proposal weights from a ptree. */
      static std::vector<double> GetWeights(boost::property_tree::ptree const& pt);

      /** Constructs the mixture proposals based on properties in the ptree. */
      static std::vector<std::shared_ptr<MCMCProposal>> GetProposals(boost::property_tree::ptree              const& pt,
                                                                     std::shared_ptr<AbstractSamplingProblem> const& prob);

      /// The proposal distribution
      std::vector<std::shared_ptr<MCMCProposal>> proposals;
      std::vector<double> weights;

      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> const& currentState) override;

      virtual double LogDensity(std::shared_ptr<SamplingState> const& currState,
                                std::shared_ptr<SamplingState> const& propState) override;


    };

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
