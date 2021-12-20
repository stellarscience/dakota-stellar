#ifndef MIDUMMYKERNEL_H_
#define MIDUMMYKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"
#include "MUQ/SamplingAlgorithms/MIInterpolation.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief Dummy kernel for Multiindex MC methods
        @details This kernel simply accepts any proposal. It can be used to turn a MLMCMC/MIMCMC method into a
        Multilevel/Multiindex Monte Carlo method with the proposal density as its parameter distribution.
     */
    class MIDummyKernel : public TransitionKernel {
    public:

      MIDummyKernel(boost::property_tree::ptree const& pt,
                    std::shared_ptr<AbstractSamplingProblem> problem,
                    std::shared_ptr<MCMCProposal> proposal,
                    std::shared_ptr<MCMCProposal> coarse_proposal,
                    std::shared_ptr<MIInterpolation> proposalInterpolation,
                    std::shared_ptr<SingleChainMCMC> coarse_chain);

      ~MIDummyKernel();

      virtual std::shared_ptr<MCMCProposal> Proposal(){return proposal;};

      virtual void PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state) override;

      virtual void PrintStatus(std::string prefix) const override;

      virtual std::vector<std::shared_ptr<SamplingState>> Step(unsigned int const t, std::shared_ptr<SamplingState> prevState) override;

    protected:
      std::shared_ptr<MCMCProposal> proposal;
      std::shared_ptr<MCMCProposal> coarse_proposal;
      std::shared_ptr<MIInterpolation> proposalInterpolation;
      std::shared_ptr<SingleChainMCMC> coarse_chain;

      unsigned int numCalls = 0;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
