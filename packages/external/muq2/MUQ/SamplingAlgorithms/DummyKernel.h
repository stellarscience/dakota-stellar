#ifndef DUMMYKERNEL_H_
#define DUMMYKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief MCMC kernel that always accepts proposals.
        @details This kernel simply accepts any proposal. It can be used to turn a MCMC method into a simple
        Monte Carlo method with the proposal density as its parameter distribution.
     */
    class DummyKernel : public TransitionKernel {
    public:

      DummyKernel(boost::property_tree::ptree const& pt,
                  std::shared_ptr<AbstractSamplingProblem> problem,
                  std::shared_ptr<MCMCProposal> proposal);

      ~DummyKernel();

      virtual std::shared_ptr<MCMCProposal> Proposal(){return proposal;};

      virtual void PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state) override;

      virtual void PrintStatus(std::string prefix) const override;

      virtual std::vector<std::shared_ptr<SamplingState>> Step(unsigned int const t, std::shared_ptr<SamplingState> prevState) override;

    protected:
      std::shared_ptr<MCMCProposal> proposal;

      unsigned int numCalls = 0;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
