#ifndef MHKERNEL_H_
#define MHKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
      @ingroup MCMCKernels
      @class MHKernel
      @brief An implementation of the standard Metropolis-Hastings transition kernel.
      @details <B>Configuration Parameters:</B>

      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "Proposal"  | String | - | A string pointing to a block of proposal options. |

     */
    class MHKernel : public TransitionKernel {
    public:

      MHKernel(boost::property_tree::ptree const& pt,
               std::shared_ptr<AbstractSamplingProblem> problem);

      MHKernel(boost::property_tree::ptree const& pt,
               std::shared_ptr<AbstractSamplingProblem> problem,
               std::shared_ptr<MCMCProposal> proposalIn);

      virtual ~MHKernel() = default;

      virtual inline std::shared_ptr<MCMCProposal> Proposal() {return proposal;};

      virtual void PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state) override;

      virtual std::vector<std::shared_ptr<SamplingState>> Step(unsigned int const t, std::shared_ptr<SamplingState> prevState) override;

      virtual void PrintStatus(std::string prefix) const override;


      virtual inline double AcceptanceRate() const {return double(numAccepts)/double(numCalls);};


    protected:
      std::shared_ptr<MCMCProposal> proposal;

      unsigned int numCalls = 0;
      unsigned int numAccepts = 0;

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
