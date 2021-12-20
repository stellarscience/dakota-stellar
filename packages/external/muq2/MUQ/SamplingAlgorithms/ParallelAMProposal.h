#ifndef PARALLELAMPROPOSAL_H_
#define PARALLELAMPROPOSAL_H_

#include "MUQ/SamplingAlgorithms/AMProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    class ParallelAMProposal : public AMProposal {
    public:

      ParallelAMProposal(boost::property_tree::ptree                  pt,
                         std::shared_ptr<AbstractSamplingProblem>     problem);

      ParallelAMProposal(boost::property_tree::ptree                  pt,
                         std::shared_ptr<AbstractSamplingProblem>     problem,
                         std::shared_ptr<parcer::Communicator> const& newcomm);

      virtual ~ParallelAMProposal() = default;

      /// Adapt the proposal after each step
      /**
	 Adapt the proposal covariance.
	 @param[in] t The current step
	 @param[in] state The current state
       */
      virtual void Adapt(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& states) override;

    private:

      /// The total number of samples used for the adaption
      unsigned int totSamps = 0;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
