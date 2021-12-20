#ifndef CONCATENATINGINTERPOLATION_H_
#define CONCATENATINGINTERPOLATION_H_

#include "MUQ/SamplingAlgorithms/MIInterpolation.h"
#include "MUQ/Utilities/MultiIndices/MultiIndex.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
     * @brief A simple implementation concatenating coarse and fine sample vectors.
     *
     * @details This interpolation takes the coarse sample vector and appends additional fine components from the fine sample.
     * This behaviour matches the one assumed by the theoretical MLMCMC literature. Should be sufficient for most applications.
     */
    class ConcatenatingInterpolation : public MIInterpolation {
    public:
      ConcatenatingInterpolation(std::shared_ptr<muq::Utilities::MultiIndex> const& index) : index(index) {
    	}

    	virtual std::shared_ptr<SamplingState> Interpolate (std::shared_ptr<SamplingState> const& coarseProposal, std::shared_ptr<SamplingState> const& fineProposal) override {
    		int fine_part_size = fineProposal->state[0].size() - coarseProposal->state[0].size();

    		Eigen::VectorXd interpolatedState(fineProposal->state[0].size());
    		interpolatedState << coarseProposal->state[0], fineProposal->state[0].tail(fine_part_size);

    		return std::make_shared<SamplingState>(interpolatedState);
    	}

    private:
      std::shared_ptr<muq::Utilities::MultiIndex> index;
    };

  }
}


#endif
