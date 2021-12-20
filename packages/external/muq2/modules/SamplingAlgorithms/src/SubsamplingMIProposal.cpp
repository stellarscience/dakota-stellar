#include "MUQ/SamplingAlgorithms/SubsamplingMIProposal.h"
namespace muq {
  namespace SamplingAlgorithms {

    SubsamplingMIProposal::SubsamplingMIProposal (pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> prob, std::shared_ptr<SingleChainMCMC> coarseChain)
     : MCMCProposal(pt,prob), coarseChain(coarseChain),
       subsampling(pt.get<int>("Subsampling"))
    {}

    std::shared_ptr<SamplingState> SubsamplingMIProposal::Sample(std::shared_ptr<SamplingState> const& currentState) {

      sampleID += subsampling+1;
      while (coarseChain->GetSamples()->size() <= sampleID) {
        coarseChain->AddNumSamps(1);
        coarseChain->Run();
      }

      return coarseChain->GetSamples()->at(sampleID);
    }

    double SubsamplingMIProposal::LogDensity(std::shared_ptr<SamplingState> const& currState,
                                             std::shared_ptr<SamplingState> const& propState) {
      return 0;
    }


  }
}
