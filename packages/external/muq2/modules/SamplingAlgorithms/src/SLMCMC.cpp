#include "MUQ/SamplingAlgorithms/SLMCMC.h"

namespace muq {
  namespace SamplingAlgorithms {

    SLMCMC::SLMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory, std::shared_ptr<MultiIndex> index)
     : SamplingAlgorithm(std::shared_ptr<SampleCollection>(), std::shared_ptr<SampleCollection>()),
       componentFactory(componentFactory)
    {
      auto finestIndex = componentFactory->FinestIndex(); 
      
      assert(index->GetLength() == finestIndex->GetLength());
      assert(*index <= *(componentFactory->FinestIndex()));

      pt::ptree ptBlockID;
      ptBlockID.put("BlockIndex",0);
      
      auto problem = componentFactory->SamplingProblem(index);
      auto proposal = componentFactory->Proposal(index, problem);

      std::vector<std::shared_ptr<TransitionKernel>> kernels(1);
      kernels[0] = std::make_shared<MHKernel>(ptBlockID,problem,proposal);
      
      Eigen::VectorXd startingPoint = componentFactory->StartingPoint(index);

      single_chain = std::make_shared<SingleChainMCMC>(pt,kernels);
      single_chain->SetState(startingPoint);
    }
    
    SLMCMC::SLMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory)
     : SLMCMC(pt,componentFactory, componentFactory->FinestIndex()) { }

    std::shared_ptr<SampleCollection> SLMCMC::GetSamples() const {
      return nullptr;
    }
    std::shared_ptr<SampleCollection> SLMCMC::GetQOIs() const {
      return nullptr;
    }

    std::shared_ptr<SampleCollection> SLMCMC::RunImpl(std::vector<Eigen::VectorXd> const& x0) {
      return single_chain->Run();
    }

    Eigen::VectorXd SLMCMC::MeanQOI() {
        return single_chain->GetQOIs()->Mean();
    }
    
    Eigen::VectorXd SLMCMC::MeanParameter() {
        auto samps = single_chain->GetSamples();
        return samps->Mean();
    }
    
    void SLMCMC::WriteToFile(std::string filename){
        auto samps = single_chain->GetSamples();
        auto QOI = single_chain->GetQOIs();
        if(QOI != nullptr)
          QOI->WriteToFile(filename,"/qois");
        samps->WriteToFile(filename,"/samples");
    }
    
  }
}
