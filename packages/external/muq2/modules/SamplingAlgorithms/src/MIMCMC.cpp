#include "MUQ/SamplingAlgorithms/MIMCMC.h"

namespace muq {
  namespace SamplingAlgorithms {

    MIMCMC::MIMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory)
    : SamplingAlgorithm(std::shared_ptr<SampleCollection>(), std::shared_ptr<SampleCollection>()),
      pt(pt),
      componentFactory(componentFactory)
    {
      gridIndices = MultiIndexFactory::CreateFullTensor(componentFactory->FinestIndex()->GetVector());

      for (int i = 0; i < gridIndices->Size(); i++) {
        std::shared_ptr<MultiIndex> boxHighestIndex = (*gridIndices)[i];
        auto box = std::make_shared<MIMCMCBox>(componentFactory, boxHighestIndex);
        boxes.push_back(box);
      }
    }

    std::shared_ptr<MIMCMCBox> MIMCMC::GetBox(std::shared_ptr<MultiIndex> index) {
      for (std::shared_ptr<MIMCMCBox> box : boxes) {
        if (*(box->GetHighestIndex()) == *index)
          return box;
      }
      return nullptr;
    }

    std::shared_ptr<SampleCollection> MIMCMC::GetSamples() const {
      return nullptr;
    }
    std::shared_ptr<SampleCollection> MIMCMC::GetQOIs() const {
      return nullptr;
    }

    std::shared_ptr<SampleCollection> MIMCMC::RunImpl(std::vector<Eigen::VectorXd> const& x0) {
      for (auto box : boxes) {
        assert(box);
        int numSamples = pt.get<int>("NumSamples" + multiindexToConfigString(box->GetHighestIndex()));
        for (int samp = 0; samp < numSamples; samp++) {
          box->Sample();
        }
      }

      return nullptr;
    }

    Eigen::VectorXd MIMCMC::MeanQOI() {
      // Compute full QOI estimate
      Eigen::VectorXd MImean(boxes[0]->GetFinestProblem()->blockSizesQOI.sum());
      MImean.setZero();

      for (auto box : boxes) {
        Eigen::VectorXd sampMean = box->MeanQOI();

        MImean += sampMean;
      }

      return MImean;
    }

    Eigen::VectorXd MIMCMC::MeanParam() {
      Eigen::VectorXd MImean(boxes[0]->GetFinestProblem()->blockSizes.sum());
      MImean.setZero();

      for (auto box : boxes) {
        Eigen::VectorXd sampMean = box->MeanParam();

        MImean += sampMean;
      }

      return MImean;
    }

    std::shared_ptr<MIMCMCBox> MIMCMC::GetMIMCMCBox(std::shared_ptr<MultiIndex> index) {
      for (auto box : boxes) {
        if (box->GetHighestIndex() == index)
          return box;
      }
      return nullptr;
    }

    void MIMCMC::WriteToFile(std::string filename) {
      for (auto box : boxes) {
        box->WriteToFile(filename);
      }
    }


    std::shared_ptr<MultiIndexSet> MIMCMC::GetIndices() {
      return gridIndices;
    }


    std::string MIMCMC::multiindexToConfigString (std::shared_ptr<MultiIndex> index) {
      std::stringstream strs;
      for (int i = 0; i < index->GetLength(); i++) {
        strs << "_" << index->GetValue(i);
      }
      return strs.str();
    }

    void MIMCMC::Draw(bool drawSamples) {
      std::ofstream graphfile;
      graphfile.open ("graph");
      graphfile << "digraph {" << std::endl;
      graphfile << "nodesep=1.2;" << std::endl;
      graphfile << "splines=false;" << std::endl;
      for (auto box : boxes) {
        box->Draw(graphfile, drawSamples);
      }
      graphfile << "}" << std::endl;
      graphfile.close();
    }

  }
}
