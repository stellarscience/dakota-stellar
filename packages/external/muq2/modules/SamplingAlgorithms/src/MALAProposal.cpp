#include "MUQ/SamplingAlgorithms/MALAProposal.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/Utilities/AnyHelpers.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

REGISTER_MCMC_PROPOSAL(MALAProposal)

MALAProposal::MALAProposal(pt::ptree                                      pt,
                          std::shared_ptr<AbstractSamplingProblem> const& probIn) :
              MCMCProposal(pt,probIn)
{

  unsigned int problemDim = prob->blockSizes(blockInd);

  stepSize = pt.get("StepSize", 1.0);
  assert(stepSize>0);

  // compute the (diagonal) covariance for the proposal
  const Eigen::VectorXd cov = Eigen::VectorXd::Ones(problemDim);

  // created a Gaussian with scaled identity covariance
  zDist = std::make_shared<Gaussian>(Eigen::VectorXd::Zero(problemDim), cov);
}

MALAProposal::MALAProposal(pt::ptree                                       pt,
                           std::shared_ptr<AbstractSamplingProblem> const& probIn,
                           std::shared_ptr<GaussianBase>                   zDistIn) :

            MCMCProposal(pt,probIn),
            zDist(zDistIn)
{
  stepSize = pt.get("StepSize", 1.0);
  assert(stepSize>0);
}

std::shared_ptr<SamplingState> MALAProposal::Sample(std::shared_ptr<SamplingState> const& currentState) {
  assert(currentState->state.size()>blockInd);

  // the mean of the proposal is the current point
  std::vector<Eigen::VectorXd> props = currentState->state;
  assert(props.size()>blockInd);
  Eigen::VectorXd const& xc = currentState->state.at(blockInd);

  // Get the gradient from the previous step
  Eigen::VectorXd sigmaGrad; // vector holding the covariance of the proposal times the gradient

  std::stringstream blockId;
  blockId << "_" << blockInd;

  if(!currentState->HasMeta("MALA_SigmaGrad" + blockId.str())){
    Eigen::VectorXd grad; // vector holding the gradient of the log density

    if(!currentState->HasMeta("GradLogDensity" + blockId.str()))
      currentState->meta["GradLogDensity" + blockId.str()] = prob->GradLogDensity(currentState, blockInd);

    grad = AnyCast( currentState->meta["GradLogDensity" + blockId.str()] );

    currentState->meta["MALA_SigmaGrad" + blockId.str()] = zDist->ApplyCovariance(grad).col(0).eval();
  }

  sigmaGrad = AnyCast(currentState->meta["MALA_SigmaGrad" + blockId.str()]);

  // Draw a sample
  Eigen::VectorXd prop = stepSize*sigmaGrad + std::sqrt(2.0*stepSize)*zDist->Sample();
  props.at(blockInd) = xc + prop;

  // store the new state in the output
  return std::make_shared<SamplingState>(props, 1.0);
}

double MALAProposal::LogDensity(std::shared_ptr<SamplingState> const& currState,
                                std::shared_ptr<SamplingState> const& propState) {

  // Get the gradient from the previous step
  Eigen::VectorXd grad; // vector holding the gradient of the log density
  Eigen::VectorXd sigmaGrad; // vector holding the covariance of the proposal times the gradient

  std::stringstream blockId;
  blockId << "_" << blockInd;

  if(!currState->HasMeta("GradLogDensity" + blockId.str())){
    try{
      currState->meta["GradLogDensity" + blockId.str()] = prob->GradLogDensity(currState, blockInd);
    }catch(std::runtime_error &e){
      currState->meta["GradLogDensity" + blockId.str()] = Eigen::VectorXd::Zero(currState->state.at(blockInd).size()).eval();
    }
  }

  grad = AnyCast( currState->meta["GradLogDensity" + blockId.str()] );

  if(!currState->HasMeta("MALA_SigmaGrad" + blockId.str()))
    currState->meta["MALA_SigmaGrad" + blockId.str()] = zDist->ApplyCovariance(grad).col(0).eval();

  sigmaGrad = AnyCast(currState->meta["MALA_SigmaGrad" + blockId.str()]);

  Eigen::VectorXd z = (propState->state.at(blockInd) - currState->state.at(blockInd) - stepSize*sigmaGrad) / std::sqrt(2.0*stepSize);

  return zDist->LogDensity(z);
}
