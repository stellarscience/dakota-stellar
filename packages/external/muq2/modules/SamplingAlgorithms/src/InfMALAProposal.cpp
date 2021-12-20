#include "MUQ/SamplingAlgorithms/InfMALAProposal.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Utilities/AnyHelpers.h"

using namespace muq::SamplingAlgorithms;
using namespace muq::Modeling;
using namespace muq::Utilities;

REGISTER_MCMC_PROPOSAL(InfMALAProposal)

InfMALAProposal::InfMALAProposal(boost::property_tree::ptree       const& pt,
                                 std::shared_ptr<AbstractSamplingProblem> prob,
                                 std::shared_ptr<GaussianBase>            priorIn) : CrankNicolsonProposal(pt,prob,priorIn),
                                                                                     stepSize(pt.get("StepSize",1.0))
{
}

InfMALAProposal::InfMALAProposal(boost::property_tree::ptree       const& pt,
                                 std::shared_ptr<AbstractSamplingProblem> prob) : CrankNicolsonProposal(pt,prob),
                                                                                  stepSize(pt.get("StepSize",1.0))
{
}


std::shared_ptr<SamplingState> InfMALAProposal::Sample(std::shared_ptr<SamplingState> const& currentState)
{
  // the mean of the proposal is the current point
  std::vector<Eigen::VectorXd> props = currentState->state;

  std::vector<Eigen::VectorXd> hypers = GetPriorInputs(currentState->state);

   // vector holding the covariance of the proposal times the gradient
  Eigen::VectorXd sigmaGrad = GetSigmaGrad(currentState);

  Eigen::VectorXd priorSamp = priorDist->Sample(hypers);

  props.at(blockInd) = priorDist->GetMean() + sqrt(1.0-beta*beta)*(currentState->state.at(blockInd)-priorDist->GetMean()) + beta*(priorSamp - priorDist->GetMean() + stepSize*sigmaGrad);

  // store the new state in the output
  return std::make_shared<SamplingState>(props, 1.0);
}

double InfMALAProposal::LogDensity(std::shared_ptr<SamplingState> const& currState,
                                   std::shared_ptr<SamplingState> const& propState)
{
  std::vector<Eigen::VectorXd> hypers = GetPriorInputs(currState->state);
  if(hypers.size()>0)
    priorDist->ResetHyperparameters(WorkPiece::ToRefVector(hypers));

  Eigen::VectorXd sigmaGrad = GetSigmaGrad(currState);
  Eigen::VectorXd mean = priorDist->GetMean() + sqrt(1.0-beta*beta)*(currState->state.at(blockInd)-priorDist->GetMean()) + beta*stepSize*sigmaGrad;
  Eigen::VectorXd diff = (propState->state.at(blockInd) - mean)/beta;

  hypers.insert(hypers.begin(), (diff + priorDist->GetMean()).eval());
  return priorDist->LogDensity(hypers);
}

Eigen::VectorXd InfMALAProposal::GetSigmaGrad(std::shared_ptr<SamplingState> const& state) const
{
  std::stringstream blockId;
  blockId << "_" << blockInd;

  if(!state->HasMeta("MALA_SigmaGrad" + blockId.str())){
    Eigen::VectorXd grad; // vector holding the gradient of the log density

    if(!state->HasMeta("GradLogDensity" + blockId.str()))
      state->meta["GradLogDensity" + blockId.str()] = prob->GradLogDensity(state, blockInd);

    grad = AnyCast( state->meta["GradLogDensity" + blockId.str()] );

    state->meta["MALA_SigmaGrad" + blockId.str()] = priorDist->ApplyCovariance(grad).col(0).eval();
  }

  Eigen::VectorXd sigmaGrad = AnyCast(state->meta["MALA_SigmaGrad" + blockId.str()]);
  return sigmaGrad;
}
