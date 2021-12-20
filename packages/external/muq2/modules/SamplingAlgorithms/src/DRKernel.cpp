#include "MUQ/SamplingAlgorithms/DRKernel.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/StringUtilities.h"

#include <iomanip>

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(DRKernel)


DRKernel::DRKernel(pt::ptree const& pt,
                   std::shared_ptr<AbstractSamplingProblem> problem) : DRKernel(pt,
                                                                                problem,
                                                                                CreateProposals(pt, problem),
                                                                                CreateScales(pt))
{



}

DRKernel::DRKernel(pt::ptree                           const& pt,
                   std::shared_ptr<AbstractSamplingProblem>   problem,
                   std::vector<std::shared_ptr<MCMCProposal>> proposalsIn,
                   std::vector<double>                        scales) : TransitionKernel(pt, problem),
                                                                        proposals(proposalsIn),
                                                                        propScales(scales)
{
  // Update the unique set of proposals
  for(auto& proposal : proposals)
    uniqueProps.insert(proposal);

  numProposalCalls = Eigen::VectorXi::Zero(proposals.size());
  numProposalAccepts = Eigen::VectorXi::Zero(proposals.size());

  for(unsigned int i=0; i<propScales.size(); ++i){
    if(std::abs(propScales.at(i)-1.0)>std::numeric_limits<double>::epsilon()){
      isScaled = true;
    }
  }

}

void DRKernel::PostStep(unsigned int const t,
                        std::vector<std::shared_ptr<SamplingState>> const& state){

  for(auto& proposal : uniqueProps)
    proposal->Adapt(t,state);
}


std::shared_ptr<SamplingState> DRKernel::SampleProposal(unsigned int                          stage,
                                                        std::shared_ptr<SamplingState> const& state) const
{
  std::shared_ptr<SamplingState> prop = proposals.at(stage)->Sample(state);
  if(isScaled){
    prop->state.at(blockInd) = (prop->state.at(blockInd) - state->state.at(blockInd))*propScales.at(stage) + state->state.at(blockInd);
  }
  return prop;
}

double DRKernel::EvaluateProposal(unsigned int                          stage,
                                  std::shared_ptr<SamplingState> const& x,
                                  std::shared_ptr<SamplingState> const& y) const
{

  if(isScaled){
    const double dim = x->state.at(blockInd).size();

    std::shared_ptr<SamplingState> yCopy = std::make_shared<SamplingState>(*y);
    yCopy->state.at(blockInd) = x->state.at(blockInd) + (yCopy->state.at(blockInd) - x->state.at(blockInd))/propScales.at(stage);
    return proposals.at(stage)->LogDensity(x, yCopy) - dim*std::log(propScales.at(stage));
  }else{
    return proposals.at(stage)->LogDensity(x,y);
  }
}



std::vector<std::shared_ptr<SamplingState>> DRKernel::Step(unsigned int const t,
                                                           std::shared_ptr<SamplingState> prevState){

  const unsigned int numStages = proposals.size();

  std::vector<double> logDensities; // vector previous density evaluations
  logDensities.reserve(numStages);
  std::vector<std::shared_ptr<SamplingState>> proposedPoints; // vector of previously proposed point
  proposedPoints.reserve(numStages);


  // If the previous state does not have a density value or it needs to be recomputed, evaluate the density
  double currentTarget;
  if( prevState->HasMeta("LogTarget") && !reeval ){
    currentTarget = AnyCast( prevState->meta["LogTarget"]);
  }else{
    currentTarget = problem->LogDensity(prevState);
    if (problem->numBlocksQOI > 0) {
      prevState->meta["QOI"] = problem->QOI();
    }
    prevState->meta["LogTarget"] = currentTarget;
  }

  logDensities.push_back(  currentTarget );
  proposedPoints.push_back( prevState );

  std::shared_ptr<SamplingState> proposedState;
  Eigen::VectorXd alphas = Eigen::VectorXd::Zero(numStages);

  // loop over delayed rejection stages, will break after acceptance
  for (int stage = 0; stage < numStages; ++stage) {

    // build and store proposal at next step:
    numProposalCalls(stage)++;

    proposedState = SampleProposal(stage, prevState);

    // The following metadata is needed by the expensive sampling problem
    if(prevState->HasMeta("iteration"))
      proposedState->meta["iteration"] = proposedState->meta["iteration"];
    proposedState->meta["IsProposal"] = true;

    // Evaluate the density at the proposal
    double propTarget = problem->LogDensity(proposedState);
    proposedState->meta["LogTarget"] = propTarget;

    // add proposed point to lists
    logDensities.push_back( propTarget );
    proposedPoints.push_back( proposedState );

    // call recursive routine to evaluate acceptance probability
    alphas(stage) = Alpha(logDensities, proposedPoints);


    // If accepted, return the accepted point
    if (RandomGenerator::GetUniform() < alphas(stage)) {

      if (problem->numBlocksQOI > 0) {
        proposedState->meta["QOI"] = problem->QOI();
      }

      numProposalAccepts(stage)++;
      proposedState->meta["IsProposal"] = false;

      return std::vector<std::shared_ptr<SamplingState>>(1, proposedState);
    }

  } // end of loop over stages

  return std::vector<std::shared_ptr<SamplingState>>(1,prevState); //could only be here if we've rejected at every stage
}

std::vector<double> DRKernel::CreateScales(boost::property_tree::ptree const& pt)
{
  std::string proposalList = pt.get<std::string>("Proposal");
  std::vector<std::string> proposalNames = StringUtilities::Split(proposalList);

  std::vector<double> output;

  // If the proposals are manually set, then set all the scales to one
  if(proposalNames.size()>1){
    output.resize(proposalNames.size(), 1.0);

  }else{

    int numStages = pt.get<int>("NumStages");
    assert(numStages>0);

    std::string scaleType = pt.get<std::string>("ScaleFunction","Power");
    double scale = pt.get("Scale", 2.0);

    output.resize(numStages);

    if(scaleType=="Power"){
      for(int i=0; i<numStages; ++i)
        output.at(i) = scale / std::pow(2.0, double(i));

    } else if(scaleType=="Linear") {
      for(int i=0; i<numStages; ++i)
        output.at(i) = scale / (double(i)+1.0);

    }else{
      std::string msg = "ERROR: In DRKernel::CreateScales, invalid scale function \"" + scaleType + "\" specified in options.   Valid options are \"Power\" or \"Linear\"";
      throw std::invalid_argument(msg);
    }
  }

  return output;

}

void DRKernel::PrintStatus(const std::string prefix) const
{
  std::stringstream msg;
  msg << std::setprecision(2);
  msg << prefix << "DR: Number of calls = " << numProposalCalls.transpose() << "\n";
  msg << prefix << "DR: Cumulative Accept Rate = ";
  double numCalls = static_cast<double>(numProposalCalls(0));
  double rate     = 100.0 * static_cast<double>(numProposalAccepts(0));
  msg << std::setw(4) << std::fixed << std::setprecision(1) << rate / numCalls << "%";

  for (int i = 1; i < numProposalAccepts.size(); ++i) {
    rate += 100.0 * static_cast<double>(numProposalAccepts(i));
    msg << ", " << std::setw(4) << std::fixed << std::setprecision(1) << rate / numCalls << "%";
  }

  std::cout << msg.str() << std::endl;
}

std::vector<std::shared_ptr<MCMCProposal>> DRKernel::CreateProposals(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> const& problem)
{
  // Extract the proposal parts from the ptree
  std::string proposalList = pt.get<std::string>("Proposal");
  std::vector<std::string> proposalNames = StringUtilities::Split(proposalList);
  std::vector<std::shared_ptr<MCMCProposal>> proposals;

  // Generate all of the proposals
  for(auto& proposalName : proposalNames){
    boost::property_tree::ptree subTree = pt.get_child(proposalName);
    subTree.put("BlockIndex", pt.get("BlockIndex",0));

    // Construct the proposal
    auto proposal = MCMCProposal::Construct(subTree, problem);
    assert(proposal);

    proposals.push_back(proposal);
  }

  if(proposals.size()==1){
    int numStages = pt.get("NumStages",-1);
    if(numStages>1){
      auto prop = proposals.at(0);
      proposals.resize(numStages, prop);
    }
  }

  return proposals;
}
