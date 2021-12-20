#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/RandomGenerator.h"

#include <iomanip>

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(MHKernel)

MHKernel::MHKernel(pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem) : TransitionKernel(pt, problem)
{

  // Extract the proposal parts from the ptree
  std::string proposalName = pt.get<std::string>("Proposal");

  boost::property_tree::ptree subTree = pt.get_child(proposalName);
  subTree.put("BlockIndex", blockInd);

  // Construct the proposal
  proposal = MCMCProposal::Construct(subTree, problem);
  assert(proposal);
}

MHKernel::MHKernel(pt::ptree const& pt,
                   std::shared_ptr<AbstractSamplingProblem> problem,
                   std::shared_ptr<MCMCProposal> proposalIn) : TransitionKernel(pt, problem),
                                                               proposal(proposalIn) {}


void MHKernel::PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state){
  proposal->Adapt(t,state);
}

std::vector<std::shared_ptr<SamplingState>> MHKernel::Step(unsigned int const t, std::shared_ptr<SamplingState> prevState){

  assert(proposal);

  // propose a new point
  std::shared_ptr<SamplingState> prop = proposal->Sample(prevState);

  // The following metadata is needed by the expensive sampling problem
  if(prevState->HasMeta("iteration"))
    prop->meta["iteration"] = prevState->meta["iteration"];
  prop->meta["IsProposal"] = true;

  // compute acceptance probability
  double propTarget;
  double currentTarget;

  if( prevState->HasMeta("LogTarget") && (prevState->HasMeta("QOI") || problem->numBlocksQOI == 0) && !reeval ){
    currentTarget = AnyCast( prevState->meta["LogTarget"]);
  }else{
    currentTarget = problem->LogDensity(prevState);
    if (problem->numBlocksQOI > 0) {
      prevState->meta["QOI"] = problem->QOI();
    }
    prevState->meta["LogTarget"] = currentTarget;
  }

  try{
    propTarget = problem->LogDensity(prop);

  }catch(std::runtime_error &e){

    // If we couldn't compute the log density, reject the proposal
    return std::vector<std::shared_ptr<SamplingState>>(1, prevState);
  }

  prop->meta["LogTarget"] = propTarget;

  // Aceptance probability
  const double forwardPropDens = proposal->LogDensity(prevState, prop);
  const double backPropDens = proposal->LogDensity(prop, prevState);
  const double alpha = std::exp(propTarget - forwardPropDens - currentTarget + backPropDens);

  // accept/reject
  numCalls++;
  if( RandomGenerator::GetUniform()<alpha ) {
    if (problem->numBlocksQOI > 0) {
      prop->meta["QOI"] = problem->QOI();
    }
    numAccepts++;

    prop->meta["IsProposal"] = false;
    return std::vector<std::shared_ptr<SamplingState>>(1, prop);
  } else {
    return std::vector<std::shared_ptr<SamplingState>>(1, prevState);
  }
}

void MHKernel::PrintStatus(const std::string prefix) const
{
  std::stringstream msg;
  msg << std::setprecision(2);
  msg << prefix << "MHKernel acceptance Rate = "  << 100.0*double(numAccepts)/double(numCalls) << "%";

  std::cout << msg.str() << std::endl;
}
