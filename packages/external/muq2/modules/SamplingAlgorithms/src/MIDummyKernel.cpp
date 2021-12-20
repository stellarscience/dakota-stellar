#include "MUQ/SamplingAlgorithms/MIDummyKernel.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/Utilities/AnyHelpers.h"

#include <iomanip>

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

MIDummyKernel::MIDummyKernel(boost::property_tree::ptree const& pt,
                   std::shared_ptr<AbstractSamplingProblem> problem,
                   std::shared_ptr<MCMCProposal> proposal,
                   std::shared_ptr<MCMCProposal> coarse_proposal,
                   std::shared_ptr<MIInterpolation> proposalInterpolation,
                   std::shared_ptr<SingleChainMCMC> coarse_chain)
  : TransitionKernel(pt, problem),
    proposal(proposal),
    coarse_proposal(coarse_proposal),
    proposalInterpolation(proposalInterpolation),
    coarse_chain(coarse_chain)
  {}


MIDummyKernel::~MIDummyKernel() {}

void MIDummyKernel::PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state){
  proposal->Adapt(t,state);
}

std::vector<std::shared_ptr<SamplingState>> MIDummyKernel::Step(unsigned int const t, std::shared_ptr<SamplingState> prevState){

  assert(proposal);

  numCalls++;

  // If no coarse sample is specified, assume it's the first one
  std::shared_ptr<SamplingState> coarsePrevState;
  if(prevState->HasMeta("coarseSample")) {
    coarsePrevState = AnyCast(prevState->meta["coarseSample"]);
  } else {
    coarsePrevState = coarse_chain->GetSamples()->at(0);
  }

  // New fine proposal
  std::shared_ptr<SamplingState> fineProp = proposal->Sample(prevState);
  std::shared_ptr<SamplingState> coarseProp = coarse_proposal->Sample(coarsePrevState);
  std::shared_ptr<SamplingState> prop = proposalInterpolation->Interpolate (coarseProp, fineProp);

  // Let's retrieve density and QOI as usual, but just not do anything with it
  problem->LogDensity(prop);
  if (problem->numBlocksQOI > 0) {
    prop->meta["QOI"] = problem->QOI();
  }
  return std::vector<std::shared_ptr<SamplingState>>(1, prop);
}

void MIDummyKernel::PrintStatus(const std::string prefix) const
{
  std::stringstream msg;
  msg << prefix << "MI dummy kernel was called "  << numCalls << " times";

  std::cout << msg.str() << std::endl;
}
