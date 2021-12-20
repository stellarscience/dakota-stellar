#include "MUQ/SamplingAlgorithms/DummyKernel.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/Utilities/AnyHelpers.h"

#include <iomanip>

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

DummyKernel::DummyKernel(boost::property_tree::ptree const& pt,
                   std::shared_ptr<AbstractSamplingProblem> problem,
                   std::shared_ptr<MCMCProposal> proposal)
  : TransitionKernel(pt, problem),
    proposal(proposal)
  {}


DummyKernel::~DummyKernel() {}

void DummyKernel::PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state){
  proposal->Adapt(t,state);
}

std::vector<std::shared_ptr<SamplingState>> DummyKernel::Step(unsigned int const t, std::shared_ptr<SamplingState> prevState){

  assert(proposal);

  numCalls++;

  std::shared_ptr<SamplingState> prop = proposal->Sample(prevState);
  // Let's retrieve density and QOI as usual, but just not do anything with it
  problem->LogDensity(prop);
  if (problem->numBlocksQOI > 0) {
    prop->meta["QOI"] = problem->QOI();
  }
  return std::vector<std::shared_ptr<SamplingState>>(1, prop);
}

void DummyKernel::PrintStatus(const std::string prefix) const
{
  std::stringstream msg;
  msg << prefix << "Dummy kernel was called "  << numCalls << " times";

  std::cout << msg.str() << std::endl;
}
