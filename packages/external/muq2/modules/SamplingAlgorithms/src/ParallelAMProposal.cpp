#include "MUQ/SamplingAlgorithms/ParallelAMProposal.h"
#include "MUQ/Utilities/Cereal/BoostAnySerializer.h"

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;

REGISTER_MCMC_PROPOSAL(ParallelAMProposal)

ParallelAMProposal::ParallelAMProposal(boost::property_tree::ptree                  pt,
                                       std::shared_ptr<AbstractSamplingProblem>     problem) : ParallelAMProposal(pt, problem, std::make_shared<parcer::Communicator>()){}

ParallelAMProposal::ParallelAMProposal(pt::ptree                                    pt ,
                                       std::shared_ptr<AbstractSamplingProblem>     problem,
                                       std::shared_ptr<parcer::Communicator> const& newcomm) : AMProposal(pt, problem)
{
  SetCommunicator(newcomm);
}

void ParallelAMProposal::Adapt(unsigned int const t, std::vector<std::shared_ptr<SamplingState> > const& states) {
  assert(comm);

  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
    // get the samples from the other processors (or send your own to them)
    std::vector<std::shared_ptr<SamplingState> > otherStates;
    if( i==comm->GetRank() ) { otherStates = states; }
    comm->Bcast(otherStates, i);

    // adapt
    totSamps += otherStates.size();
    AMProposal::Adapt(totSamps, otherStates);
  }
}
