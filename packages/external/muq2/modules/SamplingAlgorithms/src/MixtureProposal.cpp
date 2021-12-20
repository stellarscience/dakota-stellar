#include "MUQ/SamplingAlgorithms/MixtureProposal.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/StringUtilities.h"
#include "MUQ/Utilities/RandomGenerator.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

REGISTER_MCMC_PROPOSAL(MixtureProposal)


MixtureProposal::MixtureProposal(boost::property_tree::ptree                     pt,
                                 std::shared_ptr<AbstractSamplingProblem> const& prob) : MixtureProposal(pt,
                                                                                                         prob,
                                                                                                         GetProposals(pt,prob),
                                                                                                         GetWeights(pt))
{
}

MixtureProposal::MixtureProposal(boost::property_tree::ptree                       pt,
                                 std::shared_ptr<AbstractSamplingProblem>   const& prob,
                                 std::vector<std::shared_ptr<MCMCProposal>> const& proposalsIn,
                                 std::vector<double>                        const& weightsIn) : MCMCProposal(pt,prob),
                                                                                                proposals(proposalsIn),
                                                                                                weights(weightsIn)
{
  if(weights.size()==0)
    weights.resize(proposals.size(), 1.0);

  assert(proposals.size()==weights.size());

  // Make sure the weights are positive and normalize them
  double wtSum = 0.0;
  for(int i=0; i<weights.size(); ++i){
    assert(weights.at(i)>0.0);
    wtSum += weights.at(i);
  }
  for(int i=0; i<weights.size();++i)
    weights.at(i) /= wtSum;
}


std::vector<std::shared_ptr<MCMCProposal>> MixtureProposal::GetProposals(boost::property_tree::ptree              const& pt,
                                                                         std::shared_ptr<AbstractSamplingProblem> const& prob)
{
  std::vector<std::string> propStrings = StringUtilities::Split(pt.get<std::string>("Components"));
  assert(propStrings.size()>0);

  std::vector<std::shared_ptr<MCMCProposal>> props(propStrings.size());

  for(int i=0; i<props.size(); ++i){
    boost::property_tree::ptree subTree = pt.get_child(propStrings.at(i));
    subTree.put("BlockIndex", pt.get("BlockIndex",0));
    props.at(i) = MCMCProposal::Construct(subTree, prob);
    assert(props.at(i));
  }

  return props;
}

std::vector<double> MixtureProposal::GetWeights(boost::property_tree::ptree const& pt)
{
  std::vector<std::string> wtStrings = StringUtilities::Split(pt.get("Weights",""));

  if(wtStrings.size()==0){
    return std::vector<double>();

  }else{

    std::vector<double> weights(wtStrings.size());
    double wtSum = 0.0;

    for(int i=0; i<weights.size(); ++i)
      weights.at(i) = std::stof(wtStrings.at(i));

    return weights;
  }
}


std::shared_ptr<SamplingState> MixtureProposal::Sample(std::shared_ptr<SamplingState> const& currentState)
{

  // Pick a component at random
  int randInd = RandomGenerator::GetDiscrete(weights);

  // Return a sample of the random index
  return proposals.at(randInd)->Sample(currentState);
}

double MixtureProposal::LogDensity(std::shared_ptr<SamplingState> const& currState,
                                   std::shared_ptr<SamplingState> const& propState)
{
  double density = 0.0;
  for(int i=0; i<proposals.size(); ++i)
    density += weights.at(i) * std::exp(proposals.at(i)->LogDensity(currState, propState));

  return std::log(density);
}
