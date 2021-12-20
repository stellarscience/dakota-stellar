#include "MUQ/Utilities/AnyHelpers.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/SamplingAlgorithms/AMProposal.h"
#include "MUQ/Utilities/RandomGenerator.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;
using namespace muq::Modeling;


REGISTER_MCMC_PROPOSAL(AMProposal)

AMProposal::AMProposal(pt::ptree                                       pt,
                       std::shared_ptr<AbstractSamplingProblem> const& prob,
                       Eigen::MatrixXd                          const& initialCov) : MCMCProposal(pt,prob),
                                                                                     propCov(initialCov),
                                                                                     adaptSteps(pt.get<unsigned int>("AdaptSteps")),
                                                                                     adaptStart(pt.get<unsigned int>("AdaptStart")),
                                                                                     adaptEnd(pt.get<unsigned int>("AdaptEnd",std::numeric_limits<unsigned int>::max())),
                                                                                     adaptScale(pt.get("AdaptScale",2.4*2.4/initialCov.rows()))
{
  propChol = propCov.selfadjointView<Eigen::Lower>().llt();
  newSamps.resize(propCov.rows(), adaptSteps);
  numAdaptSamps = 1;
}

AMProposal::AMProposal(pt::ptree                                       pt ,
                       std::shared_ptr<AbstractSamplingProblem> const& prob) : AMProposal(pt,prob, ConstructCovariance(pt,prob))

{}

Eigen::MatrixXd AMProposal::ConstructCovariance(pt::ptree                                const& pt ,
                                                std::shared_ptr<AbstractSamplingProblem> const& prob)
{
  int dim = prob->blockSizes(pt.get("BlockIndex",0));
  double propVar = pt.get("InitialVariance",1.0);
  return propVar * Eigen::MatrixXd::Identity(dim,dim);
}


void AMProposal::Adapt(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& states) {


  if((t>=adaptStart)&&(t<adaptEnd)){

    if((numAdaptSamps==1)&&(numNewSamps==0)){
      mean = states.at(0)->state.at(blockInd);
    }

    // Copy the states into the matrix of new samples
    for(unsigned int i=0; (i<states.size())&&(numNewSamps<adaptSteps); i++){
      newSamps.col(numNewSamps) = states.at(i)->state.at(blockInd);
      numNewSamps++;
    }

    if((t%adaptSteps==0)&&(numNewSamps>0)){

      Eigen::VectorXd oldMean;
      std::swap(mean,oldMean);

      // Update the mean, covariance, and Cholesky factorization
      mean = (oldMean*numAdaptSamps + newSamps.leftCols(numNewSamps).rowwise().mean()*numNewSamps)/(numAdaptSamps+numNewSamps);

      propChol.rankUpdate(oldMean, double(numAdaptSamps));
      propChol.rankUpdate(mean, -1.0*double(numAdaptSamps+numNewSamps));

      for(unsigned int i=0; i<numNewSamps; ++i)
        propChol.rankUpdate(newSamps.col(i), 1.0);

      numAdaptSamps += numNewSamps;
      numNewSamps = 0;
    }
  }
}


Eigen::MatrixXd AMProposal::ProposalCovariance() const
{
  return adaptScale*propChol.reconstructedMatrix()/double(numAdaptSamps);
}

std::shared_ptr<SamplingState> AMProposal::Sample(std::shared_ptr<SamplingState> const& currentState)
{

  // the mean of the proposal is the current point
  Eigen::VectorXd const& xc = currentState->state.at(blockInd);

  std::vector<Eigen::VectorXd> props = currentState->state;
  props.at(blockInd) = xc + std::sqrt(adaptScale/double(numAdaptSamps))*(propChol.matrixL() * RandomGenerator::GetNormal(xc.size())).eval();

  // store the new state in the output
  return std::make_shared<SamplingState>(props, 1.0);
}

double AMProposal::LogDensity(std::shared_ptr<SamplingState> const& currState,
                              std::shared_ptr<SamplingState> const& propState)
{
  Eigen::VectorXd diff = (propState->state.at(blockInd) - currState->state.at(blockInd))/std::sqrt(adaptScale/double(numAdaptSamps));
  return -0.5*diff.dot(propChol.solve(diff));
}
