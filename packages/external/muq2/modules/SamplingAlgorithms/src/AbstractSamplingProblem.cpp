#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"

#include <iostream>

using namespace muq;
using namespace SamplingAlgorithms;

AbstractSamplingProblem::AbstractSamplingProblem(Eigen::VectorXi const& blockSizesIn,
                                                 Eigen::VectorXi const& blockSizesQOIIn) :
                         numBlocks(blockSizesIn.size()),
                         blockSizes(blockSizesIn),
                         numBlocksQOI(blockSizesQOIIn.size()),
                         blockSizesQOI(blockSizesQOIIn)
{
  assert(blockSizes.size()==numBlocks);
  assert(blockSizesQOI.size()==numBlocksQOI);
}

AbstractSamplingProblem::AbstractSamplingProblem(Eigen::VectorXi const& blockSizesIn) :
                         AbstractSamplingProblem(blockSizesIn, Eigen::VectorXi::Zero(0))
{

}

std::shared_ptr<SamplingState> AbstractSamplingProblem::QOI() {
  return nullptr;
}

Eigen::VectorXd AbstractSamplingProblem::GradLogDensity(std::shared_ptr<SamplingState> const& state,
                                                        unsigned                       const  blockWrt)
{
    std::cerr << "ERROR: AbstractSamplingProblem::GradLogDensity is not yet implemented!" << std::endl;
    assert(false);
}
