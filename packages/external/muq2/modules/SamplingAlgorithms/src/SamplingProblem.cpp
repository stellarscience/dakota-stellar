#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SamplingState.h"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

SamplingProblem::SamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& targetIn) : AbstractSamplingProblem(targetIn->inputSizes),
                                                                                             target(targetIn),
                                                                                             qoi(nullptr) {}

SamplingProblem::SamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& targetIn,
                                 std::shared_ptr<muq::Modeling::ModPiece> const& qoiIn) : AbstractSamplingProblem(targetIn->inputSizes,qoiIn->outputSizes),
                                                                                          target(targetIn),
                                                                                          qoi(qoiIn) {}


double SamplingProblem::LogDensity(std::shared_ptr<SamplingState> const& state) {
  assert(target);

  lastState = state;

  return target->Evaluate(state->state).at(0)(0);
}

std::shared_ptr<SamplingState> SamplingProblem::QOI() {
  assert(target);
  assert(lastState);
  if (qoi == nullptr)
    return nullptr;

  return std::make_shared<SamplingState>(qoi->Evaluate(lastState->state));
}

Eigen::VectorXd SamplingProblem::GradLogDensity(std::shared_ptr<SamplingState> const& state,
                                                unsigned                       const  blockWrt)
{
  return target->Gradient(0,blockWrt, state->state, Eigen::VectorXd::Ones(1).eval());
}

std::vector<int> SamplingProblem::GetBlockSizes(std::shared_ptr<ModPiece> const& target)
{
  int numBlocks = target->inputSizes.size();

  std::vector<int> output(numBlocks);
  for(int i=0; i<numBlocks; ++i)
    output.at(i) = target->inputSizes(i);

  return output;
}
