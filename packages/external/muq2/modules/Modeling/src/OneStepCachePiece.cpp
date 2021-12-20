#include "MUQ/Modeling/OneStepCachePiece.h"

using namespace muq::Modeling;

OneStepCachePiece::OneStepCachePiece(std::shared_ptr<ModPiece> baseModPiece, const double& prec)
  : ModPiece(baseModPiece->inputSizes, baseModPiece->outputSizes),
    baseModPiece(baseModPiece),
    prec(prec)
{}

void OneStepCachePiece::EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) {

  // Find out whether we have a cache hit, i.e. input equals last call
  bool cacheHit = true;

  if (firstEvaluation) {
    firstEvaluation = false;
    cacheHit = false;
    lastInput.resize(input.size());
  } else {
    for (int i = 0; i < input.size(); i++) {
      if (!input[i].get().isApprox(lastInput[i], prec)) {
        cacheHit = false;
        break;
      }
    }
  }

  // Return last output if hit; otherwise, return new evaluation
  if (cacheHit) {
    hits++;
    outputs = lastOutputs;
    return;
  }

  misses++;

  outputs = baseModPiece->Evaluate(input);
  lastOutputs = outputs;
  for (int i = 0; i < input.size(); i++) {
    lastInput[i] = input[i].get();
  }
}

double OneStepCachePiece::HitRatio() {
  return (double)hits / (hits + misses);
}