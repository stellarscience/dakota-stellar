#include "MUQ/Modeling/JacobianPiece.h"

using namespace muq::Modeling;


JacobianPiece::JacobianPiece(std::shared_ptr<ModPiece> const& basePieceIn,
                             unsigned int              const  outWrtIn,
                             unsigned int              const  inWrtIn) : ModPiece(GetInputSizes(basePieceIn, inWrtIn),
                                                                                  GetOutputSizes(basePieceIn, outWrtIn)),
                                                                         basePiece(basePieceIn),
                                                                         outWrt(outWrtIn),
                                                                         inWrt(inWrtIn)
{
}

Eigen::VectorXi JacobianPiece::GetInputSizes(std::shared_ptr<ModPiece> const& basePiece,
                                             unsigned int              const inWrt)
{
  assert(inWrt < basePiece->inputSizes.size());

  unsigned int numInputs = basePiece->inputSizes.size() + 1;
  Eigen::VectorXi inputSizes(numInputs);
  inputSizes.head(numInputs-1) = basePiece->inputSizes;
  inputSizes(numInputs-1) = basePiece->inputSizes(inWrt);

  return inputSizes;
}

Eigen::VectorXi JacobianPiece::GetOutputSizes(std::shared_ptr<ModPiece> const& basePiece,
                                              unsigned int              const outWrt)
{
  Eigen::VectorXi outputSizes(1);
  outputSizes(0) = basePiece->outputSizes(outWrt);
  return outputSizes;
}


void JacobianPiece::EvaluateImpl(ref_vector<Eigen::VectorXd> const& input)
{
  ref_vector<Eigen::VectorXd> baseInputs(input.begin(), input.end()-1);

  outputs.resize(1);
  outputs.at(0) = basePiece->ApplyJacobian(outWrt, inWrt, baseInputs, input.at(input.size()-1).get());
}
