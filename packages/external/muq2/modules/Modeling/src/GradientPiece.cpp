#include "MUQ/Modeling/GradientPiece.h"

using namespace muq::Modeling;




GradientPiece::GradientPiece(std::shared_ptr<ModPiece> const& basePieceIn,
                             unsigned int              const  outWrtIn,
                             unsigned int              const  inWrtIn) : ModPiece(GetInputSizes(basePieceIn,outWrtIn),
                                                                                  GetOutputSizes(basePieceIn, inWrtIn)),
                                                                        basePiece(basePieceIn),
                                                                        outWrt(outWrtIn),
                                                                        inWrt(inWrtIn)
{
}


void GradientPiece::EvaluateImpl(ref_vector<Eigen::VectorXd> const& input)
{
  ref_vector<Eigen::VectorXd> baseInputs(input.begin(), input.end()-1);

  outputs.resize(1);
  outputs.at(0) = basePiece->Gradient(outWrt, inWrt, baseInputs, input.at(input.size()-1).get());
}

void GradientPiece::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                      unsigned int                const  inputDimWrt,
                                      ref_vector<Eigen::VectorXd> const& input,
                                      Eigen::VectorXd             const& vec)
{
  ref_vector<Eigen::VectorXd> baseInputs(input.begin(), input.end()-1);
  jacobianAction = basePiece->ApplyHessian(outputDimWrt, inWrt, inputDimWrt, baseInputs, input.at(input.size()-1).get(), vec);
}


Eigen::VectorXi GradientPiece::GetInputSizes(std::shared_ptr<ModPiece> const& basePiece,
                                             unsigned int              const outWrt)
{
  unsigned int numInputs = basePiece->inputSizes.size()+1;

  Eigen::VectorXi inputSizes(numInputs);
  inputSizes.head(numInputs-1) = basePiece->inputSizes;
  inputSizes(numInputs-1) = basePiece->outputSizes(outWrt); // extra input for adjoint variable

  return inputSizes;
}

Eigen::VectorXi GradientPiece::GetOutputSizes(std::shared_ptr<ModPiece> const& basePiece,
                                              unsigned int              const inWrt)
{
  Eigen::VectorXi outputSizes(1);
  outputSizes(0) = basePiece->inputSizes(inWrt);
  return outputSizes;
}
