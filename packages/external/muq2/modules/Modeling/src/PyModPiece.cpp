#include "MUQ/Modeling/PyModPiece.h"

using namespace muq::Modeling;

PyModPiece::PyModPiece(Eigen::VectorXi const& inputSizes,
                       Eigen::VectorXi const& outputSizes)
                       : ModPiece(inputSizes, outputSizes) {}

void PyModPiece::EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) {
  this->EvaluateImpl(ToStdVec(input));
}

void PyModPiece::GradientImpl(unsigned int                const  outputDimWrt,
                              unsigned int                const  inputDimWrt,
                              ref_vector<Eigen::VectorXd> const& input,
                              Eigen::VectorXd             const& sensitivity) {
  this->GradientImpl(outputDimWrt, inputDimWrt, ToStdVec(input), sensitivity);
}

void PyModPiece::GradientImpl(unsigned int                 const  outputDimWrt,
                          unsigned int                 const  inputDimWrt,
                          std::vector<Eigen::VectorXd> const& input,
                          Eigen::VectorXd              const& sensitivity)
{
    gradient = GradientByFD(outputDimWrt, inputDimWrt, input,sensitivity);
}

void PyModPiece::JacobianImpl(unsigned int                const  outputDimWrt,
                              unsigned int                const  inputDimWrt,
                              ref_vector<Eigen::VectorXd> const& input) {
  this->JacobianImpl(outputDimWrt, inputDimWrt, ToStdVec(input));
}

void PyModPiece::JacobianImpl(unsigned int                 const  outputDimWrt,
                          unsigned int                 const  inputDimWrt,
                          std::vector<Eigen::VectorXd> const& input)
{
  jacobian = JacobianByFD(outputDimWrt, inputDimWrt, input);
}


void PyModPiece::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                   unsigned int                const  inputDimWrt,
                                   ref_vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd             const& vec) {
  this->ApplyJacobianImpl(outputDimWrt, inputDimWrt, ToStdVec(input), vec);
}

void PyModPiece::ApplyJacobianImpl(unsigned int                 const  outputDimWrt,
                                   unsigned int                 const  inputDimWrt,
                                   std::vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd              const& vec)
{
  jacobianAction = ApplyJacobianByFD(outputDimWrt, inputDimWrt, input, vec);
}

void PyModPiece::ApplyHessianImpl(unsigned int                 const  outputDimWrt,
                                  unsigned int                 const  inputDimWrt1,
                                  unsigned int                 const  inputDimWrt2,
                                  ref_vector<Eigen::VectorXd> const& input,
                                  Eigen::VectorXd              const& sens,
                                  Eigen::VectorXd              const& vec)
{
  this->ApplyHessianImpl(outputDimWrt, inputDimWrt1, inputDimWrt2, ToStdVec(input), sens, vec);
}

void PyModPiece::ApplyHessianImpl(unsigned int                 const  outputDimWrt,
                                   unsigned int                 const  inputDimWrt1,
                                   unsigned int                 const  inputDimWrt2,
                                   std::vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd              const& sens,
                                   Eigen::VectorXd              const& vec)
{
  hessAction = ApplyHessianByFD(outputDimWrt, inputDimWrt1, inputDimWrt2, input, sens, vec);
}
