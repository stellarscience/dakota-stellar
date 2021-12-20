#include "MUQ/Modeling/LinearAlgebra/HessianOperator.h"

using namespace muq::Modeling;



HessianOperator::HessianOperator(std::shared_ptr<ModPiece>    const& pieceIn,
                                 std::vector<Eigen::VectorXd> const& inputsIn,
                                 unsigned int                        outWrtIn,
                                 unsigned int                        inWrt1In,
                                 unsigned int                        inWrt2In,
                                 Eigen::VectorXd              const& sensIn,
                                 double                              scaleIn,
                                 double                              nuggetIn) : LinearOperator(pieceIn->inputSizes(inWrt1In), pieceIn->inputSizes(inWrt2In)),
                                                                                 basePiece(pieceIn),
                                                                                 inputs(inputsIn),
                                                                                 outWrt(outWrtIn),
                                                                                 inWrt1(inWrt1In),
                                                                                 inWrt2(inWrt2In),
                                                                                 sens(sensIn),
                                                                                 scale(scaleIn),
                                                                                 nugget(nuggetIn)
{
  assert(basePiece);
  assert(basePiece->inputSizes.size()>inWrt1In);
  assert(basePiece->inputSizes.size()>inWrt2In);
  assert(basePiece->outputSizes.size()>outWrt);
  assert(sens.size()==basePiece->outputSizes(outWrt));
  assert(nugget>=0.0);
}

Eigen::MatrixXd HessianOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  Eigen::MatrixXd output(basePiece->inputSizes(inWrt2),x.cols());
  for(unsigned int i=0; i<x.cols(); ++i)
    output.col(i) = nugget*x.col(i) + scale*basePiece->ApplyHessian(outWrt,inWrt1,inWrt2, inputs, sens, x.col(i));
  return output;
}

Eigen::MatrixXd HessianOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  assert(inWrt1==inWrt2);
  return Apply(x);// Assume symmetry based on Scharwz theorem
}
