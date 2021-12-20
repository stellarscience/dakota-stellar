#include "MUQ/Modeling/LinearAlgebra/GaussNewtonOperator.h"

#include "MUQ/Modeling/Distributions/Density.h"

using namespace muq::Modeling;


GaussNewtonOperator::GaussNewtonOperator(std::shared_ptr<ModPiece>     const& forwardModelIn,
                                         std::shared_ptr<ModPiece>     const& noiseModelIn,
                                         std::vector<Eigen::VectorXd>  const& inputsIn,
                                         unsigned int                         inWrtIn,
                                         double                               scaleIn,
                                         double                               nuggetIn) : LinearOperator(forwardModelIn->inputSizes(inWrtIn), forwardModelIn->inputSizes(inWrtIn)),
                                                                                           forwardModel(forwardModelIn),
                                                                                           noiseModel(noiseModelIn),
                                                                                           inputs(inputsIn),
                                                                                           noiseInputs(forwardModelIn->Evaluate(inputsIn)),
                                                                                           inWrt(inWrtIn),
                                                                                           scale(scaleIn),
                                                                                           nugget(nuggetIn)
{
  assert(noiseModelIn->inputSizes.size()==1);
  assert(noiseModelIn->outputSizes.size()==1);
  assert(noiseModelIn->outputSizes(0)==1);

  assert(forwardModel->outputSizes.size()==1);
  assert(forwardModel->outputSizes(0)==noiseModelIn->inputSizes(0));
  assert(forwardModel->inputSizes.size()>inWrtIn);
}


Eigen::MatrixXd GaussNewtonOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  Eigen::MatrixXd output(rows(),x.cols());

  for(unsigned int i=0; i<x.cols(); ++i){
    Eigen::VectorXd temp = forwardModel->ApplyJacobian(0,inWrt,inputs,x.col(i).eval());
    temp = noiseModel->ApplyHessian(0, 0, 0, noiseInputs, Eigen::VectorXd::Ones(1), temp);
    output.col(i) = forwardModel->Gradient(0,inWrt,inputs,temp);
    output.col(i) *= scale;
    output.col(i) += nugget*x.col(i);
  }
  return output;
}

Eigen::MatrixXd GaussNewtonOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  return Apply(x);
}
