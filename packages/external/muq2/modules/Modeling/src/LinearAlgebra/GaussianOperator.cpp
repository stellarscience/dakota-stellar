#include "MUQ/Modeling/LinearAlgebra/GaussianOperator.h"

using namespace muq::Modeling;


GaussianOperator::GaussianOperator(std::shared_ptr<GaussianBase> const& gaussIn,
                                     Gaussian::Mode                       precOrCovIn) : LinearOperator(gaussIn->Dimension(), gaussIn->Dimension()),
                                                                                         gauss(gaussIn),
                                                                                         precOrCov(precOrCovIn)
{
}


Eigen::MatrixXd GaussianOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  if(precOrCov==Gaussian::Precision){
    return gauss->ApplyPrecision(x);
  }else{
    return gauss->ApplyCovariance(x);
  }
}


Eigen::MatrixXd GaussianOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  return Apply(x);
}
