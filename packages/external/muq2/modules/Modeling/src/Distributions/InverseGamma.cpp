#include "MUQ/Modeling/Distributions/InverseGamma.h"
#include "MUQ/Utilities/RandomGenerator.h"
#include <Eigen/Core>

using namespace muq::Modeling;
using namespace muq::Utilities;

InverseGamma::InverseGamma(Eigen::VectorXd const& alphaIn,
                           Eigen::VectorXd const& betaIn) : Distribution(alphaIn.size()),
                                                  alpha(alphaIn),
                                                  beta(betaIn),
                                                  logConst(ComputeConstant(alphaIn, betaIn))
                                                  {};

InverseGamma::InverseGamma(double       alphaIn,
                           double       betaIn) : InverseGamma(alphaIn*Eigen::VectorXd::Ones(1),
                                                               betaIn*Eigen::VectorXd::Ones(1)){};

double InverseGamma::ComputeConstant(Eigen::VectorXd const& alphaIn,
                                     Eigen::VectorXd const& betaIn)
{
    double logConst = 0;
    for(int i=0; i<alphaIn.size(); ++i)
      logConst += alphaIn(i)*std::log(betaIn(i)) - std::lgamma(alphaIn(i));

    return logConst;
}

double InverseGamma::LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  Eigen::VectorXd const& x = inputs.at(0).get();

  if(x.minCoeff()<std::numeric_limits<double>::epsilon())
    return -1.0*std::numeric_limits<double>::infinity();

  return logConst + ((-alpha.array()-1.0)*x.array().log() - beta.array() / x.array()).sum();
}


Eigen::VectorXd InverseGamma::SampleImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  Eigen::VectorXd output(alpha.size());
  for(int i=0; i<alpha.size(); ++i)
    output(i) = 1.0/RandomGenerator::GetGamma(alpha(i),1.0/beta(i));

  return output;
}
