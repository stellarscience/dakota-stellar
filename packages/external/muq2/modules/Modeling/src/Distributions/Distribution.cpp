#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/RandomVariable.h"

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;

std::shared_ptr<Density> Distribution::AsDensity()
{
  return std::make_shared<Density>(shared_from_this());
}

std::shared_ptr<RandomVariable> Distribution::AsVariable()
{
  return std::make_shared<RandomVariable>(shared_from_this());
}

ref_vector<const Eigen::VectorXd> Distribution::ToRefVector(std::vector<Eigen::VectorXd> const& vec) const {

  ref_vector<const Eigen::VectorXd> refs;
  refs.reserve(vec.size());

  // populate the input vector
  for(int i=0; i<vec.size(); ++i)
    refs.push_back(std::cref(vec.at(i)));

  return refs;
}

double Distribution::LogDensity(ref_vector<Eigen::VectorXd> const& inputs) {
  return LogDensityImpl(inputs);
}


Eigen::VectorXd Distribution::Sample(ref_vector<Eigen::VectorXd> const& inputs) {
  return SampleImpl(inputs);
}

Eigen::VectorXd Distribution::Sample() {
  return Sample(ref_vector<Eigen::VectorXd>());
}

Eigen::VectorXd Distribution::GradLogDensity(unsigned int wrt, ref_vector<Eigen::VectorXd> const& inputs)
{
  assert(wrt<inputs.size());
  return GradLogDensityImpl(wrt, inputs);
}

Eigen::VectorXd Distribution::GradLogDensityImpl(unsigned int                       wrt,
                                                 ref_vector<Eigen::VectorXd> const& inputs)
{
  // Default to finite difference
  ref_vector<Eigen::VectorXd> newInputs = inputs;
  Eigen::VectorXd newIn = newInputs.at(wrt).get();
  newInputs.at(wrt) = std::cref(newIn);

  const double f0 = LogDensity(newInputs);

  const int dim = inputs.at(wrt).get().size();

  Eigen::VectorXd output(dim);
  for(int i=0; i<dim; ++i){

    double eps = std::max(1e-8, 1e-10*std::abs(inputs.at(wrt).get()(i)));
    newIn(i) += eps;

    double newF = LogDensity(newInputs);
    output(i) = (newF-f0)/eps;
    newIn(i) = inputs.at(wrt).get()(i);
  }

  return output;
}

Eigen::VectorXd Distribution::ApplyLogDensityHessian(unsigned int                const  inWrt1,
                                                     unsigned int                const  inWrt2,
                                                     ref_vector<Eigen::VectorXd> const& input,
                                                     Eigen::VectorXd             const& vec)
{
  assert(inWrt1<hyperSizes.size()+1);
  assert(inWrt2<hyperSizes.size()+1);
  assert(input.size() == hyperSizes.size()+1);
  if(inWrt2==0){
    assert(vec.size()==varSize);
  }else{
    assert(vec.size()==hyperSizes(inWrt2-1));
  }

  return ApplyLogDensityHessianImpl(inWrt1,inWrt2,input,vec);
}

Eigen::VectorXd Distribution::ApplyLogDensityHessianImpl(unsigned int                const  inWrt1,
                                                         unsigned int                const  inWrt2,
                                                         ref_vector<Eigen::VectorXd> const& input,
                                                         Eigen::VectorXd             const& vec)
{
  const double stepSize = 1e-8 / vec.norm();
  Eigen::VectorXd grad1 = GradLogDensity(inWrt1, input);
  Eigen::VectorXd grad2;

  ref_vector<Eigen::VectorXd> input2 = input;
  Eigen::VectorXd x2 = input.at(inWrt2).get() + stepSize * vec;
  input2.at(inWrt2) = std::cref(x2);
  grad2 = GradLogDensity(inWrt1, input2);

  return (grad2 - grad1)/stepSize;
}
