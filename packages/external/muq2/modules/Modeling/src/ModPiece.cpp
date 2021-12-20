#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Utilities/Exceptions.h"
#include "MUQ/Utilities/Demangler.h"

#include <chrono>

using namespace muq::Modeling;

ModPiece::ModPiece(Eigen::VectorXi const& inputSizesIn,
                   Eigen::VectorXi const& outputSizesIn) :
  WorkPiece(std::vector<std::string>(inputSizesIn.size(), typeid(Eigen::VectorXd).name()),
            std::vector<std::string>(outputSizesIn.size(), typeid(Eigen::VectorXd).name())),
  inputSizes(inputSizesIn),
  outputSizes(outputSizesIn){};


std::vector<Eigen::VectorXd> const& ModPiece::Evaluate(std::vector<Eigen::VectorXd> const& input)
{
  return Evaluate(ToRefVector(input));
}

bool ModPiece::ExistsInCache(ref_vector<Eigen::VectorXd> const& input) const
{
  // Check to see if the input is the same as what's cached
  if(input.size()!=cacheInput.size()){
    return false;
  }

  // Check the size of each input vector
  for(int i=0; i<input.size(); ++i){
    if(input.at(i).get().size()!=cacheInput.at(i).size()){
      return false;
    }
  }

  // Check the contents of each input vector
  for(int i=0; i<input.size(); ++i){
    Eigen::VectorXd const& inVec = input.at(i).get();
    for(int j=0; j<cacheInput[i].size(); ++j){
      if(std::abs(inVec(j)-cacheInput[i](j))>std::numeric_limits<double>::epsilon()){
        return false;
      }
    }
  }

  return true;
}

std::vector<Eigen::VectorXd> const& ModPiece::Evaluate(ref_vector<Eigen::VectorXd> const& input)
{
  CheckInputs(input,"Evaluate");

  // If we're using the one-step cache, check to see if the inputs are the same as the previous evaluation
  if(cacheEnabled){
    if(ExistsInCache(input)){
      return outputs;

    }else{
      // Copy the contents
      cacheInput.resize(input.size());
      for(int i=0; i<input.size(); ++i)
        cacheInput.at(i) = input.at(i);
    }
  }

  // Otherwise, evaluate the model
  numEvalCalls++;
  auto start_time = std::chrono::high_resolution_clock::now();

  EvaluateImpl(input);

  auto end_time = std::chrono::high_resolution_clock::now();
  evalTime += static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());

  return outputs;
}

Eigen::VectorXd const& ModPiece::Gradient(unsigned int                 const  outputDimWrt,
                                          unsigned int                 const  inputDimWrt,
                                          std::vector<Eigen::VectorXd> const& input,
                                          Eigen::VectorXd              const& sensitivity)
{
  return Gradient(outputDimWrt, inputDimWrt, ToRefVector(input), sensitivity);
}


Eigen::VectorXd const& ModPiece::Gradient(unsigned int                const  outputDimWrt,
                                          unsigned int                const  inputDimWrt,
                                          ref_vector<Eigen::VectorXd> const& input,
                                          Eigen::VectorXd             const& sensitivity)
{
  CheckInputs(input,"Gradient");

  numGradCalls++;
  auto start_time = std::chrono::high_resolution_clock::now();

  GradientImpl(outputDimWrt, inputDimWrt, input, sensitivity);

  auto end_time = std::chrono::high_resolution_clock::now();
  gradTime += static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());

  return gradient;
}


Eigen::VectorXd ModPiece::GradientByFD(unsigned int                 const  outputDimWrt,
                                       unsigned int                 const  inputDimWrt,
                                       std::vector<Eigen::VectorXd> const& input,
                                       Eigen::VectorXd              const& sensitivity)
{
  return GradientByFD(outputDimWrt, inputDimWrt, ToRefVector(input), sensitivity);
}

Eigen::VectorXd ModPiece::GradientByFD(unsigned int                const  outputDimWrt,
                                       unsigned int                const  inputDimWrt,
                                       ref_vector<Eigen::VectorXd> const& input,
                                       Eigen::VectorXd             const& sensitivity)
{
  numGradFDCalls++;
  return JacobianByFD(outputDimWrt,inputDimWrt, input).transpose()*sensitivity;
}

Eigen::MatrixXd const& ModPiece::Jacobian(unsigned int                 const  outputDimWrt,
                                          unsigned int                 const  inputDimWrt,
                                          std::vector<Eigen::VectorXd> const& input)
{
  return Jacobian(outputDimWrt, inputDimWrt, ToRefVector(input));
}

Eigen::MatrixXd const& ModPiece::Jacobian(unsigned int                const  outputDimWrt,
                                          unsigned int                const  inputDimWrt,
                                          ref_vector<Eigen::VectorXd> const& input)
{
  CheckInputs(input,"Jacobian");

  numJacCalls++;
  auto start_time = std::chrono::high_resolution_clock::now();

  JacobianImpl(outputDimWrt, inputDimWrt, input);

  auto end_time = std::chrono::high_resolution_clock::now();
  jacTime += static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());

  return jacobian;
}


Eigen::VectorXd const& ModPiece::ApplyJacobian(unsigned int                 const  outputDimWrt,
                                               unsigned int                 const  inputDimWrt,
                                               std::vector<Eigen::VectorXd> const& input,
                                               Eigen::VectorXd              const& vec)
{
  return ApplyJacobian(outputDimWrt, inputDimWrt, ToRefVector(input), vec);
}

Eigen::VectorXd const& ModPiece::ApplyJacobian(unsigned int                const  outputDimWrt,
                                               unsigned int                const  inputDimWrt,
                                               ref_vector<Eigen::VectorXd> const& input,
                                               Eigen::VectorXd             const& vec)
{
  CheckInputs(input,"ApplyJacobian");

  numJacActCalls++;
  auto start_time = std::chrono::high_resolution_clock::now();

  ApplyJacobianImpl(outputDimWrt, inputDimWrt, input, vec);

  auto end_time = std::chrono::high_resolution_clock::now();
  jacActTime += static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());

  return jacobianAction;
}

void ModPiece::EvaluateImpl(ref_vector<boost::any> const& inputs){
  ref_vector<Eigen::VectorXd> eigenInputs;
  for(int i=0; i<inputs.size(); ++i) {
    eigenInputs.push_back( std::cref(boost::any_cast<Eigen::VectorXd const&>(inputs.at(i))) );
  }

  EvaluateImpl(eigenInputs);

  WorkPiece::outputs.resize(outputs.size());
  for(int i=0; i<outputs.size(); ++i)
    WorkPiece::outputs.at(i) = boost::any(outputs.at(i));
}


// Eigen::VectorXd ModPiece::ApplyHessian(int                         const  outputDimWrt,
//                                        int                         const  inputDimWrt1,
//                                        int                         const  inputDimWrt2,
//                                        ref_vector<Eigen::VectorXd> const& input,
//                                        Eigen::VectorXd             const& sensitivity,
//                                        Eigen::VectorXd             const& vec)
// {
//  CheckInputs(input);
//
//  numHessCalls++;
//  auto start_time = std::chrono::high_resolution_clock::now();
//
//  ApplyHessianImpl(outputDimWrt, inputDimWrt1, inputDimWrt2, input, sensitivity, vec);
//
//  auto end_time = std::chrono::high_resolution_clock::now();
//  hessTime += 1e6*static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());
// }

void ModPiece::CheckInputs(ref_vector<Eigen::VectorXd> const& input, std::string const& funcName)
{
  bool errorOccured = false;
  std::string msg;

  if(input.size() != inputSizes.size()){
    msg += "  - Wrong number of input arguments.  Expected " + std::to_string(inputSizes.size()) + " inputs, but " + std::to_string(input.size()) + " were given.\n";
    errorOccured=true;
  }

  for(int i=0; i<std::min<int>(input.size(), inputSizes.size()); ++i){
    if(input.at(i).get().size() != inputSizes(i)){
      msg += "  - Input " + std::to_string(i) + " has the wrong size.  Expected size " + std::to_string(inputSizes(i)) + ", but given input with size " + std::to_string(input.at(i).get().size()) + "\n";
      errorOccured = true;
    }
  }

  if(errorOccured){
    std::string className = muq::Utilities::demangle(typeid(*this).name());
    msg = "\nError evaluating " + className + "::" + funcName + ":\n" + msg;

    throw muq::WrongSizeError(msg);
  }
}


void ModPiece::GradientImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input,
                            Eigen::VectorXd             const& sensitivity)
{
  gradient = Jacobian(outputDimWrt, inputDimWrt, input).transpose() * sensitivity;
}

void ModPiece::JacobianImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input)
{
  jacobian = JacobianByFD(outputDimWrt, inputDimWrt, input);
}

void ModPiece::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                 unsigned int                const  inputDimWrt,
                                 ref_vector<Eigen::VectorXd> const& input,
                                 Eigen::VectorXd             const& vec)
{
  jacobianAction = ApplyJacobianByFD(outputDimWrt, inputDimWrt, input, vec);
}

// void ModPiece::ApplyHessianImpl(int                         const  outputDimWrt,
//                                 int                         const  inputDimWrt1,
//                                 int                         const  inputDimWrt2,
//                                 ref_vector<Eigen::VectorXd> const& input,
//                                 Eigen::VectorXd             const& sensitivity
//                                 Eigen::VectorXd             const& vec)
// {
//   hessianAction = ApplyHessianByFD(outputDimWrt, inputDimWrt1, inputDimWrt2, input, sensitvity, vec);
// }

Eigen::MatrixXd ModPiece::JacobianByFD(unsigned int                 const  outputDimWrt,
                                       unsigned int                 const  inputDimWrt,
                                       std::vector<Eigen::VectorXd> const& input)
{
  return JacobianByFD(outputDimWrt, inputDimWrt, ToRefVector(input));
}

Eigen::MatrixXd ModPiece::JacobianByFD(unsigned int                const  outputDimWrt,
                                       unsigned int                const  inputDimWrt,
                                       ref_vector<Eigen::VectorXd> const& input)
{
  numJacFDCalls++;

  Eigen::VectorXd f0 = Evaluate(input).at(outputDimWrt);
  Eigen::VectorXd f;

  double eps;
  Eigen::VectorXd newInput(input.at(inputDimWrt).get());
  ref_vector<Eigen::VectorXd> newInputVec = input;

  Eigen::MatrixXd jac(outputSizes(outputDimWrt), inputSizes(inputDimWrt));
  for(int i=0; i<inputSizes(inputDimWrt); ++i){
    eps  = std::max(1e-8, 1e-10*std::abs(input.at(inputDimWrt)(i)));

    newInput(i) = input.at(inputDimWrt)(i) + eps;
    newInputVec.at(inputDimWrt) = std::cref(newInput);

    f = Evaluate(newInputVec).at(outputDimWrt);

    jac.col(i) = (f-f0)/eps;

    newInput(i) = input.at(inputDimWrt)(i);
  }

  return jac;
}

Eigen::VectorXd ModPiece::ApplyJacobianByFD(unsigned int                 const  outputDimWrt,
                                            unsigned int                 const  inputDimWrt,
                                            std::vector<Eigen::VectorXd> const& input,
                                            Eigen::VectorXd              const& vec)
{
  return ApplyJacobianByFD(outputDimWrt, inputDimWrt, ToRefVector(input), vec);
}
Eigen::VectorXd ModPiece::ApplyJacobianByFD(unsigned int                const  outputDimWrt,
                                            unsigned int                const  inputDimWrt,
                                            ref_vector<Eigen::VectorXd> const& input,
                                            Eigen::VectorXd             const& vec)
{
  numJacActFDCalls++;

  const double eps = 1e-4;
  double vecNorm = vec.norm();
  const Eigen::VectorXd stepDir = vec/vecNorm;

  ref_vector<Eigen::VectorXd> newInputVec = input;
  Eigen::VectorXd newInput = input.at(inputDimWrt).get() - 0.5*eps*stepDir;
  newInputVec.at(inputDimWrt) = std::cref(newInput);
  Eigen::VectorXd f0 = Evaluate(newInputVec).at(outputDimWrt);

  newInput = input.at(inputDimWrt).get() + 0.5*eps*stepDir;
  newInputVec.at(inputDimWrt) = std::cref(newInput);

  Eigen::VectorXd f  = Evaluate(newInputVec).at(outputDimWrt);

  return vecNorm*(f-f0)/eps;
}

Eigen::VectorXd ModPiece::ApplyHessian(unsigned int                const  outWrt,
                                       unsigned int                const  inWrt1,
                                       unsigned int                const  inWrt2,
                                       std::vector<Eigen::VectorXd> const& input,
                                       Eigen::VectorXd             const& sens,
                                       Eigen::VectorXd             const& vec)
{
  return ApplyHessian(outWrt,inWrt1, inWrt2, ToRefVector(input), sens, vec);
}

Eigen::VectorXd ModPiece::ApplyHessian(unsigned int                const  outWrt,
                                       unsigned int                const  inWrt1,
                                       unsigned int                const  inWrt2,
                                       ref_vector<Eigen::VectorXd> const& input,
                                       Eigen::VectorXd             const& sens,
                                       Eigen::VectorXd             const& vec)
{
  assert(inWrt2<inputSizes.size()+1);
  assert(outWrt<sens.size());
  assert(outputSizes(outWrt)==sens.size());
  if(inWrt2<inputSizes.size()){
    assert(inputSizes(inWrt2)==vec.size());
  }else{
    assert(vec.size()==outputSizes(outWrt));
  }

  numHessActCalls++;
  auto start_time = std::chrono::high_resolution_clock::now();

  ApplyHessianImpl(outWrt, inWrt1, inWrt2, input, sens, vec);

  auto end_time = std::chrono::high_resolution_clock::now();
  hessActTime += static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());

  return hessAction;
}

void ModPiece::ApplyHessianImpl(unsigned int                const  outWrt,
                                unsigned int                const  inWrt1,
                                unsigned int                const  inWrt2,
                                ref_vector<Eigen::VectorXd> const& input,
                                Eigen::VectorXd             const& sens,
                                Eigen::VectorXd             const& vec)
{
  //std::cout << "Just before FD..." << std::endl;
  hessAction = ApplyHessianByFD(outWrt, inWrt1, inWrt2, input, sens, vec);
  //std::cout << "After FD, hessAction = " << hessAction << std::endl;
}

Eigen::VectorXd ModPiece::ApplyHessianByFD(unsigned int                const  outWrt,
                                       unsigned int                const  inWrt1,
                                       unsigned int                const  inWrt2,
                                       std::vector<Eigen::VectorXd> const& input,
                                       Eigen::VectorXd             const& sens,
                                       Eigen::VectorXd             const& vec)
{
  return ApplyHessianByFD(outWrt,inWrt1, inWrt2, ToRefVector(input), sens, vec);
}

Eigen::VectorXd ModPiece::ApplyHessianByFD(unsigned int                const  outWrt,
                                           unsigned int                const  inWrt1,
                                           unsigned int                const  inWrt2,
                                           ref_vector<Eigen::VectorXd> const& input,
                                           Eigen::VectorXd             const& sens,
                                           Eigen::VectorXd             const& vec)
{
  numHessActFDCalls++;

  const double stepSize = std::max(1e-4, 1e-8*vec.norm());

  Eigen::VectorXd grad1, grad2;

  // If the Hessian is wrt to one of the inputs, not the sensitivity vector
  if(inWrt2<inputSizes.size()){

    ref_vector<Eigen::VectorXd> input2 = input;
    Eigen::VectorXd x2 = input.at(inWrt2).get() - 0.5*stepSize * vec;
    input2.at(inWrt2) = std::cref(x2);

    grad1 = Gradient(outWrt, inWrt1, input2, sens);

    x2 = input.at(inWrt2).get() + 0.5*stepSize * vec;
    input2.at(inWrt2) = std::cref(x2);
    grad2 = Gradient(outWrt, inWrt1, input2, sens);

  // Otherwise, we want the Jacobian of the Gradient piece wrt to the sensitivity vector
  }else{
    Eigen::VectorXd sens2 = sens - 0.5*stepSize * vec;
    grad1 = Gradient(outWrt, inWrt1, input, sens2);

    sens2 = sens + 0.5*stepSize * vec;
    grad2 = Gradient(outWrt, inWrt1, input, sens2);
  }

  return (grad2 - grad1)/stepSize;
}


double ModPiece::GetRunTime(const std::string& method) const
{
  const double toMilli = 1e-6;

  if (method.compare("Evaluate") == 0) {
    return (numEvalCalls == 0) ? -1.0 : toMilli *static_cast<double>(evalTime) / static_cast<double>(numEvalCalls);
  } else if (method.compare("Gradient") == 0) {
    return (numGradCalls == 0) ? -1.0 : toMilli *static_cast<double>(gradTime) / static_cast<double>(numGradCalls);
  } else if (method.compare("Jacobian") == 0) {
    return (numJacCalls == 0) ? -1.0 : toMilli *static_cast<double>(jacTime) / static_cast<double>(numJacCalls);
  } else if (method.compare("JacobianAction") == 0) {
    return (numJacActCalls ==
            0) ? -1.0 : toMilli *static_cast<double>(jacActTime) / static_cast<double>(numJacActCalls);
  } else if (method.compare("HessianAction") == 0) {
    return (numHessActCalls == 0) ? -1.0 : toMilli *static_cast<double>(hessActTime) / static_cast<double>(numHessActCalls);
  } else {
    assert(method.compare("Evaluate") == 0 || method.compare("Gradient") == 0 || method.compare(
             "Jacobian") == 0 || method.compare("JacobianAction") == 0 || method.compare("HessianAction") == 0);
    return -999.0;
  }
}


void ModPiece::ResetCallTime()
{
  numEvalCalls    = 0;
  numGradCalls    = 0;
  numJacCalls     = 0;
  numJacActCalls  = 0;
  numHessActCalls = 0;

  evalTime    = 0;
  gradTime    = 0;
  jacTime     = 0;
  jacActTime  = 0;
  hessActTime = 0;
}

unsigned long int ModPiece::GetNumCalls(const std::string& method) const
{
  if (method.compare("Evaluate") == 0) {
    return numEvalCalls;
  } else if (method.compare("Gradient") == 0) {
    return numGradCalls;
  } else if (method.compare("Jacobian") == 0) {
    return numJacCalls;
  } else if (method.compare("JacobianAction") == 0) {
    return numJacActCalls;
  } else if (method.compare("HessianAction") == 0) {
    return numHessActCalls;
  } else if (method.compare("GradientFD") == 0) {
    return numGradFDCalls;
  } else if (method.compare("JacobianFD") == 0) {
    return numJacFDCalls;
  } else if (method.compare("JacobianActionFD") == 0) {
    return numJacActFDCalls;
  } else if (method.compare("HessianActionFD") == 0) {
    return numHessActFDCalls;
  } else {
    assert( (method.compare("Evaluate") == 0) || (method.compare("Gradient") == 0)
          || (method.compare("Jacobian") == 0) || (method.compare("JacobianAction") == 0)
          || (method.compare("HessianAction") == 0) || (method.compare("GradientFD")==0)
          || (method.compare("JacobianFD")==0) || (method.compare("JacobianActionFD")==0)
          || (method.compare("HessianActionFD")==0));
    return -999;
  }
}
// Eigen::VectorXd ModPiece::ApplyHessianByFD(int                         const  outputDimWrt,
//                                            int                         const  inputDimWrt1,
//                                            int                         const  inputDimWrt2,
//                                            ref_vector<Eigen::VectorXd> const& input,
//                                            Eigen::VectorXd             const& sensitivity
//                                            Eigen::VectorXd             const& vec)
// {
//   Eigen::MatrixXd hess(inputSizes(inputDimWrt1), inputSizes(inputDimWrt2));
//
//   const double eps = std::max(1e-8, 1e-10*vec.norm());
//
//   Eigen::VectorXd newInput(input.at(inputDimWrt2));
//   ref_vector<Eigen::VectorXd> newInputVec = input;
//
//   Eigen::VectorXd newInput = input.at(inputDimWrt2) + eps*vec;
//
//   Eigen::VectorXd f0 = Gradient(outputDimWrt, inputDimWrt1, input, sensitvity);
//   Eigen::VectorXd f  = Gradient(outputDimWrt, inputDimWrt1, newInput, sensitvity);
//
//   return (f-f0)/eps;
// }

std::vector<Eigen::VectorXd> ModPiece::ToStdVec(ref_vector<Eigen::VectorXd> const& input) {
  std::vector<Eigen::VectorXd> newIns(input.size());

  for (int i=0; i<input.size(); ++i)
    newIns.at(i) = input.at(i).get();

  return newIns;
}
