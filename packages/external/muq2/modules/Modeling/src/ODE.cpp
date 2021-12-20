#include "MUQ/Modeling/ODE.h"
#include "MUQ/Modeling/ODEData.h"

#include <cvodes/cvodes.h>               /* access to CVODE                 */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sunnonlinsol/sunnonlinsol_newton.h>     /* access to the newton SUNNonlinearSolver      */
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h> /* access to the fixed point SUNNonlinearSolver */

#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>
#include <sunlinsol/sunlinsol_pcg.h>

//#include <sundials/sundials_types.h> // definition of type
//#include <sundials/sundials_math.h>  // contains the macros ABS, SQR, and EXP

using namespace muq::Modeling;



ODE::ODE(std::shared_ptr<ModPiece>   const& rhsIn,
                   Eigen::VectorXd             const& evalTimesIn,
                   bool                               isAuto,
                   boost::property_tree::ptree        opts) : ModPiece(GetInputSizes(rhsIn, isAuto),
                                                                       GetOutputSizes(rhsIn, evalTimesIn)),
                                                              isAutonomous(isAuto),
                                                              startTime(GetStartTime(evalTimesIn, opts)),
                                                              evalTimes(evalTimesIn),
                                                              rhs(rhsIn)
{
  // Set the nonlinear solver options if possible
  if(opts.count("NonlinearSolver")>0)
    nonlinOpts.SetOptions(opts.get_child("NonlinearSolver"));

  // Set the linear solver options if possible
  if(opts.count("LinearSolver")>0)
    linOpts.SetOptions(opts.get_child("LinearSolver"));

  // Set the integrator options if possible
  if(opts.count("Integrator")>0)
    intOpts.SetOptions(opts.get_child("Integrator"));

}

double ODE::GetStartTime(Eigen::VectorXd             const& evalTimesIn,
                              boost::property_tree::ptree const& opts)
{
  // Make sure the evalTimes are in ascending order
  for(int i=0; i<evalTimesIn.size()-1; ++i){
    if(evalTimesIn(i+1)<=evalTimesIn(i)){
      std::stringstream msg;
      msg << "ERROR: Evaluation times are not strictly increasing: t_{" << i+1 << "}-t_{" << i << "}=" << evalTimesIn(i+1)-evalTimesIn(i);
      throw std::runtime_error(msg.str());
    }
  }

  // Check to see if a start time was set in the options
  double startTime;
  if(opts.count("StartTime")>0){
    startTime = opts.get<double>("StartTime");

    // Throw an error if the specified start time is AFTER the first evaluation time
    if(startTime>evalTimesIn(0) + std::numeric_limits<double>::epsilon()){
      throw std::runtime_error("ERROR: Start time specified in options to ODE is after first evaluation time.");
    }

  }else{
    startTime = evalTimesIn(0);
  }

  return startTime;
}

Eigen::VectorXi ODE::GetInputSizes(std::shared_ptr<ModPiece> const& rhsIn,
                                        bool                             isAuto)
{
  if(isAuto){
    return rhsIn->inputSizes;
  }else{
    return rhsIn->inputSizes.tail(rhsIn->inputSizes.size()-1);
  }
}

Eigen::VectorXi ODE::GetOutputSizes(std::shared_ptr<ModPiece> const& rhsIn,
                                         Eigen::VectorXd           const& evalTimes)
{
  int numTimes = evalTimes.size();
  int stateDim = rhsIn->outputSizes(0);

  Eigen::VectorXi output(1);
  output(0) = numTimes*stateDim;
  return output;
}





void ODE::NonlinearSolverOptions::SetOptions(boost::property_tree::ptree const& opts)
{
  method = opts.get("Method", "Newton");
  maxIters = opts.get("MaximumIterations", -1);
}

void ODE::LinearSolverOptions::SetOptions(boost::property_tree::ptree const& opts)
{
  method = opts.get("Method", "Dense");
  maxIters = opts.get("MaximumIterations", -1);
}

void ODE::IntegratorOptions::SetOptions(boost::property_tree::ptree const& opts)
{
  method = opts.get("Method", "BDF");
  reltol = opts.get("RelativeTolerance", 1e-5);
  abstol = opts.get("AbsoluteTolerance", 1e-5);

  maxStepSize = opts.get("MaximumStepSize", -1);
  maxNumSteps = opts.get("MaximumSteps", -1);

  checkPtGap = opts.get("CheckPointSteps", 100);
}


std::pair<SUNLinearSolver, SUNMatrix> ODE::CreateLinearSolver(void *cvode_mem, N_Vector const& state)
{
  int flag;
  int stateSize = NV_LENGTH_S(state);

  SUNMatrix rhsJac = NULL; // jacobian matrix of RHS
  SUNLinearSolver LS = NULL;

  // If the nonlinear solver is Newton's method, then we need a linear solver too.
  if(nonlinOpts.method=="Newton"){

    if(linOpts.method=="Dense")
      rhsJac = SUNDenseMatrix(stateSize, stateSize);

    if(linOpts.method=="Dense"){
      LS = SUNLinSol_Dense(state, rhsJac);
      CheckFlag((void *)LS, "SUNLinSol_Dense", 0);

    }else if(linOpts.method=="SPGMR"){
      LS = SUNLinSol_SPGMR(state, 0, linOpts.maxIters);
      CheckFlag((void *)LS, "SUNLinSol_SPGMR", 0);
    }else if(linOpts.method=="SPFGMR"){
      LS = SUNLinSol_SPFGMR(state, 0, linOpts.maxIters);
      CheckFlag((void *)LS, "SUNLinSol_SPFGMR", 0);
    }else if(linOpts.method=="SPBCGS"){
      LS = SUNLinSol_SPBCGS(state, 0, linOpts.maxIters);
      CheckFlag((void *)LS, "SUNLinSol_SPBCGS", 0);
    }else if(linOpts.method=="SPTFQMR"){
      LS = SUNLinSol_SPTFQMR(state, 0, linOpts.maxIters);
      CheckFlag((void *)LS, "SUNLinSol_SPTFQMR", 0);
    }else if(linOpts.method=="PCG"){
      LS = SUNLinSol_PCG(state, 0, linOpts.maxIters);
      CheckFlag((void *)LS, "SUNLinSol_PGC", 0);
    }else{
      std::stringstream msg;
      msg << "ERROR: Invalid linear solver method \"" << linOpts.method << "\".  Valid options are \"Dense\", \"SPGMR\", \"SPFGMR\", \"SPBCGS\", \"SPTFQMR\", or \"PCG\"\n";
      throw std::runtime_error(msg.str());
    }

    flag = CVodeSetLinearSolver(cvode_mem, LS, rhsJac);
    CheckFlag(&flag, "CVodeSetLinearSolver", 1);

    /* Set a user-supplied Jacobian function */
    if(linOpts.method=="Dense"){
      flag = CVodeSetJacFn(cvode_mem, Interface::RHSJacobian);
      CheckFlag(&flag, "CVodeSetJacFn", 1);
    }else{
      flag = CVodeSetJacTimes(cvode_mem, NULL, Interface::RHSJacobianAction);
      CheckFlag(&flag, "CVodeSetJacTimes", 1);
    }

  }

  return std::make_pair(LS, rhsJac);
}

SUNNonlinearSolver ODE::CreateNonlinearSolver(void *cvode_mem, N_Vector const& state)
{
 int flag;

 SUNNonlinearSolver NLS = NULL;
 if(nonlinOpts.method=="Newton"){
   NLS = SUNNonlinSol_Newton(state);
   CheckFlag( (void *)NLS, "SUNNonlinSol_Newton", 0);

 }else if(nonlinOpts.method=="Iter"){

   /* create fixed point nonlinear solver object */
   NLS = SUNNonlinSol_FixedPoint(state, 0);
   CheckFlag((void *)NLS, "SUNNonlinSol_FixedPoint", 0);
 }

 // Set the maximum nubmer of iterations on the nonlinear solver
 if(nonlinOpts.maxIters>0){
   flag = SUNNonlinSolSetMaxIters(NLS, nonlinOpts.maxIters);
   CheckFlag(&flag, "SUNNonlinSolSetMaxIters", 1);
 }

 flag = CVodeSetNonlinearSolver(cvode_mem, NLS);
 CheckFlag(&flag, "CVodeSetNonlinearSolver", 1);

 return NLS;
}


std::pair<SUNLinearSolver, SUNMatrix> ODE::CreateAdjointLinearSolver(void *cvode_mem, N_Vector const& state, int indexB)
{
  int flag;
  int stateSize = NV_LENGTH_S(state);

  SUNMatrix rhsJac = NULL; // jacobian matrix of RHS
  SUNLinearSolver LS = NULL;

  // If the nonlinear solver is Newton's method, then we need a linear solver too.
  if(nonlinOpts.method=="Newton"){

    if(linOpts.method=="Dense")
      rhsJac = SUNDenseMatrix(stateSize, stateSize);

    if(linOpts.method=="Dense"){
      LS = SUNLinSol_Dense(state, rhsJac);
      CheckFlag((void *)LS, "SUNLinSol_Dense", 0);

    }else if(linOpts.method=="SPGMR"){
      LS = SUNLinSol_SPGMR(state, 0, linOpts.maxIters);
      CheckFlag((void *)LS, "SUNLinSol_SPGMR", 0);
    }else if(linOpts.method=="SPFGMR"){
      LS = SUNLinSol_SPFGMR(state, 0, linOpts.maxIters);
      CheckFlag((void *)LS, "SUNLinSol_SPFGMR", 0);
    }else if(linOpts.method=="SPBCGS"){
      LS = SUNLinSol_SPBCGS(state, 0, linOpts.maxIters);
      CheckFlag((void *)LS, "SUNLinSol_SPBCGS", 0);
    }else if(linOpts.method=="SPTFQMR"){
      LS = SUNLinSol_SPTFQMR(state, 0, linOpts.maxIters);
      CheckFlag((void *)LS, "SUNLinSol_SPTFQMR", 0);
    }else if(linOpts.method=="PCG"){
      LS = SUNLinSol_PCG(state, 0, linOpts.maxIters);
      CheckFlag((void *)LS, "SUNLinSol_PGC", 0);
    }else{
      std::stringstream msg;
      msg << "ERROR: Invalid linear solver method \"" << linOpts.method << "\".  Valid options are \"Dense\", \"SPGMR\", \"SPFGMR\", \"SPBCGS\", \"SPTFQMR\", or \"PCG\"\n";
      throw std::runtime_error(msg.str());
    }

    flag = CVodeSetLinearSolverB(cvode_mem, indexB, LS, rhsJac);
    CheckFlag(&flag, "CVodeSetLinearSolverB", 1);

    /* Set a user-supplied Jacobian function */
    if(linOpts.method=="Dense"){
      flag = CVodeSetJacFnB(cvode_mem, indexB, &Interface::AdjointJacobian);
      CheckFlag(&flag, "CVodeSetJacFnB", 1);
    }else{
      flag = CVodeSetJacTimesB(cvode_mem, indexB, NULL, &Interface::AdjointJacobianAction);
      CheckFlag(&flag, "CVodeSetJacTimesB", 1);
    }
  }

  return std::make_pair(LS, rhsJac);
}

SUNNonlinearSolver ODE::CreateAdjointNonlinearSolver(void *cvode_mem, N_Vector const& state, int indexB)
{
 int flag;

 SUNNonlinearSolver NLS = NULL;
 if(nonlinOpts.method=="Newton"){
   NLS = SUNNonlinSol_Newton(state);
   CheckFlag( (void *)NLS, "SUNNonlinSol_Newton", 0);

 }else if(nonlinOpts.method=="Iter"){

   /* create fixed point nonlinear solver object */
   NLS = SUNNonlinSol_FixedPoint(state, 0);
   CheckFlag((void *)NLS, "SUNNonlinSol_FixedPoint", 0);
 }

 // Set the maximum nubmer of iterations on the nonlinear solver
 if(nonlinOpts.maxIters>0){
   flag = SUNNonlinSolSetMaxIters(NLS, nonlinOpts.maxIters);
   CheckFlag(&flag, "SUNNonlinSolSetMaxIters", 1);
 }

 flag = CVodeSetNonlinearSolverB(cvode_mem, indexB, NLS);
 CheckFlag(&flag, "CVodeSetNonlinearSolverB", 1);

 return NLS;
}


void ODE::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  // Used for monitoring the success of CVODE calls
  int flag;

  // The initial condition
  Eigen::VectorXd const& x0  = inputs.at(0).get();
  const unsigned int stateSize = x0.size();

  auto rhsData = std::make_shared<ODEData>(rhs,
                                           ref_vector<Eigen::VectorXd>(inputs.begin(),  inputs.begin()+(isAutonomous ? rhs->numInputs : rhs->numInputs-1)),
                                           isAutonomous,
                                           -1);

  // Create an NVector version of the state that CVODES understands.
  // Also map the memory in that vector to an object Eigen can understand.
  N_Vector state = N_VNew_Serial(stateSize);
  Eigen::Map<Eigen::VectorXd> stateMap(NV_DATA_S(state), stateSize);

  // Set the initial condition
  stateMap = x0;

  // Set up for forward integration
  void *cvode_mem = NULL;
  if(intOpts.method=="BDF"){
    cvode_mem = CVodeCreate(CV_BDF);
  }else if(intOpts.method=="Adams"){
    cvode_mem = CVodeCreate(CV_ADAMS);
  }else{
    std::stringstream msg;
    msg << "ERROR: Invalid integrator method \"" << intOpts.method << "\".  Valid options are \"BDF\" or \"ADAMS\"\n";
    throw std::runtime_error(msg.str());
  }

  // Initialize CVODES with the RHS function, initial time, and initial state.
  flag = CVodeInit(cvode_mem, Interface::RHS, startTime, state);
  CheckFlag(&flag, "CVodeSetUserData", 1);

  flag = CVodeSStolerances(cvode_mem, intOpts.reltol, intOpts.abstol);
  CheckFlag(&flag, "CVodeSStolerances", 1);

  flag = CVodeSetUserData(cvode_mem, rhsData.get());
  CheckFlag(&flag, "CVodeSetUserData", 1);

  /* -------- Linear Solver and Jacobian ----------- */
  SUNMatrix rhsJac = NULL;
  SUNLinearSolver LS;
  std::tie(LS,rhsJac) = CreateLinearSolver(cvode_mem, state);

  /* -------- Nonlinear solver ----------- */
  SUNNonlinearSolver NLS = NULL;
  NLS = CreateNonlinearSolver(cvode_mem, state);

  /* Set max steps between outputs */
  if(intOpts.maxNumSteps>0){
    flag = CVodeSetMaxNumSteps(cvode_mem, intOpts.maxNumSteps);
    CheckFlag(&flag, "CVodeSetMaxNumSteps", 1);
  }

  // Resize the output vector
  outputs.resize(1);
  outputs.at(0).resize(outputSizes(0));

  // The current time
  double t = startTime;
  unsigned int tind=0;
  // Check to see if the starting time is equal to the first evalTime
  if(std::fabs(evalTimes(0)-startTime)<1e-14){
    // copy the initial state into the output
    outputs.at(0).segment(tind*stateSize, stateSize) = x0;
    // increment the time index
    ++tind;
  }

  // loop through all the times where we want evaluations
  for(; tind<evalTimes.size(); ++tind ) {
    flag = CVode(cvode_mem, evalTimes(tind), state, &t, CV_NORMAL);
    CheckFlag(&flag, "CVode", 1);

    outputs.at(0).segment(tind*stateSize, stateSize) = stateMap;
  }

  N_VDestroy(state);
  if(rhsJac != NULL)
    SUNMatDestroy(rhsJac);

  CVodeFree(&cvode_mem);

  if(LS != NULL)
    SUNLinSolFree(LS);

  SUNNonlinSolFree(NLS);
}

void ODE::GradientImpl(unsigned int outWrt,
                            unsigned int inWrt,
                            ref_vector<Eigen::VectorXd> const& inputs,
                            Eigen::VectorXd const&             sens)
{
  // Used for monitoring the success of CVODE calls
  int flag;

  // The initial condition
  Eigen::VectorXd const& x0  = inputs.at(0).get();
  const unsigned int stateSize = x0.size();

  auto rhsData = std::make_shared<ODEData>(rhs,
                                           ref_vector<Eigen::VectorXd>(inputs.begin(),  inputs.begin()+(isAutonomous ? rhs->numInputs : rhs->numInputs-1)),
                                           isAutonomous,
                                           inWrt);

  // Create an NVector version of the state that CVODES understands.
  // Also map the memory in that vector to an object Eigen can understand.
  N_Vector state = N_VNew_Serial(stateSize);
  Eigen::Map<Eigen::VectorXd> stateMap(NV_DATA_S(state), stateSize);

  // Set the initial condition
  stateMap = x0;

  // Set up for forward integration
  void *cvode_mem = NULL;
  if(intOpts.method=="BDF"){
    cvode_mem = CVodeCreate(CV_BDF);
  }else if(intOpts.method=="Adams"){
    cvode_mem = CVodeCreate(CV_ADAMS);
  }else{
    std::stringstream msg;
    msg << "ERROR: Invalid integrator method \"" << intOpts.method << "\".  Valid options are \"BDF\" or \"ADAMS\"\n";
    throw std::runtime_error(msg.str());
  }

  // Initialize CVODES with the RHS function, initial time, and initial state.
  flag = CVodeInit(cvode_mem, Interface::RHS, startTime, state);
  CheckFlag(&flag, "CVodeSetUserData", 1);

  flag = CVodeSStolerances(cvode_mem, intOpts.reltol, intOpts.abstol);
  CheckFlag(&flag, "CVodeSStolerances", 1);

  flag = CVodeSetUserData(cvode_mem, rhsData.get());
  CheckFlag(&flag, "CVodeSetUserData", 1);

  /* -------- Linear Solver and Jacobian ----------- */
  SUNMatrix rhsJac = NULL;
  SUNLinearSolver LS;
  std::tie(LS,rhsJac) = CreateLinearSolver(cvode_mem, state);

  /* -------- Nonlinear solver ----------- */
  SUNNonlinearSolver NLS = NULL;
  NLS = CreateNonlinearSolver(cvode_mem, state);

  /* Set max steps between outputs */
  if(intOpts.maxNumSteps>0){
    flag = CVodeSetMaxNumSteps(cvode_mem, intOpts.maxNumSteps);
    CheckFlag(&flag, "CVodeSetMaxNumSteps", 1);
  }

  // Initialize the adjoint problem and specify the number of steps between checkpoints
  flag = CVodeAdjInit(cvode_mem, intOpts.checkPtGap, CV_HERMITE);
  CheckFlag(&flag, "CVodeAdjInit", 1);

  /* --------- Forward Integration ---------------- */

  // Resize the output vector
  outputs.resize(1);
  outputs.at(0).resize(outputSizes(0));

  // The current time
  double t = startTime;
  int tind=0;
  int ncheck = 0;

  // Check to see if the starting time is equal to the first evalTime
  if(std::fabs(evalTimes(0)-startTime)<1e-14){
    // copy the initial state into the output
    outputs.at(0).segment(tind*stateSize, stateSize) = x0;
    // increment the time index
    ++tind;
  }

  // loop through all the times where we want evaluations
  for(; tind<evalTimes.size(); ++tind ) {
    flag = CVodeF(cvode_mem, evalTimes(tind), state, &t, CV_NORMAL, &ncheck);
    CheckFlag(&flag, "CVodeF", 1);

    outputs.at(0).segment(tind*stateSize, stateSize) = stateMap;
  }


  /* --------- Set up for adjoint solve ---------------- */
  const unsigned int paramSize = inputSizes(inWrt);

  N_Vector lambda = N_VNew_Serial(stateSize);
  CheckFlag((void *)lambda, "N_VNew", 0);

  N_Vector nvGrad = N_VNew_Serial(paramSize);
  CheckFlag((void *)nvGrad, "N_VNew", 0);

  Eigen::Map<Eigen::VectorXd> lambdaMap(NV_DATA_S(lambda), stateSize);
  Eigen::Map<Eigen::VectorXd> nvGradMap(NV_DATA_S(nvGrad), paramSize);

  lambdaMap = -1.0*sens.tail(stateSize);
  nvGradMap = Eigen::VectorXd::Zero(paramSize);

  // Create the backward integrator
  int indexB;
  if(intOpts.method=="BDF"){
    flag = CVodeCreateB(cvode_mem, CV_BDF, &indexB);
  }else if(intOpts.method=="Adams"){
    flag = CVodeCreateB(cvode_mem, CV_ADAMS, &indexB);
  }else{
    std::stringstream msg;
    msg << "ERROR: Invalid integrator method \"" << intOpts.method << "\".  Valid options are \"BDF\" or \"ADAMS\"\n";
    throw std::runtime_error(msg.str());
  }
  CheckFlag(&flag, "CVodeCreateB", 1);

  // Initialize the backward integrator
  double finalTime = evalTimes(evalTimes.size()-1);
  flag = CVodeInitB(cvode_mem, indexB, Interface::AdjointRHS, finalTime, lambda);
  CheckFlag(&flag, "CVodeInitB", 1);

  // Set the scalar relative and absolute tolerance
  flag = CVodeSStolerancesB(cvode_mem, indexB, intOpts.reltol, intOpts.abstol);
  CheckFlag(&flag, "CVodeSStolerancesB", 1);

  // Attach user data
  flag = CVodeSetUserDataB(cvode_mem, indexB, rhsData.get());
  CheckFlag(&flag, "CVodeSetUserDataB", 1);


  // Set up the linear and nonlinear solvers for the adjoint problem.
  SUNMatrix adjRhsJac = NULL;
  SUNLinearSolver adjLS;
  std::tie(adjLS,adjRhsJac) = CreateAdjointLinearSolver(cvode_mem, lambda, indexB);

  SUNNonlinearSolver adjNLS = NULL;
  adjNLS = CreateAdjointNonlinearSolver(cvode_mem, lambda, indexB);

  // Initialize backward quadrature to accumulate sensitivity wrt parameters
  flag = CVodeQuadInitB(cvode_mem, indexB, &Interface::AdjointQuad, nvGrad);
  CheckFlag(&flag, "CVodeQuadInitB", 1);

  // Tell CVODES whether or not the quadrature should be included in the error control
  flag = CVodeSetQuadErrConB(cvode_mem, indexB, true);
  CheckFlag(&flag, "CVodeSetQuadErrConB", 1);

  // Set the quadrature tolerances
  flag = CVodeQuadSStolerancesB(cvode_mem, indexB, intOpts.reltol, intOpts.abstol);
  CheckFlag(&flag, "CVodeQuadSStolerancesB", 1);

  /* --------- Adjoint (Backward) Integration ---------------- */

  for(tind=evalTimes.size()-2; tind>=0; --tind ) {

    // integrate the adjoint variable back in time
    flag = CVodeB(cvode_mem, evalTimes(tind), CV_NORMAL);
    CheckFlag(&flag, "CVodeB", 1);

    // get the adjoint variable
    flag = CVodeGetB(cvode_mem, indexB, &t, lambda);
    CheckFlag(&flag, "CVodeGetB", 1);

    // get the quadrature
    flag = CVodeGetQuadB(cvode_mem, indexB, &t, nvGrad);
    CheckFlag(&flag, "CVodeGetQuadB", 1);

    // discontinuous change in variable
    lambdaMap -= sens.segment(tind*stateSize, stateSize);

    // reinitialize the adjoint integrator because of the discontinuity in the state
    flag = CVodeReInitB(cvode_mem, indexB, evalTimes(tind), lambda);
    CheckFlag(&flag, "CVodeReInitB", 1);
    flag = CVodeQuadReInitB(cvode_mem, indexB, nvGrad);
    CheckFlag(&flag, "CVodeQuadReInitB", 1);
  }

  // Make sure we integrate back to the starting time
  if( std::fabs(t-startTime)>1.0e-14 ) {
    // integrate the adjoint variable back in time
    flag = CVodeB(cvode_mem, startTime, CV_NORMAL);
    CheckFlag(&flag, "CVodeB", 1);

    // get the adjoint variable
    flag = CVodeGetB(cvode_mem, indexB, &t, lambda);
    CheckFlag(&flag, "CVodeGetB", 1);
  }

  // Extract the gradient
  flag = CVodeGetQuadB(cvode_mem, indexB, &t, nvGrad);
  CheckFlag(&flag, "CVodeGetQuadB", 1);

  gradient = nvGradMap;

  // if wrt to the initial conditions add identity*sens
  if( inWrt==0 ) {
    for(tind=0; tind<evalTimes.size(); ++tind ) {
      gradient += sens.segment(tind*stateSize, paramSize);
    }
  }

  N_VDestroy(state);
  N_VDestroy(lambda);
  N_VDestroy(nvGrad);

  CVodeFree(&cvode_mem);

  SUNNonlinSolFree(NLS);
  SUNNonlinSolFree(adjNLS);

  if(LS != NULL)
    SUNLinSolFree(LS);
  if(adjLS != NULL)
    SUNLinSolFree(adjLS);

  if(rhsJac != NULL)
    SUNMatDestroy(rhsJac);
  if(adjRhsJac != NULL)
    SUNMatDestroy(adjRhsJac);

  // Make sure to update the one-step cache since we've updated the outputs vector
  if(cacheEnabled){
    cacheInput.resize(inputs.size());
    for(int i=0; i<inputs.size(); ++i)
      cacheInput.at(i) = inputs.at(i);
  }
}



void ODE::JacobianImpl(unsigned int outWrt,
                            unsigned int inWrt,
                            ref_vector<Eigen::VectorXd> const& inputs)
{
  // Used for monitoring the success of CVODE calls
  int flag;

  // The initial condition
  Eigen::VectorXd const& x0  = inputs.at(0).get();
  const unsigned int stateSize = x0.size();

  auto rhsData = std::make_shared<ODEData>(rhs,
                                           ref_vector<Eigen::VectorXd>(inputs.begin(),  inputs.begin()+(isAutonomous ? rhs->numInputs : rhs->numInputs-1)),
                                           isAutonomous,
                                           inWrt);

  // Create an NVector version of the state that CVODES understands.
  // Also map the memory in that vector to an object Eigen can understand.
  N_Vector state = N_VNew_Serial(stateSize);
  Eigen::Map<Eigen::VectorXd> stateMap(NV_DATA_S(state), stateSize);

  // Set the initial condition
  stateMap = x0;

  // Set up for forward integration
  void *cvode_mem = NULL;
  if(intOpts.method=="BDF"){
    cvode_mem = CVodeCreate(CV_BDF);
  }else if(intOpts.method=="Adams"){
    cvode_mem = CVodeCreate(CV_ADAMS);
  }else{
    std::stringstream msg;
    msg << "ERROR: Invalid integrator method \"" << intOpts.method << "\".  Valid options are \"BDF\" or \"ADAMS\"\n";
    throw std::runtime_error(msg.str());
  }

  // Initialize CVODES with the RHS function, initial time, and initial state.
  flag = CVodeInit(cvode_mem, Interface::RHS, startTime, state);
  CheckFlag(&flag, "CVodeSetUserData", 1);

  flag = CVodeSStolerances(cvode_mem, intOpts.reltol, intOpts.abstol);
  CheckFlag(&flag, "CVodeSStolerances", 1);

  flag = CVodeSetUserData(cvode_mem, rhsData.get());
  CheckFlag(&flag, "CVodeSetUserData", 1);

  /* Set max steps between outputs */
  if(intOpts.maxNumSteps>0){
    flag = CVodeSetMaxNumSteps(cvode_mem, intOpts.maxNumSteps);
    CheckFlag(&flag, "CVodeSetMaxNumSteps", 1);
  }

  /* -------- Linear Solver and Jacobian ----------- */
  SUNMatrix rhsJac = NULL;
  SUNLinearSolver LS;
  std::tie(LS,rhsJac) = CreateLinearSolver(cvode_mem, state);

  /* -------- Nonlinear solver ----------- */
  SUNNonlinearSolver NLS = NULL;
  NLS = CreateNonlinearSolver(cvode_mem, state);

  /* -------- Sensitivity Setup ----------- */
  const unsigned int paramSize = inputSizes(inWrt);

  N_Vector *sensState = nullptr;
  sensState = N_VCloneVectorArray_Serial(paramSize, state);
  CheckFlag((void *)sensState, "N_VCloneVectorArray_Serial", 0);

  // initialize the sensitivies.
  for( int is=0; is<paramSize; ++is )
    N_VConst(0.0, sensState[is]);

  if(inWrt==0){ // If wrt the initial state, the initial sensitivities form the identity
    for(int is=0; is<paramSize; ++is){
      NV_Ith_S(sensState[is],is) = 1.0;
    }
  }

  flag = CVodeSensInit(cvode_mem, paramSize, CV_SIMULTANEOUS, Interface::SensitivityRHS, sensState);
  CheckFlag(&flag, "CVodeSensInit", 1);

  SUNNonlinearSolver sensNLS = NULL;
  if(nonlinOpts.method=="Newton"){
    sensNLS = SUNNonlinSol_NewtonSens(paramSize+1,state);
  }else{
    sensNLS = SUNNonlinSol_FixedPointSens(paramSize+1, state, 0);
  }

  flag = CVodeSetNonlinearSolverSensSim(cvode_mem, sensNLS);
  CheckFlag((void *)sensState, "CVodeSetNonlinearSolverSensSim", 0);

  // set sensitivity tolerances
  flag = CVodeSensEEtolerances(cvode_mem);
  CheckFlag(&flag, "CVodeSensEEtolerances", 1);

  //Eigen::VectorXd absTolVec = Eigen::VectorXd::Constant(paramSize, intOpts.abstol);
  // flag = CVodeSensSStolerances(cvode_mem, intOpts.reltol, absTolVec.data());
  // CheckFlag(&flag, "CVodeSensSStolerances", 1);

  // // error control strategy should test the sensitivity variables
  // flag = CVodeSetSensErrCon(cvode_mem, true);
  // CheckFlag(&flag, "CVodeSetSensErrCon", 1);

  // Resize the output vector
  outputs.resize(1);
  outputs.at(0).resize(outputSizes(0));

  // Resize the jacobian vector
  jacobian = Eigen::MatrixXd::Zero(outputSizes(0), paramSize);

  // The current time
  double t = startTime;
  unsigned int tind=0;

  // Check to see if the starting time is equal to the first evalTime
  if(std::fabs(evalTimes(0)-startTime)<1e-14){

    // copy the initial state into the output
    outputs.at(0).segment(tind*stateSize, stateSize) = x0;

    // If Jacobian is wrt state, then initial jac will be the identity, otherwise initial jac will be zero
    if(inWrt==0){
      jacobian.block(tind*stateSize, 0, stateSize, paramSize) = Eigen::MatrixXd::Identity(stateSize, paramSize);
    }

    // increment the time index
    ++tind;
  }

  // loop through all the times where we want evaluations
  for(; tind<evalTimes.size(); ++tind ) {

    // Integrate the ODE and the sensitivity equations forward in time
    flag = CVode(cvode_mem, evalTimes(tind), state, &t, CV_NORMAL);
    CheckFlag(&flag, "CVode", 1);

    // Extract the ODE state
    outputs.at(0).segment(tind*stateSize, stateSize) = stateMap;

    // Extract the sensitivities
    int flag = CVodeGetSens(cvode_mem, &t, sensState);
    CheckFlag(&flag, "CVodeGetSense", 1);

    for( unsigned int i=0; i<stateSize; ++i ) {
      for( unsigned int j=0; j<paramSize; ++j ) {
        jacobian(tind*stateSize+i, j) += NV_Ith_S(sensState[j], i);
      }
    }

  } // end of time loop


  N_VDestroy(state);
  N_VDestroyVectorArray_Serial(sensState, paramSize);

  if(rhsJac != NULL)
    SUNMatDestroy(rhsJac);

  CVodeFree(&cvode_mem);

  if(LS != NULL)
    SUNLinSolFree(LS);

  SUNNonlinSolFree(NLS);
  SUNNonlinSolFree(sensNLS);

  // Make sure to update the one-step cache since we've updated the outputs vector
  if(cacheEnabled){
    cacheInput.resize(inputs.size());
    for(int i=0; i<inputs.size(); ++i)
      cacheInput.at(i) = inputs.at(i);
  }
}



void ODE::ApplyJacobianImpl(unsigned int outWrt,
                                  unsigned int inWrt,
                                  ref_vector<Eigen::VectorXd> const& inputs,
                                  Eigen::VectorXd const& vec)
{
  // Used for monitoring the success of CVODE calls
  int flag;

  // The initial condition
  Eigen::VectorXd const& x0  = inputs.at(0).get();
  const unsigned int stateSize = x0.size();

  auto rhsData = std::make_shared<ODEData>(rhs,
                                           ref_vector<Eigen::VectorXd>(inputs.begin(),  inputs.begin()+(isAutonomous ? rhs->numInputs : rhs->numInputs-1)),
                                           isAutonomous,
                                           inWrt,
                                           vec);

  // Create an NVector version of the state that CVODES understands.
  // Also map the memory in that vector to an object Eigen can understand.
  N_Vector state = N_VNew_Serial(stateSize);
  Eigen::Map<Eigen::VectorXd> stateMap(NV_DATA_S(state), stateSize);

  // Set the initial condition
  stateMap = x0;

  // Set up for forward integration
  void *cvode_mem = NULL;
  if(intOpts.method=="BDF"){
    cvode_mem = CVodeCreate(CV_BDF);
  }else if(intOpts.method=="Adams"){
    cvode_mem = CVodeCreate(CV_ADAMS);
  }else{
    std::stringstream msg;
    msg << "ERROR: Invalid integrator method \"" << intOpts.method << "\".  Valid options are \"BDF\" or \"ADAMS\"\n";
    throw std::runtime_error(msg.str());
  }

  // Initialize CVODES with the RHS function, initial time, and initial state.
  flag = CVodeInit(cvode_mem, Interface::RHS, startTime, state);
  CheckFlag(&flag, "CVodeSetUserData", 1);

  flag = CVodeSStolerances(cvode_mem, intOpts.reltol, intOpts.abstol);
  CheckFlag(&flag, "CVodeSStolerances", 1);

  flag = CVodeSetUserData(cvode_mem, rhsData.get());
  CheckFlag(&flag, "CVodeSetUserData", 1);

  /* Set max steps between outputs */
  if(intOpts.maxNumSteps>0){
    flag = CVodeSetMaxNumSteps(cvode_mem, intOpts.maxNumSteps);
    CheckFlag(&flag, "CVodeSetMaxNumSteps", 1);
  }

  /* -------- Linear Solver and Jacobian ----------- */
  SUNMatrix rhsJac = NULL;
  SUNLinearSolver LS;
  std::tie(LS,rhsJac) = CreateLinearSolver(cvode_mem, state);

  /* -------- Nonlinear solver ----------- */
  SUNNonlinearSolver NLS = NULL;
  NLS = CreateNonlinearSolver(cvode_mem, state);

  /* -------- Sensitivity Setup ----------- */
  const unsigned int paramSize = 1;

  N_Vector *sensState = nullptr;
  sensState = N_VCloneVectorArray_Serial(paramSize, state);
  CheckFlag((void *)sensState, "N_VCloneVectorArray_Serial", 0);

  // initialize the sensitivies.
  N_VConst(0.0, sensState[0]);

  if(inWrt==0){ // If wrt the initial state, the initial sensitivity is the vector
    for(int ind=0; ind<stateSize; ++ind){
      NV_Ith_S(sensState[0],ind) = vec(ind);
    }
  }

  flag = CVodeSensInit(cvode_mem, paramSize, CV_SIMULTANEOUS, Interface::SensitivityRHS, sensState);
  CheckFlag(&flag, "CVodeSensInit", 1);

  SUNNonlinearSolver sensNLS = NULL;
  if(nonlinOpts.method=="Newton"){
    sensNLS = SUNNonlinSol_NewtonSens(paramSize+1,state);
  }else{
    sensNLS = SUNNonlinSol_FixedPointSens(paramSize+1, state, 0);
  }

  flag = CVodeSetNonlinearSolverSensSim(cvode_mem, sensNLS);
  CheckFlag((void *)sensState, "CVodeSetNonlinearSolverSensSim", 0);

  // set sensitivity tolerances
  flag = CVodeSensEEtolerances(cvode_mem);
  CheckFlag(&flag, "CVodeSensEEtolerances", 1);

  //Eigen::VectorXd absTolVec = Eigen::VectorXd::Constant(paramSize, intOpts.abstol);
  // flag = CVodeSensSStolerances(cvode_mem, intOpts.reltol, absTolVec.data());
  // CheckFlag(&flag, "CVodeSensSStolerances", 1);

  // // error control strategy should test the sensitivity variables
  // flag = CVodeSetSensErrCon(cvode_mem, true);
  // CheckFlag(&flag, "CVodeSetSensErrCon", 1);

  // Resize the output vector
  outputs.resize(1);
  outputs.at(0).resize(evalTimes.size()*stateSize);

  // Resize the jacobian vector
  jacobianAction = Eigen::VectorXd::Zero(evalTimes.size()*stateSize);

  // The current time
  double t = startTime;
  unsigned int tind=0;

  // Check to see if the starting time is equal to the first evalTime
  if(std::fabs(evalTimes(0)-startTime)<1e-14){

    // copy the initial state into the output
    outputs.at(0).segment(tind*stateSize, stateSize) = x0;

    // If Jacobian is wrt state, then initial jac will be the identity, otherwise initial jac will be zero
    if(inWrt==0){
      jacobianAction.segment(tind*stateSize,stateSize) = vec;
    }

    // increment the time index
    ++tind;
  }

  // loop through all the times where we want evaluations
  for(; tind<evalTimes.size(); ++tind ) {

    // Integrate the ODE and the sensitivity equations forward in time
    flag = CVode(cvode_mem, evalTimes(tind), state, &t, CV_NORMAL);
    CheckFlag(&flag, "CVode", 1);

    // Extract the ODE state
    outputs.at(0).segment(tind*stateSize, stateSize) = stateMap;

    // Extract the sensitivities
    int flag = CVodeGetSens(cvode_mem, &t, sensState);
    CheckFlag(&flag, "CVodeGetSense", 1);

    for( unsigned int i=0; i<stateSize; ++i ) {
      jacobianAction(tind*stateSize + i) += NV_Ith_S(sensState[0], i);
    }

  } // end of time loop


  N_VDestroy(state);
  N_VDestroyVectorArray_Serial(sensState, paramSize);

  if(rhsJac != NULL)
    SUNMatDestroy(rhsJac);

  CVodeFree(&cvode_mem);

  if(LS != NULL)
    SUNLinSolFree(LS);

  SUNNonlinSolFree(NLS);
  SUNNonlinSolFree(sensNLS);

  // Make sure to update the one-step cache since we've updated the outputs vector
  if(cacheEnabled){
    cacheInput.resize(inputs.size());
    for(int i=0; i<inputs.size(); ++i)
      cacheInput.at(i) = inputs.at(i);
  }
}


int ODE::Interface::RHS(realtype time, N_Vector state, N_Vector timeDeriv, void *userData)
{
  // get the data type
  ODEData* data = (ODEData*)userData;
  assert(data);
  assert(data->rhs);

  // set the state input
  Eigen::Map<Eigen::VectorXd> stateMap(NV_DATA_S(state), NV_LENGTH_S(state));
  data->UpdateInputs(stateMap, time);

  // evaluate the rhs
  Eigen::Map<Eigen::VectorXd> timeDerivMap(NV_DATA_S(timeDeriv), NV_LENGTH_S(timeDeriv));
  timeDerivMap = data->rhs->Evaluate(data->inputs) [0];

  return 0;
}


int ODE::Interface::RHSJacobianAction(N_Vector v,  // The vector we want to apply the Jacobian to
                                           N_Vector Jv, // The Jacobian-vector product (output of this function)
                                           double time, // The current time
                                           N_Vector state, // The current state
                                           N_Vector timeDeriv, // The current time derivative
                                           void *userData,
                                           N_Vector tmp) {

  // get the data type
  ODEData* data = (ODEData*) userData;
  assert(data);
  assert(data->rhs);

  // set the state input
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
  data->UpdateInputs(stateref, time);

  Eigen::VectorXd vEig = Eigen::Map<Eigen::VectorXd>(NV_DATA_S(v), NV_LENGTH_S(v));

  // Compute the action of the jacobian (wrt the state) on a vector
  Eigen::Map<Eigen::VectorXd> JvMap(NV_DATA_S(Jv), NV_LENGTH_S(Jv));
  JvMap = data->rhs->ApplyJacobian(0, data->autonomous? 0 : 1, data->inputs, vEig);

  return 0;
}

int ODE::Interface::RHSJacobian(double time,
                                     N_Vector state,
                                     N_Vector rhs,
                                     SUNMatrix jac,
                                     void *userData,
                                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  // get the data type
  ODEData* data = (ODEData*)userData;
  assert(data);
  assert(data->rhs);

  // set the state input
  Eigen::Map<Eigen::VectorXd> stateMap(NV_DATA_S(state), NV_LENGTH_S(state));
  data->UpdateInputs(stateMap, time);

  // evaluate the jacobian
  Eigen::Map<Eigen::MatrixXd> jacMap(SUNDenseMatrix_Data(jac),
                                     SUNDenseMatrix_Rows(jac),
                                     SUNDenseMatrix_Columns(jac));

  jacMap = data->rhs->Jacobian(0, data->autonomous? 0 : 1, data->inputs);

  return 0;
}


int ODE::Interface::SensitivityRHS(int Ns, realtype time, N_Vector y, N_Vector ydot, N_Vector *ys, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {

  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);
  unsigned int stateWrt = data->autonomous? 0 : 1; // Which input of the rhs holds the state?

  // set the state input
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(y), NV_LENGTH_S(y));
  data->UpdateInputs(stateref, time);

  for( unsigned int i=0; i<Ns; ++i ) {
    Eigen::Map<Eigen::VectorXd> rhsMap(NV_DATA_S(ySdot[i]), stateref.size());
    Eigen::Map<Eigen::VectorXd> sensMap(NV_DATA_S(ys[i]), stateref.size());

    // Fill in (df/dy)*s_i
    rhsMap = data->rhs->ApplyJacobian(0, stateWrt, data->inputs, (Eigen::VectorXd)sensMap);
  }

  /* Check if the derivative is with respect to the initial condition or wrt
     to a parameter.  If wrt to the initial conditions, then df/dp = 0 here.
  */
  if(data->wrtIn>0){

    // the derivative of the rhs wrt the parameter
    unsigned int rhsWrt = data->autonomous? data->wrtIn : data->wrtIn+1;

    if(data->isAction){

      Eigen::Map<Eigen::VectorXd> rhsMap(NV_DATA_S(ySdot[0]), stateref.size());
      //Eigen::Map<Eigen::VectorXd> sensMap(NV_DATA_S(ys[0]), stateref.size());
      rhsMap += data->rhs->ApplyJacobian(0,rhsWrt, data->inputs, data->actionVec);

    }else{
      const Eigen::MatrixXd& dfdp = data->rhs->Jacobian(0, rhsWrt, data->inputs);

      // add the df/dp_i component for each index of the parameter p_i
      for( unsigned int i=0; i<Ns; ++i ) {
        Eigen::Map<Eigen::VectorXd> rhsMap(NV_DATA_S(ySdot[i]), stateref.size());
        //Eigen::Map<Eigen::VectorXd> sensMap(NV_DATA_S(ys[i]), stateref.size());

        rhsMap += dfdp.col(i);
      }
    }
  }

  return 0;
}



int ODE::Interface::AdjointRHS(realtype time, N_Vector state, N_Vector lambda, N_Vector deriv, void *user_data) {

  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
  data->UpdateInputs(stateref, time);

  Eigen::Map<Eigen::VectorXd> lambdaMap(NV_DATA_S(lambda), NV_LENGTH_S(lambda));
  Eigen::VectorXd lam = lambdaMap;

  Eigen::Map<Eigen::VectorXd> derivMap(NV_DATA_S(deriv), NV_LENGTH_S(deriv));

  derivMap = -1.0*data->rhs->Gradient(0, data->autonomous? 0 : 1, data->inputs, lam);

  return 0;
}


int ODE::Interface::AdjointJacobianAction(N_Vector target, N_Vector output, realtype time, N_Vector state, N_Vector lambda, N_Vector adjRhs, void *user_data, N_Vector tmp) {

  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
  data->UpdateInputs(stateref, time);

  Eigen::Map<Eigen::VectorXd> targetMap(NV_DATA_S(target), NV_LENGTH_S(target));
  Eigen::VectorXd targ = targetMap;

  Eigen::Map<Eigen::VectorXd> outputMap(NV_DATA_S(output), NV_LENGTH_S(output));
  outputMap = -1.0*data->rhs->Gradient(0, data->autonomous? 0 : 1, data->inputs, targ);

  return 0;
}



int ODE::Interface::AdjointJacobian(double time, N_Vector state, N_Vector lambda, N_Vector rhs, SUNMatrix jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
  data->UpdateInputs(stateref, time);

  // evaluate the jacobian
  Eigen::Map<Eigen::MatrixXd> jacMap(SUNDenseMatrix_Data(jac),
                                     SUNDenseMatrix_Rows(jac),
                                     SUNDenseMatrix_Columns(jac));

  jacMap = -1.0*data->rhs->Jacobian(0, data->autonomous? 0 : 1, data->inputs).transpose();

  return 0;
}

int ODE::Interface::AdjointQuad(realtype time, N_Vector state, N_Vector lambda, N_Vector quadRhs, void *user_data) {

  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
  Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
  data->UpdateInputs(stateref, time);

  Eigen::Map<Eigen::VectorXd> lambdaMap(NV_DATA_S(lambda), NV_LENGTH_S(lambda));
  Eigen::VectorXd lam = lambdaMap;

  Eigen::Map<Eigen::VectorXd> nvQuadMap(NV_DATA_S(quadRhs), NV_LENGTH_S(quadRhs));

  nvQuadMap = data->rhs->Gradient(0, data->autonomous? data->wrtIn : data->wrtIn+1, data->inputs, lam);

  return 0;
}




void ODE::CheckFlag(void *returnvalue, const char *funcname, int opt)
{
  int *retval;
  std::stringstream msg;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    msg << "SUNDIALS_ERROR: " << funcname << " failed - returned NULL pointer";
    throw std::runtime_error(msg.str());

  /* Check if retval < 0 */
  }else if (opt == 1) {
    retval = (int *) returnvalue;

    if (*retval < 0){
      msg << "SUNDIALS_ERROR: " << funcname << " failed with retval = " << *retval;
      throw std::runtime_error(msg.str());
    }

  /* Check if function returned NULL pointer - no memory allocated */
  }else if (opt == 2 && returnvalue == NULL) {
    msg << "MEMORY_ERROR: " << funcname << " failed - returned NULL pointer";
    throw std::runtime_error(msg.str());
  }

}
