#include "MUQ/Modeling/LinearAlgebra/LOBPCG.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include <boost/property_tree/ptree.hpp>

using namespace muq::Modeling;
using namespace muq::Utilities;

LOBPCG::LOBPCG(int    numEigsIn,
               double eigRelTolIn,
               double eigAbsTolIn,
               int    blockSizeIn,
               double solverTolIn,
               int    maxItsIn,
               int    verbosityIn) : numEigs(numEigsIn),
                                     blockSize(blockSizeIn),
                                     solverTol(solverTolIn),
                                     eigRelTol(eigRelTolIn),
                                     eigAbsTol(eigAbsTolIn),
                                     maxIts(maxItsIn),
                                     verbosity(verbosityIn)
{
  assert(numEigs>0);
  assert(blockSize>0);
  assert(eigRelTol>=0.0);
}


LOBPCG::LOBPCG(boost::property_tree::ptree const& opts) : LOBPCG(opts.get<int>("NumEigs"),
                                                                 opts.get("RelativeTolerance",0.0),
                                                                 opts.get("AbsoluteTolerance",0.0),
                                                                 opts.get("BlockSize", 1),
                                                                 opts.get("SolverTolerance", -1.0),
                                                                 opts.get("MaxIts",-1),
                                                                 opts.get("Verbosity",0))
{}

Eigen::MatrixXd LOBPCG::Orthonormalizer::Compute(Eigen::Ref<const Eigen::MatrixXd> const& V)
{
  Eigen::MatrixXd output(V);
  if(B!=nullptr){
    ComputeInPlace(output,B->Apply(V));
  }else{
    ComputeInPlace(output,V);
  }
  return output;
}

Eigen::MatrixXd LOBPCG::Orthonormalizer::Compute(Eigen::Ref<const Eigen::MatrixXd> const& V, Eigen::Ref<const Eigen::MatrixXd> const& BVin)
{
  Eigen::MatrixXd output(V);
  ComputeInPlace(output, BVin);
  return output;
}

bool LOBPCG::Orthonormalizer::ComputeInPlace(Eigen::Ref<Eigen::MatrixXd> V)
{
  if(B!=nullptr){
    return ComputeInPlace(V, B->Apply(V));
  }else{
    return ComputeInPlace(V,V);
  }
}

bool LOBPCG::Orthonormalizer::ComputeInPlace(Eigen::Ref<Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXd> const& BVin)
{
  vDim = V.cols();
  BV = BVin;

  auto cholFact = (V.transpose()*BV).eval().selfadjointView<Eigen::Lower>().llt();
  if(cholFact.info() != Eigen::Success)
    return false;

  VBV_Chol = cholFact.matrixL();
  V = VBV_Chol.triangularView<Eigen::Lower>().solve(V.transpose()).transpose();

  if(B!=nullptr){
    BV = VBV_Chol.triangularView<Eigen::Lower>().solve(BV.transpose()).transpose();
  }else{
    BV = V;
  }

  return true;
}

Eigen::MatrixXd LOBPCG::Orthonormalizer::InverseVBV() const
{
  return VBV_Chol.triangularView<Eigen::Lower>().solve(Eigen::MatrixXd::Identity(vDim,vDim));
}


LOBPCG::Constraints::Constraints(std::shared_ptr<LinearOperator>   const& B,
                                 Eigen::Ref<const Eigen::MatrixXd> const& constMat) : Y(constMat)
{
  if(B!=nullptr){
    BY = B->Apply(constMat);
  }else{
    BY = constMat;
  }

  YBY_llt = (constMat.transpose() * BY).eval().selfadjointView<Eigen::Lower>().llt();
}

void LOBPCG::Constraints::ApplyInPlace(Eigen::Ref<Eigen::MatrixXd> x)
{
    Eigen::MatrixXd YBX = BY.transpose()*x;
    x -= Y*YBY_llt.solve(YBX);
}

LOBPCG& LOBPCG::compute(std::shared_ptr<LinearOperator> const& A,
                        Eigen::MatrixXd                 const& constMat,
                        std::shared_ptr<LinearOperator>        B,
                        std::shared_ptr<LinearOperator>        M)
{

  const int dim = A->rows();
  if(numEigs>int(std::floor(0.2*dim))){
    std::cerr << "ERROR: LOBPCG can fail when more than 20\% of the eigenvalues are computed." << std::endl;
    assert(numEigs>int(std::ceil(0.2*dim)));
  }

  // Initial message
  if(verbosity>0){
    std::cout << "Solving ";
    if(B==nullptr){
      std::cout << "standard ";
    }else{
      std::cout << "generalized ";
    }
    std::cout << "eigenvalue problem." << std::endl;

    if(verbosity>1)
      std::cout << "  matrix size = " << A->rows() << " x " << A->cols() << std::endl;
  }

  // If this is the first call to compute, initialize everything
  reset(dim);

  // Initialize the eigenvalue and eigenvector matrices if they haven't been initialized yet
  if(eigVecs.rows()==0){
    eigVecs.resize(dim, numEigs);
    eigVecs = RandomGenerator::GetUniform(dim,numEigs);
    eigVals.resize(numEigs);

  }else if(eigVecs.cols()>numEigs){
    numEigs = eigVecs.cols();

  }else if(eigVecs.cols()<numEigs){
    unsigned int oldSize = eigVecs.cols();
    eigVecs.conservativeResize(dim, numEigs);
    eigVecs.rightCols(numEigs - oldSize) = RandomGenerator::GetUniform(dim,numEigs-oldSize);
  }

  // Get an estimate of matrix norms
  { // TODO: Move this into a function
    const int numNormTest = std::min(numEigs,5);
    Eigen::MatrixXd gaussMat = RandomGenerator::GetNormal(dim, numNormTest);
    double gaussMatNorm = gaussMat.norm();
    Anorm = A->Apply(gaussMat).norm() / gaussMatNorm;

    if(B!=nullptr){
      Bnorm = B->Apply(gaussMat).norm() / gaussMatNorm;
    }else{
      Bnorm = 1.0;
    }
  }

  Eigen::VectorXd blockVals;
  Eigen::MatrixXd blockVecs;

  Eigen::MatrixXd constraints;
  unsigned int maxBlocks = std::ceil(float(numEigs)/float(blockSize));
  if(constMat.cols()+numEigs-blockSize > 0)
    constraints.resize(dim, constMat.cols()+maxBlocks*blockSize-blockSize);
  if(constMat.cols()>0)
    constraints.leftCols(constMat.cols()) = constMat;

  // Loop over the blocks needed to compute the eigenvalues
  int numComputed;
  for(numComputed=0; numComputed<numEigs; numComputed+=blockSize){

    if(numComputed>0)
      constraints.block(0,numComputed+constMat.cols()-blockSize,dim,blockSize) = eigVecs.block(0,numComputed-blockSize,dim,blockSize);

    // Make sure the eigenvalue and eigenvector matrices have enough space.  Intialize any new vectors to random values
    if(numComputed + blockSize > eigVecs.cols()){
      int numNew = numComputed + blockSize - eigVecs.cols();
      eigVals.conservativeResize(numComputed + blockSize);
      eigVecs.conservativeResize(dim, numComputed + blockSize);
      eigVecs.rightCols(numNew) = RandomGenerator::GetUniform(dim, numNew);
    }

    std::tie(blockVals, blockVecs) = ComputeBlock(A,
                                                  eigVecs.block(0,numComputed,dim,blockSize), // initial guess
                                                  constraints.leftCols(numComputed+constMat.cols()), // Make sure this block produces vectors that are perpendicular to the previously computed vectors
                                                  B,  // RHS matrix
                                                  M); // preconditioner

    eigVals.segment(numComputed, blockSize) = blockVals;
    eigVecs.block(0,numComputed, dim, blockSize) = blockVecs;

    // Check to see if we've capture the dominant eigenvalues
    if(numComputed>0){
      if( (eigVals(numComputed-1)< eigRelTol*eigVals(0))||(eigVals(numComputed-1)<eigAbsTol)){
        eigVals.conservativeResize(numComputed);
        eigVecs = eigVecs.leftCols(numComputed).eval();
        break;
      }
    }

  }

  int numToKeep;
  for(numToKeep=1; numToKeep<numEigs; ++numToKeep){
    if((eigVals(numToKeep)< eigRelTol*eigVals(0))||(eigVals(numToKeep)<eigAbsTol))
      break;
  }
  if(numToKeep<eigVecs.cols()){
    eigVals = eigVals.head(numToKeep).eval();
    eigVecs = eigVecs.leftCols(numToKeep).eval();
  }

  return *this;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> LOBPCG::ComputeBlock(std::shared_ptr<LinearOperator>   const& A,
                                                                 Eigen::Ref<const Eigen::MatrixXd> const& X0,
                                                                 Eigen::Ref<const Eigen::MatrixXd> const& constMat,
                                                                 std::shared_ptr<LinearOperator>          B,
                                                                 std::shared_ptr<LinearOperator>          M)
{
  // Make sure X0 has the right shape
  assert(X0.cols()==blockSize);
  assert(X0.rows()==A->rows());
  assert(A->rows()==A->cols());

  if(B!=nullptr){
    assert(B->rows()==A->rows());
    assert(B->rows()==B->cols());
  }

  if(M!=nullptr){
    assert(M->rows()==A->rows());
    assert(M->rows()==M->cols());
  }

  if(constMat.rows()>0){
    assert(constMat.rows()==A->rows());
  }

  Eigen::MatrixXd X = X0;

  // To start with, all of the indices are active
  unsigned int numActive = blockSize;

  // A vector full of booleans measuring if each eigenvalue is fixed or not
  std::vector<bool> isActive(blockSize,true);

  // Get ready to apply constraints (e.g., build Y, BY, )
  std::shared_ptr<Constraints> consts;
  if(constMat.rows()>0){
    consts = std::make_shared<Constraints>(B, constMat);
    consts->ApplyInPlace(X);
  }

  // Make the current vectors orthogonal wrt B
  Orthonormalizer bOrtho(B);
  bOrtho.ComputeInPlace(X);

  Eigen::MatrixXd BX;
  if(B != nullptr){
    BX = bOrtho.GetBV();
  }else{
    BX = X;
  }

  // Compute the initial Ritz vectors by solving the reduced eigenproblem.
  Eigen::MatrixXd AX = A->Apply(X);

  Eigen::MatrixXd XAX = X.transpose() * AX;
  XAX = 0.5*(XAX+XAX.transpose()).eval();

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ritzSolver(XAX);
  Eigen::VectorXd subEigVals = ritzSolver.eigenvalues();
  Eigen::MatrixXd subEigVecs = ritzSolver.eigenvectors();

  auto swaps = GetSortSwaps(subEigVals);
  SortVec(swaps,subEigVals);
  SortCols(swaps, subEigVecs);


  X = (X*subEigVecs).eval();
  AX = (AX*subEigVecs).eval();
  if(B!=nullptr){
    BX = (BX*subEigVecs).eval();
  }else{
    BX = X;
  }

  //
  Eigen::MatrixXd resids(X.rows(),X.cols());
  Eigen::VectorXd residNorms;
  Eigen::MatrixXd AR, P, BP, AP;
  AR.resize(X.rows(),X.cols());

  Eigen::VectorXd convCrit, oldConvCrit;
  convCrit = Eigen::VectorXd::Constant(X.cols(), std::numeric_limits<double>::infinity());

  for(unsigned it=0; it<maxIts; ++it){

    bool useP = it!=0;

    if(consts!=nullptr)
      consts->ApplyInPlace(X);

    if(verbosity>0)
      std::cout << "Iteration " << it << std::endl;

    // Compute the residuals
    resids.leftCols(numActive) = AX.leftCols(numActive) - BX.leftCols(numActive)*subEigVals.head(numActive).asDiagonal();

    // Compute the norm of the residuals
    residNorms = resids.colwise().norm();

    // Compute the convergence criteria
    oldConvCrit = convCrit;
    convCrit = residNorms.array() / ((Anorm + Bnorm*subEigVals.array())*X.colwise().norm().transpose().array());

    // Figure out how many indices are active
    numActive = 0;
    for(unsigned int i=0; i<residNorms.size(); ++i){

      // Fix a vector if we've converged or haven't improved
      if((convCrit(i)<solverTol)||(std::abs(convCrit(i)-oldConvCrit(i))<1e-14)){
        isActive.at(i) = false;
      }
      numActive += int(isActive.at(i));
    }

    if(verbosity>2)
      std::cout << "  numActive = " << numActive << std::endl;

    if(verbosity>2)
      std::cout << "  eigenvalues = " << subEigVals.transpose() << std::endl;

    if(verbosity>1){
      std::cout << "  residNorms = " << residNorms.transpose() << std::endl;
      std::cout << "  convergence criteria = " << convCrit.transpose() << std::endl;
    }

    // If we've converged for all the vectors, break
    if(numActive==0){
      auto swaps = GetSortSwaps(subEigVals);
      SortVec(swaps, subEigVals);//.head(numActive));
      SortCols(swaps, X);//.leftCols(numActive));

      return std::make_pair(subEigVals, X);
    }


    // Sort everything so the first columns corresponds to the largest residuals
    auto swaps = GetSortSwaps(residNorms,isActive);

    SortVec(swaps, residNorms);
    SortVec(swaps, isActive);
    SortVec(swaps, subEigVals);//.head(numActive));
    SortCols(swaps, X);//.leftCols(numActive));
    SortCols(swaps, AX);//.leftCols(numActive));
    SortCols(swaps, BX);//.leftCols(numActive));
    SortCols(swaps, resids);

    if(it!=0){
      SortCols(swaps, P);//.leftCols(numActive));
      SortCols(swaps, AP);//.leftCols(numActive));
      SortCols(swaps, BP);//.leftCols(numActive));
    }

    // Apply the preconditioner
    if(M!=nullptr)
      resids.leftCols(numActive) = M->Apply(resids.leftCols(numActive));

    // Apply the constraints
    if(consts!=nullptr)
      consts->ApplyInPlace(resids.leftCols(numActive));

    resids.leftCols(numActive) -= (X * BX.transpose() * resids.leftCols(numActive)).eval();

    // Orthonormalize residuals wrt B-inner product
    bOrtho.ComputeInPlace(resids.leftCols(numActive));

    Eigen::MatrixXd BR;
    if(B != nullptr){
      BR = bOrtho.GetBV();
    } else {
      BR = resids.leftCols(numActive);
    }

    AR.leftCols(numActive) = A->Apply(resids.leftCols(numActive));

    if(useP){

      // Orthonormalize P
      bool success = bOrtho.ComputeInPlace(P.leftCols(numActive),BP.leftCols(numActive));
      if(!success){
        useP = false;
      }else{

        AP.leftCols(numActive) = bOrtho.VBV_Chol.triangularView<Eigen::Lower>().solve(AP.leftCols(numActive).transpose()).transpose();

        if(B!=nullptr){
          BP = B->Apply(P);
        }else{
          BP = P;
        }
      } // else if(!success)
    }

    Eigen::MatrixXd XAR = X.transpose()*AR.leftCols(numActive);
    Eigen::MatrixXd RAR = resids.leftCols(numActive).transpose()*AR.leftCols(numActive);
    RAR = 0.5*(RAR+RAR.transpose()).eval();
    Eigen::MatrixXd XAX = X.transpose()*AX;
    XAX = 0.5*(XAX+XAX.transpose()).eval();

    Eigen::MatrixXd XBX, RBR, XBR;

    XBX = BX.transpose() * X;
    XBX = 0.5*(XBX+XBX.transpose()).eval();
    RBR = BR.transpose() * resids.leftCols(numActive);
    XAX = 0.5*(RBR+RBR.transpose()).eval();
    XBR = X.transpose() * BR.leftCols(numActive);

    Eigen::MatrixXd gramA, gramB;

    if(useP){
      gramA.resize(blockSize+2*numActive,blockSize+2*numActive);
      gramB.resize(blockSize+2*numActive,blockSize+2*numActive);

      Eigen::MatrixXd XAP = X.transpose() * AP.leftCols(numActive);
      Eigen::MatrixXd RAP = resids.leftCols(numActive).transpose()*AP.leftCols(numActive);
      Eigen::MatrixXd PAP = P.leftCols(numActive).transpose()*AP.leftCols(numActive);
      PAP = 0.5*(PAP+PAP.transpose()).eval();
      Eigen::MatrixXd XBP = X.transpose()*BP.leftCols(numActive);
      Eigen::MatrixXd RBP = resids.leftCols(numActive).transpose()*BP.leftCols(numActive);

      gramA.block(0,0,blockSize,blockSize) = subEigVals.asDiagonal();
      gramA.block(0,blockSize,blockSize,numActive) = XAR;
      gramA.block(0,blockSize+numActive,blockSize,numActive) = XAP;
      gramA.block(blockSize,0,numActive,blockSize) = XAR.transpose();
      gramA.block(blockSize,blockSize,numActive,numActive) = RAR;
      gramA.block(blockSize,blockSize+numActive,numActive,numActive) = RAP;
      gramA.block(blockSize+numActive,0,numActive,blockSize) = XAP.transpose();
      gramA.block(blockSize+numActive,blockSize,numActive,numActive) = RAP.transpose();
      gramA.block(blockSize+numActive,blockSize+numActive,numActive,numActive) = PAP;

      gramB.block(0,0,blockSize,blockSize) = Eigen::MatrixXd::Identity(blockSize,blockSize);
      gramB.block(0,blockSize,blockSize,numActive) = XBR;
      gramB.block(0,blockSize+numActive,blockSize,numActive) = XBP;
      gramB.block(blockSize,0,numActive,blockSize) = XBR.transpose();
      gramB.block(blockSize,blockSize,numActive,numActive) = Eigen::MatrixXd::Identity(numActive,numActive);
      gramB.block(blockSize,blockSize+numActive,numActive,numActive) = RBP;
      gramB.block(blockSize+numActive,0,numActive,blockSize) = XBP.transpose();
      gramB.block(blockSize+numActive,blockSize,numActive,numActive) = RBP.transpose();
      gramB.block(blockSize+numActive,blockSize+numActive,numActive,numActive) = Eigen::MatrixXd::Identity(numActive,numActive);

    }else{

      gramA.resize(blockSize+numActive,blockSize+numActive);
      gramB.resize(blockSize+numActive,blockSize+numActive);

      gramA.block(0,0,blockSize,blockSize) = subEigVals.asDiagonal();
      gramA.block(0,blockSize,blockSize,numActive) = XAR;
      gramA.block(blockSize,0,numActive,blockSize) = XAR.transpose();
      gramA.block(blockSize,blockSize,numActive,numActive) = RAR;

      gramB.block(0,0,blockSize,blockSize) = Eigen::MatrixXd::Identity(blockSize,blockSize);
      gramB.block(0,blockSize,blockSize,numActive) = XBR;
      gramB.block(blockSize,0,numActive,blockSize) = XBR.transpose();
      gramB.block(blockSize,blockSize,numActive,numActive) = Eigen::MatrixXd::Identity(numActive,numActive);
    }

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ritzSolver(gramA, gramB);
    assert(ritzSolver.info()==Eigen::Success);

    // Get and sort the ritz values
    subEigVals = ritzSolver.eigenvalues().tail(blockSize);
    if(subEigVals.size()>1)
      assert(subEigVals(0)<subEigVals(1));

    subEigVecs = ritzSolver.eigenvectors().rightCols(blockSize);

    swaps = GetSortSwaps(subEigVals);
    SortVec(swaps,subEigVals);
    SortCols(swaps, subEigVecs);

    Eigen::Ref<const Eigen::MatrixXd> eigX = subEigVecs.topRows(blockSize);
    Eigen::Ref<const Eigen::MatrixXd> eigR = subEigVecs.block(blockSize,0,numActive,subEigVecs.cols());

    Eigen::MatrixXd pp  = resids.leftCols(numActive) * eigR;
    Eigen::MatrixXd app = AR.leftCols(numActive)*eigR;
    Eigen::MatrixXd bpp = BR.leftCols(numActive)*eigR;

    if(useP){
      Eigen::Ref<const Eigen::MatrixXd> eigP = subEigVecs.bottomRows(numActive);

      pp  += P.leftCols(numActive)*eigP;
      app += AP.leftCols(numActive)*eigP;
      bpp += BP.leftCols(numActive)*eigP;
    }

    X = (X*eigX + pp).eval();
    AX = (AX*eigX + app).eval();
    BX = (BX*eigX + bpp).eval();

    P = pp;
    AP = app;
    BP = bpp;

  }

  std::cerr << "WARNING: LOBPCG reached maximum iterations." << std::endl;

  // Sort the results so that the eigenvalues are descending
  swaps = GetSortSwaps(subEigVals);
  SortVec(swaps, subEigVals);//.head(numActive));
  SortCols(swaps, X);//.leftCols(numActive));

  return std::make_pair(subEigVals, X);
}


void LOBPCG::InitializeVectors(Eigen::MatrixXd const& X0)
{
  if(eigVecs.rows()>0)
    assert(X0.rows()==eigVecs.rows());

  eigVals.resize(X0.cols());
  eigVecs = X0;
}

LOBPCG& LOBPCG::reset(int dim)
{
  if(solverTol<0.0)
    solverTol = dim * std::sqrt(std::numeric_limits<double>::epsilon());

  if(maxIts<0)
    maxIts = 500;

  return *this;
}
