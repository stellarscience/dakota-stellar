#include "MUQ/Modeling/LinearAlgebra/StochasticEigenSolver.h"

#include "MUQ/Utilities/RandomGenerator.h"
#include <boost/property_tree/ptree.hpp>

using namespace muq::Modeling;
using namespace muq::Utilities;

StochasticEigenSolver::StochasticEigenSolver(int    numEigsIn,
                                             double eigRelTolIn,
                                             double eigAbsTolIn,
                                             int    expectedRankIn,
                                             int    samplingFactorIn,
                                             int    blockSizeIn,
                                             int    verbosityIn) : numEigs(numEigsIn),
                                                                   eigRelTol(eigRelTolIn),
                                                                   eigAbsTol(eigAbsTolIn),
                                                                   expectedRank((expectedRankIn<0)?numEigsIn:expectedRankIn),
                                                                   samplingFactor((samplingFactorIn<0)?(0.1*numEigsIn):samplingFactorIn),
                                                                   blockSize(blockSizeIn),
                                                                   verbosity(verbosityIn)
{
  assert(numEigs>0);
  assert(eigRelTol>=0.0);
}

StochasticEigenSolver::StochasticEigenSolver(boost::property_tree::ptree const& opts) : StochasticEigenSolver(opts.get<int>("NumEigs"),
                                                                                                              opts.get("RelativeTolerance",0.0),
                                                                                                              opts.get("AbsoluteTolerance",0.0),
                                                                                                              opts.get("ExpectedRank", -1),
                                                                                                              opts.get("OversamplingFactor", -1),
                                                                                                              opts.get("BlockSize",10),
                                                                                                              opts.get("Verbosity",0))
{}


StochasticEigenSolver& StochasticEigenSolver::compute(std::shared_ptr<LinearOperator> const& A,
                                                      std::shared_ptr<LinearOperator>        B,
                                                      std::shared_ptr<LinearOperator>        Binv)
{
    assert((B!=nullptr)==(Binv!=nullptr));

    const int dim = A->cols();

    Eigen::MatrixXd randMat; // <- Will hold random draws from standard normal
    Eigen::MatrixXd Y; // <- will hold B^{-1}*A*randMat
    bool hasConverged = false;

    randMat = RandomGenerator::GetNormal(dim, expectedRank+samplingFactor);

    if(B!=nullptr){
      Y = Binv->Apply(A->Apply(randMat));
    }else{
      Y =  A->Apply(randMat);
    }

    int it=0;
    while(!hasConverged){

      Eigen::MatrixXd Q,R;
      std::tie(Q,R) = CholeskyQR(Y, B);

      Eigen::MatrixXd T = Q.transpose()*A->Apply(Q);

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolver(T);
      eigVals = eigSolver.eigenvalues();

      auto swaps = GetSortSwaps(eigVals);
      SortVec(swaps,eigVals);

      // Check for convergence
      double smallestVal = eigVals(eigVals.size()-1);
      double largestVal = eigVals(0);

      // Convergence because smallestVal is less than absolute tolerance
      if(smallestVal<eigAbsTol){
        if(verbosity>0){
          std::cout << "Converged: Minimum eigenvalue (" << smallestVal << ") satisfies absolute tolerance (" << eigAbsTol << ")." << std::endl;
        }
        hasConverged = true;
      }

      // Convergence because smallestVal is less than relative tolerance
      if(smallestVal<eigRelTol*largestVal){
        if(verbosity>0){
          std::cout << "Converged: Minimum eigenvalue (" << smallestVal << ") satisfies relative tolerance (" << eigRelTol*largestVal << ")." << std::endl;
        }
        hasConverged = true;
      }

      // Convergence because all eigenvalues of operator have been found
      if(Q.cols()<Y.cols()){
        if(verbosity>0){
          std::cout << "Converged: All nonzero eigenvalues have been found (likely) or samples are degenerate (unlikely)." << std::endl;
        }
        hasConverged = true;
      }

      // Converge because the maximum number of eigenvalues has been found
      if(eigVals.size()>=numEigs){
        if(verbosity>0){
          std::cout << "Converged: Reached maximum number of eigenvalues." << std::endl;
        }
        hasConverged = true;
      }

      if(verbosity>1){
        std::cout << "After iteration " << it << ", " << eigVals.size() << " eigenvalues in [" << smallestVal << "," << largestVal << "]" << std::endl;
        std::cout << "  Y.shape = " << Y.rows() << " x " << Y.cols() << std::endl;
      }
      it++;

      // If we haven't found enough eigenvalues yet, add to the random matrix
      if(!hasConverged){
        randMat.conservativeResize(Eigen::NoChange, randMat.cols()+blockSize);
        Y.conservativeResize(Eigen::NoChange, Y.cols()+blockSize);

        randMat.rightCols(blockSize) = RandomGenerator::GetNormal(dim,blockSize);

        if(B!=nullptr){
          Y.rightCols(blockSize) = Binv->Apply(A->Apply(randMat.rightCols(blockSize)));
        }else{
          Y.rightCols(blockSize) =  A->Apply(randMat.rightCols(blockSize));
        }
      }

      // If we've converged, set the eigenvectors
      if(hasConverged){
        eigVecs = Q*eigSolver.eigenvectors();
        SortCols(swaps, eigVecs);
      }

    }

    return *this;
}


std::pair<Eigen::MatrixXd, Eigen::MatrixXd> StochasticEigenSolver::CholeskyQR(Eigen::MatrixXd                 const& Y,
                                                                              std::shared_ptr<LinearOperator> const& B) const
{

  auto qr = Y.colPivHouseholderQr();
  unsigned int rank = qr.rank();

  Eigen::MatrixXd Z = qr.householderQ().setLength(qr.nonzeroPivots()) * Eigen::MatrixXd::Identity(Y.rows(),rank);
  Eigen::MatrixXd Ry = qr.matrixR().topLeftCorner(rank, rank).template triangularView<Eigen::Upper>();

  if(B==nullptr){
    return std::make_pair(Z,Ry);
  }

  auto chol = (Z.transpose()*B->Apply(Z)).llt();

  Eigen::MatrixXd Q = chol.matrixL().solve(Z.transpose()).transpose();
  Eigen::MatrixXd R = chol.matrixL()*Ry;

  return std::make_pair(Q,R);
}
