#include "MUQ/Modeling/LinearAlgebra/GeneralizedEigenSolver.h"

using namespace muq::Modeling;



void GeneralizedEigenSolver::SortVec(std::vector<std::pair<int,int>> const& swapInds,
                                     Eigen::Ref<Eigen::VectorXd>            matrix)
{
  for(auto& swap : swapInds)
    std::swap(matrix(swap.first),matrix(swap.second));
}

void GeneralizedEigenSolver::SortVec(std::vector<std::pair<int,int>> const& swapInds,
                                     std::vector<bool>                    & vec)
{
  for(auto& swap : swapInds)
    std::swap(vec[swap.first],vec[swap.second]);
}

void GeneralizedEigenSolver::SortCols(std::vector<std::pair<int,int>> const& swapInds,
                                      Eigen::Ref<Eigen::MatrixXd>            matrix)
{
  for(auto& swap : swapInds)
    matrix.col(swap.first).swap(matrix.col(swap.second));
}

std::vector<std::pair<int,int>> GeneralizedEigenSolver::GetSortSwaps(Eigen::Ref<const Eigen::VectorXd> const& residNorms)
{
  std::vector<bool> allActive(residNorms.size(), true);
  return GetSortSwaps(residNorms, allActive);
}

std::vector<std::pair<int,int>> GeneralizedEigenSolver::GetSortSwaps(Eigen::Ref<const Eigen::VectorXd> const& resids,
                                                                     std::vector<bool>                 const& isActive)
{
  Eigen::VectorXd newResids = resids;
  std::vector<bool> newIsActive = isActive;

  const unsigned int size = resids.size();

  std::vector<std::pair<int,int>> swaps;
  unsigned int maxInd;

  // First, put all the active indices at the begining
  int firstInactiveInd = 0;
  while((firstInactiveInd<isActive.size())&&(newIsActive.at(firstInactiveInd)))
    firstInactiveInd++;

  int lastActiveInd = isActive.size()-1;
  while((lastActiveInd>=0)&&(!newIsActive.at(lastActiveInd)))
    lastActiveInd--;

  while(firstInactiveInd<lastActiveInd){

    swaps.push_back(std::make_pair(firstInactiveInd,lastActiveInd));
    std::swap(newIsActive.at(firstInactiveInd), newIsActive.at(lastActiveInd));
    std::swap(newResids(firstInactiveInd), newResids(lastActiveInd));

    while((firstInactiveInd<isActive.size())&&(newIsActive.at(firstInactiveInd)))
      firstInactiveInd++;

    while((lastActiveInd>=0)&&(!newIsActive.at(lastActiveInd)))
      lastActiveInd--;
  }

  unsigned int numActive = lastActiveInd+1;

  // Now, sort all of the active residuals according to magnitude
  for(unsigned int i=0; i<numActive; ++i)
  {
    // Find the maximum index
    maxInd = std::distance(newResids.data(), std::max_element(&newResids(i), newResids.data()+numActive));

    // Swap indices if needed
    if(maxInd!=i){
      swaps.push_back( std::make_pair(i, maxInd) );
      std::swap(newResids(i), newResids(maxInd));
    }
  }

  return swaps;
}
