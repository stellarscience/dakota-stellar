#ifndef GENERALIZEDEIGENSOLVER_H
#define GENERALIZEDEIGENSOLVER_H


#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

#include <boost/property_tree/ptree_fwd.hpp>
#include <Eigen/Core>

namespace muq{
namespace Modeling{

  /** @class GeneralizedEigenSolver
      @brief Abstract base class for operator based generalized eigenvalue solvers.
      @details Base class for algorithms that solve problems of the form
      \f[
      Av = \lambda B v
      \f]
      @seealso LOBPCG
      @seealso StochasticEigenSolver
  */
  class GeneralizedEigenSolver {

  public:

    virtual ~GeneralizedEigenSolver() = default;


    // /**
    // Compute the generalized eigenvalues and eigenvectors of a symmetric system \f$Av = \lambda Bv\f$ .
    // @param[in] A A shared pointer to a muq::Modeling::LinearOperator instance specifying the matrix A.
    // @param[in] B A shared pointer to a muq::Modeling::LinearOperator instance specifying the matrix B.  If not specified or set to nullptr, B will be set to the identity.
    // */
    // virtual GeneralizedEigenSolver& compute(std::shared_ptr<LinearOperator> const& A,
    //                                         std::shared_ptr<LinearOperator>        B = nullptr) = 0;

    /** Return a reference to the computed vector of eigenvalues.  The vector
        will only be valid after calling compute.
    */
    Eigen::VectorXd const& eigenvalues() const{return eigVals;}

    /** Return a matrix whose columns contain the computed eigenvectors.  The
        matrix will only be valid after calling compute.
    */
    Eigen::MatrixXd const& eigenvectors() const{return eigVecs;};

  protected:

    /**
    Sorts the columns of the matrix using precomputed swaps from the GetSortSwaps function.
    */
    static void SortCols(std::vector<std::pair<int,int>> const& swapInds,
                         Eigen::Ref<Eigen::MatrixXd>            matrix);

    static void SortVec(std::vector<std::pair<int,int>> const& swapInds,
                             Eigen::Ref<Eigen::VectorXd>       matrix);

    static void SortVec(std::vector<std::pair<int,int>> const& swapInds,
                        std::vector<bool>                    & vec);
    /**
    Returns a vector of swaps needed to sort the provided matrix using a selection sort.
    */
    static std::vector<std::pair<int,int>> GetSortSwaps(Eigen::Ref<const Eigen::VectorXd> const& residNorms,
                                                        std::vector<bool>                 const& isActive);

    static std::vector<std::pair<int,int>> GetSortSwaps(Eigen::Ref<const Eigen::VectorXd> const& residNorms);
    
    Eigen::VectorXd eigVals;
    Eigen::MatrixXd eigVecs;

    std::shared_ptr<LinearOperator> A, B, M;
  }; // class GeneralizedEigenSolver

}
}

#endif // #ifndef GENERALIZEDEIGENSOLVER_H
