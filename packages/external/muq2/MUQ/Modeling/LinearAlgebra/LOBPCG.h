#ifndef LOBPCG_H
#define LOBPCG_H


#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/GeneralizedEigenSolver.h"

#include <boost/property_tree/ptree_fwd.hpp>
#include <Eigen/Core>

namespace muq{
namespace Modeling{

  /** @class LOBPCG
      @ingroup LinearOperators
      @brief The Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG) method for matrix-free computation of eigenvalues and eigenvectors.
      @details This class solves generalized eigenvalue problems of the form \f$Av = \lambda Bv\f$ when the matrix \f$A\f$ is symmetric and the matrix \f$B\f$ is symmetric positive definite.  It uses the The Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG) method described in
      "TOWARD THE OPTIMAL PRECONDITIONED EIGENSOLVER: LOCALLY OPTIMAL BLOCK PRECONDITIONED CONJUGATE GRADIENT METHOD" by ANDREW V. KNYAZEV.
   */
  class LOBPCG : public GeneralizedEigenSolver
  {
  public:

    /**
    @param[in] numEigsIn The maximum number of eigenvalues to compute used as a stopping criteria.
    @param[in] blockSizeIn The number of eigenvalues to compute at once
    @param[in] eigTolIn Fraction of the largest eigenvalue used as stopping criteria.
    @param[in] solverTolIn Solver tolerance
    @param[in] maxItsIn The maximum number of iterations the solver is allowed to take.
    @param[in] largestIn A boolean specifying if the largest (true) or smallest (false) eigenvalues should be computed.
    @param[in] verbosityIn An integer {0,1,2,3} specifying how much information the solver should print.  When verbosityIn=0 (default), nothing is printed.
    */
    LOBPCG(int    numEigsIn,
           double eigRelTolIn=0.0,
           double eigAbsTolIn=0.0,
           int    blockSizeIn=-1,
           double solverTolIn=-1,
           int    maxItsIn=-1,
           int    verbosityIn=0);

    /**
      <B>Configuration Parameters:</B>
      Parameter Key  | Type | Default Value | Description |
      -------------  | ------------- | ------------- | ------------- |
      "NumEigs"      | integer       | -             | The maximum number of eigenvalues to compute. Used as stopping criteria. |
      "RelativeTolerance" | double        | 0.0           | Fraction of the largest eigenvalue used as stopping criteria. |
      "AbsoluteTolerance" | double        | 0.0           | Value of smallest eigenvalue to compute.  Used as stopping criteria. |
      "BlockSize"    | integer       | 1       | How many eigenvalues and eigenvectors should be computed simultaneously. |
      "SolverTolerance" | double        | \f$d\sqrt{\epsilon}\f$ | Tolerance used for stopping criteria.  |
      "MaxIts"      | integer       | min(d,20)   | Maximum number of iterations to take. |
      "Verbosity"   | integer       | 0  |  Specifies how much is printed to the screen.  Valid values are 0 (print nothing), 1, 2, or 3       |
    */
    LOBPCG(boost::property_tree::ptree const& options);

    virtual ~LOBPCG() = default;

    /**
    Compute the generalized eigenvalues and eigenvectors of a symmetric system \f$Av = \lambda Bv\f$.  If compute has been previously called, the eigenvalues and eigenvectors from the previous call will be used as an initial guess for this call to the solver.  This can speed up the solver if A changes slightly and the eigenvalues need to be recomputed (as in DILI MCMC).

    @param[in] A A shared pointer to a muq::Modeling::LinearOperator instance specifying the matrix A.
    @param[in] B A shared pointer to a muq::Modeling::LinearOperator instance specifying the matrix B.  If not specified or set to nullptr, B will be set to the identity.
    @param[in] M An optional preconditioner.  If specified, the LinearOperator M should approximate the inverse of A.
    */
    LOBPCG& compute(std::shared_ptr<LinearOperator> const& A,
                    std::shared_ptr<LinearOperator>        B=nullptr,
                    std::shared_ptr<LinearOperator>        M=nullptr){return compute(A,Eigen::MatrixXd(),B,M);};

    /**
    Compute the generalized eigenvalues and eigenvectors of a symmetric system \f$Av = \lambda Bv\f$ that are in the B-orthogonal complement to some constraints \f$Y\f$.  More precisely, the computed eigenvectors \f$v_i\v$ will satisfy \f$Y^T B v_i=0\f$.
    @param[in] A A shared pointer to a muq::Modeling::LinearOperator instance specifying the matrix A.
    @param[in] constMat A matrix defining the constraints.  This matrix should have the same number of rows as A.
    @param[in] B A shared pointer to a muq::Modeling::LinearOperator instance specifying the matrix B.  If not specified or set to nullptr, B will be set to the identity.
    @param[in] M An optional preconditioner.  If specified, the LinearOperator M should approximate the inverse of A.
    */
    LOBPCG& compute(std::shared_ptr<LinearOperator> const& A,
                    Eigen::MatrixXd                 const& constMat,
                    std::shared_ptr<LinearOperator>        B = nullptr,
                    std::shared_ptr<LinearOperator>        M = nullptr);


    /** Initializes the eigenvectors to a specific value.
    */
    void InitializeVectors(Eigen::MatrixXd const& vecs);

    /** Resets the current eigenvalues and eigenvectors so that a random initial
        guess is used during the next call to compute instead of reusing the
        previously computed eigenvalues and eigenvectors as initial guesses.
        @param[in] dim The dimension of the problem
     */
    LOBPCG& reset(int dim);

  private:

    /** Computes a single block of eigenvalues and eigenvectors that are orthogonal to the columns of constmat.
    */
    std::pair<Eigen::VectorXd, Eigen::MatrixXd> ComputeBlock(std::shared_ptr<LinearOperator>   const& A,
                                                             Eigen::Ref<const Eigen::MatrixXd> const& X0,
                                                             Eigen::Ref<const Eigen::MatrixXd> const& constMat,
                                                             std::shared_ptr<LinearOperator>          B,
                                                             std::shared_ptr<LinearOperator>          M);

    /** Makes the columns of a matrix V orthonormal wrt the B inner product \f$v^T B v\f$. */
    class Orthonormalizer{
    public:
      Orthonormalizer(std::shared_ptr<LinearOperator> const& Bin) : B(Bin){};

      bool ComputeInPlace(Eigen::Ref<Eigen::MatrixXd> V);
      bool ComputeInPlace(Eigen::Ref<Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXd> const& BVin);

      Eigen::MatrixXd Compute(Eigen::Ref<const Eigen::MatrixXd> const& V);
      Eigen::MatrixXd Compute(Eigen::Ref<const Eigen::MatrixXd> const& V, Eigen::Ref<const Eigen::MatrixXd> const& BVin);

      Eigen::MatrixXd InverseVBV() const;
      Eigen::MatrixXd const& GetBV() const{return BV;};
      Eigen::MatrixXd& GetBV(){return BV;};

      int vDim;
      std::shared_ptr<LinearOperator> B;
      Eigen::MatrixXd BV;
      Eigen::MatrixXd VBV_Chol;
    };

    class Constraints{
    public:
      Constraints(std::shared_ptr<LinearOperator>   const& B,
                  Eigen::Ref<const Eigen::MatrixXd> const& constVec);

      void ApplyInPlace(Eigen::Ref<Eigen::MatrixXd> x);

      int size() const{return BY.cols();};

    private:
      Eigen::MatrixXd BY;
      Eigen::Ref<const Eigen::MatrixXd> const& Y;
      Eigen::LLT<Eigen::MatrixXd> YBY_llt;

    }; // class LOBPCG::Constraints

    /// The maximum number of eigenvalues to compute
    int numEigs;

    /// Number of eigenvalues to compute simultaneously
    int blockSize;

    /// Solver tolerance
    double solverTol;

    /// Fraction of largest eigenvalue used to terminate
    double eigRelTol, eigAbsTol;

    /// Maximum number of iterations
    int maxIts;

    /// Controls how much information we want to print
    int verbosity;

    //Eigen::VectorXd eigVals;
    //Eigen::MatrixXd eigVecs;

    //std::shared_ptr<LinearOperator> A, B, M;
    double Anorm, Bnorm;

  }; // class LOBPCG
}
}



#endif // LOBPCG_H
