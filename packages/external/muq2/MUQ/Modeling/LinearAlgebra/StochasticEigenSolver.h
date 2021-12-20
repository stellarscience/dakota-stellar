#ifndef STOCHASTICEIGENSOLVER_H
#define STOCHASTICEIGENSOLVER_H

#include "MUQ/Modeling/LinearAlgebra/GeneralizedEigenSolver.h"
#include <boost/property_tree/ptree_fwd.hpp>

namespace muq{
namespace Modeling{

/** @class StochasticEigenSolver
    @brief Two-pass stochastic algorithm for computing generalized eigenvalues from matrix products.
    @details This class implements the two-pass stochastic algorithm outline in
     [Saibaba et al., 2010] and [Villa et al., 2019].  While the standard two-pass
     algorithm computes a fixed number of eigenvalues and eigenvectors $r$, this
     implementation will continue to add new samples and recompute the decomposition
     until the computed eigenvalues satisfy one of three stopping criteria:
     1. The number of computed eigenvalues is greater than some threshold "NumEigs"
     2. The smallest eigenvalue is smaller than some fraction 'RelativeTolerance' of the maximum eigenvalue.
     3. The smallesst eigenvalues is smaller than some absolute tolerance 'AbsoluteTolerance'.

     The values of these stopping criteria are passed to the constructor of this
     class in a boost property_tree.

     REFERENCES:
     - Saibaba, Kitanidis, Lee, (2010) "Randomized algorithms for Generalized Hermitian Eigenvalue Problems with application to computing Karhunen-Loeve expansion." NUMERICAL LINEAR ALGEBRA WITH APPLICATIONS
     - Villa, Petra, Ghattas, (2019) "hIPPYlib:  AN EXTENSIBLE SOFTWARE FRAMEWORK FORLARGE-SCALE INVERSE PROBLEMS GOVERNED BY PDES;PART I: DETERMINISTIC INVERSION AND LINEARIZEDBAYESIAN INFERENCE." ACM TOMS 2019
*/
class StochasticEigenSolver : public GeneralizedEigenSolver {

public:

  StochasticEigenSolver(int    numEigsIn,
                        double eigRelTolIn=0.0,
                        double eigAbsTolIn=0.0,
                        int    expectedRankIn=-1,
                        int    samplingFactorIn=-1,
                        int    blockSize=10,
                        int    verbosityIn=0);
  /**
  <B>Configuration Parameters:</B>
  Parameter Key  | Type | Default Value | Description |
  -------------  | ------------- | ------------- | ------------- |
  "NumEigs"      | integer       | -             | The maximum number of eigenvalues to compute. Used as stopping criteria. |
  "RelativeTolerance" | double        | 0.0           | Fraction of the largest eigenvalue used as stopping criteria. |
  "AbsoluteTolerance" | double        | 0.0           | Value of smallest eigenvalue to compute.  Used as stopping criteria. |
  "ExpectedRank"    | integer       | 0.1*dim       | The expected number of eigenvalues that are larger than the the tolerances.  This is used to decide on the initial number of samples to compute. |
  "OversamplingFactor" | integer        | 0.5*expected rank | When computing a rank $r$ eigenvalue decomposition, $r+l$ matvecs are used, where $l$ is the oversampling factor.  |
  "BlockSize"   | integer       | 10 | The number of new random matvecs to compute when the current eigenvalues do not satisfy the stopping criteria. |
  "Verbosity"   | integer       | 0  |  Specifies how much is printed to the screen.  Valid values are 0 (print nothing), 1, 2, or 3       |
  */
  StochasticEigenSolver(boost::property_tree::ptree const& options);

  /** Runs the two pass algorithm to compute the (generalized) eigenvalues and eigenvectors
      of the matrix \f$A\f$.
  */
  virtual StochasticEigenSolver& compute(std::shared_ptr<LinearOperator> const& A,
                                         std::shared_ptr<LinearOperator>        B = nullptr,
                                         std::shared_ptr<LinearOperator>        Binv = nullptr);


private:

  /** Implements the `PreCholQR` algorithm of [Villa et al., 2019], which computes
      a QR decompsotion of the matrix \f$Y\f$ such that the columns of \f$Q\f$ are
      B-orthogonal, i.e., \f$Q^T B Q=I\f$.
  */
  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> CholeskyQR(Eigen::MatrixXd                 const& Y,
                                                         std::shared_ptr<LinearOperator> const& B) const;

  int numEigs;
  double eigRelTol;
  double eigAbsTol;
  int    expectedRank;
  int    samplingFactor;
  int    blockSize;
  int    verbosity;
};

} // namespace Modeling
} // namespace muq

#endif // #ifndef STOCHASTICEIGENSOLVER_H
