#ifndef GAUSSIANOPERATOR_H
#define GAUSSIANOPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/Distributions/GaussianBase.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"

#include <memory>

namespace muq
{
namespace Modeling
{


/** @class GaussianOperator
 *  @ingroup LinearAlgebra
 *  @brief Creates a linear operator for the action of the covariance or precision
    matrix of a Gaussian distribution.
    @seealso GaussNewtonOperator, HessianOperator
 */
class GaussianOperator : public LinearOperator {
public:

  /**
  @params[in] gaussIn A shared_ptr to the GaussianBase instance defining the
                      precision or covariance that we want to apply.
  @params[in] precOrCovIn A flag deciding whether the precision or covariance of
                          the Gaussian should be used.  Valid options are
                          Gaussian::Precision or Gaussian::Covariance.
  */
  GaussianOperator(std::shared_ptr<GaussianBase> const& gaussIn,
                   Gaussian::Mode                       precOrCovIn);

  virtual ~GaussianOperator() = default;

  /** Apply the linear operator to a vector */
  virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

  /** Apply the transpose of the linear operator to a vector. */
  virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

protected:
  std::shared_ptr<GaussianBase> gauss;
  Gaussian::Mode precOrCov;

};

} // namespace Modeling
} // namespace MUQ



#endif
