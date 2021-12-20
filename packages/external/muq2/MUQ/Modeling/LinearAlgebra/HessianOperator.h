#ifndef HESSIANOPERATOR_H
#define HESSIANOPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

#include <memory>

namespace muq
{
namespace Modeling
{


/** @class HessianOperator
 *  @ingroup LinearAlgebra
 *  @brief Creates a linear operator for the action of the Hessian of a ModPiece
           on a vector.  Useful for computing the Hessian spectrum with iteratives solvers like LOBPCG.
    @seealso GaussNewtonOperator
 */
class HessianOperator : public LinearOperator {
public:

  HessianOperator(std::shared_ptr<ModPiece>    const& pieceIn,
                  std::vector<Eigen::VectorXd> const& inputsIn,
                  unsigned int                        outWrtIn,
                  unsigned int                        inWrt1In,
                  unsigned int                        inWrt2In,
                  Eigen::VectorXd              const& sensIn,
                  double                              scaleIn=1.0,
                  double                              nuggetIn=0.0);

  virtual ~HessianOperator() = default;

  /** Apply the linear operator to a vector */
  virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

  /** Apply the transpose of the linear operator to a vector. */
  virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

protected:
  std::shared_ptr<ModPiece> basePiece;

  const std::vector<Eigen::VectorXd> inputs;

  const unsigned int outWrt;
  const unsigned int inWrt1;
  const unsigned int inWrt2;
  const Eigen::VectorXd sens;
  const double scale;
  const double nugget;

};

} // namespace Modeling
} // namespace MUQ



#endif
