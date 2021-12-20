#ifndef JACOBIANPIECE_H_
#define JACOBIANPIECE_H_

#include "MUQ/Modeling/ModPiece.h"

namespace muq{
namespace Modeling{

  /**
  @class JacobianPiece
  @ingroup Modeling
  @brief A wrapper around another ModPiece that evaluates the action of the other
  piece's Jacobian on a vector.
  */
  class JacobianPiece : public ModPiece {

  public:
    JacobianPiece(std::shared_ptr<ModPiece> const& basePieceIn,
                  unsigned int              const  outWrt,
                  unsigned int              const  inWrt);

    virtual ~JacobianPiece() = default;

    std::shared_ptr<ModPiece> BasePiece(){return basePiece;};

  protected:

    virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override;

    static Eigen::VectorXi GetInputSizes(std::shared_ptr<ModPiece> const& basePiece,
                                         unsigned int              const outWrt);

    static Eigen::VectorXi GetOutputSizes(std::shared_ptr<ModPiece> const& basePiece,
                                          unsigned int              const inWrt);
                                          
    std::shared_ptr<ModPiece> basePiece;

    const unsigned int outWrt, inWrt;

  }; // class JacobianPiece

} // namespace Modeling
} // namespace muq

#endif
