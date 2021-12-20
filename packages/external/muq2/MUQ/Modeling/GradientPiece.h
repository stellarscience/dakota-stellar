#ifndef GRADIENTPIECE_H_
#define GRADIENTPIECE_H_

#include "MUQ/Modeling/ModPiece.h"

namespace muq{
namespace Modeling{

  class GradientPiece : public ModPiece {

  public:
    GradientPiece(std::shared_ptr<ModPiece> const& basePieceIn,
                  unsigned int              const  outWrt,
                  unsigned int              const  inWrt);

    virtual ~GradientPiece() = default;

    std::shared_ptr<ModPiece> BasePiece(){return basePiece;};

  protected:

    virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override;

    virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                   unsigned int                const  inputDimWrt,
                                   ref_vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd             const& vec) override;

    static Eigen::VectorXi GetInputSizes(std::shared_ptr<ModPiece> const& basePiece,
                                         unsigned int              const outWrt);

    static Eigen::VectorXi GetOutputSizes(std::shared_ptr<ModPiece> const& basePiece,
                                          unsigned int              const inWrt);

    std::shared_ptr<ModPiece> basePiece;

    const unsigned int outWrt, inWrt;

  }; // class GradientPiece

} // namespace Modeling
} // namespace muq



#endif
