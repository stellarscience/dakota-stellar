#ifndef ONESTEPCACHEPIECE_H
#define ONESTEPCACHEPIECE_H

#include "MUQ/Modeling/ModPiece.h"


namespace muq{
  namespace Modeling{
    class OneStepCachePiece : public ModPiece {
    public:
      OneStepCachePiece(std::shared_ptr<ModPiece> baseModPiece, const double& prec = Eigen::NumTraits<double>::dummy_precision());

      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override;

      double HitRatio();

    private:
      unsigned int hits = 0;
      unsigned int misses = 0;
      bool firstEvaluation = true;
      std::vector<Eigen::VectorXd> lastInput;
      std::vector<Eigen::VectorXd> lastOutputs;

      std::shared_ptr<ModPiece> baseModPiece;
      const double& prec;
    };
  }
}


#endif