#ifndef DENSITY_H
#define DENSITY_H

#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/ModPiece.h"

namespace muq{
  namespace Modeling{

    class DensityBase : public Distribution, public ModPiece{
    public:
      DensityBase(Eigen::VectorXi const& inputSizes);

      virtual ~DensityBase() = default;

    protected:

      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void GradientImpl(unsigned int                const  outputDimWrt,
                                unsigned int                const  inputDimWrt,
                                ref_vector<Eigen::VectorXd> const& input,
                                Eigen::VectorXd             const& sensitivity) override;

      virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                                unsigned int                const  inputDimWrt,
                                ref_vector<Eigen::VectorXd> const& input) override;

      virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                     unsigned int                const  inputDimWrt,
                                     ref_vector<Eigen::VectorXd> const& input,
                                     Eigen::VectorXd             const& vec) override;

     virtual void ApplyHessianImpl(unsigned int                const  outWrt,
                                   unsigned int                const  inWrt1,
                                   unsigned int                const  inWrt2,
                                   ref_vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd             const& sens,
                                   Eigen::VectorXd             const& vec) override;

    };


    class Density : public DensityBase{

    public:
      Density(std::shared_ptr<Distribution> distIn);

      virtual ~Density() = default;

      // Return the distribution this density is built from
      virtual std::shared_ptr<Distribution> GetDistribution(){ return dist;};

    protected:
      std::shared_ptr<Distribution> dist;

      virtual double LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual Eigen::VectorXd GradLogDensityImpl(unsigned int wrt,
                                                 ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual Eigen::VectorXd ApplyLogDensityHessianImpl(unsigned int                const  inWrt1,
                                                         unsigned int                const  inWrt2,
                                                         ref_vector<Eigen::VectorXd> const& input,
                                                         Eigen::VectorXd             const& vec) override;

      virtual Eigen::VectorXd SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      static Eigen::VectorXi GetInputSizes(std::shared_ptr<Distribution> distIn);

    }; // class Density

  } // namespace Modeling
} // namespace muq



#endif // #ifndef DENSITY_H
