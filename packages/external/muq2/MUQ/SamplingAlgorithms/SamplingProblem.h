#ifndef SAMPLINGPROBLEM_H_
#define SAMPLINGPROBLEM_H_

// include Density and not ModPiece so that if a SamplingProblem is constructed with a Density the compiler knows it is a child of ModPiece
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
    @ingroup SamplingAlgorithms
    @class SamplingProblem
    @brief Class for sampling problems based purely on a density function.
    */
    class SamplingProblem : public AbstractSamplingProblem{
    public:

      /**
	     @param[in] target The target distribution
       */
      SamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& targetIn);

      /**
	     @param[in] target The target distribution
	     @param[in] qoi Quantity of interest associated with model
       */
      SamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& targetIn,
                      std::shared_ptr<muq::Modeling::ModPiece> const& qoiIn);

      virtual ~SamplingProblem() = default;


      virtual double LogDensity(std::shared_ptr<SamplingState> const& state) override;

      virtual Eigen::VectorXd GradLogDensity(std::shared_ptr<SamplingState> const& state,
                                             unsigned                       const  blockWrt) override;

      std::shared_ptr<muq::Modeling::ModPiece> GetDistribution(){return target;};

      virtual std::shared_ptr<SamplingState> QOI() override;

    protected:

      /// The target distribution (the prior in the inference case)
      std::shared_ptr<muq::Modeling::ModPiece> target;

      std::shared_ptr<muq::Modeling::ModPiece> qoi;

    private:

      static unsigned GetNumBlocks(std::shared_ptr<muq::Modeling::ModPiece> const& target);
      static std::vector<int> GetBlockSizes(std::shared_ptr<muq::Modeling::ModPiece> const& target);

      std::shared_ptr<SamplingState> lastState;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
