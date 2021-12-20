#ifndef MIMCMC_H
#define MIMCMC_H

#include <boost/property_tree/ptree.hpp>

#include "MUQ/SamplingAlgorithms/MIMCMCBox.h"
#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

namespace pt = boost::property_tree;

/** @defgroup MIMCMC
    @ingroup MCMC
    @brief Tools for defininig and running multilevel and multiindex MCMC algorithms.
    @details Multiindex MCMC methods are built on an entire arbitriry-dimensional
    grid of sampling problems, in contrast to classical MCMC only sampling from
    a single distribution.

    In order to be effective, MIMCMC methods require distributions closely related to each
    other, where evaluating the coarser ones should be significantly computationally cheaper.
    A typical example are models based on the Finite Element method, where strongly varying mesh
    resolutions lead to a hierarchy of models. Then, coarser models are a
    good approximation of finer ones while being far cheaper to compute.

    In analogy to the mathematical definition, running a MIMCMC method requires defining
    a grid of models (as well as other needed components) via MIComponentFactory.
    Most importantly, how proposals will be drawn from coarser chains has to be defined
    as well as how to combine them (via MIInterpolation) with a proposal for the current
    fine chain.

    Refer to the MCMC examples for complete code examples.
*/

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief Multiindex MCMC method.
        @details A basic MIMCMC method based on a fixed
        number of samples for all model indices.
     */
    class MIMCMC : public SamplingAlgorithm {
    public:
      MIMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory);

      virtual std::shared_ptr<SampleCollection> GetSamples() const override;
      virtual std::shared_ptr<SampleCollection> GetQOIs() const override;

      Eigen::VectorXd MeanQOI();

      /**
       * @brief Access an MIMCMCBox
       *
       * @param index Model index for which to retrieve the box.
       * @return std::shared_ptr<MIMCMCBox> The MIMCMCBox representing the Multiindex telescoping sum component
       * associated with that model.
       */
      std::shared_ptr<MIMCMCBox> GetBox(std::shared_ptr<MultiIndex> index);

      /**
       * @brief Evaluate parameter mean via MI telescoping sum.
       */
      Eigen::VectorXd MeanParam();

      /**
       * @brief Draw MI structure (mostly debugging purposes)
       */
      void Draw(bool drawSamples = true);

      std::shared_ptr<MIMCMCBox> GetMIMCMCBox(std::shared_ptr<MultiIndex> index);

      /**
       * @brief Get set of indices of boxes set up by the method.
       */
      std::shared_ptr<MultiIndexSet> GetIndices();

      virtual std::shared_ptr<SampleCollection> RunImpl(std::vector<Eigen::VectorXd> const& x0) override;

      /**
       * @brief Write HDF5 output for the entire MIMCMC method
       */
      void WriteToFile(std::string filename);

    private:
      pt::ptree pt;
      std::shared_ptr<MultiIndexSet> gridIndices;
      std::shared_ptr<MIComponentFactory> componentFactory;
      std::vector<std::shared_ptr<MIMCMCBox>> boxes;

      std::string multiindexToConfigString (std::shared_ptr<MultiIndex> index);
    };

  }
}

#endif
