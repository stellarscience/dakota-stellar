#ifndef PARALLELIZABLEMICOMPONENTFACTORY_H_
#define PARALLELIZABLEMICOMPONENTFACTORY_H_

#if MUQ_HAS_MPI

#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"

#if !MUQ_HAS_PARCER
#error
#endif


namespace muq {
  namespace SamplingAlgorithms {

    /**
     * @brief MIComponentFactory subclass allowing a communicator to be passed.
     * @details This allows passing a Communicator to processes responsible for
     * setting up parallel SamplingProblems. Needed for parallel MIMCMC and
     * related methods. Implementations will typically pass the given communicator
     * to the model implementations it provides.
     */
    class ParallelizableMIComponentFactory : public MIComponentFactory {
    public:
      virtual void SetComm(std::shared_ptr<parcer::Communicator> const& comm) = 0;
    };

  }
}

#endif

#endif
