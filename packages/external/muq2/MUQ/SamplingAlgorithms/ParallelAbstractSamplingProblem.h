#ifndef PARALLELABSTRACTSAMPLINGPROBLEM_H_
#define PARALLELABSTRACTSAMPLINGPROBLEM_H_

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif

#include <parcer/Communicator.h>
#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/ParallelFlags.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
     * @brief Parallelization layer supporting parallel models underlying a SamplingProblem.
     *
     * @details This class passes a LogDensity call to all other ranks of the given communicator,
     * initiated by rank zero.
     *
     * Each other rank receives a copy of the parameter vector to work with. This allows computing
     * a model (e.g. a PDE) underlying the LogDensity in parallel across all these ranks.
     * It is assumed that rank 0 returns the exact LogDensity value. The actual model parallelization
     * of course lies in the responsiblity of the user, since parallelization strategies vary extremely between
     * different models.
     *
     * This class is used by the controller rank responsible for handling the MCMC chains in
     * a distributed model MCMC setting. To the controller, the ParallelAbstractSamplingProblem looks
     * like a sequential SamplingProblem. The controller rank handles MCMC proposals, MCMC kernels and sample
     * storage while the other ranks stand by, only joining the LogDensity computations.
     */
    class ParallelAbstractSamplingProblem : public AbstractSamplingProblem
    {
    public:

      ParallelAbstractSamplingProblem(std::shared_ptr<parcer::Communicator> comm, int id, std::shared_ptr<AbstractSamplingProblem> abstractSamplingProblem)
       : AbstractSamplingProblem(abstractSamplingProblem->blockSizes, abstractSamplingProblem->blockSizesQOI)
      {
        this->comm = comm;
        this->id = id;
        this->abstractSamplingProblem = abstractSamplingProblem;
      }

      virtual ~ParallelAbstractSamplingProblem() = default;

      virtual double LogDensity(std::shared_ptr<SamplingState> const& state) override {

        if (comm->GetRank() == 0) {
          for (int dest = 1; dest < comm->GetSize(); dest++) {
            comm->Send(ControlFlag::LOGDENSITY, dest, WorkgroupTag);
            comm->Send(id, dest, WorkgroupTag);
            comm->Send(state->state[0], dest, WorkgroupTag);
          }
        }

        return abstractSamplingProblem->LogDensity(state);
      }

      virtual std::shared_ptr<SamplingState> QOI() override {
        return abstractSamplingProblem->QOI();
      }

      std::shared_ptr<AbstractSamplingProblem> GetSequentialProblem() {
        return abstractSamplingProblem;
      }

    private:
      int id;
      std::shared_ptr<parcer::Communicator> comm;
      std::shared_ptr<AbstractSamplingProblem> abstractSamplingProblem;
    };

  }
}

#endif

#endif
