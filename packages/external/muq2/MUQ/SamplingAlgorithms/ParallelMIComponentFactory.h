#ifndef PARALLELMICOMPONENTFACTORY_H_
#define PARALLELMICOMPONENTFACTORY_H_

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif

#include "spdlog/spdlog.h"
#include "MUQ/SamplingAlgorithms/ParallelAbstractSamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
     * @brief Wrapper for MIComponentFactory supporting parallel model setup
     * @details This MIComponentFactory provides parallelized model setup
     * initiated by rank zero of the given communicator and makes those models
     * available to be controlled by the root rank. This allows rank zero
     * to set up SaplingProblem's through this factory exactly as if it was sequential;
     * this class ensures that the other worker ranks perform the same SamplingProblem
     * construction at the same time. In particular, parallel models are guaranteed
     * to be set up in sync, so most external model codes should just work once the appropriate
     * communicator is passed to them. Further, once SamplingProblem's are set up, worker processes
     * listen for LogDensity requests from the respective ParallelAbstractSamplingProblem.
     */
    class ParallelMIComponentFactory : public MIComponentFactory {

    public:

      ParallelMIComponentFactory (std::shared_ptr<parcer::Communicator> comm, std::shared_ptr<parcer::Communicator> global_comm, std::shared_ptr<MIComponentFactory> componentFactory)
       : comm(comm), global_comm(global_comm), componentFactory(componentFactory)
      {

        if (comm->GetRank() != 0) {
          while (true) {
            spdlog::trace("Parallel factory rank {} waiting...", comm->GetRank());
            ControlFlag command = comm->Recv<ControlFlag>(0, WorkgroupTag);
            spdlog::trace("Parallel factory rank {} received command", comm->GetRank(), command);
            if (command == ControlFlag::FINALIZE) {
              samplingProblems.clear(); // Tear down models synchronously
              comm->Barrier();
              spdlog::trace("Parallel factory rank {} passed finalize barrier", comm->GetRank());
              break;
            }
            if (command == ControlFlag::INIT_PROBLEM) {
              auto index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(0, WorkgroupTag));
              int id = comm->Recv<int>(0, WorkgroupTag);
              spdlog::trace("Parallel factory rank {} building model index {}", comm->GetRank(), *index);
              samplingProblems[id] = componentFactory->SamplingProblem(index);//std::make_shared<MySamplingProblem>(index, comm, id, measurements);
            }
            else if (command == ControlFlag::LOGDENSITY) {
              int id = comm->Recv<int>(0, WorkgroupTag);
              auto state = std::make_shared<SamplingState>(comm->Recv<Eigen::VectorXd>(0, WorkgroupTag));
              samplingProblems[id]->LogDensity(state);
            }
            else if (command == ControlFlag::TEST) {
              int id = comm->Recv<int>(0, WorkgroupTag);
              auto state = std::make_shared<SamplingState>(comm->Recv<Eigen::VectorXd>(0, WorkgroupTag));

              double density;
              std::cerr << "Not implemented!!!" << std::endl;
            }
            else if (command == ControlFlag::QOI) {
              int id = comm->Recv<int>(0, WorkgroupTag);

              samplingProblems[id]->QOI();
            } else {
              std::cerr << "Unexpected command!" << std::endl;
              exit(43);
            }
          }
        }
      }
      virtual bool IsInverseProblem() override {
        return componentFactory->IsInverseProblem();
      }

      /**
       * @brief Stops worker processes command loop, freeing them for other tasks.
       *
       */
      void finalize() {
        if (comm->GetRank() != 0)
          return;
        spdlog::trace("Parallel factory sending finalize");
        for (int dest = 1; dest < comm->GetSize(); dest++)
          comm->Send(ControlFlag::FINALIZE, dest, WorkgroupTag);
        samplingProblems.clear(); // Tear down models synchronously
        comm->Barrier();
        spdlog::trace("Parallel factory finalized");
      }

      virtual std::shared_ptr<MCMCProposal> Proposal (std::shared_ptr<MultiIndex> const& index, std::shared_ptr<AbstractSamplingProblem> const& samplingProblem) override {
        return componentFactory->Proposal(index, samplingProblem);
      }

      virtual std::shared_ptr<MultiIndex> FinestIndex() override {
        return componentFactory->FinestIndex();
      }

      virtual std::shared_ptr<MCMCProposal> CoarseProposal (std::shared_ptr<MultiIndex> const& index,
                                                            std::shared_ptr<AbstractSamplingProblem> const& coarseProblem,
                                                            std::shared_ptr<SingleChainMCMC> const& coarseChain) override {
        return componentFactory->CoarseProposal(index, coarseProblem, coarseChain);
      }

      virtual std::shared_ptr<AbstractSamplingProblem> SamplingProblem (std::shared_ptr<MultiIndex> const& index) override {
        idcnt++;
        if (comm->GetRank() == 0) {
          spdlog::debug("Rank {} requesting model {} from parallel factory", comm->GetRank(), *index);
          for (int dest = 1; dest < comm->GetSize(); dest++) {
            comm->Send(ControlFlag::INIT_PROBLEM, dest, WorkgroupTag);
            comm->Send(*index, dest, WorkgroupTag);
            comm->Send(idcnt, dest, WorkgroupTag);
          }
        }
        samplingProblems[idcnt] = std::make_shared<ParallelAbstractSamplingProblem>(comm, idcnt, componentFactory->SamplingProblem(index));
        return samplingProblems[idcnt];
      }

      virtual std::shared_ptr<MIInterpolation> Interpolation (std::shared_ptr<MultiIndex> const& index) override {
        return componentFactory->Interpolation(index);
      }

      virtual Eigen::VectorXd StartingPoint (std::shared_ptr<MultiIndex> const& index) override {
        return componentFactory->StartingPoint(index);
      }

    private:
      int idcnt = 0;
      std::shared_ptr<parcer::Communicator> comm;
      std::shared_ptr<parcer::Communicator> global_comm;
      std::shared_ptr<MIComponentFactory> componentFactory;

      std::map<int, std::shared_ptr<AbstractSamplingProblem>> samplingProblems;

    };

  }
}

#endif

#endif
