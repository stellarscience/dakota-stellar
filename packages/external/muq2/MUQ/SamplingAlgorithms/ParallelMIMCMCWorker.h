#ifndef PARALLELMIMCMCWORKER_H_
#define PARALLELMIMCMCWORKER_H_

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif

#include <chrono>
#include <list>
#include <thread>
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include "MUQ/SamplingAlgorithms/MarkovChain.h"
#include "MUQ/SamplingAlgorithms/DistributedCollection.h"
#include "MUQ/SamplingAlgorithms/ParallelFlags.h"
#include "MUQ/SamplingAlgorithms/ParallelizableMIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/ParallelMIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/ParallelMIMCMCBox.h"
#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/Cereal/MultiIndexSerializer.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
     * @brief High-level communication wrapper for controlling SampleCollectors.
     * @details This takes care about assigning workers to a set of collectors for a specific
     * model index/level, sending commands to them and finally unassigning them again.
     */
    class CollectorClient {
    public:
      CollectorClient(std::shared_ptr<parcer::Communicator> comm, std::vector<int> subgroup, std::shared_ptr<MultiIndex> modelindex)
       : comm(comm), subgroup(subgroup), boxHighestIndex(modelindex) {
         for (int dest : subgroup) {
           comm->Send(ControlFlag::ASSIGN_COLLECTOR, dest, ControlTag);
           comm->Send(subgroup, dest, ControlTag);
           comm->Send(*boxHighestIndex, dest, ControlTag);
         }

         // Set up Multiindex box
         boxLowestIndex = MultiIndex::Copy(boxHighestIndex);
         --(*boxLowestIndex);
         std::shared_ptr<MultiIndex> boxSize = std::make_shared<MultiIndex>(*boxHighestIndex - *boxLowestIndex);
         boxIndices = MultiIndexFactory::CreateFullTensor(boxSize->GetVector());
      }

      std::shared_ptr<MultiIndex> GetModelIndex() const {
        return MultiIndex::Copy(boxHighestIndex);
      }

      void CollectSamples (int numSamples) {
        spdlog::debug("Kick off collection of {} samples", numSamples);
        sampling = true;
        int samplesAssigned = 0;
        for (int i = 0; i < subgroup.size(); i++) {
          int dest = subgroup[i];
          comm->Send(ControlFlag::SAMPLE_BOX, dest, ControlTag);
          int assigning = (numSamples - samplesAssigned) / (subgroup.size() - i); // Ensure mostly even assignment of samples
          comm->Send(assigning, dest, ControlTag);
          samplesAssigned += assigning;
        }
        spdlog::debug("Collection kick off done");
      }

      void ComputeMeans() {
        computingMeans = true;
        for (int dest : subgroup) {
          comm->Send(ControlFlag::MEANS, dest, ControlTag);
        }
      }

      void WriteToFile(std::string filename) {
        for (int dest : subgroup) {
          comm->Send(ControlFlag::WRITE_TO_FILE, dest, ControlTag);
          comm->Send(filename, dest, ControlTag);
          bool sync = comm->Recv<bool>(dest, ControlTag);
        }
      }

      bool Receive (ControlFlag command, const MPI_Status& status) {
        if (status.MPI_SOURCE != subgroup[0])
          return false;

        if (command == ControlFlag::SAMPLE_BOX_DONE) {
          sampling = false;
        } else if (command == ControlFlag::MEANS_DONE) {

          for (uint i = 0; i < boxIndices->Size(); i++) {
            std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];

            Eigen::VectorXd chainSampleMean = comm->Recv<Eigen::VectorXd>(status.MPI_SOURCE, ControlTag);
            Eigen::VectorXd chainQOIMean = comm->Recv<Eigen::VectorXd>(status.MPI_SOURCE, ControlTag);

            std::shared_ptr<MultiIndex> index = std::make_shared<MultiIndex>(*boxLowestIndex + *boxIndex);
            auto indexDiffFromTop = std::make_shared<MultiIndex>(*boxHighestIndex - *index);

            //std::cout << "Contribution of model " << *index << ": " << chainQOIMean << std::endl;
            if (i == 0) {
              if (indexDiffFromTop->Sum() % 2 == 0) {
                boxQOIMean = chainQOIMean;
              } else {
                boxQOIMean = -chainQOIMean;
              }
            } else {
              if (indexDiffFromTop->Sum() % 2 == 0) {
                boxQOIMean += chainQOIMean;
              } else {
                boxQOIMean -= chainQOIMean;
              }
            }
          }
          computingMeans = false;
        } else {
          std::cerr << "Unexpected command!" << std::endl;
          exit(43);
        }

        return true;
      }

      Eigen::VectorXd GetQOIMean() {
        return boxQOIMean;
      }

      void Unassign() {
        for (int dest : subgroup) {
          comm->Send(ControlFlag::UNASSIGN, dest, ControlTag);
        }
      }

      bool IsSampling() {
        return sampling;
      }
      bool IsComputingMeans() {
        return computingMeans;
      }

    private:

      std::shared_ptr<parcer::Communicator> comm;
      std::vector<int> subgroup;

      bool sampling = false;
      bool computingMeans = false;
      std::shared_ptr<MultiIndexSet> boxIndices;
      Eigen::VectorXd boxQOIMean;

      std::shared_ptr<MultiIndex> boxHighestIndex;
      std::shared_ptr<MultiIndex> boxLowestIndex;

      std::map<std::shared_ptr<MultiIndex>, Eigen::VectorXd, MultiPtrComp> means;

    };

    /**
     * @brief High-level communication wrapper for controlling worker processes.
     * @details This takes care about assigning workers to a worker group for a specific
     * model index/level, sending commands to them and finally unassigning them again.
     */
    class WorkerClient {
    public:
      WorkerClient(std::shared_ptr<parcer::Communicator> comm, std::shared_ptr<PhonebookClient> phonebookClient, int RootRank)
       : comm(comm), phonebookClient(phonebookClient) {
      }


      void assignGroup (std::vector<int> subgroup, std::shared_ptr<MultiIndex> modelindex) {
        for (int dest : subgroup) {
          comm->Send(ControlFlag::ASSIGN, dest, ControlTag);
          comm->Send(subgroup, dest, ControlTag);
          comm->Send(*modelindex, dest, ControlTag);
        }
        phonebookClient->Register(modelindex, subgroup[0]);
      }

      std::vector<int> UnassignGroup (std::shared_ptr<MultiIndex> modelIndex, int groupRootRank) {
        spdlog::trace("UnRegister {}", groupRootRank);
        phonebookClient->UnRegister(modelIndex, groupRootRank);
        spdlog::trace("Sending unassign to {}", groupRootRank);
        comm->Ssend(ControlFlag::UNASSIGN, groupRootRank, ControlTag);
        std::vector<int> groupMembers = comm->Recv<std::vector<int>>(groupRootRank, ControlTag);
        return groupMembers;
      }

      void UnassignAll() {
        std::shared_ptr<MultiIndex> largest = nullptr;
        do {
          largest = phonebookClient->LargestIndex();
          spdlog::trace("Unassigning model {}", *largest);
          std::vector<int> ranks = phonebookClient->GetWorkgroups(largest);
          for (int rank : ranks) {
            UnassignGroup(largest, rank);
          }
        } while (largest->Max() != 0);
      }

      void Finalize() {
        for (int dest = 1; dest < comm->GetSize(); dest++)
          comm->Send(ControlFlag::QUIT, dest, ControlTag);
      }

    private:
      std::shared_ptr<parcer::Communicator> comm;
      std::shared_ptr<PhonebookClient> phonebookClient;
    };

    /**
     * @brief Implements the actual sampling / collecting logic for parallel MIMCMC.
     * @details Workers will, in a loop until finalized, wait for instructions to
     * join a set of collectors or a worker group for sampling. They then listen
     * for commands to execute these tasks until unassigned from that task again.
     * As a collector, workers will request MCMC samples (more specifically,
     * samples and their coarser ancesters in analogy to the sequential MIMCMCBox)
     * from sampling worker groups. As part of a sampling worker group, they will
     * compute MCMC samples and provide them to other processes via the phonebook.
     */
    class WorkerServer {
    public:
      WorkerServer(boost::property_tree::ptree const& pt, std::shared_ptr<parcer::Communicator> comm, std::shared_ptr<PhonebookClient> phonebookClient, int RootRank, std::shared_ptr<ParallelizableMIComponentFactory> componentFactory, std::shared_ptr<muq::Utilities::OTF2TracerBase> tracer) {

        while (true) {
          ControlFlag command = comm->Recv<ControlFlag>(RootRank, ControlTag);
          if (command == ControlFlag::ASSIGN) {
            std::vector<int> subgroup_proc = comm->Recv<std::vector<int>>(0, ControlTag);
            auto samplingProblemIndex = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(0, ControlTag));

            spdlog::trace("Setting up MPI subcommunicator for group");
            MPI_Group world_group;
            MPI_Comm_group (MPI_COMM_WORLD, &world_group);
            MPI_Group subgroup;
            MPI_Group_incl (world_group, subgroup_proc.size(), &subgroup_proc[0], &subgroup);

            MPI_Comm subcomm_raw;
            MPI_Comm_create_group(MPI_COMM_WORLD, subgroup, ControlTag, &subcomm_raw);
            auto subcomm = std::make_shared<parcer::Communicator>(subcomm_raw);

            componentFactory->SetComm(subcomm);
            spdlog::trace("Setting up ParallelMIComponentFactory");
            auto parallelComponentFactory = std::make_shared<ParallelMIComponentFactory>(subcomm, comm, componentFactory);

            if (subcomm->GetRank() == 0) {
              std::cout << "Subgroup root is global " << comm->GetRank() << std::endl;
              auto finestProblem = parallelComponentFactory->SamplingProblem(parallelComponentFactory->FinestIndex());

              spdlog::trace("Setting up ParallelMIMCMCBox");
              auto box = std::make_shared<ParallelMIMCMCBox>(pt, parallelComponentFactory, samplingProblemIndex, comm, phonebookClient, tracer);

              spdlog::debug("Rank {} begins sampling", comm->GetRank());
              const int subsampling = pt.get<int>("MLMCMC.Subsampling");
              tracer->enterRegion(TracerRegions::Sampling);
              for (int i = 0; i < 2 + subsampling; i++) // TODO: Really subsampling on every level? Maybe subsample when requesting samples?
                box->Sample();
              tracer->leaveRegion(TracerRegions::Sampling);
              phonebookClient->SetWorkerReady(samplingProblemIndex, comm->GetRank());
              spdlog::trace("Awaiting instructions");

              //Dune::Timer timer_idle;
              //Dune::Timer timer_full;
              while (true) {
                MPI_Status status;
                //timer_idle.start();
                command = comm->Recv<ControlFlag>(MPI_ANY_SOURCE, ControlTag, &status);
                //timer_idle.stop();
                if (command == ControlFlag::UNASSIGN) {
                  comm->Send<std::vector<int>>(subgroup_proc, status.MPI_SOURCE, ControlTag);
                  break;
                } else if (command == ControlFlag::SAMPLE) {
                  spdlog::trace("Send sample from {} to rank {}", comm->GetRank(), status.MPI_SOURCE);

                  auto sampleCollection = box->FinestChain()->GetSamples(); // TODO: last() function for collection? // TODO: Do not store chains here
                  auto latestSample = sampleCollection->at(sampleCollection->size()-1);
                  // TODO: Send "full" sample via parcer?
                  comm->Send<Eigen::VectorXd>(latestSample->state[0], status.MPI_SOURCE, ControlTag);
                  comm->Send<double>(AnyCast(latestSample->meta["LogTarget"]), status.MPI_SOURCE, ControlTag);
                  if (latestSample->HasMeta("QOI")) {
                    std::shared_ptr<SamplingState> qoi = AnyCast(latestSample->meta["QOI"]);
                    comm->Send<Eigen::VectorXd>(qoi->state[0], status.MPI_SOURCE, ControlTag);
                  } else {
                    std::cerr << "No QOI!" << sampleCollection->size() << std::endl;
                    exit(47);
                  }


                  spdlog::trace("Sampling");
                  tracer->enterRegion(TracerRegions::Sampling);
                  for (int i = 0; i < 1 + subsampling; i++) // TODO: Really subsampling on every level? Maybe subsample when requesting samples?
                    box->Sample();
                  tracer->leaveRegion(TracerRegions::Sampling);
                  phonebookClient->SetWorkerReady(samplingProblemIndex, comm->GetRank());
                } else if (command == ControlFlag::SAMPLE_BOX) {
                  spdlog::trace("Send box from {} to rank {}", comm->GetRank(), status.MPI_SOURCE);

                  for (int i = 0; i < box->NumChains(); i++) {
                    auto sampleCollection = box->GetChain(i)->GetSamples(); // TODO: last() function for collection? // TODO: Do not store chains here
                    auto latestSample = sampleCollection->at(sampleCollection->size()-1);
                    // TODO: Send "full" sample via parcer?
                    comm->Send<Eigen::VectorXd>(latestSample->state[0], status.MPI_SOURCE, ControlTag);
                    if (latestSample->HasMeta("QOI")) {
                      std::shared_ptr<SamplingState> qoi = AnyCast(latestSample->meta["QOI"]);
                      comm->Send<Eigen::VectorXd>(qoi->state[0], status.MPI_SOURCE, ControlTag);
                    } else {
                      std::cerr << "No QOI!" << sampleCollection->size() << std::endl;
                      exit(47);
                    }
                  }

                  tracer->enterRegion(TracerRegions::Sampling);
                  for (int i = 0; i < 5; i++)
                    box->Sample();
                  tracer->leaveRegion(TracerRegions::Sampling);
                  phonebookClient->SetWorkerReady(samplingProblemIndex, comm->GetRank());
                } else {
                  std::cerr << "Unexpected command!" << std::endl;
                  exit(43);
                }

              }
              //std::cout << "Worker Controller " << comm->GetRank() << " idle time:\t" << timer_idle.elapsed() << " of:\t" << timer_full.elapsed() << std::endl;

            }

            parallelComponentFactory->finalize();
            spdlog::trace("Rank {} finalized", comm->GetRank());
          } else if (command == ControlFlag::ASSIGN_COLLECTOR) {
            std::vector<int> subgroup_proc = comm->Recv<std::vector<int>>(RootRank, ControlTag);
            auto boxHighestIndex = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(RootRank, ControlTag));

            // Set up subcommunicator
            MPI_Group world_group;
            MPI_Comm_group (MPI_COMM_WORLD, &world_group);
            MPI_Group subgroup;
            MPI_Group_incl (world_group, subgroup_proc.size(), &subgroup_proc[0], &subgroup);

            MPI_Comm subcomm_raw;
            MPI_Comm_create_group(MPI_COMM_WORLD, subgroup, ControlTag, &subcomm_raw);
            auto subcomm = std::make_shared<parcer::Communicator>(subcomm_raw);


            // Set up Multiindex box
            auto boxLowestIndex = MultiIndex::Copy(boxHighestIndex);
            --(*boxLowestIndex);
            std::shared_ptr<MultiIndex> boxSize = std::make_shared<MultiIndex>(*boxHighestIndex - *boxLowestIndex);
            std::shared_ptr<MultiIndexSet> boxIndices = MultiIndexFactory::CreateFullTensor(boxSize->GetVector());

            std::vector<std::shared_ptr<DistributedCollection>> sampleCollections(boxIndices->Size());
            std::vector<std::shared_ptr<DistributedCollection>> qoiCollections(boxIndices->Size());
            for (uint i = 0; i < boxIndices->Size(); i++) {
              auto sampleCollection = std::make_shared<MarkovChain>();
              sampleCollections[i] = std::make_shared<DistributedCollection>(sampleCollection, subcomm);
              auto qoiCollection = std::make_shared<MarkovChain>();
              qoiCollections[i] = std::make_shared<DistributedCollection>(qoiCollection, subcomm);
            }


            //Dune::Timer timer_idle;
            //Dune::Timer timer_full;
            while (true) {
              MPI_Status status;
              //timer_idle.start();
              command = comm->Recv<ControlFlag>(MPI_ANY_SOURCE, ControlTag, &status);
              //timer_idle.stop();
              if (command == ControlFlag::UNASSIGN)
                break;
              tracer->enterRegion(TracerRegions::CollectorBusy);
              if (command == ControlFlag::SAMPLE_BOX) {

                int numSamples = comm->Recv<int>(0, ControlTag);
                spdlog::debug("Collecting {} samples for box {}", numSamples, *boxHighestIndex);

            		//std::cout << "Rank " << comm->GetRank() << " requesting sample from rank " << remoteRank << std::endl;

                for (int i = 0; i < numSamples; i++) {
                  spdlog::trace("Requesting sample box for model {}", *boxHighestIndex);
                  if (i % (numSamples / 10) == 0)
                    spdlog::debug("Collected {} out of {} samples for model {}", i, numSamples, *boxHighestIndex);
                  int remoteRank = phonebookClient->Query(boxHighestIndex, boxHighestIndex, false);
                  comm->Send(ControlFlag::SAMPLE_BOX, remoteRank, ControlTag); // TODO: Receive sample in one piece?
                  for (uint i = 0; i < boxIndices->Size(); i++) {
                    //std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];
                    sampleCollections[i]->Add(std::make_shared<SamplingState>(comm->Recv<Eigen::VectorXd>(remoteRank, ControlTag)));
                    qoiCollections[i]->Add(std::make_shared<SamplingState>(comm->Recv<Eigen::VectorXd>(remoteRank, ControlTag)));
                		//Eigen::VectorXd remoteQOI = comm->Recv<Eigen::VectorXd>(remoteRank, ControlTag);
                  }
                }

                if (subcomm->GetRank() == 0)
                  comm->Send(ControlFlag::SAMPLE_BOX_DONE, RootRank, ControlTag); // TODO: Receive sample in one piece?
                //  box->Sample();
              } else if (command == ControlFlag::MEANS) {
                std::list<Eigen::VectorXd> sampleMeans;
                std::list<Eigen::VectorXd> qoiMeans;
                for (uint i = 0; i < boxIndices->Size(); i++) {
                  sampleMeans.push_back(sampleCollections[i]->GlobalMean());
                  qoiMeans.push_back(qoiCollections[i]->GlobalMean());
                }
                if (subcomm->GetRank() == 0) {
                  comm->Send(ControlFlag::MEANS_DONE, RootRank, ControlTag);
                  auto qoiMean = qoiMeans.begin();
                  for (auto sampleMean = sampleMeans.begin(); sampleMean != sampleMeans.end(); sampleMean++) {
                    comm->Send(*sampleMean, RootRank, ControlTag);
                    comm->Send(*qoiMean, RootRank, ControlTag);
                    qoiMean++;
                  }
                }
              } else if (command == ControlFlag::WRITE_TO_FILE) {
                std::string filename = comm->Recv<std::string>(status.MPI_SOURCE, ControlTag);
                for (uint i = 0; i < boxIndices->Size(); i++) {
                  std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];
                  sampleCollections[i]->WriteToFile(filename, "/Collector_model" + boxHighestIndex->ToString() + "_subchain_" + boxIndex->ToString() + "_samples_rank_" + std::to_string(subcomm->GetRank()));
                  qoiCollections[i]->WriteToFile(filename, "/Collector_model" + boxHighestIndex->ToString() + "_subchain_" + boxIndex->ToString() + "_qois_rank_" + std::to_string(subcomm->GetRank()));
                }
                comm->Send(true, status.MPI_SOURCE, ControlTag);
              } else {
                std::cerr << "Unexpected command!" << std::endl;
                exit(43);
              }
              tracer->leaveRegion(TracerRegions::CollectorBusy);

            }
            //std::cout << "Collector " << comm->GetRank() << " idle time:\t" << timer_idle.elapsed() << " of:\t" << timer_full.elapsed() << std::endl;


          } else if (command == ControlFlag::QUIT) {
            spdlog::trace("Rank {} quit", comm->GetRank());
            break;
          } else {
            std::cerr << "Unexpected command!" << std::endl;
            exit(42);
          }
        }

      }

    };

  }
}

#endif

#endif
