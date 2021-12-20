#ifndef PHONEBOOK_H_
#define PHONEBOOK_H_

#if MUQ_HAS_MPI


#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string>
#include <map>
#include <mpi.h>
#include "spdlog/fmt/ostr.h"
#include "spdlog/spdlog.h"
#include <parcer/Communicator.h>
#include "MUQ/SamplingAlgorithms/ParallelFlags.h"
#include "MUQ/Utilities/Cereal/MultiIndexSerializer.h"
#include <deque>
#include <chrono>
#include "MUQ/Utilities/OTF2Tracer.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
     * @brief Phonebook implementation facilitating communication between different roles in a parallel MIMCMC type
     * method. Also is responsible for dynamic scheduling.
     * @details The phonebook process tracks which worker processes are assigned to which models. When samples are
     * computed, the phonebook is informed by the respective worker group controller. When samples are requested
     * by other controllers (in order to be used as proposals for finer chains) or by collectors, the phonebook
     * sends to them the rank of a controller process that has such a sample available.
     *
     * Since the phonebook can also track which model indices/levels are busy and which are not, it can also
     * take care of initiating the dynamic rescheduling of worker processes as needed.
     */
    class PhonebookServer {
    public:

      PhonebookServer(std::shared_ptr<parcer::Communicator> comm,
                      bool scheduling_active = true,
                      std::shared_ptr<muq::Utilities::OTF2TracerBase> tracer = std::make_shared<muq::Utilities::OTF2TracerDummy>())
        : comm(comm), scheduling_active(scheduling_active), tracer(tracer)
      {
      }

      void Run() {
        //Dune::Timer timer_idle;
        //Dune::Timer timer_full;


        while (true) {
          MPI_Status status;
          //timer_idle.start();


          if (scheduling_active && !rescheduling_in_progress) {
            for ( auto &indexWorkerListPair : phonebook ) {
              std::shared_ptr<MultiIndex> index = indexWorkerListPair.first;
              WorkerList& workerList = indexWorkerListPair.second;
              bool work_run_out = workerList.NumWorkersReady() == workerList.NumWorkers() && getNumQueuedTasksForIndex(index, workerList) == 0.0;

              if (workerList.NormalizedRegisteredReadyCounter() > .5 || work_run_out || workerList.recheck) {

                if (workerList.NormalizedRegisteredReadyCounter() > .5) {
                  workerList.ResetTimer();
                }

                if (workerList.NumWorkers() <= 1) // Only reschedule if we have more than 1 worker on this model
                  continue;

                if (workerList.NumWorkersReady() == 0) {
                  workerList.recheck = true;
                  continue;
                }
                workerList.recheck = false;

                spdlog::debug("Timer triggered for {}, idle fraction {}, {} workers on that model, work run out {}", *indexWorkerListPair.first, workerList.GetIdleFraction(), workerList.NumWorkers(), work_run_out);

                //if (index->GetValue(0) == 1) // FIXME: Fix rescheduling model 1
                //  continue;

                double my_load = getLoadFactor(index, workerList);
                double largest_others_load = .0;
                std::shared_ptr<MultiIndex> most_loaded_index = nullptr;
                for ( auto &indexWorkerListPair : phonebook ) {
                  std::shared_ptr<MultiIndex> index = indexWorkerListPair.first;
                  WorkerList& workerList = indexWorkerListPair.second;
                  double others_load = getLoadFactor(index, workerList);
                  spdlog::debug("Load on model {}: {}", *index, others_load);
                  if (others_load > largest_others_load) {
                    largest_others_load = others_load;
                    most_loaded_index = index;
                  }
                }
                assert (most_loaded_index != nullptr);

                if (*index == *most_loaded_index)
                  continue;

                if (work_run_out || my_load + .5 / (double)workerList.NumWorkers() + .01 < largest_others_load) {

                  spdlog::debug("Reassigning from model {} to model {}, {} ready here", *index, *most_loaded_index, workerList.NumWorkersReady());

                  const int RootNode = 0; // Send to root

                  int rescheduleRank = workerList.GetWorkersReady()[0];
                  UnRegister(index, rescheduleRank); // Already unregister this rank so it won't be offered for sampling anymore while being rescheduled!

                  comm->Send(ControlFlag::SCHEDULING_NEEDED, RootNode, ControlTag);
                  comm->Send(*index, RootNode, ControlTag);
                  comm->Send(rescheduleRank, RootNode, ControlTag);
                  comm->Send(*most_loaded_index, RootNode, ControlTag);

                  rescheduling_in_progress = true;
                  break; // Phonebook needs to be ready for further communication after rescheduling a process, so don't reschedule another one right now
                }

              }

            }
          }

          ControlFlag command = comm->Recv<ControlFlag>(MPI_ANY_SOURCE, ControlTag, &status);
          //timer_idle.stop();

          for ( auto &indexWorkerListPair : phonebook ) {
            WorkerList& workerList = indexWorkerListPair.second;
            workerList.tick();
          }



          tracer->enterRegion(TracerRegions::PhonebookBusy);

          if (command == ControlFlag::GET_WORKGROUP) {
            auto requested_sample_index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(status.MPI_SOURCE, ControlTag));
            auto request_source_index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(status.MPI_SOURCE, ControlTag));
            bool high_priority = comm->Recv<bool>(status.MPI_SOURCE, ControlTag);
            if (high_priority)
              requests.push_front(SampleRequest{.requestedSampleIndex = requested_sample_index, .sourceModelIndex = request_source_index, .sourceMPIRank = status.MPI_SOURCE, .highPriority = high_priority});
            else
              requests.push_front(SampleRequest{.requestedSampleIndex = requested_sample_index, .sourceModelIndex = request_source_index, .sourceMPIRank = status.MPI_SOURCE, .highPriority = high_priority});
          } else if (command == ControlFlag::SCHEDULING_DONE) {
            rescheduling_in_progress = false;
          } else if (command == ControlFlag::SCHEDULING_STOP) {
            rescheduling_in_progress = true;
          } else if (command == ControlFlag::SET_WORKGROUP) {
            auto index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(status.MPI_SOURCE, ControlTag));
            int rank = comm->Recv<int>(status.MPI_SOURCE, ControlTag, &status);

            if (!phonebook.count(index)) {
              phonebook[index] = WorkerList();
            }
            phonebook[index].AddWorker(rank);
            //comm->Send(ControlFlag::HANDSHAKE, status.MPI_SOURCE, ControlTag);
            spdlog::trace("Phonebook entry for {} set", *index);
          } else if (command == ControlFlag::UNSET_WORKGROUP) {
            auto index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(status.MPI_SOURCE, ControlTag));
            int rank = comm->Recv<int>(status.MPI_SOURCE, ControlTag, &status);

            UnRegister(index, rank);
          } else if (command == ControlFlag::GET_WORKGROUPS) {
            auto index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(status.MPI_SOURCE, ControlTag));

            if (!phonebook.count(index)) {
              std::cerr << "getting workers for nonexistent model!" << std::endl;
            }
            spdlog::trace("Getting workers from phonebook map");
            std::vector<int> sendvec = phonebook[index].GetWorkers();
            spdlog::trace("Sending {} workgroups", sendvec.size());
            comm->Send(sendvec, status.MPI_SOURCE, ControlTag);
          } else if (command == ControlFlag::GET_LARGEST_INDEX) {
            if (phonebook.empty()) {
              comm->Send(-1, status.MPI_SOURCE, ControlTag);
              spdlog::trace("Sent empty largest index");
            } else {
              comm->Send(*phonebook.rbegin()->first, status.MPI_SOURCE, ControlTag);
              spdlog::trace("Sent largest index {} unset", *phonebook.rbegin()->first);
            }
          } else if (command == ControlFlag::SET_WORKER_READY) {

            auto index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(status.MPI_SOURCE, ControlTag));
            int rank = comm->Recv<int>(status.MPI_SOURCE, ControlTag, &status);
            if (!phonebook.count(index)) {
              std::cerr << "setting ready for nonexistent model!" << std::endl;
              continue;
            }
            phonebook[index].SetWorkerReady(rank);
          } else if (command == ControlFlag::QUIT) {
            tracer->leaveRegion(TracerRegions::PhonebookBusy);
            spdlog::trace("Rank {} quit", comm->GetRank());
            break;
          }

          for (auto request_iter = requests.begin(); request_iter != requests.end();) {
            if (!phonebook.count(request_iter->requestedSampleIndex)) {
              std::cerr << "checking request for nonexistent model!" << std::endl;
            }
            int workerRank = phonebook[request_iter->requestedSampleIndex].NextWorker();
            if (workerRank == -1) {
              request_iter++;
            } else {
              comm->Send(workerRank, request_iter->sourceMPIRank, ControlTag);
              request_iter = requests.erase(request_iter);
            }
          }



          /*for (auto request_iter = requests.begin(); request_iter != requests.end();) {
            std::shared_ptr<MultiIndex> index = std::get<0>(*request_iter);
            int sender = std::get<1>(*request_iter);
            bool high_priority = std::get<2>(*request_iter);
            if (high_priority) {
              request_iter++;
              continue;
            }

            if (!phonebook.count(index)) {
              std::cerr << "checking request for nonexistent model!" << std::endl;
            }
            int workerRank = phonebook[index].NextWorker();
            if (workerRank == -1) {
              request_iter++;
            } else {
              comm->Send(workerRank, sender, ControlTag);
              request_iter = requests.erase(request_iter);
            }
          }*/

          tracer->leaveRegion(TracerRegions::PhonebookBusy);
        }
        //std::cout << "Phonebook " << comm->GetRank() << " idle time:\t" << timer_idle.elapsed() << " of:\t" << timer_full.elapsed() << std::endl;

      }

    private:

      void UnRegister(std::shared_ptr<MultiIndex> modelIndex, int rank) {
        /*if (!phonebook.count(modelIndex)) {
          std::cerr << "unsetting nonexistent entry!" << std::endl;
        }*/
        phonebook[modelIndex].RemoveWorker(rank);
        if (phonebook[modelIndex].NumWorkers() == 0) {
          spdlog::debug("Phonebook erasing modelIndex", *modelIndex);
          phonebook.erase(modelIndex);
        }
        //comm->Send(ControlFlag::HANDSHAKE, status.MPI_SOURCE, ControlTag);
        spdlog::trace("Phonebook entry for {} unset", *modelIndex);
      }

      class WorkerList {
      public:

        WorkerList() {
          ResetTimer();
        }

        int NextWorker() {
          //assert(workers.size() > 0);
          if (workers.size() == 0)
            return -1;

          if (workers_ready.size() == 0)
            return -1;

          //tick();
          int worker = workers_ready.front();
          workers_ready.pop_front();
          return worker;
        }
        void AddWorker(int worker) {
          //tick();
          workers.push_back(worker);
        }
        void RemoveWorker(int worker) {
          //tick();
          workers_ready.erase(std::remove(workers_ready.begin(), workers_ready.end(), worker), workers_ready.end());
          workers.erase(std::remove(workers.begin(), workers.end(), worker), workers.end());
          spdlog::trace("Removed worker {}, {} remaining", worker, NumWorkers());
        }
        int NumWorkers() const {
          return workers.size();
        }
        int NumWorkersReady() const {
          return workers_ready.size();
        }
        void SetWorkerReady(int worker) {
          //tick();
          registeredReadyCounter++;

          workers_ready.push_back(worker);
        }
        std::vector<int> GetWorkers() {
          return workers;
        }
        std::deque<int> GetWorkersReady() {
          return workers_ready;
        }

        double GetIdleFraction() {
          if (total_time == std::chrono::nanoseconds::zero())
            return .0;
          return ((double)idle_time.count()) / (double)(total_time.count());
        }
        double NormalizedRegisteredReadyCounter() {
          if (workers.size() == 0)
            return -1;
          return (double)registeredReadyCounter / (double)workers.size();
        }
        void ResetTimer() {
          registeredReadyCounter = 0;
          idle_time = std::chrono::nanoseconds::zero();
          total_time = std::chrono::nanoseconds::zero();
          begin_tick = std::chrono::high_resolution_clock::now();
          last_tick = begin_tick;
        }
        void tick() {
          std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();

          std::chrono::nanoseconds measurement_period = now - last_tick;
          idle_time += measurement_period * workers_ready.size();
          total_time += measurement_period * workers.size();

          last_tick = now;
        }
      public:
        bool recheck = false;

      private:
        int registeredReadyCounter;
        std::chrono::nanoseconds idle_time;
        std::chrono::nanoseconds total_time;
        std::chrono::high_resolution_clock::time_point begin_tick, last_tick;

        std::vector<int> workers;
        std::deque<int> workers_ready;
      };


      double getNumQueuedTasksForIndex (std::shared_ptr<MultiIndex> index, WorkerList& worker_list) {
        double in_queue = 0;
        for (auto request_iter = requests.begin(); request_iter != requests.end(); request_iter++) {
          if (*index == *(request_iter->requestedSampleIndex))
            in_queue += request_iter->highPriority ? 1.0 : .1;
        }
        return in_queue;
      }

      int getNumTasksQueuedTasksFromIndex (std::shared_ptr<MultiIndex> index) {
        int numQueuedTasks = 0;
        for (SampleRequest& request : requests) {
          if (request.highPriority)
            numQueuedTasks += *(request.sourceModelIndex) == *index;
        }
        return numQueuedTasks;
      }


      double getLoadFactor (std::shared_ptr<MultiIndex> index, WorkerList& worker_list) {
        double load_factor = 1.0 - worker_list.GetIdleFraction() - (double)getNumTasksQueuedTasksFromIndex(index) / worker_list.NumWorkers();

        load_factor += getNumQueuedTasksForIndex(index, worker_list) / worker_list.NumWorkers();
        return load_factor;
      }

      struct SampleRequest {
        std::shared_ptr<MultiIndex> requestedSampleIndex;
        std::shared_ptr<MultiIndex> sourceModelIndex;
        int sourceMPIRank;
        bool highPriority;
      };


      std::deque< SampleRequest > requests;


      std::map<std::shared_ptr<MultiIndex>, WorkerList, MultiPtrComp> phonebook;
      std::shared_ptr<parcer::Communicator> comm;
      bool scheduling_active;
      std::shared_ptr<muq::Utilities::OTF2TracerBase> tracer;
      bool rescheduling_in_progress = false;
    };

    /**
     * @brief High-level wrapper for communicating with the phonebook process.
     * @details This wraps communication with the phonebook in an easy to use
     * interface. Needed for example by sampling worker processes to announce when
     * samples are available to be collected.
     */
    class PhonebookClient {

    public:

      PhonebookClient(std::shared_ptr<parcer::Communicator> comm, int phonebookRank)
        : comm(comm), phonebookRank(phonebookRank)
      {
      }

      int Query(std::shared_ptr<MultiIndex> remoteIndex, std::shared_ptr<MultiIndex> sourceIndex, bool high_priority) {
        comm->Send(ControlFlag::GET_WORKGROUP, phonebookRank, ControlTag);
        comm->Send<MultiIndex>(*remoteIndex, phonebookRank, ControlTag);
        comm->Send<MultiIndex>(*sourceIndex, phonebookRank, ControlTag);
        comm->Send(high_priority, phonebookRank, ControlTag);
        return comm->Recv<int>(phonebookRank, ControlTag);
      }

      std::vector<int> GetWorkgroups(std::shared_ptr<MultiIndex> modelIndex) {
        spdlog::debug("GetWorkgroups call for model {}", *modelIndex);
        comm->Send(ControlFlag::GET_WORKGROUPS, phonebookRank, ControlTag);
        comm->Send(*modelIndex, phonebookRank, ControlTag);
        spdlog::debug("GetWorkgroups call for model {}, retrieving", *modelIndex);
        std::vector<int> ret = comm->Recv<std::vector<int>>(phonebookRank, ControlTag);
        spdlog::debug("GetWorkgroups call for model {}, returning", *modelIndex);
        return ret;
      }

      std::shared_ptr<MultiIndex> LargestIndex() {
        comm->Send(ControlFlag::GET_LARGEST_INDEX, phonebookRank, ControlTag);
        return std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(phonebookRank, ControlTag));
      }

      void Register(std::shared_ptr<MultiIndex> modelIndex, int rank) {
        comm->Send(ControlFlag::SET_WORKGROUP, phonebookRank, ControlTag);
        comm->Send(*modelIndex, phonebookRank, ControlTag);
        comm->Ssend(rank, phonebookRank, ControlTag);
      }

      void UnRegister(std::shared_ptr<MultiIndex> modelIndex, int rank) {
        comm->Send(ControlFlag::UNSET_WORKGROUP, phonebookRank, ControlTag);
        comm->Send(*modelIndex, phonebookRank, ControlTag);
        comm->Ssend(rank, phonebookRank, ControlTag);
        //if (comm->Recv<ControlFlag>(phonebookRank, ControlTag) != ControlFlag::HANDSHAKE)
        //  std::cerr << "Failed handshake in UnRegister()!" << std::endl;
      }

      void SetWorkerReady(std::shared_ptr<MultiIndex> modelIndex, int rank) {
        comm->Send(ControlFlag::SET_WORKER_READY, phonebookRank, ControlTag);
        comm->Send(*modelIndex, phonebookRank, ControlTag);
        comm->Send(rank, phonebookRank, ControlTag);
      }

      void SchedulingDone() {
        comm->Send(ControlFlag::SCHEDULING_DONE, phonebookRank, ControlTag);
      }
      void SchedulingStop() {
        comm->Ssend(ControlFlag::SCHEDULING_STOP, phonebookRank, ControlTag);
      }

    private:
      std::shared_ptr<parcer::Communicator> comm;
      int phonebookRank;
    };
  }
}

#endif

#endif
