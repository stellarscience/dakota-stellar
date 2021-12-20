#ifndef PARALLELFLAGS_H_
#define PARALLELFLAGS_H_

#if MUQ_HAS_MPI

namespace muq {
  namespace SamplingAlgorithms {

    /**
     * @brief Tags separating MPI communication between
     * control level processes and internal work group communication.
     */
    const int ControlTag = 1;
    const int WorkgroupTag = 2;

    /**
     * @brief Flags used by parallel MCMC/MIMCMC type methods.
     * @details Commands sent between the ranks of a parallel MCMC/MIMCMC
     * method are distinguished by these flags.
     */
    enum ControlFlag : const int {

      // SamplingProblem
      FINALIZE,
      INIT_PROBLEM,
      LOGDENSITY,
      TEST,
      QOI,

      // Assignment
      ASSIGN,
      ASSIGN_COLLECTOR,
      UNASSIGN,
      SAMPLE,
      SAMPLE_BOX,
      QUIT,

      SAMPLE_BOX_DONE,

      // Collector
      MEANS,
      MEANS_DONE,
      WRITE_TO_FILE,

      // Workgroup phonebook
      GET_WORKGROUP,
      SET_WORKGROUP,
      GET_WORKGROUPS,
      UNSET_WORKGROUP,
      GET_LARGEST_INDEX,
      SET_WORKER_READY,
      SCHEDULING_NEEDED,
      SCHEDULING_DONE,
      SCHEDULING_STOP,

      HANDSHAKE
    };

  }
}

#endif

#endif
