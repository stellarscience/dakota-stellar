# USAGE:
# CreateCompileGroup(
#     <group name>
#     <brief description>
#     <library_name>
#     <other compile group dependencies>
#     <required dependencies>
#     <optional dependencies>
#     <source file 1>
#     <source file 2>
#     ...
#     <source file N>
# )

CreateCompileGroup(
    SAMPLING_PROBLEMS
    "Implement sampling sampling problems for MCMC."
    muqSamplingAlgorithms
    "MODELING_CORE;APPROXIMATION_REGRESSION"
    ""
    ""

    SamplingProblem.cpp
)

CreateCompileGroup(
    FANCY_SAMPLING_PROBLEMS
    "Implement more advanced sampling problem."
    muqSamplingAlgorithms
    "MODELING_CORE;APPROXIMATION_REGRESSION"
    ""
    ""

    ExpensiveSamplingProblem.cpp
)


CreateCompileGroup(
    SAMPLING_ALGORITHM
    "Implement sampling algorithms such as Monte Carlo, importance sampling, and Markov chain Monte Carlo"
    muqSamplingAlgorithms
    "UTILITIES_CORE;UTILITIES_MULTIINDEX;MODELING_CORE;MODELING_DISTRIBUTIONS;UTILITIES_CORE;UTILITIES_HDF5;SAMPLING_PROBLEMS"
    ""
    ""
    AbstractSamplingProblem.cpp
    SamplingAlgorithm.cpp
    #MonteCarlo.cpp
    ImportanceSampling.cpp
    SingleChainMCMC.cpp

    SamplingState.cpp
    SampleCollection.cpp
    #DistributedCollection.cpp
    MarkovChain.cpp

    TransitionKernel.cpp
    DummyKernel.cpp
    MIDummyKernel.cpp
    #MCKernel.cpp
    #ISKernel.cpp
    MHKernel.cpp
    DRKernel.cpp
    DILIKernel.cpp
    MIKernel.cpp
    GreedyMLMCMC.cpp
    MIMCMC.cpp
    MIMCMCBox.cpp
    SLMCMC.cpp
    SubsamplingMIProposal.cpp
    GMHKernel.cpp

    MCMCFactory.cpp

    ThinScheduler.cpp

    MCMCProposal.cpp
    MHProposal.cpp
    AMProposal.cpp
    MALAProposal.cpp
    InfMALAProposal.cpp
    MixtureProposal.cpp

    InverseGammaProposal.cpp
    CrankNicolsonProposal.cpp
)

if(MUQ_USE_MPI)

  CreateCompileGroup(
      PARALLEL_SAMPLING_ALGORITHM
      "Tools for sampling that require MPI distributed parallelism."
      muqSamplingAlgorithms
      "MODELING_CORE;UTILITIES_CORE;SAMPLING_ALGORITHM"
      "PARCER;SPDLOG;OTF2"
      ""
      ParallelAMProposal.cpp
      DistributedCollection.cpp
      ParallelAbstractSamplingProblem.cpp
      ParallelMIComponentFactory.cpp
  )

endif()
