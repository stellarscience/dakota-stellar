#ifndef RUNTIME_ENVIRONMENT_HPP
#define RUNTIME_ENVIRONMENT_HPP

//#include "MPIFunctions.hpp"

#define RUNTIME_MPI_COMM_WORLD   0
#define RUNTIME_MPI_COMM_NULL    0

namespace Pecos {

typedef int MPI_Comm;

class ParallelObject;
class RuntimeEnvironment;

class MPICommunicator
{
protected:

  /// MPI Communicator
  MPI_Comm MPICommunicator_;
  
public:  

  MPICommunicator( MPI_Comm mpi_comm )
  {
    MPICommunicator_ =  mpi_comm;
  }

  ~MPICommunicator(){};

  friend class RuntimeEnvironment;
  friend class ParallelObject;
};

/**
 * \class Iterator 
 * \brief The base class for all parallel iterators.
 */
class RuntimeEnvironment
{

protected:

  /// MPI Communicator on which Heat is running
  MPICommunicator MPICommunicator_;
  
  /// The unqiue identifier of the processor in MPICommunicator_
  int environmentProcessorId_;

  /// The number of processors in MPICommunicator_
  int numProcessors_;

  /// Flag specifiying ownership of MPI_Init/MPI_Finalize
  bool ownMPI_;
  
  /// Defines the verbosity level of i/o
  int verbosity_;

public:

  /**
   *  \brief Default library mode constructor (assumes MPI_COMM_WORLD)
   */
  RuntimeEnvironment()
    : MPICommunicator_( RUNTIME_MPI_COMM_WORLD ), environmentProcessorId_( 0 ),
      numProcessors_( 1 ), ownMPI_( false ), verbosity_( 0 )
  { 
    initialise();
  };

  /**
   * \brief Stand-alone mode constructor. This constructor is the one 
   * used by main.cpp
   */
  RuntimeEnvironment( int argc, char **argv )
    : MPICommunicator_( RUNTIME_MPI_COMM_WORLD ), environmentProcessorId_( 0 ),
      numProcessors_( 1 ), ownMPI_( true ), verbosity_( 0 )
  {
#ifdef ENABLE_LIBHEAT_MPI
    MPI_Init( &argc, &argv );
#endif
    initialise();
  };

  /**
   * \brief Library mode constructor accepting communicator
   */
  RuntimeEnvironment( MPI_Comm mpi_comm )
    : MPICommunicator_( mpi_comm ), environmentProcessorId_( 0 ), 
      numProcessors_( 1 ), ownMPI_( false ), verbosity_( 0 )
  {
    initialise();
  };

  /// Deconstructor. Call MPI finalise if running in parallel mode and ownMPI_
  ~RuntimeEnvironment()
  {
#ifdef ENABLE_LIBHEAT_MPI
    if ( ownMPI_ )
      MPI_Finalize();
#endif
  };
  
  void initialise()
  {

#ifdef ENABLE_LIBHEAT_MPI
    // Get total number of processes
    MPI_Comm_size( MPICommunicator_.MPICommunicator_, &numProcessors_ ); 

    // Get process number for this process
    MPI_Comm_rank( MPICommunicator_.MPICommunicator_, &environmentProcessorId_ );
#endif

    if ( verbosity_ > 0 )
      std::cout << "Running MPI Executable on " << numProcessors_ << 
	" processors\n";
  };

  /// Return the unqiue identifier of the current processor. 
  int environment_processor_id() 
  {
    return environmentProcessorId_;
  };

  /// Return the number of processors in the environement
  int num_processors()
  {
    return numProcessors_;
  };

  /// Set the verbosity level of i/o
  void verbosity( int verbosity_in )
  {
    verbosity_ = verbosity_in;
  }

  MPICommunicator* get_mpi_communicator()
  {
    return &MPICommunicator_;
  }

};

class ParallelObject
{
 
protected:
  /// The unqiue identifier of the processor
  int processorId_;

  /// The unqiue identifier of the master processor
  int masterProcessorId_;

  /// The number of processors in MPICommunicator_
  int numProcessors_;

  /// MPI Communicator on which ParallelObject is running
  MPI_Comm MPICommunicator_;

  /// Defines the verbosity level of i/o
  int verbosity_;

public:
  /// Default constructor
  ParallelObject() :
    processorId_( 0 ), masterProcessorId_( 0 ), numProcessors_( 1 ),
    MPICommunicator_( RUNTIME_MPI_COMM_NULL )
  {};

 ~ParallelObject()
  {};

  /// Set the MPI communicator
  void mpi_communicator( MPICommunicator* mpi_comm )
  {
    MPICommunicator_ =  mpi_comm->MPICommunicator_;
#ifdef ENABLE_LIBHEAT_MPI
    // Get total number of processes
    MPI_Comm_size( MPICommunicator_, &numProcessors_ ); 
    
    // Get process number for this process
    MPI_Comm_rank( MPICommunicator_, &processorId_ );
#endif
  };

  int num_processors()
  {
    return numProcessors_;
  };

  /// Set the verbosity level of i/o
  void verbosity( int verbosity_in )
  {
    verbosity_ = verbosity_in;
  }
  
  /// Determine if the current processor is the master
  bool is_master()
  {
    return ( processorId_ == masterProcessorId_ );
  };

  /// Return the id of the master processor
  int master_processor_id()
  {
    return masterProcessorId_;
  };

  /// Return the unqiue identifier of the current processor. 
  int processor_id() 
  {
    return processorId_;
  };

};

} // namespace Pecos

#endif //RUNTIME_ENVIRONMENT_HPP
