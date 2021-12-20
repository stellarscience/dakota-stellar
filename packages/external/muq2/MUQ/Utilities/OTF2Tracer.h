#ifndef OTF2TRACER_H
#define OTF2TRACER_H

#if MUQ_HAS_OTF2==1
#include <otf2/otf2.h>
#include <otf2/OTF2_MPI_Collectives.h>
#endif

namespace muq {
  namespace Utilities{

    enum TracerRegions {
      Setup,
      Finalize,
      PhonebookBusy,
      BurnIn,
      Sampling,
      CollectorBusy,
      RetrievingProposal,
      FetchingProposal
    };

    /**
     * @brief Base interface for OTF2 tracer implemetations
     *
     */
    class OTF2TracerBase {
    public:

      virtual void enterRegion(TracerRegions region) = 0;

      virtual void leaveRegion(TracerRegions region) = 0;

      virtual void write() = 0;
    };

    /**
     * @brief Fallback dummy implementation not doing anything; Does not require libotf2
     *
     */
    class OTF2TracerDummy : public OTF2TracerBase {
    public:

      OTF2TracerDummy() {}
      OTF2TracerDummy(std::string archive_path, std::string archive_name) {}

      void enterRegion(TracerRegions region) override {}

      void leaveRegion(TracerRegions region) override {}

      void write() override {
        spdlog::warn("write() has been called on an OTF2TracerDummy; no OTF2 trace file will be written! You may need to install libotf2.");
      }
    };


#if MUQ_HAS_OTF2==1

    static OTF2_TimeStamp
    get_time( void )
    {
      double t = MPI_Wtime() * 1e9;
      return ( uint64_t )t;
    }
    static OTF2_FlushType
    pre_flush( void*            userData,
                OTF2_FileType    fileType,
                OTF2_LocationRef location,
                void*            callerData,
                bool             final )
    {
      return OTF2_FLUSH;
    }
    static OTF2_TimeStamp
    post_flush( void*            userData,
                OTF2_FileType    fileType,
                OTF2_LocationRef location )
    {
      return get_time();
    }
    static OTF2_FlushCallbacks flush_callbacks =
    {
      .otf2_pre_flush  = pre_flush,
      .otf2_post_flush = post_flush
    };

    /**
     * @brief Tracer implementation writing to OTF2 via libotf2
     * The result can be viewed by several programs, one
     * example is "vite" developed at INRIA.
     */
    class OTF2Tracer : public OTF2TracerBase {


    public:

      OTF2Tracer(std::string archive_path, std::string archive_name) {
        MPI_Comm_size( MPI_COMM_WORLD, &size );
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        archive = OTF2_Archive_Open( archive_path.c_str(),
                                    archive_name.c_str(),
                                    OTF2_FILEMODE_WRITE,
                                    1024 * 1024 /* event chunk size */,
                                    4 * 1024 * 1024 /* def chunk size */,
                                    OTF2_SUBSTRATE_POSIX,
                                    OTF2_COMPRESSION_NONE );
        OTF2_Archive_SetFlushCallbacks( archive, &flush_callbacks, NULL );
        OTF2_MPI_Archive_SetCollectiveCallbacks( archive,
                                                MPI_COMM_WORLD,
                                                MPI_COMM_NULL );
        OTF2_Archive_OpenEvtFiles( archive );

        evt_writer = OTF2_Archive_GetEvtWriter( archive, rank );

        epoch_start = get_time();

        setRegionName(TracerRegions::BurnIn, "Burnin");
        setRegionName(TracerRegions::CollectorBusy, "CollectorBusy");
        setRegionName(TracerRegions::Finalize, "Finalize");
        setRegionName(TracerRegions::PhonebookBusy, "PhonebookBusy");
        setRegionName(TracerRegions::RetrievingProposal, "RetrievingProposal");
        setRegionName(TracerRegions::Sampling, "Sampling");
        setRegionName(TracerRegions::Setup, "Setup");
        setRegionName(TracerRegions::FetchingProposal, "FetchingProposal");
      }

      /**
       * @brief Call this to mark that a certain tracer region has been entered.
       */
      void enterRegion(TracerRegions region) override {
        ensureRegionName(region);
        OTF2_EvtWriter_Enter( evt_writer,
                              NULL,
                              get_time(),
                              region );
      }

      /**
       * @brief Call this to mark that a certain tracer region has been left.
       */
      void leaveRegion(TracerRegions region) override {
        ensureRegionName(region);
        OTF2_EvtWriter_Leave( evt_writer,
                              NULL,
                              get_time(),
                              region );
      }

    private:
      /**
       * @brief Makes sure a region name has been defined, otherwise generates a default region name
       */
      void ensureRegionName(TracerRegions region) {
        if (regionNames.count(region) == 0)
          regionNames[region] = "Unnamed region " + std::to_string(region);
      }

      void setRegionName(TracerRegions region, std::string name) {
        regionNames[region] = name;
      }

    public:

      void write() override {
        uint64_t epoch_end = get_time();
        OTF2_Archive_CloseEvtWriter( archive, evt_writer );
        OTF2_Archive_CloseEvtFiles( archive );
        OTF2_Archive_OpenDefFiles( archive );
        OTF2_DefWriter* def_writer = OTF2_Archive_GetDefWriter( archive,
                                                                rank );
        OTF2_Archive_CloseDefWriter( archive, def_writer );
        OTF2_Archive_CloseDefFiles( archive );
        uint64_t global_epoch_start;
        MPI_Reduce( &epoch_start,
                    &global_epoch_start,
                    1, OTF2_MPI_UINT64_T, MPI_MIN,
                    0, MPI_COMM_WORLD );
        uint64_t global_epoch_end;
        MPI_Reduce( &epoch_end,
                    &global_epoch_end,
                    1, OTF2_MPI_UINT64_T, MPI_MAX,
                    0, MPI_COMM_WORLD );
        if ( 0 == rank )
        {
          OTF2_GlobalDefWriter* global_def_writer = OTF2_Archive_GetGlobalDefWriter( archive );
          OTF2_GlobalDefWriter_WriteClockProperties( global_def_writer,
                                                    1000000000,
                                                    global_epoch_start,
                                                    global_epoch_end - global_epoch_start + 1 );
          OTF2_GlobalDefWriter_WriteString( global_def_writer, 0, "" );
          OTF2_GlobalDefWriter_WriteString( global_def_writer, 1, "Master Thread" );
          OTF2_GlobalDefWriter_WriteString( global_def_writer, 2, "MPI_Barrier" );
          OTF2_GlobalDefWriter_WriteString( global_def_writer, 3, "PMPI_Barrier" );
          OTF2_GlobalDefWriter_WriteString( global_def_writer, 4, "barrier" );
          OTF2_GlobalDefWriter_WriteString( global_def_writer, 5, "MyHost" );
          OTF2_GlobalDefWriter_WriteString( global_def_writer, 6, "node" );
          OTF2_GlobalDefWriter_WriteString( global_def_writer, 7, "MPI" );
          OTF2_GlobalDefWriter_WriteString( global_def_writer, 8, "MPI_COMM_WORLD" );
          int num_strings = 9;

          for ( const auto &regionNamePair : regionNames ) {

            OTF2_GlobalDefWriter_WriteString( global_def_writer, num_strings, regionNamePair.second.c_str() );
            OTF2_GlobalDefWriter_WriteRegion( global_def_writer,
                                              regionNamePair.first /* id */,
                                              num_strings /* region name  */,
                                              num_strings /* alternative name */,
                                              num_strings /* description */,
                                              OTF2_REGION_ROLE_CODE,
                                              OTF2_PARADIGM_MPI,
                                              OTF2_REGION_FLAG_NONE,
                                              7 /* source file */,
                                              0 /* begin lno */,
                                              0 /* end lno */ );
            num_strings++;
          }
          OTF2_GlobalDefWriter_WriteSystemTreeNode( global_def_writer,
                                                    0 /* id */,
                                                    5 /* name */,
                                                    6 /* class */,
                                                    OTF2_UNDEFINED_SYSTEM_TREE_NODE /* parent */ );

          for ( int r = 0; r < size; r++ )
          {
            char process_name[ 32 ];
            sprintf( process_name, "MPI Rank %d", r );
            OTF2_GlobalDefWriter_WriteString( global_def_writer,
                                              num_strings + r,
                                              process_name );
            OTF2_GlobalDefWriter_WriteLocationGroup( global_def_writer,
                                                    r /* id */,
                                                    num_strings + r /* name */,
                                                    OTF2_LOCATION_GROUP_TYPE_PROCESS,
                                                    0 /* system tree */ );
            OTF2_GlobalDefWriter_WriteLocation( global_def_writer,
                                                r /* id */,
                                                1 /* name */,
                                                OTF2_LOCATION_TYPE_CPU_THREAD,
                                                4 /* # events */,
                                                r /* location group */ );
          }
          uint64_t comm_locations[ size ];
          for ( int r = 0; r < size; r++ )
          {
            comm_locations[ r ] = r;
          }
          OTF2_GlobalDefWriter_WriteGroup( global_def_writer,
                                          0 /* id */,
                                          7 /* name */,
                                          OTF2_GROUP_TYPE_COMM_LOCATIONS,
                                          OTF2_PARADIGM_MPI,
                                          OTF2_GROUP_FLAG_NONE,
                                          size,
                                          comm_locations );
          OTF2_GlobalDefWriter_WriteGroup( global_def_writer,
                                          1 /* id */,
                                          0 /* name */,
                                          OTF2_GROUP_TYPE_COMM_GROUP,
                                          OTF2_PARADIGM_MPI,
                                          OTF2_GROUP_FLAG_NONE,
                                          size,
                                          comm_locations );
          OTF2_GlobalDefWriter_WriteComm( global_def_writer,
                                          0 /* id */,
                                          8 /* name */,
                                          1 /* group */,
                                          OTF2_UNDEFINED_COMM /* parent */ );
          OTF2_Archive_CloseGlobalDefWriter( archive,
                                            global_def_writer );
        }

        OTF2_Archive_Close( archive );
      }

    private:
      int size;
      int rank;
      uint64_t epoch_start;

      OTF2_Archive* archive;
      OTF2_EvtWriter* evt_writer;

      std::map<int, std::string> regionNames;
    };
#else

    using OTF2Tracer = OTF2TracerDummy; // Fall back to dummy implementation if OTF2 library is not found.

#endif

  }
}

#endif