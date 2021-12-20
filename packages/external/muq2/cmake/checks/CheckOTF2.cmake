set(CMAKE_REQUIRED_LIBRARIES ${OTF2_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${OTF2_INCLUDE_DIRS})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
  "
 // This test code is from otf2_writer_example.c provided by libotf2. See included libotf2 package for its BSD style license.
 #include <otf2/otf2.h>
 #include <stdlib.h>
 static OTF2_TimeStamp
 get_time( void )
 {
     static uint64_t sequence;
     return sequence++;
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
 int
 main( int    argc,
       char** argv )
 {
     OTF2_Archive* archive = OTF2_Archive_Open( \"ArchivePath\",
                                                \"ArchiveName\",
                                                OTF2_FILEMODE_WRITE,
                                                1024 * 1024 /* event chunk size */,
                                                4 * 1024 * 1024 /* def chunk size */,
                                                OTF2_SUBSTRATE_POSIX,
                                                OTF2_COMPRESSION_NONE );
     OTF2_Archive_SetFlushCallbacks( archive, &flush_callbacks, NULL );
     OTF2_Archive_SetSerialCollectiveCallbacks( archive );
     OTF2_Archive_OpenEvtFiles( archive );
     OTF2_EvtWriter* evt_writer = OTF2_Archive_GetEvtWriter( archive, 0 );
     OTF2_EvtWriter_Enter( evt_writer,
                           NULL,
                           get_time(),
                           0 /* region */ );
     OTF2_EvtWriter_Leave( evt_writer,
                           NULL,
                           get_time(),
                           0 /* region */ );
     OTF2_Archive_CloseEvtWriter( archive, evt_writer );
     OTF2_Archive_CloseEvtFiles( archive );
     OTF2_Archive_OpenDefFiles( archive );
     OTF2_DefWriter* def_writer = OTF2_Archive_GetDefWriter( archive, 0 );
     OTF2_Archive_CloseDefWriter( archive, def_writer );
     OTF2_Archive_CloseDefFiles( archive );
     OTF2_GlobalDefWriter* global_def_writer = OTF2_Archive_GetGlobalDefWriter( archive );
     OTF2_GlobalDefWriter_WriteClockProperties( global_def_writer,
                                                1 /* 1 tick per second */,
                                                0 /* epoch */,
                                                2 /* length */ );
     OTF2_GlobalDefWriter_WriteString( global_def_writer, 0, \"\" );
     OTF2_GlobalDefWriter_WriteString( global_def_writer, 1, \"Master Process\" );
     OTF2_GlobalDefWriter_WriteString( global_def_writer, 2, \"Main Thread\" );
     OTF2_GlobalDefWriter_WriteString( global_def_writer, 3, \"MyFunction\" );
     OTF2_GlobalDefWriter_WriteString( global_def_writer, 4, \"Alternative function name (e.g. mangled one)\" );
     OTF2_GlobalDefWriter_WriteString( global_def_writer, 5, \"Computes something\" );
     OTF2_GlobalDefWriter_WriteString( global_def_writer, 6, \"MyHost\" );
     OTF2_GlobalDefWriter_WriteString( global_def_writer, 7, \"node\" );
     OTF2_GlobalDefWriter_WriteRegion( global_def_writer,
                                       0 /* id */,
                                       3 /* region name  */,
                                       4 /* alternative name */,
                                       5 /* description */,
                                       OTF2_REGION_ROLE_FUNCTION,
                                       OTF2_PARADIGM_USER,
                                       OTF2_REGION_FLAG_NONE,
                                       0 /* source file */,
                                       0 /* begin lno */,
                                       0 /* end lno */ );
     OTF2_GlobalDefWriter_WriteSystemTreeNode( global_def_writer,
                                               0 /* id */,
                                               6 /* name */,
                                               7 /* class */,
                                               OTF2_UNDEFINED_SYSTEM_TREE_NODE /* parent */ );
     OTF2_GlobalDefWriter_WriteLocationGroup( global_def_writer,
                                              0 /* id */,
                                              1 /* name */,
                                              OTF2_LOCATION_GROUP_TYPE_PROCESS,
                                              0 /* system tree */ );
     OTF2_GlobalDefWriter_WriteLocation( global_def_writer,
                                         0 /* id */,
                                         2 /* name */,
                                         OTF2_LOCATION_TYPE_CPU_THREAD,
                                         2 /* # events */,
                                         0 /* location group */ );
     OTF2_Archive_Close( archive );
     return EXIT_SUCCESS;
 }
   "
  OTF2_COMPILES)

if (NOT SPDLOG_COMPILES)
  set(SPDLOG_TEST_FAIL 1)
else ()
  set(SPDLOG_TEST_FAIL 0)
endif()
