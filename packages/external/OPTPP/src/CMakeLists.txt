set(local_library optpp)

## Top-level library sources.
set(lib${local_library}_sources globals.C)

# Build a single library from all subdirectories
set(local_src_subdirs Base Constraints GSS Newton PDS Utils)

# Aggregate all the sources into a single library for installation
foreach(sd_it ${local_src_subdirs})

  add_subdirectory(${sd_it})

  # grab the trailing part of any paths
  string(REGEX MATCH "[^/]+$" last_path_name ${sd_it})

  # extract the variables ${last_path_name}_sources from the subdir's env
  get_directory_property(sd_sources DIRECTORY ${sd_it} 
    DEFINITION ${last_path_name}_sources)

  # add the sources to the library's list and includes
  foreach ( src_file ${sd_sources} )
    list ( APPEND lib${local_library}_sources "${sd_it}/${src_file}" )
  endforeach()

  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${sd_it})

endforeach()

add_library(${local_library} ${lib${local_library}_sources})
target_link_libraries(${local_library} ${Teuchos_LIBRARIES} ${TEUCHOS_LIBRARY} ${BLAS_LIBS} ${LAPACK_LIBS})

# TODO: is this needed?
set_target_properties(${local_library} PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${optpp_BINARY_DIR}/lib
  ARCHIVE_OUTPUT_DIRECTORY ${optpp_BINARY_DIR}/lib)

install(TARGETS ${local_library} EXPORT ${ExportTarget} DESTINATION lib)
