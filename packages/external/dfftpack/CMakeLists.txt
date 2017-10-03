cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(DFFTPACK C Fortran)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${DFFTPACK_SOURCE_DIR}")
  include(FortranCInterface)
elseif(NOT COMMAND FortranCInterface_HEADER)
  message(FATAL_ERROR "The top-level project must include(FortranCInterface)")
endif()

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")

option(BUILD_SHARED_LIBS "Build DFFTPACK with shared libraries?" ON)

## dfftpack library sources.
SET(dfftpack_SOURCES
  zfftb.f		
  cfftb1.f	
  zfftf.f		
  cfftf1.f	
  zffti.f		
  cffti1.f	
  dcosqb.f	
  cosqb1.f	
  dcosqf.f	
  cosqf1.f	
  dcosqi.f	
  dcost.f		
  dcosti.f	
  ezfft1.f	
  dzfftb.f	
  dzfftf.f	
  dzffti.f	
  passb.f		
  passb2.f	
  passb3.f	
  passb4.f	
  passb5.f	
  passf.f		
  passf2.f	
  passf3.f	
  passf4.f	
  passf5.f	
  radb2.f		
  radb3.f		
  radb4.f		
  radb5.f		
  radbg.f		
  radf2.f		
  radf3.f		
  radf4.f		
  radf5.f		
  radfg.f		
  dfftb.f		
  rfftb1.f	
  dfftf.f		
  rfftf1.f	
  dffti.f		
  rffti1.f	
  dsinqb.f	
  dsinqf.f	
  dsinqi.f	
  dsint.f		
  sint1.f		
  dsinti.f
)

if(BUILD_SHARED_LIBS AND WIN32 AND NOT CYGWIN)
  # Compute public interface symbol names after Fortran mangling.
  string(TO${FortranCInterface_GLOBAL_CASE}
    "${FortranCInterface_GLOBAL_PREFIX}zfftb${FortranCInterface_GLOBAL_SUFFIX}" ZFFTB)
  string(TO${FortranCInterface_GLOBAL_CASE}
    "${FortranCInterface_GLOBAL_PREFIX}zffti${FortranCInterface_GLOBAL_SUFFIX}" ZFFTI)

  # Configure a module definition file to export the public interface.
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/dfftpack.def.in
                 ${CMAKE_CURRENT_BINARY_DIR}/dfftpack.def @ONLY)
  list(APPEND dfftpack_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/dfftpack.def)

  # This platform information is missing in CMake <= 2.8.4.
  if(NOT CMAKE_LINK_DEF_FILE_FLAG
      AND "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_LINK_DEF_FILE_FLAG "-Wl,")
  endif()
endif()

if(${CMAKE_Fortran_COMPILER_ID} MATCHES "XL")
  set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "-r")
endif()

## Set the dfftpack library name.
add_library(dfftpack ${dfftpack_SOURCES})
#set_target_properties(dfftpack PROPERTIES LINKER_LANGUAGE Fortran)
install(TARGETS dfftpack EXPORT ${ExportTarget} DESTINATION lib)
