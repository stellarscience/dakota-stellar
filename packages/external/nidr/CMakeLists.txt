# NIDR build process
#
# * Always build libnidr (all target)
#   from sources (which may be checked into the repo)
#   nidr-scanner.c nidr-parser.c avltree.c nidr.c
#
# * If NIDR_NIDRGEN, build nidrgen and Dt binaries
#   not currently installed, but perhaps should be for JAGUAR benefit
#
# * If NIDR_REGENERATE_YACCLEX and have lex/yacc
#   build src-scanner: nidr-scanner.c nidr-parser.h nidr-parser.c nidrgen.c
#   NOTE: target "nidr-generated-files" must be explictly invoked

cmake_minimum_required(VERSION 2.6)

project(NIDR)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

function(file_replace input output find replace)
  ADD_CUSTOM_COMMAND(
    OUTPUT "${output}"
    DEPENDS "${input}"
    COMMAND ${CMAKE_COMMAND}
    ARGS
      -Dinfile:FILEPATH=${input}
      -Doutfile:FILEPATH=${output}
      -Dfind:STRING=${find}
      -Dreplace:STRING=${replace}
      -P ${NIDR_SOURCE_DIR}/cmake/replace.cmake
    VERBATIM)
endfunction()

OPTION(NIDR_REGENERATE_YACCLEX "Regenerate YACC and LEX files" OFF)
IF(NIDR_REGENERATE_YACCLEX)
  FIND_PROGRAM(YACC_EXECUTABLE
    NAMES yacc bison
    PATHS /usr/bin
    DOC "Yacc or Bison executable")
  FIND_PROGRAM(FLEX_EXECUTABLE
    NAMES flex
    PATHS /usr/bin
    DOC "Flex executable")
  MARK_AS_ADVANCED(YACC_EXECUTABLE FLEX_EXECUTABLE)
  IF(YACC_EXECUTABLE)
    SET(BISON_FLAGS)
    IF(YACC_EXECUTABLE MATCHES "bison")
      SET(BISON_FLAGS "--yacc")
    ENDIF(YACC_EXECUTABLE MATCHES "bison")
    SET(src "${CMAKE_CURRENT_SOURCE_DIR}/nidrparse.y")
    SET(hdr "${CMAKE_CURRENT_BINARY_DIR}/nidrparse.1.h")
    ADD_CUSTOM_COMMAND(
      OUTPUT "${hdr}"
      "${CMAKE_CURRENT_BINARY_DIR}/nidrparse.1.c"
      #OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/nidrparse.1.c"
      #       "${hdr}"
      DEPENDS "${src}"
      COMMAND
      ${YACC_EXECUTABLE}
      ARGS
        ${BISON_FLAGS} --defines=${hdr}
        -o${CMAKE_CURRENT_BINARY_DIR}/nidrparse.1.c ${src}
      )

    #sed rules for nidr-parser.c
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrparse.1.c" "${CMAKE_CURRENT_BINARY_DIR}/nidrparse.2.c" "yy" "nidr") 
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrparse.2.c" "${CMAKE_CURRENT_BINARY_DIR}/nidrparse.3.c" "reset_lex_state" "reset_nidrlex_state") 
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrparse.3.c" "${CMAKE_CURRENT_BINARY_DIR}/nidrparse.4.c" "(#line [0-9]+ \"[^\"]+\")" "/* \\1 */") 
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrparse.4.c" "${CMAKE_CURRENT_BINARY_DIR}/nidrparse.5.c" "${CMAKE_CURRENT_BINARY_DIR}" " ") 
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrparse.5.c" "${CMAKE_CURRENT_SOURCE_DIR}/nidr-parser.c" "${CMAKE_CURRENT_SOURCE_DIR}" " ") 
    
    #sed rules for nidr-parser.h
    file_replace("${hdr}" "${CMAKE_CURRENT_BINARY_DIR}/nidrparse.2.h" "yy" "nidr") 
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrparse.2.h" "${CMAKE_CURRENT_BINARY_DIR}/nidrparse.3.h" "${CMAKE_CURRENT_BINARY_DIR}/" " ") 
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrparse.3.h" "${CMAKE_CURRENT_SOURCE_DIR}/nidr-parser.h" "${CMAKE_CURRENT_SOURCE_DIR}/" " ") 

    ADD_CUSTOM_TARGET(RerunYacc DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/nidr-parser.c")
  ENDIF(YACC_EXECUTABLE)
  IF(FLEX_EXECUTABLE)
    SET(LexInputFiles nidrgen nidrscan)
    SET(lex_files)
    FOREACH(name ${LexInputFiles})
      SET(src "${CMAKE_CURRENT_SOURCE_DIR}/${name}.l")
      SET(dst "${CMAKE_CURRENT_BINARY_DIR}/${name}.1.c")
      SET(hdr "${CMAKE_CURRENT_BINARY_DIR}/${name}.h")
      ADD_CUSTOM_COMMAND(
        OUTPUT "${dst}"
        DEPENDS "${src}"
        COMMAND
        ${FLEX_EXECUTABLE}
        --header-file=${hdr} -o${dst} ${src}
        )
      SET(lex_files ${lex_files} "${dst}")
    ENDFOREACH(name)
    ADD_CUSTOM_TARGET(RerunLex DEPENDS ${lex_files})

    #sed rules for nidrgen.c
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrgen.1.c" "${CMAKE_CURRENT_BINARY_DIR}/nidrgen.2.c" "(#line [0-9]+ \"[^\"]+\")" "/* \\1 */") 
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrgen.2.c" "${CMAKE_CURRENT_BINARY_DIR}/nidrgen.3.c" "${CMAKE_CURRENT_BINARY_DIR}" " ") 
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrgen.3.c" "${CMAKE_CURRENT_SOURCE_DIR}/nidrgen.c" "${CMAKE_CURRENT_SOURCE_DIR}" " ") 

    #sed rules for nidr-scanner.c
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrscan.1.c" "${CMAKE_CURRENT_BINARY_DIR}/nidrscan.2.c" "yy" "nidr") 
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrscan.2.c" "${CMAKE_CURRENT_BINARY_DIR}/nidrscan.3.c" "reset_lex_state" "reset_nidrlex_state") 
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrscan.3.c" "${CMAKE_CURRENT_BINARY_DIR}/nidrscan.4.c" "(#line [0-9]+ \"[^\"]+\")" "/* \\1 */") 
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrscan.4.c" "${CMAKE_CURRENT_BINARY_DIR}/nidrscan.5.c" "${CMAKE_CURRENT_BINARY_DIR}" " ") 
    file_replace("${CMAKE_CURRENT_BINARY_DIR}/nidrscan.5.c" "${CMAKE_CURRENT_SOURCE_DIR}/nidr-scanner.c" "${CMAKE_CURRENT_SOURCE_DIR}" " ") 

  ENDIF(FLEX_EXECUTABLE)

  # target nidr-generated-files must be explictly invoked
  add_custom_target(nidr-generated-files
    DEPENDS nidr-parser.c
            nidr-parser.h
	    nidr-scanner.c
	    nidrgen.c
    VERBATIM
    )

ENDIF(NIDR_REGENERATE_YACCLEX)

include(CheckIncludeFile)

option(NIDR_DL_SOLVER "Toggle NIDR DL Solvers, default is disabled." OFF)
add_definitions("-DNO_NIDR_DYNLIB")

# Build the nidrgen and Dt binaries if requested
OPTION(NIDR_NIDRGEN "Build nidrgen utilities" ON)
if (NIDR_NIDRGEN)

  if(CMAKE_SYSTEM_NAME MATCHES Darwin)
    set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS
        "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
  endif()
  
  add_executable(nidrgen nidrgen.c)
  set(nidrgen-libs "")

  if(NIDR_DL_SOLVER)
    find_package(DL)

    if(DL_FOUND)
      set(nidrgen-libs ${CMAKE_DL_LIBS})
      remove_definitions("-DNO_NIDR_DYNLIB")
    else()
      message(WARNING "DL Solver requested, but NO dl* support found; disabling.")
      add_definitions("-DNO_NIDR_DYNLIB")
    endif(DL_FOUND)
  endif(NIDR_DL_SOLVER)

  find_library(libm m)
  if(libm)
    list(APPEND nidrgen-libs ${libm})
  endif()

  if(nidrgen-libs)
    target_link_libraries(nidrgen ${nidrgen-libs})
  endif()

  # Not building Dt for now due to need for localtime_r alternative on Windows
  #add_executable(Dt Dt.c)

endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# libnidr has extern symbols which aren't resolved until link time
# TODO: generalize for non-gcc compilers?
if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS
      "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
endif()

# Always build libnidr, possibly (most likely) from sources checked into the repo
if(NIDR_DL_SOLVER)
  find_package(DL)

  if(DL_FOUND)
    message(STATUS "NIDR DL Solver enabled, LibDL = ${CMAKE_DL_LIBS}")
    remove_definitions("-DNO_NIDR_DYNLIB")
  else()
    message(WARNING "DL Solver requested, but NO dl* support found; disabling.")
    add_definitions("-DNO_NIDR_DYNLIB")
    set(NIDR_DL_SOLVER OFF CACHE BOOL 
      "NIDR DLsolver is disabled -- NO dlfcn support on this platform" FORCE)
  endif() # DL_FOUND
endif(NIDR_DL_SOLVER)

add_library(nidr nidr-scanner.c nidr-parser.c avltree.c nidr.c)

# Install headers.
install(
  FILES nidr.h
  DESTINATION include
  )

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget NIDRTargets CACHE STRING 
  "Export set name for ${CMAKE_PROJECT_NAME}")

# Install binaries.
install(
  TARGETS nidr
  EXPORT  ${ExportTarget}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  )

# Install package configuration and import files.
configure_file(${NIDR_SOURCE_DIR}/cmake/nidr-config.cmake.in
               ${NIDR_BINARY_DIR}/cmake/nidr-config.cmake @ONLY)
install(FILES ${NIDR_BINARY_DIR}/cmake/nidr-config.cmake
        DESTINATION lib/cmake/nidr)

# Workaround until we fix Dakota's overall use of exports
if(ExportTarget STREQUAL NIDRTargets)
  message(STATUS 
    "Installing NIDR exported targets to lib/cmake/NIDR")
  install(EXPORT ${ExportTarget} DESTINATION lib/cmake/NIDR)
endif()

