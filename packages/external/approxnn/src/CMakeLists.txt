# Simple CMakeLists.txt for ANN/src only added for Dakota project

cmake_minimum_required(VERSION 2.8)

project(ApproxNN CXX)

set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")

set(ApproxNN_SOURCES 
  ANN.cpp brute.cpp kd_tree.cpp kd_util.cpp kd_split.cpp
  kd_dump.cpp kd_search.cpp kd_pr_search.cpp kd_fix_rad_search.cpp
  bd_tree.cpp bd_search.cpp bd_pr_search.cpp bd_fix_rad_search.cpp
  perf.cpp
  )

set(ApproxNN_HEADERS
  kd_tree.h kd_split.h kd_util.h kd_search.h
  kd_pr_search.h kd_fix_rad_search.h pr_queue.h pr_queue_k.h
  )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

# When building the library, dllexport the symbols BMA: Would prefer
# to use target properties, but we build static by default, so we just
# add the definition.  This approach is only safe because no other
# sources in this directory or below import the header to use the
# library.
add_definitions(-DDLL_EXPORTS)
add_library(approxnn ${ApproxNN_SOURCES})
#set_target_properties(approxnn PROPERTIES DEFINE_SYMBOL DLL_EXPORTS)
#set_target_properties(approxnn PROPERTIES COMPILE_DEFINITIONS ...)

install(TARGETS approxnn EXPORT ${ExportTarget} DESTINATION lib)
install(FILES ${ApproxNN_HEADERS} DESTINATION include)
