cmake_minimum_required(VERSION 2.8)
project(PSUADE)

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")

if(Boost_DIR)
  include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
  add_definitions(-DHAVE_BOOST)
else()
  find_package(Boost REQUIRED)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
    add_definitions(-DHAVE_BOOST)
  endif()
endif()

add_library(psuade MOATAnalyzer.cpp MOATSampling.cpp DakotaPsuade.C)
install(TARGETS psuade
  EXPORT ${ExportTarget} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  )
