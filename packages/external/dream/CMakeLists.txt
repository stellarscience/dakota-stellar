cmake_minimum_required(VERSION 2.8)
PROJECT("Dream" CXX)

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")

add_library(dream dream.cpp pdflib.cpp rnglib.cpp)

option(DREAM_ENABLE_TESTING "Enable DREAM test drivers" OFF)
if (DREAM_ENABLE_TESTING)
  add_subdirectory(test)
endif()

install(TARGETS dream EXPORT ${ExportTarget} DESTINATION lib)

########### install files ###############
install(FILES 
  dream.hpp
  dream_user.hpp
  pdflib.hpp
  rnglib.hpp  
  DESTINATION include
  )
