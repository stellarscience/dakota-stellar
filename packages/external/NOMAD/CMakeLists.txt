cmake_minimum_required(VERSION 2.8)
project("NOMAD" CXX)

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")

SUBDIRS(src)
