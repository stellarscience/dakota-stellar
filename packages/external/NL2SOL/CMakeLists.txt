cmake_minimum_required(VERSION 2.6)
project("NL2SOL" C)

## NL2SOL library sources.
file(GLOB Nl2sol_src *.c)

set(Nl2sol_src ${Nl2sol_src}
)

## Set the NL2SOL library name.
add_library(cport ${Nl2sol_src})

#### WJB - ToDo: investigate implicit dependency on libf2c

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")

install(TARGETS cport EXPORT ${ExportTarget} DESTINATION lib)
install(FILES f2c.h DESTINATION include)