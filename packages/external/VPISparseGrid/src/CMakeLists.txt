set(sparsegrid_LIB_SRCS
   sandia_cubature.cpp
   sandia_rules.cpp
   sandia_sgmg.cpp
   sandia_sgmga.cpp
   sandia_sgmgg.cpp
)

add_library(sparsegrid STATIC ${sparsegrid_LIB_SRCS})

# Compile for position-independent code for use in the higher-level Pecos shared lib
if(NOT WIN32)
  if(DEFINED CMAKE_SHARED_LIBRARY_CXX_FLAGS)
    set_target_properties(sparsegrid PROPERTIES COMPILE_FLAGS ${CMAKE_SHARED_LIBRARY_CXX_FLAGS})
  elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "XL")
    set_target_properties(sparsegrid PROPERTIES COMPILE_FLAGS "-qpic")
  elseif(CMAKE_SYSTEM_NAME MATCHES "[Dd][Aa][Rr][Ww][Ii][Nn]")
    set_target_properties(sparsegrid PROPERTIES COMPILE_FLAGS "-fPIC")
  elseif(NOT DEFINED CMAKE_SHARED_LIBRARY_CXX_FLAGS AND 
         CMAKE_SYSTEM_NAME MATCHES "[Cc][Yy][Gg][Ww][Ii][Nn]")
    message("Building a STATIC, cygwin library in ${CMAKE_CURRENT_BINARY_DIR}")
  else()
    message(FATAL_ERROR "Unable to determine the proper PIC flag for compiling
            in ${CMAKE_CURRENT_BINARY_DIR}")
  endif()
endif()

#set_target_properties(sparsegrid PROPERTIES VERSION 4.2.0 SOVERSION 4)
install(TARGETS sparsegrid EXPORT ${ExportTarget} DESTINATION lib)

########### install files ###############
install(FILES sandia_cubature.hpp sandia_rules.hpp sandia_sgmg.hpp sandia_sgmgg.hpp sandia_sgmga.hpp DESTINATION include)
