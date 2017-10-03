cmake_minimum_required(VERSION 2.8)
project(LHS C Fortran)
set(CMAKE_MODULE_PATH ${LHS_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})
if("${CMAKE_SOURCE_DIR}" STREQUAL "${LHS_SOURCE_DIR}")
  set(LHS_STANDALONE 1)
else()
  set(LHS_STANDALONE 0)
endif()

if(LHS_STANDALONE)
  include(FortranCInterface)
elseif(NOT COMMAND FortranCInterface_HEADER)
  message(FATAL_ERROR "The top-level project must include(FortranCInterface)")
endif()

if(NOT LHS_STANDALONE AND CMAKE_SYSTEM_NAME MATCHES Darwin)
  set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS
     "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

# Provide mangling information to C programs.
FortranCInterface_HEADER(include/LHS.h MACRO_NAMESPACE LHS_)
install(FILES ${LHS_BINARY_DIR}/include/LHS.h DESTINATION include)

# Organize output files.  On Windows this also keeps .dll files next
# to the .exe files that need them, making tests easy to run.
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${LHS_BINARY_DIR}/bin)
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LHS_BINARY_DIR}/lib)
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LHS_BINARY_DIR}/lib)
endif()

# Workaround Intel VS plugin bug that fails to create directories.
# CMake 2.8.4.20110223 and higher already have a workaround.
if(${CMAKE_VERSION} VERSION_LESS 2.8.4.20110223 AND
    ${CMAKE_GENERATOR} MATCHES "Visual Studio")
  foreach(c ${CMAKE_CONFIGURATION_TYPES})
    file(MAKE_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/${c})
    file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${c})
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${c})
  endforeach()
endif()

option(BUILD_SHARED_LIBS "Build LHS with shared libraries?" ON)

# On Windows-GNU builds try to provide MS import libraries too.
if(BUILD_SHARED_LIBS)
  include(GNUtoMS)
endif()

if(GNUtoMS)
  set(LHS_GNUtoMS_IMPORT ${LHS_SOURCE_DIR}/CMake/LHS-GNUtoMS.cmake)
  set(LHS_GNUtoMS_INSTALL "include(\${_SELF_DIR}/LHS-GNUtoMS.cmake)")
  set(LHS_GNUtoMS_BUILD "include(\"${LHS_GNUtoMS_IMPORT}\")")
endif()

# Set the export name for install targets; parent packages likely want
# to override this to the name of their target
set(ExportTarget ${CMAKE_PROJECT_NAME}-targets CACHE STRING 
  "Name for the export target for ${CMAKE_PROJECT_NAME}")

macro(lhs_install_library lib)
  install(TARGETS ${lib} 
    EXPORT ${ExportTarget}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
  )
  if(GNUtoMS)
    install(FILES ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/lib${lib}.lib DESTINATION lib)
  endif()
endmacro()

option(LHS_Intel "Build LHS with intel settings?" OFF)

if(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel" OR
    ${CMAKE_GENERATOR} MATCHES "Visual Studio") # Intel plugin
  set(LHS_Intel 1)
endif()

# hack to be fixed when CMake does it automatically
if(LHS_Intel)
  set(FIXED_FLAG "-fixed")
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(FIXED_FLAG "-ffixed-form")
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
  set(FIXED_FLAG "-Mfixed")
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "VisualAge")
  set(FIXED_FLAG "-qfixed")
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "XL")
  set(FIXED_FLAG "-qfixed")
  #set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "-r")
  set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "-r")

  #set(CMAKE_Fortran_FLAGS "-Wl,--allow-multiple-definition")
  set(CMAKE_Fortran_FLAGS "-Wl,-zmuldefs")
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "SunPro")
  set(FIXED_FLAG "-fixed")
endif()

# For Windows DLLs with Intel Fortran preprocess the source files and
# substitute "DEC$" for LHS_EXPORT_DEC to activate dllexport markup.
if(BUILD_SHARED_LIBS AND WIN32 AND LHS_Intel)
  set(FPP_FLAG "-fpp")
  set(LHS_EXPORT_DEC "LHS_EXPORT_DEC=DEC$")
endif()

if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(LEGACY_FLAG "-std=legacy")
endif()

add_subdirectory(mods)

set(LHS_UTIL_SRC
  Banner.f90 Beta.f90 Betafn.f90 Betaic.f90 Betaln.f90 Binom.f90 
  Chkdat.f90 Chkdim.f90 Chkemp.f90 Chkstr.f90 Chkzro.f90 
  Chlsky.f90 Cmcrd.f90 Corcal.f90 Corout.f90 Cumulc.f90 
  Cumuld.f90 C_wrap.f90 Datout.f90 Datsqz.f90 Dmfsd.f90 Dminmax.f90 
  Dsinv.f90 Entrpf.f90 Entrpy.f90 Errchk.f90 Errget.f90 
  Erstgt.f90 Erxset.f90 Expon.f90 Factor.f90 Factr2.f90 
  Fileoc.f90 Findit.f90 Finvnor.f90 Gammab.f90 Gamma.f90 
  Gammam.f90 Geom.f90 Hgeom.f90 Histo.f90 Hpsrt.f90 Hstout.f90 
  Hypgeo.f90 Igaus1.f90 Igaus.f90 Igausf.f90 Imtql2.f90 
  Interp.f90 Intrpd.f90 Ljust.f90 Matinv.f90 Mix.f90 Nbinom.f90 
  Normal.f90 Outcrd.f90 Outdat.f90 Pareto.f90 Pmtrx.f90 
  Poison.f90 Posdef.f90 Ranker.f90 Rierfc1.f90 Rmcnp2.f90 
  Rmcnp.f90 Rmcnpi2.f90 Rmcnpi.f90  
  Samout.f90 SamStor.f90 Setdef.f90 Sift.f90 Sspev.f90 Table.f90 
  Tqlrat.f90 Trbak3.f90 Tred3.f90 Triang.f90 Unifrm.f90 Vif.f90 
  Weibul.f90 Wrtcrd.f90 Wrtpar.f90 Gumbel.f90 Frechet.f90 
  Lhssetseed.f90 DefaultRnum1.f90 DefaultRnum2.f90 
)

set(LHS_CORE_SRC
  Lhs_cls.f90 Lhs_cnst.f90 Lhs_corr.f90 Lhs_cout.f90 Lhs_dist.f90 
  Lhs_file.f90 Lhs_imem.f90 Lhs_init.f90 Lhs_optn.f90 Lhs_prep.f90
  Lhs_rsed.f90 Lhs_run.f90 Lhs_same.f90 Lhs_sdst.f90 Lhs_udst.f90
)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/mods)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/mods/mod)

add_library(lhs ${LHS_UTIL_SRC} ${LHS_CORE_SRC})

target_link_libraries(lhs lhs_mod lhs_mods)
set_property(TARGET lhs PROPERTY COMPILE_FLAGS
  "${FIXED_FLAG} ${FPP_FLAG} ${LEGACY_FLAG}"
)
if(LHS_EXPORT_DEC)
  set_property(TARGET lhs APPEND PROPERTY COMPILE_DEFINITIONS ${LHS_EXPORT_DEC})
endif()

lhs_install_library(lhs)
# TODO: install .mod files to include/ ?


set(LHSDRV_SRC
  fileleng.f90 flname.f90 lhs.f90 lhsdrv.f90 lread.f90
  newcrd.f90 rdpar2.f90 rdpar.f90 read.f90
)

add_executable(lhsdrv ${LHSDRV_SRC})
target_link_libraries(lhsdrv lhs)
set_property(TARGET lhsdrv PROPERTY COMPILE_FLAGS ${FIXED_FLAG})
install(TARGETS lhsdrv DESTINATION bin)

set(DATAFILES InpData.txt KEYWORD.DAT SIPRA.INI)

foreach (DATAFILE ${DATAFILES})
  add_custom_command(TARGET lhsdrv POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${LHS_SOURCE_DIR}/${DATAFILE}
    ${LHS_BINARY_DIR}/${DATAFILE})
endforeach(DATAFILE)

option(LHS_ENABLE_TESTS "Enable LHS-specific tests?" ON)
if(BUILD_TESTING AND LHS_ENABLE_TESTS)
  ## LHS test program.
  add_executable(TestRun TestRun.f90)
  target_link_libraries(TestRun lhs)
  set_property(TARGET TestRun PROPERTY COMPILE_FLAGS ${FIXED_FLAG})
  add_test(TestLHS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TestRun)
endif() # LHS Testing

if(LHS_STANDALONE)
  configure_file(${LHS_SOURCE_DIR}/CMake/LHSConfig-build.cmake.in
    ${LHS_BINARY_DIR}/LHSConfig.cmake @ONLY)
  configure_file(${LHS_SOURCE_DIR}/CMake/LHSConfig-install.cmake.in
    ${LHS_BINARY_DIR}/CMakeFiles/LHSConfig.cmake @ONLY)
  export(TARGETS lhs lhs_mods lhs_mod FILE ${ExportTarget}.cmake)
  install(FILES
    ${LHS_GNUtoMS_IMPORT}
    ${LHS_BINARY_DIR}/CMakeFiles/LHSConfig.cmake
    DESTINATION lib/cmake/LHS
    )
  install(EXPORT ${ExportTarget} DESTINATION lib/cmake/LHS)
endif()
