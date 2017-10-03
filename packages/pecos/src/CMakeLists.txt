# CMakeLists.txt for pecos source (WJB: Note the source file GLOB below)

include_directories(${Boost_INCLUDE_DIR} ${Boost_INCLUDE_DIRS}
  ${Teuchos_INCLUDE_DIRS})

remove_definitions("-DHAVE_CONFIG_H")

option(BOOST_DISABLE_ASSERTS "Toggle BoostAssert, default is AssertDisabled" ON)
if (BOOST_DISABLE_ASSERTS)
  add_definitions("-DBOOST_DISABLE_ASSERTS")
endif(BOOST_DISABLE_ASSERTS)

if (${CMAKE_CXX_COMPILER_ID} MATCHES "XL" OR
    ${CMAKE_CXX_COMPILER_ID} MATCHES "VisualAge")
  add_definitions("-DBOOST_MPL_CFG_ASSERT_USE_RELATION_NAMES")
endif()


# --- Dependency configuration based on previously defined packages ---

if(HAVE_LHS)
  include_directories(${LHS_BINARY_DIR}/include)
  list(APPEND Pecos_libs lhs lhs_mods lhs_mod)
endif()

if(HAVE_SPARSE_GRID)
  include_directories(${VPISparseGrid_SOURCE_DIR}/src)
  list(APPEND Pecos_libs sparsegrid)
endif()

if(HAVE_FFT)
  if(HAVE_FFTW)
    include_directories("${FFTW_SOURCE_DIR}/api")
    link_directories("${FFTW_BINARY_DIR}/.libs")
    list(APPEND Pecos_libs fftw3)
  endif()
  if(HAVE_DFFTPACK)
    list(APPEND Pecos_libs dfftpack)
  endif()
endif()

# Unconditional dependencies.
list(APPEND Pecos_libs ${Teuchos_LIBRARIES} ${TEUCHOS_LIBRARY}
  ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${FCLIBS}
)

## Pecos library sources.
file(GLOB Pecos_src *.cpp)
file(GLOB Pecos_hdr *.hpp)

add_library(pecos_src ${Pecos_src})
target_link_libraries(pecos_src ${Pecos_libs})
if (HAVE_FFTW)
  add_dependencies(pecos_src fftw)
endif()
install(TARGETS pecos_src EXPORT ${ExportTarget} DESTINATION lib)
install(FILES ${Pecos_hdr} DESTINATION include)
