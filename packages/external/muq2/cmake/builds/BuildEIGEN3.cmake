set(EIGEN_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/src/EIGEN3")


include(ExternalProject)
if(NOT EIGEN_EXTERNAL_SOURCE)
	set(EIGEN_EXTERNAL_SOURCE https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.bz2)
endif()

set(EIGEN3_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/muq_external/)
ExternalProject_Add(
  EIGEN3
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3
  URL ${EIGEN_EXTERNAL_SOURCE}
	BUILD_COMMAND ""
	CONFIGURE_COMMAND ""
  INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${EIGEN3_INSTALL_DIR}/include && cp -r ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/src/EIGEN3/Eigen ${EIGEN3_INSTALL_DIR}/include && cp -r ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/src/EIGEN3/unsupported ${EIGEN3_INSTALL_DIR}/include
)




set_property( TARGET EIGEN3 PROPERTY FOLDER "Externals")

set(EIGEN3_INCLUDE_DIRS ${EIGEN3_INSTALL_DIR}/include/eigen3)
message(STATUS "Adding ${EIGEN3_INSTALL_DIR}/include/eigen3 for an Eigen include directory.")
