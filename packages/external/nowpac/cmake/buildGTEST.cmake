INCLUDE(ExternalProject)

SET(NOWPAC_GTEST_URL "https://github.com/google/googletest/archive/release-1.7.0.zip")

ExternalProject_Add( GTEST
   
   PREFIX ${NOWPAC_INSTALL_PREFIX}/nowpac/external/gtest
   URL ${NOWPAC_GTEST_URL}
   INSTALL_DIR ${NOWPAC_INSTALL_PREFIX}/nowpac/external/gtest
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${NOWPAC_INSTALL_PREFIX}/nowpac/external/gtest -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   BUILD_COMMAND "make"
   INSTALL_COMMAND "make"
   LOG_DOWNLOAD 1
   LOG_CONFIGURE 1
   LOG_BUILD 1
   LOG_INSTALL 1
   BUILD_IN_SOURCE 0

)

#SET_PROPERTY( TARGET GTEST PROPERTY FOLDER "nowpac/external")
SET(NOWPAC_GTEST_INCLUDE_PATH "${NOWPAC_INSTALL_PREFIX}/nowpac/external/gtest/src/GTEST/include")
SET(NOWPAC_GTEST_LIBRARY_PATH "${NOWPAC_INSTALL_PREFIX}/nowpac/external/gtest/src/GTEST-build")

