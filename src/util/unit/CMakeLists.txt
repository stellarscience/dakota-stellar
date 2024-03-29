#  _______________________________________________________________________
#
#  DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
#  Copyright 2014-2020
#  National Technology & Engineering Solutions of Sandia, LLC (NTESS).
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Dakota directory.
#  _______________________________________________________________________

# BMA TODO: helpers to manage all these test registration steps

include(DakotaUnitTest)

# List of test files to copy into this directory
#set(dakota_unit_test_copied_files)

dakota_add_unit_test(NAME DataScalerTest
  SOURCES util_test_driver.cpp DataScalerTest.cpp
  )

# Rationale: These all include Teuchos headers and need to link to the
# unit test system
target_include_directories(DataScalerTest PRIVATE "${Teuchos_INCLUDE_DIRS}")
target_link_libraries(DataScalerTest PRIVATE dakota_util ${Teuchos_LIBRARIES})


dakota_add_unit_test(NAME MathToolsTest
  SOURCES util_test_driver.cpp MathToolsTest.cpp
  )

target_include_directories(MathToolsTest PRIVATE "${Teuchos_INCLUDE_DIRS}")
target_link_libraries(MathToolsTest PRIVATE dakota_util ${Teuchos_LIBRARIES})

dakota_add_unit_test(NAME MetricsTest
  SOURCES util_test_driver.cpp MetricsTest.cpp
  )

target_include_directories(MetricsTest PRIVATE "${Teuchos_INCLUDE_DIRS}")
target_link_libraries(MetricsTest PRIVATE dakota_util ${Teuchos_LIBRARIES})


dakota_add_unit_test(NAME TeuchosEigenConversionTest
  SOURCES util_test_driver.cpp teuchos_eigen_conversion.cpp
  )

target_include_directories(TeuchosEigenConversionTest PRIVATE "${Teuchos_INCLUDE_DIRS}")
target_link_libraries(TeuchosEigenConversionTest PRIVATE dakota_util ${Teuchos_LIBRARIES})

dakota_add_unit_test(NAME LinearSolverTest
  SOURCES util_test_driver.cpp LinearSolverTest.cpp
  )

target_include_directories(LinearSolverTest PRIVATE "${Teuchos_INCLUDE_DIRS}")
target_link_libraries(LinearSolverTest PRIVATE dakota_util ${Teuchos_LIBRARIES})

#target_include_directories(DataScalerTest PRIVATE
#  "${CMAKE_CURRENT_SOURCE_DIR}/.." "${Teuchos_INCLUDE_DIRS}")

# Add a single target for all copied files
#add_custom_target(util_unit_test_files ALL
#  DEPENDS ${util_unit_test_copied_files})

if(DAKOTA_PYTHON AND DAKOTA_PYBIND11)
  pybind11_add_module(pybind11_example pybind11_example.cpp)
  # BMA TODO: couldn't get to work using .py in source tree
  configure_file(pybind11_example.py pybind11_example.py COPYONLY)
  add_test(NAME dakota_util_pybind11_example
    COMMAND ${Python_EXECUTABLE} pybind11_example.py)
  set_tests_properties(dakota_util_pybind11_example PROPERTIES LABELS UnitTest)
endif()
