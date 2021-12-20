#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

configure_file("setup.py.cmake.in" "setup.py" @ONLY)

set(pydakota_stamp "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

# Probably need the Development component too, but wasn't probed for prior
find_package(Python REQUIRED COMPONENTS Interpreter)

add_custom_command(
  OUTPUT ${pydakota_stamp}
  COMMAND ${Python_EXECUTABLE} setup.py build
  COMMAND ${CMAKE_COMMAND} -E touch ${pydakota_stamp}
##  DEPENDS ${DEPS}
  )

add_custom_target(pydakota DEPENDS ${pydakota_stamp})

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
  "build;swigpyrun.h")

install(CODE "execute_process(COMMAND ${Python_EXECUTABLE} setup.py install)")
