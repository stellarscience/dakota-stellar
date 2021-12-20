#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

# Create Doxygen files with substituted PATH information
if( DOXYGEN_FOUND )

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    @ONLY)
  #-- Add a custom target (doc) to run Doxygen when ever the project is built
  #-- To create documentation use make doc
  add_custom_target(doc 
    COMMAND ${DOXYGEN_EXECUTABLE}  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    COMMENT "Generating API documentation with Doxygen" VERBATIM)

else()

  message( FATAL_ERROR "Doxygen is needed to build the documentation. 
			     Please install it correctly" )

endif( DOXYGEN_FOUND )

