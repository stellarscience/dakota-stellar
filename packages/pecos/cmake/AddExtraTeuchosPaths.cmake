
# Adds "extra" Teuchos directory paths for INCLUDE and LIBRARY to support the
# "autoconf-style use case" in which the Teuchos headers/library are specified
# independently rather than the standard, CMake/Trilinos_PREFIX approach; i.e.
# search non-standard "install" locations for Teuchos_config.h and libteuchos.*

include(FindPackageMessage)

function(AddExtraTeuchosPaths)
  find_file(TEUCHOS_CONFIG_H Teuchos_config.h
            ${Extra_Teuchos_INCLUDE_DIR} NO_DEFAULT_PATH)

  if(EXISTS ${TEUCHOS_CONFIG_H})
    #message("CUSTOM FIND - Teuchos include path: ${Extra_Teuchos_INCLUDE_DIR}")
    include_directories(${Extra_Teuchos_INCLUDE_DIR})
  else()
    message(FATAL_ERROR
            "Teuchos headers NOT FOUND in: ${Extra_Teuchos_INCLUDE_DIR}")
  endif() # TEUCHOS_CONFIG_H EXISTS

  find_library(TEUCHOS_LIBRARY NAMES teuchos HINTS ${Extra_Teuchos_LIBRARY_DIR})

  if(EXISTS ${TEUCHOS_LIBRARY})
    find_package_message(Teuchos "Found Teuchos: ${TEUCHOS_LIBRARY}"
      "[${Extra_Teuchos_LIBRARY_DIR}][${Extra_Teuchos_INCLUDE_DIR}]")
    set(Teuchos_FOUND 1 PARENT_SCOPE)
    set(Teuchos_DIR ${Extra_Teuchos_LIBRARY_DIR} PARENT_SCOPE)
    link_directories(${Extra_Teuchos_LIBRARY_DIR})
  else()
    message(FATAL_ERROR
            "Teuchos library NOT FOUND in: ${Extra_Teuchos_LIBRARY_DIR}")
  endif() # TEUCHOS_LIBRARY EXISTS

endfunction(AddExtraTeuchosPaths)

