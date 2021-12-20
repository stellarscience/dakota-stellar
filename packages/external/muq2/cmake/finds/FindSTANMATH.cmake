find_package(PkgConfig)
if(NOT DEFINED MUQ_STANMATH_DIR)

  pkg_check_modules(PC_STANMATH QUIET STANMATH)
  set(STANMATH_DEFINITIONS ${PC_STANMATH_CFLAGS_OTHER})

  find_path(STANMATH_INCLUDE_DIR stan/math.hpp
            HINTS ${PC_STANMATH_INCLUDEDIR} ${PC_STANMATH_INCLUDE_DIRS} /usr/local/include /usr/include
            PATH_SUFFIXES stan)

else()

  find_path(STANMATH_INCLUDE_DIR stan/math.hpp
            HINTS ${MUQ_STANMATH_DIR}
            PATH_SUFFIXES stan NO_DEFAULT_PATH)

endif()
set(STANMATH_INCLUDE_DIRS ${STANMATH_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(STANMATH DEFAULT_MSG STANMATH_INCLUDE_DIR)

mark_as_advanced(STANMATH_INCLUDE_DIR)
