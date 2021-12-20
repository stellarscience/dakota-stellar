find_package(PkgConfig)

if(NOT DEFINED MUQ_SPDLOG_DIR)

	pkg_check_modules(PC_SPDLOG QUIET spdlog)
	set(SPDLOG_DEFINITIONS ${PC_SPDLOG_CFLAGS_OTHER})

	find_path(SPDLOG_INCLUDE_DIR spdlog/spdlog.h
					HINTS ${PC_SPDLOG_INCLUDEDIR} ${PC_SPDLOG_INCLUDE_DIRS}
					PATH_SUFFIXES spdlog )

	find_library(SPDLOG_LIBRARY NAMES spdlog
             HINTS ${PC_SPDLOG_LIBDIR} ${PC_SPDLOG_LIBRARY_DIRS} )

	find_library(SPDLOG_LIBRARY_STATIC NAMES ${library_prefix}spdlog.${static_library_suffix}
             HINTS ${PC_SPDLOG_LIBDIR} ${PC_SPDLOG_LIBRARY_DIRS} )

else()
	find_path(SPDLOG_INCLUDE_DIR spdlog/spdlog.h
	          HINTS ${MUQ_SPDLOG_DIR}/include
	          PATH_SUFFIXES spdlog NO_DEFAULT_PATH)

	find_library(SPDLOG_LIBRARY NAMES spdlog
	             HINTS ${MUQ_SPDLOG_DIR}/lib NO_DEFAULT_PATH)

endif()

set(SPDLOG_LIBRARIES ${SPDLOG_LIBRARY} )
set(SPDLOG_INCLUDE_DIRS ${SPDLOG_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(SPDLOG  DEFAULT_MSG
                                  SPDLOG_LIBRARY SPDLOG_INCLUDE_DIR)

mark_as_advanced(SPDLOG_INCLUDE_DIRS SPDLOG_LIBRARIES )