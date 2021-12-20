## DEMO_TPL library sources.
set(demo_tpl_HEADERS
  demo_opt.hpp
)

set(demo_tpl_SOURCES
  demo_opt.cpp
)

## Set the DEMO_TPL library name.
add_library(demo_tpl ${demo_tpl_SOURCES})

# Define install targets for "make install"
install(TARGETS demo_tpl EXPORT ${ExportTarget} DESTINATION lib)
