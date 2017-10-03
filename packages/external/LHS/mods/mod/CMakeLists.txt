add_library(lhs_mod Parms.f90)
set_property(TARGET lhs_mod PROPERTY COMPILE_FLAGS "${FIXED_FLAG} ${FPP_FLAG}")
if(LHS_EXPORT_DEC)
  set_property(TARGET lhs_mod APPEND PROPERTY COMPILE_DEFINITIONS ${LHS_EXPORT_DEC})
endif()
lhs_install_library(lhs_mod)
# TODO: install .mod files to include/ ?
