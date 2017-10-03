
add_subdirectory(mod)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/mod)

## Define the sources for the primary LHS module.
set(LHS_MODS_SRC 
  Ccmatr.f90 Cparam.f90 Crank.f90   Csamp.f90 Cworkc.f90
  Cworkx.f90 Distnm.f90 Killfile.f90 Miscmod.f90 ModInByC.f90
)
## Set the convenience library name for the primary LHS modules.
add_library(lhs_mods ${LHS_MODS_SRC})
set_property(TARGET lhs_mods PROPERTY COMPILE_FLAGS "${FIXED_FLAG} ${FPP_FLAG}")
if(LHS_EXPORT_DEC)
  set_property(TARGET lhs_mods APPEND PROPERTY COMPILE_DEFINITIONS ${LHS_EXPORT_DEC})
endif()
target_link_libraries(lhs_mods lhs_mod)
lhs_install_library(lhs_mods)
# TODO: install .mod files to include/ ?
