set(sparsegrid_LIB_SRCS
   sandia_cubature.cpp
   sandia_rules.cpp
   sandia_sgmg.cpp
   sandia_sgmga.cpp
   sandia_sgmgg.cpp
)

add_library(sparsegrid STATIC ${sparsegrid_LIB_SRCS})
# Compile for position-independent code for use in the higher-level Pecos shared lib
set_target_properties(sparsegrid PROPERTIES POSITION_INDEPENDENT_CODE True)
set_target_properties(sparsegrid PROPERTIES VERSION 4.2.0 SOVERSION 4)
install(TARGETS sparsegrid EXPORT ${ExportTarget} DESTINATION lib)

########### install files ###############
install(FILES sandia_cubature.hpp sandia_rules.hpp sandia_sgmg.hpp sandia_sgmgg.hpp sandia_sgmga.hpp DESTINATION include)
