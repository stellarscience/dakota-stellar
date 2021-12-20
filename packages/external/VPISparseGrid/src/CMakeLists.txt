set(sparsegrid_LIB_SRCS
   sandia_cubature.cpp
   sandia_rules.cpp
   sandia_sgmg.cpp
   sandia_sgmga.cpp
   sandia_sgmgg.cpp
)

add_library(sparsegrid ${sparsegrid_LIB_SRCS})
# Compile for position-independent code for use in the higher-level Pecos shared lib
install(TARGETS sparsegrid EXPORT ${ExportTarget} DESTINATION lib)

########### install files ###############
install(FILES sandia_cubature.hpp sandia_rules.hpp sandia_sgmg.hpp sandia_sgmgg.hpp sandia_sgmga.hpp DESTINATION include)
