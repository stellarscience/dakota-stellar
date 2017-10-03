
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)

########### next target ###############

SET(sandia_cubature_prb_SRCS
   sandia_cubature_prb.cpp
)

add_executable(sandia_cubature_prb ${sandia_cubature_prb_SRCS})

target_link_libraries(sandia_cubature_prb sparsegrid)

#install(TARGETS sandia_cubature_prb DESTINATION bin)


########### next target ###############

SET(sandia_rules_prb_SRCS
   sandia_rules_prb.cpp
)

add_executable(sandia_rules_prb ${sandia_rules_prb_SRCS})

target_link_libraries(sandia_rules_prb sparsegrid)

#install(TARGETS sandia_rules_prb DESTINATION bin)


########### next target ###############

SET(sandia_rules2_prb_SRCS
   sandia_rules2_prb.cpp
   sandia_rules2.cpp
)

add_executable(sandia_rules2_prb ${sandia_rules2_prb_SRCS})

target_link_libraries(sandia_rules2_prb sparsegrid)

#install(TARGETS sandia_rules2_prb DESTINATION bin)


########### next target ###############

SET(comp_next_prb_SRCS
   comp_next_prb.cpp
)

add_executable(comp_next_prb ${comp_next_prb_SRCS})

target_link_libraries(comp_next_prb sparsegrid)

#install(TARGETS comp_next_prb DESTINATION bin)


########### next target ###############

SET(product_mixed_growth_weight_prb_SRCS
   product_mixed_growth_weight_prb.cpp
   sandia_rules2.cpp
)

add_executable(product_mixed_growth_weight_prb ${product_mixed_growth_weight_prb_SRCS})

target_link_libraries(product_mixed_growth_weight_prb sparsegrid)

#install(TARGETS product_mixed_growth_weight_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmg_index_prb_SRCS
   sandia_sgmg_index_prb.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmg_index_prb ${sandia_sgmg_index_prb_SRCS})

target_link_libraries(sandia_sgmg_index_prb sparsegrid)

#install(TARGETS sandia_sgmg_index_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmg_point_prb_SRCS
   sandia_sgmg_point_prb.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmg_point_prb ${sandia_sgmg_point_prb_SRCS})

target_link_libraries(sandia_sgmg_point_prb sparsegrid)

#install(TARGETS sandia_sgmg_point_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmg_size_prb_SRCS
   sandia_sgmg_size_prb.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmg_size_prb ${sandia_sgmg_size_prb_SRCS})

target_link_libraries(sandia_sgmg_size_prb sparsegrid)

#install(TARGETS sandia_sgmg_size_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmg_size_table_SRCS
   sandia_sgmg_size_table.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmg_size_table ${sandia_sgmg_size_table_SRCS})

target_link_libraries(sandia_sgmg_size_table sparsegrid)

#install(TARGETS sandia_sgmg_size_table DESTINATION bin)


########### next target ###############

SET(sandia_sgmg_unique_index_prb_SRCS
   sandia_sgmg_unique_index_prb.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmg_unique_index_prb ${sandia_sgmg_unique_index_prb_SRCS})

target_link_libraries(sandia_sgmg_unique_index_prb sparsegrid)

#install(TARGETS sandia_sgmg_unique_index_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmg_weight_prb_SRCS
   sandia_sgmg_weight_prb.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmg_weight_prb ${sandia_sgmg_weight_prb_SRCS})

target_link_libraries(sandia_sgmg_weight_prb sparsegrid)

#install(TARGETS sandia_sgmg_weight_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmga_aniso_normalize_prb_SRCS
   sandia_sgmga_aniso_normalize_prb.cpp
)

add_executable(sandia_sgmga_aniso_normalize_prb ${sandia_sgmga_aniso_normalize_prb_SRCS})

target_link_libraries(sandia_sgmga_aniso_normalize_prb sparsegrid)

#install(TARGETS sandia_sgmga_aniso_normalize_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmga_index_prb_SRCS
   sandia_sgmga_index_prb.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmga_index_prb ${sandia_sgmga_index_prb_SRCS})

target_link_libraries(sandia_sgmga_index_prb sparsegrid)

#install(TARGETS sandia_sgmga_index_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmga_point_prb_SRCS
   sandia_sgmga_point_prb.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmga_point_prb ${sandia_sgmga_point_prb_SRCS})

target_link_libraries(sandia_sgmga_point_prb sparsegrid)

#install(TARGETS sandia_sgmga_point_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmga_product_weight_prb_SRCS
   sandia_sgmga_product_weight_prb.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmga_product_weight_prb ${sandia_sgmga_product_weight_prb_SRCS})

target_link_libraries(sandia_sgmga_product_weight_prb sparsegrid)

#install(TARGETS sandia_sgmga_product_weight_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmga_size_prb_SRCS
   sandia_sgmga_size_prb.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmga_size_prb ${sandia_sgmga_size_prb_SRCS})

target_link_libraries(sandia_sgmga_size_prb sparsegrid)

#install(TARGETS sandia_sgmga_size_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmga_size_table_SRCS
   sandia_sgmga_size_table.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmga_size_table ${sandia_sgmga_size_table_SRCS})

target_link_libraries(sandia_sgmga_size_table sparsegrid)

#install(TARGETS sandia_sgmga_size_table DESTINATION bin)


########### next target ###############

SET(sandia_sgmga_unique_index_prb_SRCS
   sandia_sgmga_unique_index_prb.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmga_unique_index_prb ${sandia_sgmga_unique_index_prb_SRCS})

target_link_libraries(sandia_sgmga_unique_index_prb sparsegrid)

#install(TARGETS sandia_sgmga_unique_index_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmga_vcn_coef_prb_SRCS
   sandia_sgmga_vcn_coef_prb.cpp
)

add_executable(sandia_sgmga_vcn_coef_prb ${sandia_sgmga_vcn_coef_prb_SRCS})

target_link_libraries(sandia_sgmga_vcn_coef_prb sparsegrid)

#install(TARGETS sandia_sgmga_vcn_coef_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmga_vcn_prb_SRCS
   sandia_sgmga_vcn_prb.cpp
)

add_executable(sandia_sgmga_vcn_prb ${sandia_sgmga_vcn_prb_SRCS})

target_link_libraries(sandia_sgmga_vcn_prb sparsegrid)

#install(TARGETS sandia_sgmga_vcn_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmga_weight_prb_SRCS
   sandia_sgmga_weight_prb.cpp
   sandia_rules2.cpp
)

add_executable(sandia_sgmga_weight_prb ${sandia_sgmga_weight_prb_SRCS})

target_link_libraries(sandia_sgmga_weight_prb sparsegrid)

#install(TARGETS sandia_sgmga_weight_prb DESTINATION bin)


########### next target ###############

SET(sandia_sgmgg_prb_SRCS
   sandia_sgmgg_prb.cpp
)

add_executable(sandia_sgmgg_prb ${sandia_sgmgg_prb_SRCS})

target_link_libraries(sandia_sgmgg_prb sparsegrid)

#install(TARGETS sandia_sgmgg_prb DESTINATION bin)


########### install files ###############

