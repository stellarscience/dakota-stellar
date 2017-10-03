include_directories(${Dream_SOURCE_DIR})

add_executable(problem0 problem0.cpp)
target_link_libraries(problem0 dream)

# Can't test these without separating dream and rnglib functions

#add_executable(pdflib_prb pdflib_prb.cpp)
#target_link_libraries(pdflib_prb dream)

#add_executable(rnglib_prb rnglib_prb.cpp)
#target_link_libraries(rnglib_prb dream)
