add_executable(tstillum tstillum.C illum.C MAIN.c tstfcn.h)
target_link_libraries(tstillum ${test_libs})
add_test(tstillum tstillum)

add_executable(tstrosen tstrosen.C rosen.C MAIN.c tstfcn.h)
target_link_libraries(tstrosen ${test_libs})
add_test(tstrosen tstrosen)

add_executable(tstnpsollinear tstlinear.C rosen.C MAIN.c tstfcn.h)
target_link_libraries(tstnpsollinear ${test_libs})
add_test(tstnpsollinear tstnpsollinear)

add_executable(tstnpsolnonlinear tstnonlinear.C sixtyfive.C MAIN.c tstfcn.h)
target_link_libraries(tstnpsolnonlinear ${test_libs})
add_test(tstnpsolnonlinear tstnpsolnonlinear)
