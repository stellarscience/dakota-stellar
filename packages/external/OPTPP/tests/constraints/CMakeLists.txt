add_executable(tstbc tstbc.C tstfcn.C)
target_link_libraries(tstbc ${test_libs})
add_test(tstbc tstbc)

add_executable(tstcompound tstcompound.C tstfcn.C)
target_link_libraries(tstcompound ${test_libs})
add_test(tstcompound tstcompound)

add_executable(tstlinear tstlinear.C tstfcn.C)
target_link_libraries(tstlinear ${test_libs})
add_test(tstlinear tstlinear)

add_executable(tstnlp tstnlp.C tstfcn.C tstfcn.h)
target_link_libraries(tstnlp ${test_libs})
add_test(tstnlp tstnlp)

add_executable(tstnonlinear tstnonlinear.C tstfcn.C tstfcn.h)
target_link_libraries(tstnonlinear ${test_libs})
add_test(tstnonlinear tstnonlinear)
