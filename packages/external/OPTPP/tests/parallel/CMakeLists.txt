add_executable(tstfdnewtpds tstfdnewtpds.C SetupTest.C tstfcn.C tstfcn.h)
target_link_libraries(tstfdnewtpds ${test_libs})
add_test(tstfdnewtpds tstfdnewtpds)

add_executable(tstnewtpds tstnewtpds.C SetupTest.C tstfcn.C tstfcn.h)
target_link_libraries(tstnewtpds ${test_libs})
add_test(tstnewtpds tstnewtpds)

add_executable(tstpds tstpds.C SetupTest.C tstfcn.C tstfcn.h)
target_link_libraries(tstpds ${test_libs})
add_test(tstpds tstpds)

add_executable(tsttrpds tsttrpds.C SetupTest.C tstfcn.C tstfcn.h)
target_link_libraries(tsttrpds ${test_libs})
add_test(tsttrtpds tsttrpds)

add_executable(tstGSS tstGSS.C SetupTest.C tstfcn.C tstfcn.h)
target_link_libraries(tstGSS ${test_libs})
add_test(tstGSS tstGSS)
