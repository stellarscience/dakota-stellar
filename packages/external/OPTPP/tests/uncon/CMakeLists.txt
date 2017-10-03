add_executable(tstqnewton tstqnewton.C rosen.C tstfcn.h)
target_link_libraries(tstqnewton ${test_libs})
add_test(tstqnewton tstqnewton)

add_executable(tstgnewton tstgnewton.C rosen.C tstfcn.h)
target_link_libraries(tstgnewton ${test_libs})
add_test(tstgnewton tstgnewton)

add_executable(tstnewton tstnewton.C rosen.C tstfcn.h)
target_link_libraries(tstnewton ${test_libs})
add_test(tstnewton tstnewton)

add_executable(tstfdnlf1 tstfdnlf1.C rosen.C tstfcn.h)
target_link_libraries(tstfdnlf1 ${test_libs})
add_test(tstfdnlf1 tstfdnlf1)

add_executable(tstcg tstcg.C rosen.C tstfcn.h)
target_link_libraries(tstcg ${test_libs})
add_test(tstcg tstcg)

add_executable(tstLBFGS tstLBFGS.C rosen.C tstfcn.h)
target_link_libraries(tstLBFGS ${test_libs})
add_test(tstLBFGS tstLBFGS)
