set(test_libs optpp ${Teuchos_LIBRARIES} ${TEUCHOS_LIBRARY}
    ${BLAS_LIBS} ${LAPACK_LIBS}) 

if(npsol_lib)
  list(APPEND test_libs ${npsol_lib}) 
endif()
if(OPTPP_HAVE_MPI AND MPI_FOUND)
  list(APPEND test_libs ${MPI_LIBRARIES}) 
endif()


add_subdirectory(constraints)
add_subdirectory(hock)
add_subdirectory(parallel)
add_subdirectory(uncon)

if(npsol_lib)
  add_subdirectory(npsol)
endif()
