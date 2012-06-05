include_directories(${PROJECT_BINARY_DIR}/src)

project(linalg)

set(LAPACK_LIBS
    lapack
    blas
    )

add_executable(test_eig test_eig.f90)
target_link_libraries(test_eig fortran_utils ${LAPACK_LIBS})
add_test(test_eig ${PROJECT_BINARY_DIR}/test_eig)

add_executable(test_inv test_inv.f90)
target_link_libraries(test_inv fortran_utils ${LAPACK_LIBS})
add_test(test_inv ${PROJECT_BINARY_DIR}/test_inv)
