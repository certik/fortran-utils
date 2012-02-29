include_directories(${PROJECT_BINARY_DIR}/src)

project(linalg)

add_executable(test_eig test_eig.f90)
set(LAPACK_LIBS
    lapack
    blas
    )
target_link_libraries(test_eig fortran_utils ${LAPACK_LIBS})

add_test(test_eig ${PROJECT_BINARY_DIR}/test_eig)
