include_directories(${PROJECT_BINARY_DIR}/src)

project(splines)

add_executable(test_splines test_splines.f90)
set(LAPACK_LIBS
    lapack
    blas
    )
target_link_libraries(test_splines fortran_utils ${LAPACK_LIBS})

add_test(test_splines ${PROJECT_BINARY_DIR}/test_splines)
