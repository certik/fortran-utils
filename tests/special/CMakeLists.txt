include_directories(${PROJECT_BINARY_DIR}/src)

project(special)

add_executable(test_bessel_zeros test_bessel_zeros.f90)
target_link_libraries(test_bessel_zeros fortran_utils)

add_test(test_bessel_zeros ${PROJECT_BINARY_DIR}/test_bessel_zeros)
