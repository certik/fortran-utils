include_directories(${PROJECT_BINARY_DIR}/src)

project(ppm)

add_executable(test_ppm test_ppm.f90 convolution.f90)
target_link_libraries(test_ppm fortran_utils)

add_test(test_ppm ${PROJECT_BINARY_DIR}/test_ppm)
