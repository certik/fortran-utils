include_directories(${PROJECT_BINARY_DIR}/src)

project(strings)

add_executable(test_str test_str.f90)
target_link_libraries(test_str fortran_utils)

add_test(test_str ${PROJECT_BINARY_DIR}/test_str)
