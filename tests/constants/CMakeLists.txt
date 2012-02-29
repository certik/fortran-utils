include_directories(${PROJECT_BINARY_DIR}/src)

project(constants)

add_executable(test_constants test_constants.f90)
target_link_libraries(test_constants fortran_utils)

add_test(test_constants ${PROJECT_BINARY_DIR}/test_constants)
