include_directories(${PROJECT_BINARY_DIR}/src)

project(optimize)

add_executable(test_bisect test_bisect.f90)
add_executable(test_secant test_secant.f90)
target_link_libraries(test_bisect fortran_utils)
target_link_libraries(test_secant fortran_utils)
add_test(test_bisect ${PROJECT_BINARY_DIR}/test_bisect)
add_test(test_secant ${PROJECT_BINARY_DIR}/test_secant)
