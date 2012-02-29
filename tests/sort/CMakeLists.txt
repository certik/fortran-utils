include_directories(${PROJECT_BINARY_DIR}/src)

project(sort)

add_executable(test_sort test_sort.f90)
target_link_libraries(test_sort fortran_utils)

add_executable(test_argsort test_argsort.f90)
target_link_libraries(test_argsort fortran_utils)

add_test(test_sort ${PROJECT_BINARY_DIR}/test_sort)
add_test(test_argsort ${PROJECT_BINARY_DIR}/test_argsort)
