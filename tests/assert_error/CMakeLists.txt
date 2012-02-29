include_directories(${PROJECT_BINARY_DIR}/src)

project(assert_error)

add_executable(test_error test_error.f90)
target_link_libraries(test_error fortran_utils)

add_executable(test_assert_pass test_assert_pass.f90)
target_link_libraries(test_assert_pass fortran_utils)

add_executable(test_assert_fail1 test_assert_fail1.f90)
target_link_libraries(test_assert_fail1 fortran_utils)

add_executable(test_assert_fail2 test_assert_fail2.f90)
target_link_libraries(test_assert_fail2 fortran_utils)

add_test(test_error ${PROJECT_BINARY_DIR}/test_error)
set_tests_properties(test_error PROPERTIES WILL_FAIL TRUE)

add_test(test_assert_pass ${PROJECT_BINARY_DIR}/test_assert_pass)

add_test(test_assert_fail1 ${PROJECT_BINARY_DIR}/test_assert_fail1)
set_tests_properties(test_assert_fail1 PROPERTIES WILL_FAIL TRUE)

add_test(test_assert_fail2 ${PROJECT_BINARY_DIR}/test_assert_fail2)
set_tests_properties(test_assert_fail2 PROPERTIES WILL_FAIL TRUE)
