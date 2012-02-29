include_directories(${PROJECT_BINARY_DIR}/src)

project(mesh)

add_executable(test_mesh test_mesh.f90)
target_link_libraries(test_mesh fortran_utils)

add_test(test_mesh ${PROJECT_BINARY_DIR}/test_mesh)
