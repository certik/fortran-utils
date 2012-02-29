
include_directories(${PROJECT_BINARY_DIR}/src)

project(test_h5_utils)

add_executable(test_h5_utils test_h5_utils.f90)
set(HDF5_LIBS
    hdf5hl_fortran
    hdf5_hl
    hdf5_fortran
    hdf5
    )
target_link_libraries(test_h5_utils fortran_utils ${HDF5_LIBS})
add_test(test_h5_utils ${PROJECT_BINARY_DIR}/test_h5_utils)
