set(SRC
    constants.f90   sorting.f90           types.f90
    utils.f90
    mesh.f90
    ppm.f90
    optimize.f90
    special.f90
    amos.f90
)

if(WITH_LAPACK)
    set(SRC ${SRC} lapack.f90 splines.f90 linalg.f90)
endif()

if(WITH_HDF5)
    set(SRC ${SRC} h5_utils.f90)
    # TODO: make this search path more general somehow:
    # This directory contains the hdf5.mod that we need:
    include_directories($ENV{SPKG_LOCAL}/include)
endif()


add_subdirectory(legacy/amos)

add_library(fortran_utils ${SRC})
target_link_libraries(fortran_utils p::amos)
