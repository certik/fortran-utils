cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_SOURCE_DIR}/cmake/UserOverride.cmake)

enable_language(Fortran)

project(fortranutils)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(WITH_LAPACK yes
    CACHE BOOL "Build with LAPACK support")
set(WITH_HDF5 no
    CACHE BOOL "Build with HDF5 support")

# Make sure that CMAKE_BUILD_TYPE is either Debug or Release:
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug
        CACHE STRING "Build type (Debug, Release)" FORCE)
endif ()
if (NOT (CMAKE_BUILD_TYPE STREQUAL "Debug" OR
        CMAKE_BUILD_TYPE STREQUAL "Release"))
    message("${CMAKE_BUILD_TYPE}")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be one of: Debug, Release (current value: '${CMAKE_BUILD_TYPE}')")
endif ()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -dumpversion
                        OUTPUT_VARIABLE GCC_VERSION)
    if (GCC_VERSION VERSION_LESS 4.5)
        set(common "-std=f2008 -Wall -Wextra -Wimplicit-interface -fPIC -Werror -fmax-errors=1")
        set(CMAKE_Fortran_FLAGS_RELEASE "${common} -O3 -march=native -ffast-math -funroll-loops")
        set(CMAKE_Fortran_FLAGS_DEBUG   "${common} -g -fbounds-check -fcheck-array-temporaries -fbacktrace")
    endif()

    # gfortran
    # Enable this if you want to check for single/double corruption (and use
    # the Debug build):
    #set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fdefault-real-8")
endif ()

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)

message("\n")
message("Configuration results")
message("---------------------")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message("Build type: ${CMAKE_BUILD_TYPE}")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_DEBUG}")
else ()
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_RELEASE}")
endif ()
message("Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message("With LAPACK: ${WITH_LAPACK}")
message("With HDF5: ${WITH_HDF5}")
