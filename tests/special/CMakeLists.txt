include_directories(${PROJECT_BINARY_DIR}/src)

project(special)

add_executable(test_bessel_zeros test_bessel_zeros.f90)
target_link_libraries(test_bessel_zeros fortran_utils)
add_test(test_bessel_zeros ${PROJECT_BINARY_DIR}/test_bessel_zeros)

add_executable(test_besselj test_besselj.f90)
target_link_libraries(test_besselj fortran_utils)
add_test(test_besselj ${PROJECT_BINARY_DIR}/test_besselj)

add_executable(test_bessely test_bessely.f90)
target_link_libraries(test_bessely fortran_utils)
add_test(test_bessely ${PROJECT_BINARY_DIR}/test_bessely)

add_executable(test_hankel1 test_hankel1.f90)
target_link_libraries(test_hankel1 fortran_utils)
add_test(test_hankel1 ${PROJECT_BINARY_DIR}/test_hankel1)

add_executable(test_hankel2 test_hankel2.f90)
target_link_libraries(test_hankel2 fortran_utils)
add_test(test_hankel2 ${PROJECT_BINARY_DIR}/test_hankel2)

add_executable(test_besseli test_besseli.f90)
target_link_libraries(test_besseli fortran_utils)
add_test(test_besseli ${PROJECT_BINARY_DIR}/test_besseli)

add_executable(test_besselk test_besselk.f90)
target_link_libraries(test_besselk fortran_utils)
add_test(test_besselk ${PROJECT_BINARY_DIR}/test_besselk)

add_executable(test_airyai test_airyai.f90)
target_link_libraries(test_airyai fortran_utils)
add_test(test_airyai ${PROJECT_BINARY_DIR}/test_airyai)

add_executable(test_dairyai test_dairyai.f90)
target_link_libraries(test_dairyai fortran_utils)
add_test(test_dairyai ${PROJECT_BINARY_DIR}/test_dairyai)

add_executable(test_airybi test_airybi.f90)
target_link_libraries(test_airybi fortran_utils)
add_test(test_airybi ${PROJECT_BINARY_DIR}/test_airybi)

add_executable(test_dairybi test_dairybi.f90)
target_link_libraries(test_dairybi fortran_utils)
add_test(test_dairybi ${PROJECT_BINARY_DIR}/test_dairybi)

add_executable(test_spherical_bessel_yn test_spherical_bessel_yn.f90)
target_link_libraries(test_spherical_bessel_yn fortran_utils)
add_test(test_spherical_bessel_yn ${PROJECT_BINARY_DIR}/test_spherical_bessel_yn)
