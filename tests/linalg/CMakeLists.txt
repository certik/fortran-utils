include_directories(${PROJECT_BINARY_DIR}/src)

project(linalg)

set(LAPACK_LIBS
    lapack
    blas
    )

macro(ADDTEST name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} fortran_utils ${LAPACK_LIBS})
    add_test(${name} ${PROJECT_BINARY_DIR}/${name})
endmacro(ADDTEST)

ADDTEST(test_eig)
ADDTEST(test_eigh)
ADDTEST(test_eigvals)
ADDTEST(test_eigh_bigger)
ADDTEST(test_inv)
ADDTEST(test_solve)
ADDTEST(test_det)
ADDTEST(test_lstsq)
ADDTEST(test_diag)
ADDTEST(test_trace)
ADDTEST(test_svdvals)
ADDTEST(test_svd)
ADDTEST(test_cholesky)
