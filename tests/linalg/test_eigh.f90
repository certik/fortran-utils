program test_eigh
use types, only: dp
use utils, only: assert
use linalg, only: eigh, eye
use constants, only: i_
implicit none

! test eigenvalue computation for symmetric/hermitian matrices:

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(2, 2), B(2, 2), lam(2), c(2, 2), r(2), n
complex(dp) :: AC(2, 2), BC(2, 2), cc(2, 2)
complex(dp) :: rc(2), nc
integer :: i

A = reshape([1, 0, 0, -1], [2, 2])
B = reshape([1, 0, 0, 1], [2, 2])

! test generalized eigenvalue problem with real symmetric matrices:
call eigh(A, B, lam, c)
! Test eigenvalues:
call assert(all(abs(lam - [-1, 1]) < eps))
do i = 1, 2
    ! Test that c(:, i) is an eigenvector with eigenvalue lam(i):
    r = matmul(A-lam(i)*B, c(:, i))
    call assert(sqrt(dot_product(r, r)) < eps)
    ! Test that eigenvectors are properly normalized:
    n = dot_product(c(:, i), matmul(B, c(:, i)))
    call assert(abs(n - 1) < eps)
end do

! test eigenvalue problem with real symmetric matrices:
call eigh(A, lam, c)
! Test eigenvalues:
call assert(all(abs(lam - [-1, 1]) < eps))
do i = 1, 2
    ! Test that c(:, i) is an eigenvector with eigenvalue lam(i):
    r = matmul(A, c(:, i)) - lam(i) * c(:, i)
    call assert(sqrt(dot_product(r, r)) < eps)
    ! Test that eigenvectors are properly normalized:
    n = dot_product(c(:, i), c(:, i))
    call assert(abs(n - 1) < eps)
end do

! another test for the generalized real problem:
A = reshape([2, -4, -4, 2], [2, 2])
B = reshape([2, 1, 1, 2], [2, 2])
call eigh(A, B, lam, c)
! Test eigenvalues:
call assert(all(abs(lam - [-2._dp/3, 6._dp]) < eps))
do i = 1, 2
    ! Test that c(:, i) are eigenvectors:
    r = matmul(A-lam(i)*B, c(:, i))
    call assert(sqrt(dot_product(r, r)) < eps)
    ! Test that eigenvectors are properly normalized:
    n = dot_product(c(:, i), matmul(B, c(:, i)))
    call assert(abs(n - 1) < eps)
end do

! test with complex matrices:
AC = reshape([1.0_dp+0.0_dp*i_, i_, -i_, 3.0_dp+0.0_dp*i_], [2, 2])
call eigh(AC, lam, cc)
do i = 1, 2
    ! Test that c(:, i) are eigenvectors:
    rc = matmul(AC, cc(:, i)) - lam(i) * cc(:, i)
    call assert(sqrt(abs(dot_product(rc, rc))) < eps)
    ! Test that eigenvectors are properly normalized:
    nc = dot_product(cc(:, i), cc(:, i))
    call assert(abs(nc - 1) < eps)
end do

AC = reshape([1.0_dp +0*i_, 3*i_, -3*i_, 2.0_dp-0*i_], shape=[2,2])
BC = reshape([1.0_dp + 0*i_, 1*i_, -1*i_, 2.0_dp+ 0*i_], shape=[2,2])
call eigh(AC, BC, lam, cc)
! Test eigenvalues:
! TODO: add test for correctness eigenvalues
do i = 1, 2
    ! Test that cc(:, i) are eigenvectors:
    rc = matmul(AC-lam(i)*BC, cc(:, i))
    call assert(sqrt(abs(dot_product(rc, rc))) < eps)
end do
! Test that eigenvectors are properly normalized:
BC = matmul(matmul(conjg(transpose(cc)), BC), cc)  ! should be eye(2); Z^H B Z = I
call assert(maxval(abs(BC - eye(2))) < eps)

end program
