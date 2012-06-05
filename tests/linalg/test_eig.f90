program test_eig
use types, only: dp
use utils, only: assert
use linalg, only: eigh
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(2, 2), B(2, 2), lam(2), c(2, 2), r(2), n
integer :: i
A = reshape([1, 0, 0, -1], [2, 2])
B = reshape([1, 0, 0, 1], [2, 2])
call eigh(A, B, lam, c)
! Test eigenvalues:
call assert(all(abs(lam - [-1, 1]) < eps))
do i = 1, 2
    ! Test that c(:, i) are eigenvectors:
    r = matmul(A-lam(i)*B, c(:, i))
    call assert(sqrt(dot_product(r, r)) < eps)
    ! Test that eigenvectors are properly normalized:
    n = dot_product(c(:, i), matmul(B, c(:, i)))
    call assert(abs(n - 1) < eps)
end do

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

end program
