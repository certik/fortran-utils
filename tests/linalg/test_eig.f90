program test_eig
use types, only: dp
use utils, only: assert
use linalg, only: eig, eye
use constants, only : i_
implicit none

! test eigenvalue comutation for general matrices:

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(2, 2), B(5, 5)
complex(dp) :: AC(2, 2), lam(2), c(2, 2), r(2), n, lamb(5), cb(5, 5)
integer :: i

! test a matrix with complex eigenvalues/eigenvectors:
A = reshape([1, -1, 2, 1], shape=[2, 2])  ! lambda_i = 1 \pm \sqrt(2) i
call eig(A, lam, c)
! TODO: add test for correctness of eigenvalues
do i = 1, 2
    ! Test proper norm of eigenvectors:
    n = dot_product(c(:,i), c(:,i))
    call assert(abs(n - 1) < eps)
    ! Test that c(:, i) is an eigenvector with eigenvalue lam(i):
    r = matmul(A-lam(i)*eye(2), c(:, i))
    call assert(sqrt(abs(dot_product(r, r))) < eps)
end do

! test a multiple of the unit matrix:
B = 3*eye(5)
call eig(B, lamb, cb)
call assert(maxval(abs(lamb - 3.0_dp)) < eps)  ! all eigenvalues are 3
call assert(maxval(abs(cb - eye(5))) < eps)  ! eigenvectors are cartesian unit basis vectors

! test complex matrices:
AC = reshape([1.0_dp+0*i_, 2*i_, 3*i_, -4.0_dp+0*i_], shape=[2,2])
call eig(AC, lam, c)
call assert(all(abs(lam - [-1.0_dp, -2.0_dp]) < eps))
do i = 1, 2
    ! Test proper norm of eigenvectors:
    n = dot_product(c(:,i), c(:,i))
    call assert(abs(n - 1) < eps)
    ! Test that c(:, i) is an eigenvector with eigenvalue lam(i):
    r = matmul(AC-lam(i)*eye(2), c(:, i))
    call assert(sqrt(abs(dot_product(r, r))) < eps)
end do

end program
