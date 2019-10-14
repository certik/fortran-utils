program test_inv
use types, only: dp
use constants, only: i_
use utils, only: assert
use linalg, only: inv, eye
implicit none

real(dp), parameter :: eps = 1e-9_dp
complex(dp) :: A(2, 2), B(2, 2), C(2, 2)
real(dp) :: D(2,2), E(2,2), F(2,2)

! test for complex double routines:
A = reshape([0+0*i_, i_, 2+0*i_, 0+0*i_], [2, 2], order=[2, 1])
C = inv(A)
B = reshape([0+0*i_, 1+0*i_, -2*i_, 0+0*i_], [2, 2], order=[2, 1]) / 2
call assert(maxval(abs(C-B)) < eps)
call assert(maxval(abs(matmul(C, A)-eye(2))) < eps)
call assert(maxval(abs(matmul(A, C)-eye(2))) < eps)

A = reshape([1+0*i_, i_, 2+0*i_, 3*i_], [2, 2], order=[2, 1])
C = inv(A)
B = reshape([3+0*i_, -1+0*i_, 2*i_, -i_], [2, 2], order=[2, 1])
call assert(maxval(abs(C-B)) < eps)
call assert(maxval(abs(matmul(C, A)-eye(2))) < eps)
call assert(maxval(abs(matmul(A, C)-eye(2))) < eps)

A = reshape([0+0*i_, i_, -i_, 0+0*i_], [2, 2], order=[2, 1])
C = inv(A)
call assert(maxval(abs(matmul(C, A)-eye(2))) < eps)
call assert(maxval(abs(matmul(A, C)-eye(2))) < eps)

! tests for double precision routine:
D = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], shape=[2,2])
E = reshape([-2.0_dp, 1.0_dp, 1.5_dp, -0.5_dp], shape=[2,2])
F = inv(D)
call assert(maxval(abs(E - F)) < eps)
call assert(maxval(abs(matmul(F, D) - eye(2))) < eps)
call assert(maxval(abs(matmul(E, D) - eye(2))) < eps)

D = reshape([0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], shape=[2,2])  ! is its own inverse
F = inv(D)
call assert(maxval(abs(D - F)) < eps)
call assert(maxval(abs(matmul(F, D) - eye(2))) < eps)

end program
