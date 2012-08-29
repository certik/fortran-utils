program test_diag
use types, only: dp
use utils, only: assert
use linalg, only: diag, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(3,3), B(3,3)
complex(dp) :: C(3,3), D(3,3)

A = diag([1.0_dp, 1.0_dp, 1.0_dp])
call assert(maxval(abs(A-eye(3))) < eps)

A = diag([1.0_dp, 2.0_dp, 3.0_dp])
B = reshape([1,0,0,0,2,0,0,0,3], shape=[3,3])
call assert(maxval(abs(A - B)) < eps)

C = diag([1*i_, 2*i_, 3*i_])
D = i_*reshape([1,0,0,0,2,0,0,0,3], shape=[3,3])
call assert(maxval(abs(C - D)) < eps)
end program
