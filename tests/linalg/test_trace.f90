program test_trace
use types, only: dp
use utils, only: assert
use linalg, only: trace
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(3,3)
complex(dp) :: B(3,3)

A = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9], shape=[3,3], order=[2,1])
call assert(abs(trace(A) - 15.0_dp) < eps)

B = i_*reshape([1, 2, 3, 4, 5, 6, 7, 8, 9], shape=[3,3], order=[2,1])
call assert(abs(trace(B) - 15*i_) < eps)

end program
