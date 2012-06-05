program test_inv
use types, only: dp
use constants, only: i_
use utils, only: assert
use linalg, only: inv, eye
implicit none

real(dp), parameter :: eps = 1e-9_dp
complex(dp) :: A(2, 2), B(2, 2), C(2, 2)
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

end program
